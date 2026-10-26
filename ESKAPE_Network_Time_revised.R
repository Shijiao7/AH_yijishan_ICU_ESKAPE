# ============================================================
# ICU Metagenomics — Ultimate Pipeline v6 (High-Intensity Filter)
# 核心逻辑：全屋空间聚合 -> 高强度去噪过滤 -> 靶向Pearson -> Cytoscape导出
# ============================================================

# ── 0. 加载必要包 ───────────────────────────────────────────
library(SpiecEasi)
library(tidyverse)
library(igraph)
library(ggraph)
library(ggrepel)

# ── 1. 读取原始数据与病原体数据库 ───────────────────────────
cat("\n=== 1. 加载数据与病原体数据库 ===\n")
otu  <- read.csv("/home/bio/Documents/AH_ICU/Abundance_shortReads/sample_all.csv", row.names = 1, check.names = FALSE)
meta <- read_tsv("/home/bio/Documents/AH_ICU/R/metadata_AHICU_HW_match.txt", show_col_types = FALSE)

bact_db <- read.csv("/home/bio/Documents/ICU_Sample/data/pathogenDB/pathegon_bacteriaTaxa_list.csv", check.names = FALSE) %>% select(`Pathogen name`, `Biosafety level`)
fung_db <- read.csv("/home/bio/Documents/ICU_Sample/data/pathogenDB/pathegon_fungiTaxa_list.csv", check.names = FALSE) %>% select(`Pathogen name`, `Biosafety level`)
pathogen_db <- bind_rows(bact_db, fung_db) %>%
  filter(`Biosafety level` %in% c("BSL-2", "BSL-3")) %>%
  rename(Species = `Pathogen name`, BSL = `Biosafety level`) %>%
  distinct(Species, .keep_all = TRUE)

ESKAPE_species <- c(
  "Enterococcus faecium", "Staphylococcus aureus",
  "Klebsiella pneumoniae", "Acinetobacter baumannii",
  "Pseudomonas aeruginosa", "Enterobacter cloacae",
  "Enterobacter cloacae complex sp."
)

# ── 2. 空间聚合 (Spatial Pooling) ───────────────────────────
cat("\n=== 2. 构建全屋聚合超级样本 (Meta-Samples) ===\n")
time_order_map <- list(
  Room16       = c("Room16-In", "R-W1", "R-W2", "R-W3", "R-W4", "Room16-Out"),
  Room19       = c("Room19-In", "R-W3", "R-W4", "R-W5", "Room19-Out"),
  Public       = c("R-W1", "R-W2", "R-W3", "R-W4", "R-W5"),
  Healthworker = c("R-W1", "R-W2", "R-W3", "R-W4", "R-W5")
)

pooled_counts_list <- list()
for (site in names(time_order_map)) {
  t_levels <- time_order_map[[site]]
  site_meta <- meta %>% filter(Sites == site, Date %in% t_levels)
  if(nrow(site_meta) == 0) next
  
  meta_mat <- matrix(0, nrow = nrow(otu), ncol = length(t_levels), dimnames = list(rownames(otu), t_levels))
  
  for (tp in t_levels) {
    samples_in_tp <- site_meta %>% filter(Date == tp) %>% pull(Sample_ID)
    samples_exist <- intersect(samples_in_tp, colnames(otu))
    if (length(samples_exist) == 1) {
      meta_mat[, tp] <- otu[, samples_exist]
    } else if (length(samples_exist) > 1) {
      meta_mat[, tp] <- rowSums(otu[, samples_exist]) # 同一天全屋负荷相加
    }
  }
  meta_mat <- meta_mat[rowSums(meta_mat) > 0, , drop = FALSE]
  pooled_counts_list[[site]] <- meta_mat
  cat("  Site:", site, "聚合完成 | 包含", length(t_levels), "个超级时间点\n")
}

# ── 3. CLR 去组成性转换 (用于 LSA 时序分析) ─────────────────
cat("\n=== 3. 准备 LSA 连续时序数据 (CLR 转换) ===\n")
pooled_clr_list <- list()
for (site in names(pooled_counts_list)) {
  mat <- pooled_counts_list[[site]]
  ra_mat <- sweep(mat, 2, colSums(mat), "/")
  
  # 基础过滤：去除极低丰度噪音
  keep <- (rowSums(mat > 0) >= 2) & (rowMeans(ra_mat) >= 1e-4) | (rownames(mat) %in% ESKAPE_species)
  mat_filtered <- mat[keep, , drop = FALSE]
  
  mat_pseudo <- mat_filtered + 0.5
  prop <- sweep(mat_pseudo, 2, colSums(mat_pseudo), "/")
  clr_mat <- apply(prop, 2, function(x) log(x) - mean(log(x)))
  rownames(clr_mat) <- rownames(mat_filtered)
  pooled_clr_list[[site]] <- clr_mat
}

# ── 4. 靶向时序网络分析 (LSA - Pearson) ─────────────────────
cat("\n=== 4. 运行靶向时序相关性分析 (Pearson) ===\n")
cor_with_lag_pearson <- function(x, y, max_lag = 2, n_perm = 999, cor_threshold = 0.5) {
  n <- length(x); best <- list(cor = 0, lag = 0, pval = 1)
  for (lag in -max_lag:max_lag) {
    if      (lag == 0) { xi <- x;             yi <- y }
    else if (lag > 0)  { xi <- x[1:(n-lag)];  yi <- y[(lag+1):n] }
    else               { xi <- x[(-lag+1):n]; yi <- y[1:(n+lag)] }
    
    if (length(xi) < 4) next
    obs_r <- cor(xi, yi, method = "pearson") 
    if (is.na(obs_r)) next
    
    if (abs(obs_r) < cor_threshold) { pval <- 1 } else {
      perm_r <- replicate(n_perm, cor(sample(xi), yi, method = "pearson"))
      pval <- (sum(abs(perm_r) >= abs(obs_r), na.rm = TRUE) + 1) / (n_perm + 1)
    }
    if (abs(obs_r) > abs(best$cor)) best <- list(cor = obs_r, lag = lag, pval = pval)
  }
  return(best)
}

pooled_lsa_results <- list()
for (site in names(pooled_clr_list)) {
  clr_mat <- pooled_clr_list[[site]]
  sp_in <- rownames(clr_mat)
  targets <- intersect(ESKAPE_species, sp_in)
  
  results <- list(); idx <- 1
  for (t_sp in targets) {
    for (bg_sp in sp_in) {
      if (t_sp == bg_sp) next 
      res <- cor_with_lag_pearson(as.numeric(clr_mat[t_sp, ]), as.numeric(clr_mat[bg_sp, ]), cor_threshold = 0.5)
      if(abs(res$cor) >= 0.5) {
        results[[idx]] <- tibble(sp1 = t_sp, sp2 = bg_sp, cor_clr = res$cor, lag = res$lag, pval_perm = res$pval, direction = ifelse(res$cor > 0, "positive", "negative"))
        idx <- idx + 1
      }
    }
  }
  if(length(results) > 0) {
    df <- bind_rows(results) %>% arrange(desc(abs(cor_clr))) %>% filter(pval_perm < 0.05) 
    df$Site <- site
    pooled_lsa_results[[site]] <- df
    cat("  Site:", site, "LSA 找到显著宏观边数:", nrow(df), "\n")
  }
}

# ── 5. 高强度过滤版 SPIEC-EASI ──────────────────────────────
cat("\n=== 5. 运行 SPIEC-EASI (高强度过滤版，极速运行) ===\n")
spiec_results <- list()

for (site in names(pooled_counts_list)) {
  mat_raw <- pooled_counts_list[[site]]
  
  # 【核心高强度过滤逻辑】
  ra_tmp <- sweep(mat_raw, 2, colSums(mat_raw), "/")
  # 要求：在 >20% 的样本中出现，且平均相对丰度 > 0.01% (1e-4)
  keep_base <- (rowMeans(mat_raw > 0) >= 0.20) & (rowMeans(ra_tmp) >= 1e-4)
  # 强制保留 ESKAPE 靶标
  keep_force <- rownames(mat_raw) %in% ESKAPE_species
  keep_final <- keep_base | keep_force
  
  mat <- t(mat_raw[keep_final, , drop = FALSE])
  
  if (nrow(mat) < 5 || ncol(mat) < 3) next
  cat("  SPIEC-EASI | Site:", site, "→ 已过滤至 [", nrow(mat), "samples ×", ncol(mat), "species ]\n")
  
  # 降低重采样次数 (rep.num = 20) 以加速
  se <- spiec.easi(mat, method = "mb", lambda.min.ratio = 1e-2, nlambda = 20, pulsar.params = list(rep.num = 20, seed = 42))
  beta <- as.matrix(symBeta(getOptBeta(se), mode = "ave"))
  sp_names <- colnames(mat)
  rownames(beta) <- colnames(beta) <- sp_names
  
  edge_idx <- which(abs(beta) > 0 & upper.tri(beta), arr.ind = TRUE)
  if (nrow(edge_idx) > 0) {
    spiec_results[[site]] <- tibble(
      sp1 = sp_names[edge_idx[, 1]], sp2 = sp_names[edge_idx[, 2]],
      beta_coef = beta[edge_idx], direction = ifelse(beta[edge_idx] > 0, "positive", "negative"), Site = site
    )
  }
}



# ── 6. 生成网络图并导出为 Cytoscape 格式 ────────────────────
cat("\n=== 6. 生成高颜网络并导出节点/边文件 ===\n")
out_dir <- "/home/bio/Documents/AH_ICU/Abundance_shortReads/Network_Outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

build_cytoscape_network <- function(lsa_df, spiec_df, site_name, cor_cutoff = 0.6) {
  
  # --- 1. 提取边数据 (已加入 P_value) ---
  if (!is.null(lsa_df) && nrow(lsa_df) > 0) {
    edges <- lsa_df %>% filter(abs(cor_clr) >= cor_cutoff) %>%
      transmute(Source = sp1, Target = sp2, 
                Weight = cor_clr, abs_Weight = abs(cor_clr), 
                Direction = direction, Algorithm = "LSA", 
                Lag = lag, 
                P_value = pval_perm) # <--- 【新增 P_value】
  } else if (!is.null(spiec_df) && nrow(spiec_df) > 0) {
    edges <- spiec_df %>%
      transmute(Source = sp1, Target = sp2, 
                Weight = beta_coef, abs_Weight = abs(beta_coef), 
                Direction = direction, Algorithm = "SPIEC-EASI", 
                Lag = NA_integer_, 
                P_value = NA_real_)  # <--- 【SPIEC-EASI无传统P值，设为NA】
  } else { return(NULL) }
  
  edges <- edges %>% filter(Source %in% ESKAPE_species | Target %in% ESKAPE_species)
  if (nrow(edges) == 0) return(NULL)
  
  # --- 2. 构建节点属性 ---
  unique_nodes <- unique(c(edges$Source, edges$Target))
  nodes <- tibble(Id = unique_nodes, Label = unique_nodes) %>%
    left_join(pathogen_db, by = c("Id" = "Species")) %>%
    mutate(
      BSL = replace_na(BSL, "None"),
      Group = case_when(Id %in% ESKAPE_species ~ "ESKAPE", BSL == "BSL-3" ~ "Pathogen (BSL-3)", BSL == "BSL-2" ~ "Pathogen (BSL-2)", TRUE ~ "Others")
    )
  
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  
  # --- 3. 计算 ESKAPE 连接数与层级 Layer ---
  adj_list <- adjacent_vertices(g, V(g))
  nodes$eskape_connect_count <- sapply(adj_list, function(neighbors) {
    sum(V(g)$name[neighbors] %in% ESKAPE_species)
  })
  
  nodes <- nodes %>%
    mutate(Layer = case_when(
      Group != "ESKAPE" & eskape_connect_count >= 2 ~ 1, # 最内圈：多重关联枢纽
      Group == "ESKAPE" ~ 2,                             # 第二圈：ESKAPE
      Group != "ESKAPE" & eskape_connect_count == 1 ~ 3, # 第三圈：边缘关联
      TRUE ~ 4
    ))
  
  # 强制标注所有非 Others 的名字
  #nodes$Label_Show <- ifelse(nodes$Group != "Others", nodes$Id, "")
  nodes$Label_Show <- nodes$Id
  # 同步属性到图对象
  V(g)$Group <- nodes$Group
  V(g)$Degree <- degree(g)
  V(g)$Layer <- nodes$Layer
  V(g)$Label_Show <- nodes$Label_Show
  
  # --- 4. 导出 CSV (供 Cytoscape 使用) ---
  write_csv(edges, file.path(out_dir, paste0(site_name, "_Pooled_edges.csv")))
  write_csv(nodes, file.path(out_dir, paste0(site_name, "_Pooled_nodes.csv")))
  
  # --- 5. R 内部可视化调整 ---
  color_map <- c("ESKAPE" = "#E63946", "Pathogen (BSL-3)" = "#9B5DE5", "Pathogen (BSL-2)" = "#F4A261", "Others" = "#D3D3D3")
  
  set.seed(42)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = abs_Weight, color = Direction), alpha = 0.6) +
    scale_edge_color_manual(values = c("positive" = "#E24B4A", "negative" = "#378ADD")) +
    scale_edge_width(range = c(0.5, 2.5), guide = "none") +
    geom_node_point(aes(size = Degree, fill = Group), shape = 21, color = "white", stroke = 0.5) +
    scale_fill_manual(values = color_map) +
    scale_size(range = c(2, 12), guide = "none") +
    geom_node_text(aes(label = Label_Show), 
                   repel = TRUE, 
                   size = 2.8, 
                   fontface = "italic",
                   max.overlaps = Inf, # 允许无限重叠，强制显示所有标签
                   box.padding = 0.5, 
                   point.padding = 0.5) +
    labs(title = paste("Pathogen-Centric Network:", site_name), 
         subtitle = "Inner: Hub Pathogens | Mid: ESKAPE | Outer: Others") +
    theme_graph()
  
  print(p)
}

# 批量执行
for (site in names(pooled_counts_list)) {
  lsa_df <- if(site %in% names(pooled_lsa_results)) pooled_lsa_results[[site]] else NULL
  spiec_df <- if(site %in% names(spiec_results)) spiec_results[[site]] else NULL
  build_cytoscape_network(lsa_df, spiec_df, site, cor_cutoff = 0.6) 
}

cat("\n=== 🎉 终极分析已完成！请前往 Network_Outputs 文件夹提取 CSV 导入 Cytoscape！ ===\n")

#############################3
library(tidyverse)

# 1. 定义你的输出文件夹路径
out_dir <- "/home/bio/Documents/AH_ICU/Abundance_shortReads/Network_Outputs"

# 2. 自动寻找目录下所有以 _Pooled_edges.csv 结尾的文件
edge_files <- list.files(path = out_dir, pattern = "_Pooled_edges\\.csv$", full.names = TRUE)

cat("找到以下待合并的文件：\n")
print(basename(edge_files))

# 3. 批量读取、提取 Site 名称并合并
all_edges <- map_df(edge_files, function(file) {
  # 从文件名中提取 Site 名字 (例如把 "Room16_Pooled_edges.csv" 变成 "Room16")
  site_name <- basename(file) %>% str_replace("_Pooled_edges\\.csv$", "")
  
  # 读取 CSV，并强制新增一列 Site，放在第一列方便查看
  read_csv(file, show_col_types = FALSE) %>%
    mutate(Site = site_name) %>%
    relocate(Site, .before = Source)
})

# 4. 保存合并后的完整大表
merged_file <- file.path(out_dir, "AllSites_Merged_edges.csv")
write_csv(all_edges, merged_file)
cat("\n✅ 合并完成！总边数:", nrow(all_edges), "\n")
cat("💾 合并文件已保存至:", merged_file, "\n")


# ============================================================
# 💡 进阶彩蛋：一键找出在 4 个 Site 中都存在的“超级核心关系”
# ============================================================

cat("\n=== 正在分析跨 Site 普遍存在的网络边 ===\n")

core_edges <- all_edges %>%
  # 【关键技巧】：网络是无向的，Site1可能是 A-B，Site2可能是 B-A。
  # 为了准确统计，我们将两端的节点按字母顺序排序标准化为 Node1 和 Node2
  mutate(
    Node1 = pmin(Source, Target),
    Node2 = pmax(Source, Target)
  ) %>%
  group_by(Node1, Node2) %>%
  summarise(
    Site_Count = n_distinct(Site),               # 统计这条边在几个不同的 Site 出现过
    Sites_Present = paste(unique(Site), collapse = " | "), # 列出具体的 Site 名字
    Avg_Weight = mean(Weight),                   # 计算平均相关系数
    Avg_Lag = mean(Lag, na.rm = TRUE),           # 计算平均时滞
    .groups = "drop"
  ) %>%
  arrange(desc(Site_Count), desc(abs(Avg_Weight)))

# 提取在所有 4 个 Site (或者至少 3 个 Site) 中都存在的边
universal_edges <- core_edges %>% filter(Site_Count >= 3) 

if (nrow(universal_edges) > 0) {
  cat("🔥 找到了", nrow(universal_edges), "对在 >= 3 个 Site 中普遍存在的核心关系！\n")
  print(head(universal_edges, 10))
  
  # 保存核心保守网络
  write_csv(universal_edges, file.path(out_dir, "Conserved_Core_Edges.csv"))
} else {
  cat("⚠️ 没有找到在 >= 3 个 Site 中同时存在的边。这说明各个房间/区域的微生态具有极高的异质性（特异性）。\n")
}




# ============================================================
# 🎯 专项提取：Room16 与 Room19 共同保守的微生态关系
# ============================================================

cat("\n=== 正在提取 Room16 与 Room19 共同存在的核心关系 ===\n")

# 假设当前环境中已经存在 all_edges 数据框。
# 如果你重启了 R，请先运行这一行重新读取：
# all_edges <- read_csv("/home/bio/Documents/AH_ICU/Abundance_shortReads/Network_Outputs/AllSites_Merged_edges.csv", show_col_types = FALSE)

room16_19_shared <- all_edges %>%
  # 1. 仅保留这两个病房的数据
  filter(Site %in% c("Room16", "Room19")) %>%
  # 2. 标准化节点方向 (解决 A-B 和 B-A 的反向匹配问题)
  mutate(
    Node1 = pmin(Source, Target),
    Node2 = pmax(Source, Target)
  ) %>%
  # 3. 分组统计
  group_by(Node1, Node2) %>%
  summarise(
    Site_Count = n_distinct(Site),
    Avg_Weight = mean(Weight),
    Avg_abs_Weight = mean(abs(Weight)),
    # 巧妙技巧：把 Room16 和 Room19 的具体相关系数单独拆成两列，方便你直观对比
    Weight_Room16 = sum(Weight[Site == "Room16"], na.rm = TRUE),
    Weight_Room19 = sum(Weight[Site == "Room19"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # 4. 必须在两个病房都出现过 (Site_Count == 2)
  filter(Site_Count == 2) %>%
  arrange(desc(Avg_abs_Weight))

# 5. 打印结果与导出
if (nrow(room16_19_shared) > 0) {
  cat("🔥 找到了", nrow(room16_19_shared), "对在 Room16 和 Room19 中均存在的核心关系！\n")
  print(head(room16_19_shared, 10))
  
  out_file <- "/home/bio/Documents/AH_ICU/Abundance_shortReads/Network_Outputs/Room16_19_Shared_Edges.csv"
  write_csv(room16_19_shared, out_file)
  cat("\n💾 结果已保存至:", out_file, "\n")
} else {
  cat("⚠️ 没有找到同时存在于 Room16 和 Room19 的边。\n")
}

###########################################3
# ============================================================
# 绘制跨 Site 保守核心网络图 (Shared & Conserved Networks)
# 适用文件: Room16_19_Shared_Edges.csv & Conserved_Core_Edges.csv
# ============================================================
library(tidyverse)
library(igraph)
library(ggraph)
library(ggrepel)

# 确保目标文件夹路径正确
out_dir <- "/home/bio/Documents/AH_ICU/Abundance_shortReads/Network_Outputs"

# 专门为汇总表编写的绘图函数
plot_conserved_network <- function(csv_file, plot_title, plot_subtitle) {
  
  if (!file.exists(csv_file)) {
    cat("⚠️ 文件不存在，请检查路径:", csv_file, "\n")
    return(NULL)
  }
  
  # 1. 读取汇总边数据，并统一列名以便与之前的逻辑兼容
  edges_raw <- read_csv(csv_file, show_col_types = FALSE)
  
  if (nrow(edges_raw) == 0) {
    cat("⚠️ 文件中无数据:", csv_file, "\n")
    return(NULL)
  }
  
  edges <- edges_raw %>%
    rename(Source = Node1, Target = Node2, Weight = Avg_Weight) %>%
    mutate(
      abs_Weight = abs(Weight),
      Direction = ifelse(Weight > 0, "positive", "negative")
    )
  
  # 2. 依然坚持 Pathogen-Centric 原则：只看与 ESKAPE 相关的边
  edges <- edges %>% filter(Source %in% ESKAPE_species | Target %in% ESKAPE_species)
  if (nrow(edges) == 0) {
    cat("⚠️ 过滤 ESKAPE 后，", plot_title, " 无边可画。\n")
    return(NULL)
  }
  
  # 3. 构建节点属性 (匹配 BSL 数据库)
  unique_nodes <- unique(c(edges$Source, edges$Target))
  nodes <- tibble(Id = unique_nodes, Label = unique_nodes) %>%
    left_join(pathogen_db, by = c("Id" = "Species")) %>%
    mutate(
      BSL = replace_na(BSL, "None"),
      Group = case_when(
        Id %in% ESKAPE_species ~ "ESKAPE", 
        BSL == "BSL-3" ~ "Pathogen (BSL-3)", 
        BSL == "BSL-2" ~ "Pathogen (BSL-2)", 
        TRUE ~ "Others"
      )
    )
  
  g <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
  
  # 4. 计算层级 (Layer) 和向心力
  adj_list <- adjacent_vertices(g, V(g))
  nodes$eskape_connect_count <- sapply(adj_list, function(neighbors) {
    sum(V(g)$name[neighbors] %in% ESKAPE_species)
  })
  
  nodes <- nodes %>%
    mutate(Layer = case_when(
      Group != "ESKAPE" & eskape_connect_count >= 2 ~ 1, # 核心枢纽
      Group == "ESKAPE" ~ 2,                             # ESKAPE本尊
      Group != "ESKAPE" & eskape_connect_count == 1 ~ 3, # 边缘菌群
      TRUE ~ 4
    ))
  
  # 强制显示所有 BSL 和 ESKAPE 的名字
  #nodes$Label_Show <- ifelse(nodes$Group != "Others", nodes$Id, "")
  nodes$Label_Show <- nodes$Id
  # 同步属性
  V(g)$Group <- nodes$Group
  V(g)$Degree <- degree(g)
  V(g)$Layer <- nodes$Layer
  V(g)$Label_Show <- nodes$Label_Show
  
  # 5. R 内部可视化渲染 (FR 布局)
  color_map <- c("ESKAPE" = "#E63946", "Pathogen (BSL-3)" = "#9B5DE5", "Pathogen (BSL-2)" = "#F4A261", "Others" = "#D3D3D3")
  
  set.seed(42)
  p <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(width = abs_Weight, color = Direction), alpha = 0.6) +
    scale_edge_color_manual(values = c("positive" = "#E24B4A", "negative" = "#378ADD")) +
    scale_edge_width(range = c(0.8, 3.0), guide = "none") +  # 线条稍微加粗一点，因为是核心边
    geom_node_point(aes(size = Degree, fill = Group), shape = 21, color = "white", stroke = 0.5) +
    scale_fill_manual(values = color_map) +
    scale_size(range = c(3, 14), guide = "none") +
    geom_node_text(aes(label = Label_Show), 
                   repel = TRUE, 
                   size = 3.2, 
                   fontface = "italic",
                   max.overlaps = Inf, 
                   box.padding = 0.6, 
                   point.padding = 0.5) +
    labs(title = plot_title, subtitle = plot_subtitle) +
    theme_graph() +
    theme(legend.position = "right")
  
  print(p)
  
  # 如果你想再次导入 Cytoscape 做完美的同心圆，我把更新 Layer 后的数据也存下来
  prefix <- str_replace(basename(csv_file), "\\.csv", "")
  write_csv(edges, file.path(out_dir, paste0(prefix, "_CytoEdges.csv")))
  write_csv(nodes, file.path(out_dir, paste0(prefix, "_CytoNodes.csv")))
  
  return(p)
}

# ============================================================
# 执行绘图
# ============================================================

cat("\n=== 正在生成 Room16 & 19 共有保守网络 ===\n")
p1 <- plot_conserved_network(
  csv_file = file.path(out_dir, "Room16_19_Shared_Edges.csv"),
  plot_title = "Patient Rooms Shared Network (Room 16 & 19)",
  plot_subtitle = "Edges strictly conserved across both patient rooms"
)

cat("\n=== 正在生成全局核心保守网络 (>=3 Sites) ===\n")
p2 <- plot_conserved_network(
  csv_file = file.path(out_dir, "Conserved_Core_Edges.csv"),
  plot_title = "Universal ICU Core Network",
  plot_subtitle = "Highly conserved interactions present in >= 3 ICU spatial zones"
)

cat("\n✅ 出图完成！如果您需要 Cytoscape 排版，配套的 _CytoNodes.csv 已自动生成至 Network_Outputs 目录。\n")


############################################################
# ============================================================
# 🎯 修复版：ESKAPE Intra-Group Network (4-Panel Facet Plot)
# ============================================================
library(tidyverse)
library(igraph)
library(ggraph)
library(ggrepel)

ESKAPE_species <- c(
  "Enterococcus faecium", "Staphylococcus aureus",
  "Klebsiella pneumoniae", "Acinetobacter baumannii",
  "Pseudomonas aeruginosa", "Enterobacter cloacae",
  "Enterobacter cloacae complex sp."
)

cat("\n=== 正在提取 ESKAPE 病原体内部的互作网络 ===\n")

# 1. 严格筛选、去重与列位置重排
eskape_net_df <- all_edges %>%
  filter(Source %in% ESKAPE_species & Target %in% ESKAPE_species) %>%
  mutate(
    Node1 = pmin(Source, Target),
    Node2 = pmax(Source, Target)
  ) %>%
  distinct(Site, Node1, Node2, .keep_all = TRUE) %>%
  rename(from = Node1, to = Node2) %>%
  # 【核心修复】：强制将 from 和 to 移到前两列，防止 igraph 把 Site 当成节点！
  relocate(from, to)

if(nrow(eskape_net_df) == 0) {
  stop("⚠️ 没有发现 ESKAPE 物种之间的互作，请检查数据。")
}

# 2. 构建包含所有 Site 边的全局图对象
eskape_nodes <- tibble(name = unique(c(eskape_net_df$from, eskape_net_df$to)))
g <- graph_from_data_frame(eskape_net_df, directed = FALSE, vertices = eskape_nodes)

# 3. 构建包含所有 Site 边的全局图对象
eskape_nodes <- tibble(name = unique(c(eskape_net_df$from, eskape_net_df$to)))
g <- graph_from_data_frame(eskape_net_df, directed = FALSE, vertices = eskape_nodes)

# 4. ggraph 2x2 分面可视化 (修复版)
set.seed(123) 

p <- ggraph(g, layout = "circle") + 
  
  # 1. 连线：映射颜色与宽度
  geom_edge_link(
    aes(width = abs_Weight, color = Direction, label = sprintf("%.2f", abs_Weight)), 
    alpha = 0.8,
    label_size = 3.2,
    label_colour = "black",
    angle_calc = 'along',
    label_dodge = unit(3, 'mm'),
    show.legend = TRUE # 确保边颜色进入图例
  ) +
  
  # 2. 节点
  geom_node_point(size = 8, fill = "#E63946", shape = 21, color = "white", stroke = 1.2) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "italic", color = "black") +
  
  # 3. 分面：2x2 排布
  facet_edges(~ Site, nrow = 2, ncol = 2) +
  
  # 4. 颜色与粗细刻度
  scale_edge_color_manual(
    values = c("positive" = "#E24B4A", "negative" = "#378ADD"),
    labels = c("negative" = "Negative Correlation (Competition)", 
               "positive" = "Positive Correlation (Co-occurrence)"),
    name = "Ecological Interaction"
  ) +
  scale_edge_width(range = c(0.8, 3), guide = "none") + # 宽度不显示图例，避免杂乱
  
  # 5. 核心样式调整：实现物理分隔感
  theme_graph(base_family = "sans") +
  theme(
    # --- A. 分面标题栏美化 ---
    strip.text = element_text(size = 12, face = "bold", margin = margin(5,0,5,0)),
    strip.background = element_rect(fill = "gray90", color = "gray60", linewidth = 0.5), # 标题灰色背景
    
    # --- B. 子图边框线 (增加物理界限) ---
    panel.border = element_rect(colour = "gray70", fill = NA, linewidth = 0.8),
    panel.spacing = unit(1.5, "lines"), # 增加四个图之间的间距
    
    # --- C. 图例位置与样式 ---
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(1.5, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    
    plot.title = element_text(size = 16, face = "bold", margin = margin(0,0,10,0)),
    plot.subtitle = element_text(size = 11, color = "gray30")
  ) +
  
  labs(
    title = "Spatio-temporal Interactions of ESKAPE Pathogens",
    subtitle = "Conserved and site-specific patterns within the ICU environment"
  )

# 打印并导出
print(p)

ggsave("/home/bio/Documents/AH_ICU/Abundance_shortReads/Network_Outputs/ESKAPE_IntraGroup_Final_Boxed.pdf", 
       plot = p, width = 11, height = 11, bg = "white")



###############################3333
library(tidyverse)
library(igraph)
library(ggraph)
library(ggrepel)
library(tidygraph)

# ==========================================
# 1. 参数与颜色配置区
# ==========================================
# 显著性连线颜色
color_pos       <- "#378ADD"   # 正相关 (蓝色)
color_neg       <- "#E24B4A"   # 负相关 (红色)
color_ns        <- "gray80"    # 无显著相关 (灰色)

# 节点的专属颜色字典
eskape_colors <- c(
  "Enterococcus faecium"             = "#3C5488",  # 红色系 (原备注，实为深蓝，以你的代码为准)
  "Staphylococcus aureus"            = "#4DBBD5",  # 蓝色系
  "Klebsiella pneumoniae"            = "#00A087",  # 绿色系
  "Acinetobacter baumannii"          = "#E64B35",  # 深蓝/紫色系 (原备注，实为红色，以代码为准)
  "Pseudomonas aeruginosa"           = "#F39B7F",  # 珊瑚橘
  "Enterobacter cloacae"             = "#8491B4",  # 灰蓝色
  "Enterobacter cloacae complex sp." = "#91D1C2"   # 薄荷绿
)

ESKAPE_species <- names(eskape_colors)

cat("\n=== 正在构建 ESKAPE 全连接网络以对齐固定节点 ===\n")

# ==========================================
# 2. 构建全连接边数据集 
# ==========================================
sites <- unique(eskape_net_df$Site)
if(length(sites) == 0) stop("未检测到 Site 信息，请检查 eskape_net_df。")

# 生成 7 个 ESKAPE 物种的所有两两组合 (21 条边)
all_pairs <- combn(ESKAPE_species, 2, simplify = FALSE)
base_edges <- tibble(
  from = sapply(all_pairs, function(x) min(x)), 
  to   = sapply(all_pairs, function(x) max(x))
)

# 将这 21 条边扩展到每一个 Site 中
complete_edges <- expand_grid(base_edges, Site = sites)

# 将“完整结构”与“真实显著的数据”进行 Left Join 合并
plot_edges <- complete_edges %>%
  left_join(eskape_net_df, by = c("from", "to", "Site")) %>%
  mutate(
    # 根据 Direction 判定连线类型
    edge_class = case_when(
      Direction == "positive" ~ "Positive",
      Direction == "negative" ~ "Negative",
      TRUE ~ "Non-significant"
    ),
    # 无显著相关的边，标签设为空白，不显示数值
    edge_label = ifelse(is.na(Weight), "", sprintf("%.2f", Weight))
  )

# ==========================================
# 3. 循环生成并导出 4 张独立的图
# ==========================================
# 提取节点（固定顺序，保证 layout = "circle" 每次画出的位置绝对一致）
eskape_nodes <- tibble(name = ESKAPE_species)

# 设定基础输出路径
out_dir <- "/home/bio/Documents/AH_ICU/Abundance_shortReads/Network_Outputs/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (current_site in sites) {
  
  cat(sprintf("正在生成 %s 的网络图...\n", current_site))
  
  # 提取当前 Site 的边
  current_edges <- plot_edges %>% filter(Site == current_site)
  
  # 构建 tbl_graph
  g <- tbl_graph(nodes = eskape_nodes, edges = current_edges, directed = FALSE)
  
  # 开始画图
  set.seed(123) 
  p <- ggraph(g, layout = "circle") + 
    
    # 1. 连线：映射颜色、线型、宽度和标签 (取消了连续性的 abs_Weight 映射)
    geom_edge_arc(
      aes(
        color = edge_class, 
        linetype = edge_class, 
        width = edge_class,   # 宽度映射给分类变量
        label = edge_label
      ), 
      strength = 0.2,         # 弧度
      alpha = 0.8,
      label_size = 3.2,
      label_colour = "black",
      angle_calc = 'along',
      label_dodge = unit(3, 'mm'),
      show.legend = TRUE
    ) +
    
    # 2. 节点：使用 fill = name 映射到 eskape_colors
    geom_node_point(aes(fill = name), size = 9, shape = 21, color = "white", stroke = 1.2) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3.5, fontface = "italic", color = "black") +
    
    # 3. 映射比例尺设置
    # A. 边颜色
    scale_edge_color_manual(
      values = c(
        "Positive" = color_pos, 
        "Negative" = color_neg, 
        "Non-significant" = color_ns
      ),
      name = "Ecological Interaction"
    ) +
    # B. 边线型 (无显著相关为虚线)
    scale_edge_linetype_manual(
      values = c(
        "Positive" = "solid", 
        "Negative" = "solid", 
        "Non-significant" = "dashed"
      ),
      name = "Ecological Interaction"
    ) +
    # C. 边宽度 (固定粗细：显著为1.2，不显著为0.4)
    scale_edge_width_manual(
      values = c(
        "Positive" = 1.2, 
        "Negative" = 1.2, 
        "Non-significant" = 0.4
      ),
      guide = "none" # 宽度不显示图例
    ) +
    # D. 节点颜色 (调用你提供的字典)
    scale_fill_manual(
      values = eskape_colors,
      guide = "none" # 节点颜色也不需要图例，因为标签已经写在旁边了
    ) +
    
    # 4. 主题样式
    theme_graph(base_family = "sans") +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.width = unit(1.5, "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold"),
      
      plot.title = element_text(size = 16, face = "bold", margin = margin(0,0,10,0)),
      plot.subtitle = element_text(size = 12, color = "gray30", margin = margin(0,0,15,0))
    ) +
    
    labs(
      title = sprintf("Spatio-temporal Interactions: %s", current_site),
      subtitle = "ESKAPE Pathogens Intra-group Co-occurrence Network"
    )
  
  # 保存 PDF
  file_name <- paste0(out_dir, "ESKAPE_Network_", current_site, "_FixedNodes.pdf")
  ggsave(filename = file_name, plot = p, width = 8.5, height = 8.5, bg = "white")
}

cat("\n=== 绘图完成！四张独立的 PDF 已保存到 Network_Outputs 文件夹。 ===\n")


# ============================================================
# 7. PathogenOnly 图：导出对应表格并使用防重叠标签重新绘图
# ============================================================
# 该步骤独立封装在同目录的脚本中。它读取已经导出的
# *_PathogenOnly_DualCircle_edges.csv / nodes.csv，因而不会重复运行网络推断。
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  getwd()
}
pathogen_only_addon <- file.path(script_dir, "PathogenOnly_tables_and_replot.R")
if (!file.exists(pathogen_only_addon)) {
  stop("缺少 PathogenOnly_tables_and_replot.R: ", pathogen_only_addon)
}
source(pathogen_only_addon, local = FALSE)
