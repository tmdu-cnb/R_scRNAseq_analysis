# =============================================================
# 細胞構成比統計検定スクリプト（サンプル単位）
# DTAとCTLの細胞タイプ構成比の統計的有意差を検定
# 各サンプル（n=3 vs n=3）での割合を比較
# =============================================================

# 0. ライブラリ -------------------------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# =============================================================
# 1. 処理済みデータの読み込み ------------------------------------
cat("処理済みデータを読み込み中...\n")

# データ読み込み（複数のパスを試す）
data_paths <- c(
  "seu_dta_processed.rds",
  file.path(getwd(), "seu_dta_processed.rds"),
  "/Users/godakyosuke/Desktop/goda sc RNA seq/seu_dta_processed.rds"
)

dta_file <- NULL
ctl_file <- NULL

for (path in data_paths) {
  if (file.exists(path)) {
    dta_file <- path
    ctl_file <- sub("seu_dta_processed.rds", "seu_ctl_processed.rds", path)
    if (file.exists(ctl_file)) {
      break
    }
  }
}

if (is.null(dta_file) || !file.exists(dta_file)) {
  stop("エラー: seu_dta_processed.rdsが見つかりません。\n",
       "      現在の作業ディレクトリ: ", getwd(), "\n",
       "      ファイルを確認してください。")
}

if (!file.exists(ctl_file)) {
  stop("エラー: seu_ctl_processed.rdsが見つかりません。\n",
       "      ファイルを確認してください。")
}

cat("✓ データファイルを検出:\n")
cat("  DTA:", dta_file, "\n")
cat("  CTL:", ctl_file, "\n")

seu_dta <- readRDS(dta_file)
seu_ctl <- readRDS(ctl_file)

cat("✓ DTAデータ読み込み完了:", ncol(seu_dta), "cells\n")
cat("✓ CTLデータ読み込み完了:", ncol(seu_ctl), "cells\n")

# サンプル情報の確認
if (!"sample_id" %in% colnames(seu_dta@meta.data)) {
  stop("エラー: sample_idがメタデータに存在しません。\n",
       "      元のデータから再計算する必要があります。")
}

cat("\n✓ DTA群のサンプル数:", length(unique(seu_dta$sample_id)), "\n")
cat("  サンプルID:", paste(unique(seu_dta$sample_id), collapse = ", "), "\n")
cat("✓ CTL群のサンプル数:", length(unique(seu_ctl$sample_id)), "\n")
cat("  サンプルID:", paste(unique(seu_ctl$sample_id), collapse = ", "), "\n")

# celltype_10classが存在しない場合は自動分類
if (!"celltype_10class" %in% colnames(seu_dta@meta.data) || 
    !"celltype_10class" %in% colnames(seu_ctl@meta.data)) {
  cat("\ncelltype_10classがメタデータに存在しないため、自動分類を実行します...\n")
  
  # 論文ベースの小脳細胞10タイプ分類関数
  classify_cerebellar_cells <- function(seu_obj) {
    counts <- GetAssayData(seu_obj, layer = "counts")
    data_norm <- GetAssayData(seu_obj, layer = "data")
    
    available_genes <- rownames(counts)
    cell_type <- rep("Unclassified", ncol(seu_obj))
    scores <- matrix(0, nrow = ncol(seu_obj), ncol = 10)
    colnames(scores) <- c("Granule", "Purkinje", "Bergmann_glia", "Interneurons",
                           "Astrocytes", "Oligodendrocytes", "OPC", 
                           "Microglia", "Endothelial", "Fibroblasts")
    
    granule_markers <- c("Rims1", "Grm4", "Gabra6")
    granule_markers <- granule_markers[granule_markers %in% available_genes]
    if (length(granule_markers) > 0) {
      scores[, "Granule"] <- colMeans(data_norm[granule_markers, , drop = FALSE])
    }
    
    purkinje_markers <- c("Itpr1", "Calb1", "Car8")
    purkinje_markers <- purkinje_markers[purkinje_markers %in% available_genes]
    if (length(purkinje_markers) > 0) {
      scores[, "Purkinje"] <- colMeans(data_norm[purkinje_markers, , drop = FALSE])
    }
    
    bergmann_markers <- c("Tubb2b", "Aqp4")
    bergmann_markers <- bergmann_markers[bergmann_markers %in% available_genes]
    if (length(bergmann_markers) > 0) {
      scores[, "Bergmann_glia"] <- colMeans(data_norm[bergmann_markers, , drop = FALSE])
    }
    
    inter_markers <- c("Gad1", "Pvalb")
    inter_markers <- inter_markers[inter_markers %in% available_genes]
    if (length(inter_markers) > 0) {
      scores[, "Interneurons"] <- colMeans(data_norm[inter_markers, , drop = FALSE])
    }
    
    astro_markers <- c("Ttn", "Aqp1", "Gfap")
    astro_markers <- astro_markers[astro_markers %in% available_genes]
    if (length(astro_markers) > 0) {
      scores[, "Astrocytes"] <- colMeans(data_norm[astro_markers, , drop = FALSE])
    }
    
    oligo_markers <- c("Plp1", "Mbp")
    oligo_markers <- oligo_markers[oligo_markers %in% available_genes]
    if (length(oligo_markers) > 0) {
      scores[, "Oligodendrocytes"] <- colMeans(data_norm[oligo_markers, , drop = FALSE])
    }
    
    opc_markers <- c("Pdgfra", "Olig1")
    opc_markers <- opc_markers[opc_markers %in% available_genes]
    if (length(opc_markers) > 0) {
      scores[, "OPC"] <- colMeans(data_norm[opc_markers, , drop = FALSE])
    }
    
    micro_markers <- c("Cd74", "Csf1r")
    micro_markers <- micro_markers[micro_markers %in% available_genes]
    if (length(micro_markers) > 0) {
      scores[, "Microglia"] <- colMeans(data_norm[micro_markers, , drop = FALSE])
    }
    
    endo_markers <- c("Cldn5", "Vwf")
    endo_markers <- endo_markers[endo_markers %in% available_genes]
    if (length(endo_markers) > 0) {
      scores[, "Endothelial"] <- colMeans(data_norm[endo_markers, , drop = FALSE])
    }
    
    fibro_markers <- c("Dcn", "Apod")
    fibro_markers <- fibro_markers[fibro_markers %in% available_genes]
    if (length(fibro_markers) > 0) {
      scores[, "Fibroblasts"] <- colMeans(data_norm[fibro_markers, , drop = FALSE])
    }
    
    max_types <- apply(scores, 1, which.max)
    for (i in 1:length(cell_type)) {
      cell_type[i] <- colnames(scores)[max_types[i]]
    }
    
    return(list(celltype = cell_type, scores = scores))
  }
  
  if (!"celltype_10class" %in% colnames(seu_dta@meta.data)) {
    cat("  DTAデータを分類中...\n")
    result_dta <- classify_cerebellar_cells(seu_dta)
    seu_dta$celltype_10class <- result_dta$celltype
    cat("  ✓ DTA分類完了\n")
  }
  
  if (!"celltype_10class" %in% colnames(seu_ctl@meta.data)) {
    cat("  CTLデータを分類中...\n")
    result_ctl <- classify_cerebellar_cells(seu_ctl)
    seu_ctl$celltype_10class <- result_ctl$celltype
    cat("  ✓ CTL分類完了\n")
  }
  
  cat("✓ 自動分類完了\n\n")
}

# =============================================================
# 2. 各サンプルでの細胞タイプ構成比を計算 -----------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("各サンプルでの細胞タイプ構成比を計算中...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# DTA群：各サンプルでの構成比
dta_samples <- unique(seu_dta$sample_id)
dta_prop_by_sample <- list()

for (sample in dta_samples) {
  seu_sample <- subset(seu_dta, subset = sample_id == sample)
  counts <- table(seu_sample$celltype_10class)
  prop <- prop.table(counts)
  
  dta_prop_by_sample[[sample]] <- data.frame(
    Sample = sample,
    Group = "DTA",
    CellType = names(prop),
    Count = as.numeric(counts),
    Proportion = as.numeric(prop),
    TotalCells = ncol(seu_sample)
  )
}

dta_prop_df <- bind_rows(dta_prop_by_sample)

# CTL群：各サンプルでの構成比
ctl_samples <- unique(seu_ctl$sample_id)
ctl_prop_by_sample <- list()

for (sample in ctl_samples) {
  seu_sample <- subset(seu_ctl, subset = sample_id == sample)
  counts <- table(seu_sample$celltype_10class)
  prop <- prop.table(counts)
  
  ctl_prop_by_sample[[sample]] <- data.frame(
    Sample = sample,
    Group = "CTL",
    CellType = names(prop),
    Count = as.numeric(counts),
    Proportion = as.numeric(prop),
    TotalCells = ncol(seu_sample)
  )
}

ctl_prop_df <- bind_rows(ctl_prop_by_sample)

# 結合
prop_by_sample_df <- bind_rows(dta_prop_df, ctl_prop_df)

cat("各サンプルでの構成比:\n\n")
print(prop_by_sample_df)

# =============================================================
# 3. 各細胞タイプごとに統計検定（n=3 vs n=3） -------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("統計検定: 各サンプルでの割合を比較（n=3 vs n=3）\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 細胞タイプの順序
celltype_order <- c(
  "Granule", "Oligodendrocytes", "Bergmann_glia", "Astrocytes",
  "Interneurons", "Endothelial", "Purkinje", "OPC",
  "Microglia", "Fibroblasts"
)

celltype_order <- celltype_order[celltype_order %in% unique(prop_by_sample_df$CellType)]

# 各細胞タイプごとに検定
statistical_results <- data.frame(
  CellType = character(),
  DTA_mean = numeric(),
  DTA_sd = numeric(),
  CTL_mean = numeric(),
  CTL_sd = numeric(),
  t_test_p = numeric(),
  mw_test_p = numeric(),
  p_adj_fdr_t = numeric(),
  p_adj_fdr_mw = numeric(),
  stringsAsFactors = FALSE
)

for (celltype in celltype_order) {
  # DTA群の割合
  dta_props <- prop_by_sample_df$Proportion[
    prop_by_sample_df$Group == "DTA" & prop_by_sample_df$CellType == celltype
  ]
  
  # CTL群の割合
  ctl_props <- prop_by_sample_df$Proportion[
    prop_by_sample_df$Group == "CTL" & prop_by_sample_df$CellType == celltype
  ]
  
  if (length(dta_props) == 3 && length(ctl_props) == 3) {
    # t検定（正規性を仮定）
    t_test <- tryCatch(
      t.test(dta_props, ctl_props, alternative = "two.sided"),
      error = function(e) NULL
    )
    
    # Mann-Whitney U test（ノンパラメトリック）
    mw_test <- wilcox.test(dta_props, ctl_props, alternative = "two.sided")
    
    statistical_results <- rbind(
      statistical_results,
      data.frame(
        CellType = celltype,
        DTA_mean = mean(dta_props),
        DTA_sd = sd(dta_props),
        CTL_mean = mean(ctl_props),
        CTL_sd = sd(ctl_props),
        t_test_p = ifelse(!is.null(t_test), t_test$p.value, NA),
        mw_test_p = mw_test$p.value,
        p_adj_fdr_t = NA,
        p_adj_fdr_mw = NA
      )
    )
  }
}

# 多重比較補正
statistical_results$p_adj_fdr_t <- p.adjust(statistical_results$t_test_p, method = "fdr")
statistical_results$p_adj_fdr_mw <- p.adjust(statistical_results$mw_test_p, method = "fdr")

# Log2(FC)を計算
statistical_results$Log2FC <- log2(statistical_results$DTA_mean / statistical_results$CTL_mean)

# Cohen's dを計算（効果量）
# Cohen's d = (mean1 - mean2) / pooled_sd
# pooled_sd = sqrt(((n1-1)*sd1^2 + (n2-1)*sd2^2) / (n1+n2-2))
n_dta <- 3
n_ctl <- 3
statistical_results$Cohens_d <- (statistical_results$DTA_mean - statistical_results$CTL_mean) / 
  sqrt(((n_dta - 1) * statistical_results$DTA_sd^2 + (n_ctl - 1) * statistical_results$CTL_sd^2) / (n_dta + n_ctl - 2))

# 効果量の大きさを分類
statistical_results$EffectSize <- ifelse(
  abs(statistical_results$Cohens_d) < 0.2, "negligible",
  ifelse(
    abs(statistical_results$Cohens_d) < 0.5, "small",
    ifelse(
      abs(statistical_results$Cohens_d) < 0.8, "medium",
      "large"
    )
  )
)

cat("統計検定結果:\n\n")
print(statistical_results, row.names = FALSE)

cat("\n【Cohen's d効果量の解釈】\n")
cat("|d| < 0.2: negligible (実質的に差がない)\n")
cat("0.2 ≤ |d| < 0.5: small (小さな差)\n")
cat("0.5 ≤ |d| < 0.8: medium (中程度の差)\n")
cat("|d| ≥ 0.8: large (大きな差)\n\n")

cat("\n")

# =============================================================
# 4. 結果の解釈と表示 -------------------------------------------
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("有意な差がある細胞タイプ（Mann-Whitney U test, FDR補正後、q < 0.05）\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

significant_mw <- statistical_results[statistical_results$p_adj_fdr_mw < 0.05, ]
if (nrow(significant_mw) > 0) {
  significant_mw <- significant_mw[order(significant_mw$p_adj_fdr_mw), ]
  for (i in 1:nrow(significant_mw)) {
    row <- significant_mw[i, ]
    direction <- ifelse(row$Log2FC > 0, "増加", "減少")
    cat(sprintf("%-18s: p = %.4e, FDR補正後 = %.4e, Log2(FC) = %.3f (%s)\n",
                row$CellType,
                row$mw_test_p,
                row$p_adj_fdr_mw,
                row$Log2FC,
                direction))
    cat(sprintf("  DTA: %.2f%% ± %.2f%%, CTL: %.2f%% ± %.2f%%\n",
                row$DTA_mean * 100, row$DTA_sd * 100,
                row$CTL_mean * 100, row$CTL_sd * 100))
  }
} else {
  cat("  なし\n")
}

cat("\n")

# =============================================================
# 5. 結果をCSVファイルに保存 -----------------------------------
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("結果をCSVファイルに保存中...\n")

output_dir <- dirname(dta_file)
if (output_dir == ".") {
  output_dir <- getwd()
}

# サンプル単位の構成比データを保存
write.csv(prop_by_sample_df,
          file = file.path(output_dir, "統計検定_サンプル単位構成比.csv"),
          row.names = FALSE)

# 統計検定結果を保存
write.csv(statistical_results,
          file = file.path(output_dir, "統計検定_サンプル単位検定結果.csv"),
          row.names = FALSE)

cat("✓ 結果を保存しました:\n")
cat("  -", file.path(output_dir, "統計検定_サンプル単位構成比.csv"), "\n")
cat("  -", file.path(output_dir, "統計検定_サンプル単位検定結果.csv"), "\n")

# =============================================================
# 6. グラフ作成 -------------------------------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("グラフを作成中...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 1. 各サンプルでの構成比をプロット（ドットプロット）
cat("1. サンプル単位構成比ドットプロットを作成中...\n")

plot_data_dot <- prop_by_sample_df |>
  mutate(
    CellType = factor(CellType, levels = celltype_order),
    Group = factor(Group, levels = c("DTA", "CTL"))
  )

p_dot <- ggplot(plot_data_dot, aes(x = CellType, y = Proportion, color = Group)) +
  geom_point(position = position_dodge(width = 0.5), size = 3, alpha = 0.7) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 5,
    shape = 21,
    fill = "white",
    position = position_dodge(width = 0.5)
  ) +
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_color_manual(values = c("DTA" = "#E31A1C", "CTL" = "#1F78B4")) +
  labs(
    x = "Cell Type",
    y = "Proportion",
    title = "Cell Type Proportions by Sample (n=3 per group)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  filename = file.path(output_dir, "統計検定_サンプル単位構成比ドットプロット.png"),
  plot = p_dot,
  width = 14,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("  ✓ ドットプロットを保存しました\n")

# 2. ボルケーノプロット（Mann-Whitney U test）
cat("2. ボルケーノプロットを作成中...\n")

plot_data_volcano <- statistical_results |>
  mutate(
    neg_log10_p = -log10(mw_test_p),
    neg_log10_p_adj = -log10(p_adj_fdr_mw),
    Significant = ifelse(p_adj_fdr_mw < 0.05, "Yes", "No")
  )

p_volcano <- ggplot(plot_data_volcano, aes(x = Log2FC, y = neg_log10_p_adj)) +
  geom_point(
    aes(color = Significant, fill = Significant),
    size = 4,
    alpha = 0.7,
    shape = 21
  ) +
  geom_text_repel(
    data = plot_data_volcano[plot_data_volcano$p_adj_fdr_mw < 0.05, ],
    aes(label = CellType),
    size = 4,
    fontface = "bold",
    box.padding = 0.5,
    point.padding = 0.3
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(
    values = c("Yes" = "red", "No" = "gray50"),
    name = "Significant\n(FDR < 0.05)"
  ) +
  scale_fill_manual(
    values = c("Yes" = "red", "No" = "gray50"),
    name = "Significant\n(FDR < 0.05)"
  ) +
  labs(
    x = "Log2(FC) (DTA / CTL)",
    y = "-log10(FDR adjusted p-value)",
    title = "Volcano Plot: Sample-level Comparison (n=3 vs n=3)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  filename = file.path(output_dir, "統計検定_サンプル単位ボルケーノプロット.png"),
  plot = p_volcano,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("  ✓ ボルケーノプロットを保存しました\n")

# 3. Log2(FC)の棒グラフ
cat("3. Log2(FC)棒グラフを作成中...\n")

plot_data_bar <- statistical_results |>
  mutate(
    CellType = factor(CellType, levels = celltype_order),
    Significant = ifelse(p_adj_fdr_mw < 0.05, "Yes", "No"),
    fill_color = ifelse(Log2FC > 0, "#E31A1C", "#1F78B4")
  )

p_bar_log2fc <- ggplot(plot_data_bar, aes(x = CellType, y = Log2FC, fill = fill_color)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) +
  scale_fill_identity() +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(
    x = "Cell Type",
    y = "Log2(FC) (DTA / CTL)",
    title = "Log2(FC) by Cell Type (n=3 vs n=3)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  filename = file.path(output_dir, "統計検定_サンプル単位_Log2FC棒グラフ.png"),
  plot = p_bar_log2fc,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("  ✓ Log2(FC)棒グラフを保存しました\n")

# 4. 構成比の棒グラフ（DTA vs CTL）
cat("4. 構成比比較棒グラフを作成中...\n")

plot_data_prop_bar <- statistical_results |>
  select(CellType, DTA_mean, CTL_mean) |>
  pivot_longer(cols = c(DTA_mean, CTL_mean), names_to = "Group", values_to = "Proportion") |>
  mutate(
    CellType = factor(CellType, levels = celltype_order),
    Group = factor(ifelse(Group == "DTA_mean", "DTA", "CTL"), levels = c("DTA", "CTL"))
  )

p_bar_prop <- ggplot(plot_data_prop_bar, aes(x = CellType, y = Proportion, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7, alpha = 0.8) +
  scale_fill_manual(values = c("DTA" = "#E31A1C", "CTL" = "#1F78B4")) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(
    x = "Cell Type",
    y = "Proportion",
    title = "Cell Type Proportions: DTA vs CTL (n=3 per group)",
    fill = "Group"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  filename = file.path(output_dir, "統計検定_サンプル単位_構成比棒グラフ.png"),
  plot = p_bar_prop,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("  ✓ 構成比比較棒グラフを保存しました\n")

# 5. Cohen's dの棒グラフ
cat("5. Cohen's d棒グラフを作成中...\n")

plot_data_cohensd <- statistical_results |>
  mutate(
    CellType = factor(CellType, levels = celltype_order),
    fill_color = ifelse(Cohens_d > 0, "#E31A1C", "#1F78B4"),
    EffectSize_label = paste0(EffectSize, " (|d|=", sprintf("%.2f", abs(Cohens_d)), ")")
  )

p_bar_cohensd <- ggplot(plot_data_cohensd, aes(x = CellType, y = Cohens_d, fill = fill_color)) +
  geom_bar(stat = "identity", width = 0.7, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 1) +
  # 効果量の基準線を追加
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dashed", color = "gray50", linewidth = 0.5, alpha = 0.5) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "gray50", linewidth = 0.5, alpha = 0.5) +
  geom_hline(yintercept = c(-0.8, 0.8), linetype = "dashed", color = "gray50", linewidth = 0.5, alpha = 0.5) +
  scale_fill_identity() +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(
    x = "Cell Type",
    y = "Cohen's d (Effect Size)",
    title = "Cohen's d Effect Size by Cell Type (n=3 vs n=3)",
    subtitle = "Reference lines: |d| = 0.2 (small), 0.5 (medium), 0.8 (large)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  filename = file.path(output_dir, "統計検定_サンプル単位_Cohens_d棒グラフ.png"),
  plot = p_bar_cohensd,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("  ✓ Cohen's d棒グラフを保存しました\n")

# 6. Cohen's dとLog2(FC)の比較プロット
cat("6. Cohen's d vs Log2(FC)比較プロットを作成中...\n")

plot_data_comparison <- statistical_results |>
  mutate(
    CellType = factor(CellType, levels = celltype_order),
    EffectSize = factor(EffectSize, levels = c("negligible", "small", "medium", "large"))
  )

p_comparison <- ggplot(plot_data_comparison, aes(x = Log2FC, y = Cohens_d, color = EffectSize)) +
  geom_point(size = 5, alpha = 0.7) +
  geom_text_repel(aes(label = CellType), size = 4, fontface = "bold", box.padding = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(
    values = c("negligible" = "gray70", "small" = "yellow", "medium" = "orange", "large" = "red"),
    name = "Effect Size"
  ) +
  labs(
    x = "Log2(FC) (DTA / CTL)",
    y = "Cohen's d",
    title = "Effect Size (Cohen's d) vs Log2(FC)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  filename = file.path(output_dir, "統計検定_サンプル単位_Cohens_d_vs_Log2FC.png"),
  plot = p_comparison,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("  ✓ Cohen's d vs Log2(FC)比較プロットを保存しました\n")

cat("\n✓ グラフを保存しました:\n")
cat("  -", file.path(output_dir, "統計検定_サンプル単位構成比ドットプロット.png"), "\n")
cat("  -", file.path(output_dir, "統計検定_サンプル単位ボルケーノプロット.png"), "\n")
cat("  -", file.path(output_dir, "統計検定_サンプル単位_Log2FC棒グラフ.png"), "\n")
cat("  -", file.path(output_dir, "統計検定_サンプル単位_構成比棒グラフ.png"), "\n")
cat("  -", file.path(output_dir, "統計検定_サンプル単位_Cohens_d棒グラフ.png"), "\n")
cat("  -", file.path(output_dir, "統計検定_サンプル単位_Cohens_d_vs_Log2FC.png"), "\n")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("解析完了！\n")
cat("※ このスクリプトは各サンプル（n=3 vs n=3）での割合を比較しています。\n")
cat("   疑似反復（pseudoreplication）の問題を回避した正しい統計検定です。\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

