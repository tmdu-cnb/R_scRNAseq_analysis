# =============================================================
# 細胞構成比統計検定スクリプト（ベイズ解析版）
# 論文準拠: 89%信頼区間におけるlog2(FC)の事後分布を計算
# 信頼区間が0を含まない場合に有意と判断
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

# =============================================================
# 3. ベイズ統計による事後分布の推定 ----------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("ベイズ統計による事後分布の推定（論文準拠）\n")
cat("89%信頼区間におけるlog2(FC)の事後分布を計算\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 細胞タイプの順序
celltype_order <- c(
  "Granule", "Oligodendrocytes", "Bergmann_glia", "Astrocytes",
  "Interneurons", "Endothelial", "Purkinje", "OPC",
  "Microglia", "Fibroblasts"
)

celltype_order <- celltype_order[celltype_order %in% unique(prop_by_sample_df$CellType)]

# ベイズ統計による解析
# 各細胞タイプについて、ベータ分布を仮定して事後分布を推定
bayesian_results <- data.frame(
  CellType = character(),
  DTA_mean = numeric(),
  DTA_ci_lower = numeric(),
  DTA_ci_upper = numeric(),
  CTL_mean = numeric(),
  CTL_ci_lower = numeric(),
  CTL_ci_upper = numeric(),
  Log2FC_mean = numeric(),
  Log2FC_ci_lower = numeric(),
  Log2FC_ci_upper = numeric(),
  Significant = logical(),
  stringsAsFactors = FALSE
)

# 信頼区間のレベル（論文では89%）
ci_level <- 0.89
alpha <- 1 - ci_level  # 0.11

set.seed(123)
n_samples_bayes <- 10000  # 事後分布からのサンプリング数

for (celltype in celltype_order) {
  # DTA群のデータ
  dta_props <- prop_by_sample_df$Proportion[
    prop_by_sample_df$Group == "DTA" & prop_by_sample_df$CellType == celltype
  ]
  dta_counts <- prop_by_sample_df$Count[
    prop_by_sample_df$Group == "DTA" & prop_by_sample_df$CellType == celltype
  ]
  dta_totals <- prop_by_sample_df$TotalCells[
    prop_by_sample_df$Group == "DTA" & prop_by_sample_df$CellType == celltype
  ]
  
  # CTL群のデータ
  ctl_props <- prop_by_sample_df$Proportion[
    prop_by_sample_df$Group == "CTL" & prop_by_sample_df$CellType == celltype
  ]
  ctl_counts <- prop_by_sample_df$Count[
    prop_by_sample_df$Group == "CTL" & prop_by_sample_df$CellType == celltype
  ]
  ctl_totals <- prop_by_sample_df$TotalCells[
    prop_by_sample_df$Group == "CTL" & prop_by_sample_df$CellType == celltype
  ]
  
  # ベイズ統計：各サンプルでの割合の事後分布を推定
  # ベータ分布を仮定（無情報事前分布: Beta(1,1)）
  
  # DTA群の事後分布からのサンプリング
  # 各サンプルの事後分布からサンプリングして平均を取る
  dta_posterior_samples <- numeric(n_samples_bayes)
  for (i in 1:length(dta_props)) {
    # Beta(alpha + count, beta + total - count) where alpha=beta=1 (無情報事前分布)
    alpha_post <- 1 + dta_counts[i]
    beta_post <- 1 + dta_totals[i] - dta_counts[i]
    dta_posterior_samples <- dta_posterior_samples + rbeta(n_samples_bayes, alpha_post, beta_post)
  }
  dta_posterior_samples <- dta_posterior_samples / length(dta_props)
  
  # CTL群の事後分布からのサンプリング
  ctl_posterior_samples <- numeric(n_samples_bayes)
  for (i in 1:length(ctl_props)) {
    alpha_post <- 1 + ctl_counts[i]
    beta_post <- 1 + ctl_totals[i] - ctl_counts[i]
    ctl_posterior_samples <- ctl_posterior_samples + rbeta(n_samples_bayes, alpha_post, beta_post)
  }
  ctl_posterior_samples <- ctl_posterior_samples / length(ctl_props)
  
  # Log2(FC)の事後分布
  log2fc_posterior <- log2(dta_posterior_samples / ctl_posterior_samples)
  
  # 信頼区間の計算
  dta_ci <- quantile(dta_posterior_samples, probs = c(alpha/2, 1 - alpha/2))
  ctl_ci <- quantile(ctl_posterior_samples, probs = c(alpha/2, 1 - alpha/2))
  log2fc_ci <- quantile(log2fc_posterior, probs = c(alpha/2, 1 - alpha/2))
  
  # 信頼区間が0を含まない場合に有意と判断
  significant <- log2fc_ci[1] > 0 || log2fc_ci[2] < 0
  
  bayesian_results <- rbind(
    bayesian_results,
    data.frame(
      CellType = celltype,
      DTA_mean = mean(dta_posterior_samples),
      DTA_ci_lower = dta_ci[1],
      DTA_ci_upper = dta_ci[2],
      CTL_mean = mean(ctl_posterior_samples),
      CTL_ci_lower = ctl_ci[1],
      CTL_ci_upper = ctl_ci[2],
      Log2FC_mean = mean(log2fc_posterior),
      Log2FC_ci_lower = log2fc_ci[1],
      Log2FC_ci_upper = log2fc_ci[2],
      Significant = significant
    )
  )
}

cat("ベイズ統計解析結果:\n\n")
print(bayesian_results, row.names = FALSE)

cat("\n")

# =============================================================
# 4. 結果の解釈と表示 -------------------------------------------
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("有意な差がある細胞タイプ（89%信頼区間が0を含まない）\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

significant_cells <- bayesian_results[bayesian_results$Significant == TRUE, ]
if (nrow(significant_cells) > 0) {
  for (i in 1:nrow(significant_cells)) {
    row <- significant_cells[i, ]
    direction <- ifelse(row$Log2FC_mean > 0, "増加", "減少")
    cat(sprintf("%-18s: Log2FC = %.3f [%.3f, %.3f] (%s)\n",
                row$CellType,
                row$Log2FC_mean,
                row$Log2FC_ci_lower,
                row$Log2FC_ci_upper,
                direction))
    cat(sprintf("  DTA: %.2f%% [%.2f%%, %.2f%%], CTL: %.2f%% [%.2f%%, %.2f%%]\n",
                row$DTA_mean * 100, row$DTA_ci_lower * 100, row$DTA_ci_upper * 100,
                row$CTL_mean * 100, row$CTL_ci_lower * 100, row$CTL_ci_upper * 100))
  }
} else {
  cat("  なし（全ての細胞タイプで89%信頼区間が0を含む）\n")
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

# ベイズ統計解析結果を保存
write.csv(bayesian_results,
          file = file.path(output_dir, "統計検定_ベイズ解析結果.csv"),
          row.names = FALSE)

cat("✓ 結果を保存しました:\n")
cat("  -", file.path(output_dir, "統計検定_ベイズ解析結果.csv"), "\n")

# =============================================================
# 6. グラフ作成（論文のPanel Dに類似） -------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("グラフを作成中（論文のPanel Dに類似）...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Ridge plot風の図（論文のPanel D）
cat("1. Log2(FC)の事後分布Ridge Plotを作成中...\n")

# 事後分布のサンプルを保存するためのデータフレーム
posterior_samples_all <- data.frame(
  CellType = character(),
  Log2FC = numeric(),
  stringsAsFactors = FALSE
)

set.seed(123)
for (celltype in celltype_order) {
  dta_props <- prop_by_sample_df$Proportion[
    prop_by_sample_df$Group == "DTA" & prop_by_sample_df$CellType == celltype
  ]
  dta_counts <- prop_by_sample_df$Count[
    prop_by_sample_df$Group == "DTA" & prop_by_sample_df$CellType == celltype
  ]
  dta_totals <- prop_by_sample_df$TotalCells[
    prop_by_sample_df$Group == "DTA" & prop_by_sample_df$CellType == celltype
  ]
  
  ctl_props <- prop_by_sample_df$Proportion[
    prop_by_sample_df$Group == "CTL" & prop_by_sample_df$CellType == celltype
  ]
  ctl_counts <- prop_by_sample_df$Count[
    prop_by_sample_df$Group == "CTL" & prop_by_sample_df$CellType == celltype
  ]
  ctl_totals <- prop_by_sample_df$TotalCells[
    prop_by_sample_df$Group == "CTL" & prop_by_sample_df$CellType == celltype
  ]
  
  dta_posterior_samples <- numeric(n_samples_bayes)
  for (i in 1:length(dta_props)) {
    alpha_post <- 1 + dta_counts[i]
    beta_post <- 1 + dta_totals[i] - dta_counts[i]
    dta_posterior_samples <- dta_posterior_samples + rbeta(n_samples_bayes, alpha_post, beta_post)
  }
  dta_posterior_samples <- dta_posterior_samples / length(dta_props)
  
  ctl_posterior_samples <- numeric(n_samples_bayes)
  for (i in 1:length(ctl_props)) {
    alpha_post <- 1 + ctl_counts[i]
    beta_post <- 1 + ctl_totals[i] - ctl_counts[i]
    ctl_posterior_samples <- ctl_posterior_samples + rbeta(n_samples_bayes, alpha_post, beta_post)
  }
  ctl_posterior_samples <- ctl_posterior_samples / length(ctl_props)
  
  log2fc_posterior <- log2(dta_posterior_samples / ctl_posterior_samples)
  
  posterior_samples_all <- rbind(
    posterior_samples_all,
    data.frame(
      CellType = rep(celltype, length(log2fc_posterior)),
      Log2FC = log2fc_posterior
    )
  )
}

# 細胞タイプの順序を逆順に（論文では下から上へ）
posterior_samples_all$CellType <- factor(
  posterior_samples_all$CellType,
  levels = rev(celltype_order)
)

# 色の定義
celltype_colors <- c(
  "Granule" = "#E31A1C",
  "Oligodendrocytes" = "#FF7F00",
  "Bergmann_glia" = "#B2DF8A",
  "Astrocytes" = "#33A02C",
  "Interneurons" = "#1F78B4",
  "Endothelial" = "#A6CEE3",
  "Purkinje" = "#6A3D9A",
  "OPC" = "#CAB2D6",
  "Microglia" = "#FB9A99",
  "Fibroblasts" = "#E7298A"
)

# Ridge plot（ggridgesが必要）
if (!requireNamespace("ggridges", quietly = TRUE)) {
  cat("警告: ggridgesパッケージがインストールされていません。\n")
  cat("      install.packages('ggridges')を実行してください。\n")
} else {
  library(ggridges)
  
  p_ridge <- ggplot(posterior_samples_all, aes(x = Log2FC, y = CellType, fill = CellType)) +
    geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) +
    scale_fill_manual(values = celltype_colors, guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    # 89%信頼区間を表示
    geom_segment(
      data = bayesian_results,
      aes(x = Log2FC_ci_lower, xend = Log2FC_ci_upper, y = CellType, yend = CellType),
      color = ifelse(bayesian_results$Significant, "red", "gray50"),
      linewidth = 2,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = bayesian_results,
      aes(x = Log2FC_mean, y = CellType),
      color = "black",
      size = 2,
      inherit.aes = FALSE
    ) +
    scale_x_continuous(
      limits = c(-2, 2),
      breaks = seq(-2, 2, by = 0.5),
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    labs(
      x = "Log2(FC) (DTA / CTL)",
      y = "",
      title = "Posterior Distribution of Log2(FC) (89% Credible Interval)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 14, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      axis.line = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(
    filename = file.path(output_dir, "統計検定_ベイズ解析_RidgePlot.png"),
    plot = p_ridge,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  cat("  ✓ Ridge Plotを保存しました\n")
}

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("解析完了！\n")
cat("※ このスクリプトは論文準拠のベイズ統計解析を行っています。\n")
cat("   89%信頼区間が0を含まない場合に有意と判断します。\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

