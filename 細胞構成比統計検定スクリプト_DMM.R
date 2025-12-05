# =============================================================
# 細胞構成比統計検定スクリプト（ディリクレ多項式モデリング版）
# 論文準拠: DMM (Dirichlet Multinomial Modeling) using rstan
# 89%信頼区間におけるlog2(FC)の事後分布を計算
# =============================================================

# 0. ライブラリ -------------------------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

# rstanとbayestestRが必要
if (!requireNamespace("rstan", quietly = TRUE)) {
  stop("エラー: rstanパッケージがインストールされていません。\n",
       "      install.packages('rstan')を実行してください。\n",
       "      注意: rstanのインストールには時間がかかります。")
}

if (!requireNamespace("bayestestR", quietly = TRUE)) {
  stop("エラー: bayestestRパッケージがインストールされていません。\n",
       "      install.packages('bayestestR')を実行してください。")
}

library(rstan)
library(bayestestR)

# Stanのオプション設定
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

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
# 2. 各サンプルでの細胞タイプカウントマトリックスを作成 --------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("各サンプルでの細胞タイプカウントマトリックスを作成中...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 細胞タイプの順序
celltype_order <- c(
  "Granule", "Oligodendrocytes", "Bergmann_glia", "Astrocytes",
  "Interneurons", "Endothelial", "Purkinje", "OPC",
  "Microglia", "Fibroblasts"
)

celltype_order <- celltype_order[celltype_order %in% 
  unique(c(unique(seu_dta$celltype_10class), unique(seu_ctl$celltype_10class)))]

# DTA群：各サンプルでのカウント
dta_samples <- unique(seu_dta$sample_id)
dta_count_matrix <- matrix(0, nrow = length(dta_samples), ncol = length(celltype_order))
rownames(dta_count_matrix) <- dta_samples
colnames(dta_count_matrix) <- celltype_order

for (i in 1:length(dta_samples)) {
  sample <- dta_samples[i]
  seu_sample <- subset(seu_dta, subset = sample_id == sample)
  counts <- table(seu_sample$celltype_10class)
  dta_count_matrix[i, names(counts)] <- as.numeric(counts)
}

# CTL群：各サンプルでのカウント
ctl_samples <- unique(seu_ctl$sample_id)
ctl_count_matrix <- matrix(0, nrow = length(ctl_samples), ncol = length(celltype_order))
rownames(ctl_count_matrix) <- ctl_samples
colnames(ctl_count_matrix) <- celltype_order

for (i in 1:length(ctl_samples)) {
  sample <- ctl_samples[i]
  seu_sample <- subset(seu_ctl, subset = sample_id == sample)
  counts <- table(seu_sample$celltype_10class)
  ctl_count_matrix[i, names(counts)] <- as.numeric(counts)
}

cat("DTA群カウントマトリックス:\n")
print(dta_count_matrix)
cat("\nCTL群カウントマトリックス:\n")
print(ctl_count_matrix)

# =============================================================
# 3. Stanモデルの定義（ディリクレ多項式モデル） ----------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Stanモデルを定義中（ディリクレ多項式モデル）...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Stanモデルコード
stan_model_code <- "
data {
  int<lower=1> N;  // サンプル数
  int<lower=1> K;  // 細胞タイプ数
  int<lower=0> counts[N, K];  // カウントマトリックス
  real<lower=0> alpha_prior;  // ディリクレ事前分布のパラメータ（10^-7）
}

parameters {
  simplex[K] theta;  // 細胞タイプの割合（simplex制約: 合計=1）
}

model {
  // ディリクレ事前分布（均一な事前分布）
  theta ~ dirichlet(rep_vector(alpha_prior, K));
  
  // 多項式分布で各サンプルをモデル化
  for (n in 1:N) {
    counts[n] ~ multinomial(theta);
  }
}

generated quantities {
  // 事後予測分布のサンプル（オプション）
  int<lower=0> counts_pred[N, K];
  for (n in 1:N) {
    int total = sum(counts[n]);
    counts_pred[n] = multinomial_rng(theta, total);
  }
}
"

cat("Stanモデルをコンパイル中...\n")
stan_model <- stan_model(model_code = stan_model_code)
cat("✓ Stanモデルのコンパイル完了\n")

# =============================================================
# 4. DTA群の事後分布を推定 -------------------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("DTA群の事後分布を推定中（HMCサンプリング）...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# データ準備
dta_stan_data <- list(
  N = nrow(dta_count_matrix),
  K = ncol(dta_count_matrix),
  counts = dta_count_matrix,
  alpha_prior = 1e-7  # 論文通り: 10^-7
)

cat("DTA群のStanデータ:\n")
cat("  サンプル数 (N):", dta_stan_data$N, "\n")
cat("  細胞タイプ数 (K):", dta_stan_data$K, "\n")
cat("  事前分布パラメータ (alpha_prior):", dta_stan_data$alpha_prior, "\n\n")

# HMCサンプリング
cat("HMCサンプリングを実行中（時間がかかります）...\n")
dta_fit <- sampling(
  stan_model,
  data = dta_stan_data,
  iter = 4000,
  chains = 4,
  warmup = 2000,
  thin = 1,
  control = list(adapt_delta = 0.95)
)

cat("✓ DTA群のサンプリング完了\n")

# 事後分布からサンプルを抽出
dta_samples_post <- extract(dta_fit)$theta
cat("  事後分布サンプル数:", nrow(dta_samples_post), "\n")

# =============================================================
# 5. CTL群の事後分布を推定 -------------------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("CTL群の事後分布を推定中（HMCサンプリング）...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# データ準備
ctl_stan_data <- list(
  N = nrow(ctl_count_matrix),
  K = ncol(ctl_count_matrix),
  counts = ctl_count_matrix,
  alpha_prior = 1e-7
)

cat("CTL群のStanデータ:\n")
cat("  サンプル数 (N):", ctl_stan_data$N, "\n")
cat("  細胞タイプ数 (K):", ctl_stan_data$K, "\n")
cat("  事前分布パラメータ (alpha_prior):", ctl_stan_data$alpha_prior, "\n\n")

# HMCサンプリング
cat("HMCサンプリングを実行中（時間がかかります）...\n")
ctl_fit <- sampling(
  stan_model,
  data = ctl_stan_data,
  iter = 4000,
  chains = 4,
  warmup = 2000,
  thin = 1,
  control = list(adapt_delta = 0.95)
)

cat("✓ CTL群のサンプリング完了\n")

# 事後分布からサンプルを抽出
ctl_samples_post <- extract(ctl_fit)$theta
cat("  事後分布サンプル数:", nrow(ctl_samples_post), "\n")

# =============================================================
# 6. Log2(FC)の事後分布を計算 -----------------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Log2(FC)の事後分布を計算中...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 論文の方法: 対照群の事後確率分布をATから減算
# つまり、log2(DTA/CTL)の事後分布を計算

# 各細胞タイプについて、Log2(FC)の事後分布を計算
n_post_samples <- min(nrow(dta_samples_post), nrow(ctl_samples_post))
log2fc_posterior <- matrix(0, nrow = n_post_samples, ncol = length(celltype_order))
colnames(log2fc_posterior) <- celltype_order

for (i in 1:n_post_samples) {
  dta_theta <- dta_samples_post[i, ]
  ctl_theta <- ctl_samples_post[i, ]
  log2fc_posterior[i, ] <- log2(dta_theta / ctl_theta)
}

cat("✓ Log2(FC)の事後分布を計算完了\n")

# =============================================================
# 7. 89%信頼区間を計算 -----------------------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("89%信頼区間を計算中...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

ci_level <- 0.89
alpha <- 1 - ci_level  # 0.11

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

for (i in 1:length(celltype_order)) {
  celltype <- celltype_order[i]
  
  # DTA群の事後分布
  dta_post <- dta_samples_post[, i]
  dta_ci <- quantile(dta_post, probs = c(alpha/2, 1 - alpha/2))
  
  # CTL群の事後分布
  ctl_post <- ctl_samples_post[, i]
  ctl_ci <- quantile(ctl_post, probs = c(alpha/2, 1 - alpha/2))
  
  # Log2(FC)の事後分布
  log2fc_post <- log2fc_posterior[, i]
  log2fc_ci <- quantile(log2fc_post, probs = c(alpha/2, 1 - alpha/2))
  
  # 89%信頼区間が0を含まない場合に有意と判断
  significant <- log2fc_ci[1] > 0 || log2fc_ci[2] < 0
  
  bayesian_results <- rbind(
    bayesian_results,
    data.frame(
      CellType = celltype,
      DTA_mean = mean(dta_post),
      DTA_ci_lower = dta_ci[1],
      DTA_ci_upper = dta_ci[2],
      CTL_mean = mean(ctl_post),
      CTL_ci_lower = ctl_ci[1],
      CTL_ci_upper = ctl_ci[2],
      Log2FC_mean = mean(log2fc_post),
      Log2FC_ci_lower = log2fc_ci[1],
      Log2FC_ci_upper = log2fc_ci[2],
      Significant = significant
    )
  )
}

cat("ベイズ統計解析結果（DMM）:\n\n")
print(bayesian_results, row.names = FALSE)

cat("\n")

# =============================================================
# 8. 結果の解釈と表示 -------------------------------------------
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
# 9. 結果をCSVファイルに保存 -----------------------------------
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("結果をCSVファイルに保存中...\n")

output_dir <- dirname(dta_file)
if (output_dir == ".") {
  output_dir <- getwd()
}

# ベイズ統計解析結果を保存
write.csv(bayesian_results,
          file = file.path(output_dir, "統計検定_DMM解析結果.csv"),
          row.names = FALSE)

cat("✓ 結果を保存しました:\n")
cat("  -", file.path(output_dir, "統計検定_DMM解析結果.csv"), "\n")

# =============================================================
# 10. グラフ作成（論文のPanel Dに類似） ------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("グラフを作成中（論文のPanel Dに類似）...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# Ridge plot（ggridgesが必要）
if (!requireNamespace("ggridges", quietly = TRUE)) {
  cat("警告: ggridgesパッケージがインストールされていません。\n")
  cat("      install.packages('ggridges')を実行してください。\n")
} else {
  library(ggridges)
  
  # 事後分布のサンプルをデータフレームに変換
  posterior_samples_all <- data.frame(
    CellType = rep(celltype_order, each = n_post_samples),
    Log2FC = as.vector(log2fc_posterior)
  )
  
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
  
  # Ridge plot
  p_ridge <- ggplot(posterior_samples_all, aes(x = Log2FC, y = CellType, fill = CellType)) +
    geom_density_ridges(alpha = 0.7, scale = 1.2, rel_min_height = 0.01) +
    scale_fill_manual(values = celltype_colors, guide = "none") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    # 89%信頼区間を表示（有意なものは赤、そうでないものはグレー）
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
      title = "Posterior Distribution of Log2(FC) (89% Credible Interval)\nDirichlet Multinomial Modeling"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 14, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      axis.line = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(
    filename = file.path(output_dir, "統計検定_DMM_RidgePlot.png"),
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
cat("※ このスクリプトは論文準拠のDMM（ディリクレ多項式モデリング）解析を行っています。\n")
cat("   rstanを使用したHMCサンプリングにより、89%信頼区間を計算しています。\n")
cat("   信頼区間が0を含まない場合に有意と判断します。\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

