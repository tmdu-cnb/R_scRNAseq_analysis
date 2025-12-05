# =============================================================
# 細胞構成比統計検定スクリプト
# DTAとCTLの細胞タイプ構成比の統計的有意差を検定
# =============================================================

# 0. ライブラリ -------------------------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)  # ラベル重複回避用

# =============================================================
# 1. 処理済みデータの読み込み ------------------------------------
cat("処理済みデータを読み込み中...\n")

# データ読み込み（複数のパスを試す）
data_paths <- c(
  "seu_dta_processed.rds",  # 現在の作業ディレクトリ
  file.path(getwd(), "seu_dta_processed.rds"),  # 明示的に作業ディレクトリ
  "/Users/godakyosuke/Desktop/goda sc RNA seq/seu_dta_processed.rds"  # 絶対パス
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
    
    # 各細胞タイプのマーカー遺伝子でスコア計算
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
# 2. 細胞タイプの構成比を計算 -----------------------------------
cat("\n細胞タイプの構成比を計算中...\n")

# DTAとCTLの細胞タイプ数をカウント
dta_counts <- table(seu_dta$celltype_10class)
ctl_counts <- table(seu_ctl$celltype_10class)

# 構成比を計算
dta_prop <- prop.table(dta_counts)
ctl_prop <- prop.table(ctl_counts)

# データフレームに整理
prop_data <- bind_rows(
  data.frame(
    Group = "DTA",
    CellType = names(dta_prop),
    Count = as.numeric(dta_counts),
    Proportion = as.numeric(dta_prop)
  ),
  data.frame(
    Group = "CTL",
    CellType = names(ctl_prop),
    Count = as.numeric(ctl_counts),
    Proportion = as.numeric(ctl_prop)
  )
)

# 細胞タイプの順序を定義
celltype_order <- c(
  "Granule", "Oligodendrocytes", "Bergmann_glia", "Astrocytes",
  "Interneurons", "Endothelial", "Purkinje", "OPC",
  "Microglia", "Fibroblasts"
)

celltype_order <- celltype_order[celltype_order %in% unique(prop_data$CellType)]
prop_data$CellType <- factor(prop_data$CellType, levels = celltype_order)
prop_data$Group <- factor(prop_data$Group, levels = c("DTA", "CTL"))

cat("✓ 構成比計算完了\n")

# =============================================================
# 3. カイ二乗検定（全体の構成比） -----------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("【統計検定1: カイ二乗検定（全体の構成比）】\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# データを2×nの分割表に変換
dta_counts_vec <- prop_data$Count[prop_data$Group == "DTA"]
ctl_counts_vec <- prop_data$Count[prop_data$Group == "CTL"]

# 分割表を作成
contingency_table <- matrix(
  c(dta_counts_vec, ctl_counts_vec),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Group = c("DTA", "CTL"),
    CellType = celltype_order
  )
)

cat("分割表:\n")
print(contingency_table)
cat("\n")

# カイ二乗検定
chi2_test <- chisq.test(contingency_table)

cat("カイ二乗検定の結果:\n")
cat(sprintf("  カイ二乗統計量: %.4f\n", chi2_test$statistic))
cat(sprintf("  自由度: %d\n", chi2_test$parameter))
cat(sprintf("  p-value: %.4e\n", chi2_test$p.value))
cat("\n")

# 解釈
if (chi2_test$p.value < 0.001) {
  cat("  → 統計的に非常に有意な差があります（p < 0.001）\n")
} else if (chi2_test$p.value < 0.01) {
  cat("  → 統計的に有意な差があります（p < 0.01）\n")
} else if (chi2_test$p.value < 0.05) {
  cat("  → 統計的に有意な差があります（p < 0.05）\n")
} else {
  cat("  → 統計的に有意な差はありません（p >= 0.05）\n")
}

# 期待度数を確認
cat("\n期待度数:\n")
print(chi2_test$expected)

# 期待度数が5未満のセルがあるか確認
if (any(chi2_test$expected < 5)) {
  cat("\n⚠ 警告: 期待度数が5未満のセルがあります。\n")
  cat("  Fisher's exact testの使用を検討してください。\n")
  
  # Fisher's exact testも実行
  cat("\n【Fisher's exact test（参考）】\n")
  fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE, B = 10000)
  cat(sprintf("  p-value: %.4e\n", fisher_test$p.value))
}

cat("\n")

# =============================================================
# 4. 各細胞タイプごとの二項検定 + FDR補正 ----------------------
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("【統計検定2: 各細胞タイプごとの二項検定 + FDR補正】\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

dta_total <- sum(dta_counts_vec)
ctl_total <- sum(ctl_counts_vec)

cat(sprintf("DTA総細胞数: %d\n", dta_total))
cat(sprintf("CTL総細胞数: %d\n\n", ctl_total))

# 各細胞タイプごとに二項検定
binomial_results <- data.frame(
  CellType = character(),
  DTA_count = numeric(),
  DTA_prop = numeric(),
  CTL_count = numeric(),
  CTL_prop = numeric(),
  p_value = numeric(),
  p_adj_fdr = numeric(),
  p_adj_bonf = numeric(),
  stringsAsFactors = FALSE
)

for (i in 1:length(celltype_order)) {
  celltype <- celltype_order[i]
  
  # DTA群での該当細胞タイプの数
  dta_count <- prop_data$Count[prop_data$Group == "DTA" & prop_data$CellType == celltype]
  dta_prop <- prop_data$Proportion[prop_data$Group == "DTA" & prop_data$CellType == celltype]
  
  # CTL群での該当細胞タイプの数
  ctl_count <- prop_data$Count[prop_data$Group == "CTL" & prop_data$CellType == celltype]
  ctl_prop <- prop_data$Proportion[prop_data$Group == "CTL" & prop_data$CellType == celltype]
  
  # 二項検定: DTAの割合がCTLの割合と異なるか
  binom_test <- binom.test(
    x = dta_count,
    n = dta_total,
    p = ctl_prop,
    alternative = "two.sided"
  )
  
  binomial_results <- rbind(
    binomial_results,
    data.frame(
      CellType = celltype,
      DTA_count = dta_count,
      DTA_prop = dta_prop,
      CTL_count = ctl_count,
      CTL_prop = ctl_prop,
      p_value = binom_test$p.value,
      p_adj_fdr = NA,
      p_adj_bonf = NA
    )
  )
}

# 多重比較補正
binomial_results$p_adj_fdr <- p.adjust(binomial_results$p_value, method = "fdr")
binomial_results$p_adj_bonf <- p.adjust(binomial_results$p_value, method = "bonferroni")

# Log2(FC)も計算
binomial_results$Log2FC <- log2(binomial_results$DTA_prop / binomial_results$CTL_prop)

# 結果を表示
cat("各細胞タイプごとの二項検定結果:\n\n")
print(binomial_results, row.names = FALSE)

cat("\n")

# 有意な差がある細胞タイプを表示
cat("【有意な差がある細胞タイプ（FDR補正後、q < 0.05）】\n")
significant_fdr <- binomial_results[binomial_results$p_adj_fdr < 0.05, ]
if (nrow(significant_fdr) > 0) {
  significant_fdr <- significant_fdr[order(significant_fdr$p_adj_fdr), ]
  for (i in 1:nrow(significant_fdr)) {
    row <- significant_fdr[i, ]
    direction <- ifelse(row$Log2FC > 0, "増加", "減少")
    cat(sprintf("  %-18s: p = %.4e, FDR補正後 = %.4e, Log2(FC) = %.3f (%s)\n",
                row$CellType,
                row$p_value,
                row$p_adj_fdr,
                row$Log2FC,
                direction))
  }
} else {
  cat("  なし\n")
}

cat("\n")

cat("【有意な差がある細胞タイプ（Bonferroni補正後、p < 0.05）】\n")
significant_bonf <- binomial_results[binomial_results$p_adj_bonf < 0.05, ]
if (nrow(significant_bonf) > 0) {
  significant_bonf <- significant_bonf[order(significant_bonf$p_adj_bonf), ]
  for (i in 1:nrow(significant_bonf)) {
    row <- significant_bonf[i, ]
    direction <- ifelse(row$Log2FC > 0, "増加", "減少")
    cat(sprintf("  %-18s: p = %.4e, Bonferroni補正後 = %.4e, Log2(FC) = %.3f (%s)\n",
                row$CellType,
                row$p_value,
                row$p_adj_bonf,
                row$Log2FC,
                direction))
  }
} else {
  cat("  なし\n")
}

cat("\n")

# =============================================================
# 5. 結果をCSVファイルに保存 -----------------------------------
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("結果をCSVファイルに保存中...\n")

# 出力ディレクトリを設定
output_dir <- dirname(dta_file)
if (output_dir == ".") {
  output_dir <- getwd()
}

# カイ二乗検定の結果を保存
chi2_summary <- data.frame(
  Test = "Chi-square test",
  Statistic = as.numeric(chi2_test$statistic),
  df = chi2_test$parameter,
  p_value = chi2_test$p.value,
  Significant = ifelse(chi2_test$p.value < 0.05, "Yes", "No")
)

write.csv(chi2_summary, 
          file = file.path(output_dir, "統計検定_カイ二乗検定結果.csv"),
          row.names = FALSE)

# 二項検定の結果を保存
write.csv(binomial_results, 
          file = file.path(output_dir, "統計検定_二項検定結果.csv"),
          row.names = FALSE)

cat("✓ 結果を保存しました:\n")
cat("  -", file.path(output_dir, "統計検定_カイ二乗検定結果.csv"), "\n")
cat("  -", file.path(output_dir, "統計検定_二項検定結果.csv"), "\n")

# =============================================================
# 6. グラフ作成 -------------------------------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("グラフを作成中...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# グラフ用のデータ準備
plot_data <- binomial_results |>
  mutate(
    neg_log10_p = -log10(p_value),
    neg_log10_p_adj_fdr = -log10(p_adj_fdr),
    Significant = ifelse(p_adj_fdr < 0.05, "Yes", "No"),
    Direction = ifelse(Log2FC > 0, "Increased in DTA", "Decreased in DTA")
  )

# 細胞タイプの順序をLog2(FC)でソート
plot_data$CellType <- factor(
  plot_data$CellType,
  levels = plot_data$CellType[order(plot_data$Log2FC)]
)

# 1. ボルケーノプロット風の図（Log2(FC) vs -log10(p値)）
cat("1. ボルケーノプロットを作成中...\n")

p_volcano <- ggplot(plot_data, aes(x = Log2FC, y = neg_log10_p_adj_fdr)) +
  geom_point(
    aes(color = Significant, fill = Significant),
    size = 4,
    alpha = 0.7,
    shape = 21
  ) +
  geom_text_repel(
    data = plot_data[plot_data$p_adj_fdr < 0.05, ],
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
    title = "Volcano Plot: Cell Type Proportion Changes"
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
  filename = file.path(output_dir, "統計検定_ボルケーノプロット.png"),
  plot = p_volcano,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("  ✓ ボルケーノプロットを保存しました\n")

# 2. Log2(FC)のバープロット（有意なものをハイライト）
cat("2. Log2(FC)バープロットを作成中...\n")

# 細胞タイプをLog2(FC)でソート
plot_data_bar <- plot_data |>
  arrange(Log2FC) |>
  mutate(
    CellType = factor(CellType, levels = CellType),
    BarColor = ifelse(p_adj_fdr < 0.05, ifelse(Log2FC > 0, "Increased", "Decreased"), "Not significant")
  )

p_bar_log2fc <- ggplot(plot_data_bar, aes(x = CellType, y = Log2FC, fill = BarColor)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
  scale_fill_manual(
    values = c("Increased" = "#E31A1C", "Decreased" = "#1F78B4", "Not significant" = "gray70"),
    name = "Significance",
    labels = c("Increased (FDR < 0.05)", "Decreased (FDR < 0.05)", "Not significant")
  ) +
  labs(
    x = "Cell Type",
    y = "Log2(FC) (DTA / CTL)",
    title = "Log2 Fold Change of Cell Type Proportions"
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
  filename = file.path(output_dir, "統計検定_Log2FCバープロット.png"),
  plot = p_bar_log2fc,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("  ✓ Log2(FC)バープロットを保存しました\n")

# 3. p値のバープロット（FDR補正前後を比較）
cat("3. p値バープロットを作成中...\n")

plot_data_p <- plot_data_bar |>
  select(CellType, p_value, p_adj_fdr) |>
  pivot_longer(
    cols = c(p_value, p_adj_fdr),
    names_to = "Test",
    values_to = "p_value_plot"
  ) |>
  mutate(
    Test = factor(Test, levels = c("p_value", "p_adj_fdr"),
                  labels = c("p-value (raw)", "p-value (FDR adjusted)")),
    neg_log10_p = -log10(p_value_plot),
    Significant = ifelse(p_value_plot < 0.05, "Yes", "No")
  )

p_bar_pvalue <- ggplot(plot_data_p, aes(x = CellType, y = neg_log10_p, fill = Test)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(
    values = c("p-value (raw)" = "gray50", "p-value (FDR adjusted)" = "steelblue"),
    name = ""
  ) +
  labs(
    x = "Cell Type",
    y = "-log10(p-value)",
    title = "P-values: Raw vs FDR Adjusted"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  filename = file.path(output_dir, "統計検定_p値バープロット.png"),
  plot = p_bar_pvalue,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("  ✓ p値バープロットを保存しました\n")

# 4. 構成比の比較バープロット（有意なものをハイライト）
cat("4. 構成比比較バープロットを作成中...\n")

plot_data_prop <- prop_data |>
  left_join(
    binomial_results |> select(CellType, p_adj_fdr),
    by = "CellType"
  ) |>
  mutate(
    Significant = p_adj_fdr < 0.05,
    CellType = factor(CellType, levels = levels(plot_data_bar$CellType)),
    BarColor = ifelse(Significant, "Significant", "Not significant")
  )

p_bar_prop <- ggplot(plot_data_prop, aes(x = CellType, y = Proportion, fill = interaction(Group, BarColor))) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(
    values = c(
      "DTA.Significant" = "#E31A1C",
      "CTL.Significant" = "#1F78B4",
      "DTA.Not significant" = "gray70",
      "CTL.Not significant" = "gray90"
    ),
    name = "",
    labels = c("DTA (significant)", "CTL (significant)", "DTA (not significant)", "CTL (not significant)")
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = "Cell Type",
    y = "Proportion",
    title = "Cell Type Proportions: DTA vs CTL"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 14, face = "bold"),
    legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  filename = file.path(output_dir, "統計検定_構成比比較バープロット.png"),
  plot = p_bar_prop,
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("  ✓ 構成比比較バープロットを保存しました\n")

cat("\n✓ グラフを保存しました:\n")
cat("  -", file.path(output_dir, "統計検定_ボルケーノプロット.png"), "\n")
cat("  -", file.path(output_dir, "統計検定_Log2FCバープロット.png"), "\n")
cat("  -", file.path(output_dir, "統計検定_p値バープロット.png"), "\n")
cat("  -", file.path(output_dir, "統計検定_構成比比較バープロット.png"), "\n")

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("解析完了！\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

