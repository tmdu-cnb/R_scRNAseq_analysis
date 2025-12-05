# =============================================================
# 細胞構成比比較スクリプト
# CTLとDTAの細胞タイプ構成比を比較・可視化
# =============================================================

# 0. ライブラリ -------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)
library(ggridges)  # ridge plot用
library(tidyr)     # pivot_wider, unnest用

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

# 出力ディレクトリを設定（データファイルと同じディレクトリ）
output_dir <- dirname(dta_file)
if (output_dir == ".") {
  output_dir <- getwd()
}
cat("✓ 出力ディレクトリ:", output_dir, "\n")

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
    data_norm <- GetAssayData(seu_obj, layer = "data")  # 正規化データも使用
    
    # 利用可能な遺伝子リストを確認
    available_genes <- rownames(counts)
    
    # 各細胞型のマーカー発現をチェック（スコアベース）
    cell_type <- rep("Unclassified", ncol(seu_obj))
    scores <- matrix(0, nrow = ncol(seu_obj), ncol = 10)
    colnames(scores) <- c("Granule", "Purkinje", "Bergmann_glia", "Interneurons",
                           "Astrocytes", "Oligodendrocytes", "OPC", 
                           "Microglia", "Endothelial", "Fibroblasts")
    
    # 1. Granule cells: RIMS1, GRM4, GABRA6
    granule_markers <- c("Rims1", "Grm4", "Gabra6")
    granule_markers <- granule_markers[granule_markers %in% available_genes]
    if (length(granule_markers) > 0) {
      scores[, "Granule"] <- colMeans(data_norm[granule_markers, , drop = FALSE])
    }
    
    # 2. Purkinje cells: ITPR1, CALB1, CAR8
    purkinje_markers <- c("Itpr1", "Calb1", "Car8")
    purkinje_markers <- purkinje_markers[purkinje_markers %in% available_genes]
    if (length(purkinje_markers) > 0) {
      scores[, "Purkinje"] <- colMeans(data_norm[purkinje_markers, , drop = FALSE])
    }
    
    # 3. Bergmann glia: TUBB2B, AQP4
    bergmann_markers <- c("Tubb2b", "Aqp4")
    bergmann_markers <- bergmann_markers[bergmann_markers %in% available_genes]
    if (length(bergmann_markers) > 0) {
      scores[, "Bergmann_glia"] <- colMeans(data_norm[bergmann_markers, , drop = FALSE])
    }
    
    # 4. Interneurons: GAD1, PVALB
    inter_markers <- c("Gad1", "Pvalb")
    inter_markers <- inter_markers[inter_markers %in% available_genes]
    if (length(inter_markers) > 0) {
      scores[, "Interneurons"] <- colMeans(data_norm[inter_markers, , drop = FALSE])
    }
    
    # 5. Astrocytes: TTN, AQP1 (+ GFAP)
    astro_markers <- c("Ttn", "Aqp1", "Gfap")
    astro_markers <- astro_markers[astro_markers %in% available_genes]
    if (length(astro_markers) > 0) {
      scores[, "Astrocytes"] <- colMeans(data_norm[astro_markers, , drop = FALSE])
    }
    
    # 6. Oligodendrocytes: PLP1, MBP
    oligo_markers <- c("Plp1", "Mbp")
    oligo_markers <- oligo_markers[oligo_markers %in% available_genes]
    if (length(oligo_markers) > 0) {
      scores[, "Oligodendrocytes"] <- colMeans(data_norm[oligo_markers, , drop = FALSE])
    }
    
    # 7. OPC: PDGFRA, OLIG1
    opc_markers <- c("Pdgfra", "Olig1")
    opc_markers <- opc_markers[opc_markers %in% available_genes]
    if (length(opc_markers) > 0) {
      scores[, "OPC"] <- colMeans(data_norm[opc_markers, , drop = FALSE])
    }
    
    # 8. Microglia: CD74, CSF1R
    micro_markers <- c("Cd74", "Csf1r")
    micro_markers <- micro_markers[micro_markers %in% available_genes]
    if (length(micro_markers) > 0) {
      scores[, "Microglia"] <- colMeans(data_norm[micro_markers, , drop = FALSE])
    }
    
    # 9. Endothelial cells: CLDN5, VWF
    endo_markers <- c("Cldn5", "Vwf")
    endo_markers <- endo_markers[endo_markers %in% available_genes]
    if (length(endo_markers) > 0) {
      scores[, "Endothelial"] <- colMeans(data_norm[endo_markers, , drop = FALSE])
    }
    
    # 10. Fibroblasts: DCN, APOD
    fibro_markers <- c("Dcn", "Apod")
    fibro_markers <- fibro_markers[fibro_markers %in% available_genes]
    if (length(fibro_markers) > 0) {
      scores[, "Fibroblasts"] <- colMeans(data_norm[fibro_markers, , drop = FALSE])
    }
    
    # 全ての細胞を必ず分類（閾値なし）
    max_scores <- apply(scores, 1, max)
    max_types <- apply(scores, 1, which.max)
    
    # 全細胞に最もスコアが高いタイプを割り当て
    for (i in 1:length(cell_type)) {
      cell_type[i] <- colnames(scores)[max_types[i]]
    }
    
    return(list(celltype = cell_type, scores = scores))
  }
  
  # DTAとCTLに適用
  if (!"celltype_10class" %in% colnames(seu_dta@meta.data)) {
    cat("  DTAデータを分類中...\n")
    result_dta <- classify_cerebellar_cells(seu_dta)
    seu_dta$celltype_10class <- result_dta$celltype
    seu_dta$celltype_10class_score <- apply(result_dta$scores, 1, max)
    cat("  ✓ DTA分類完了\n")
  }
  
  if (!"celltype_10class" %in% colnames(seu_ctl@meta.data)) {
    cat("  CTLデータを分類中...\n")
    result_ctl <- classify_cerebellar_cells(seu_ctl)
    seu_ctl$celltype_10class <- result_ctl$celltype
    seu_ctl$celltype_10class_score <- apply(result_ctl$scores, 1, max)
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

# 細胞タイプの順序を定義（画像に合わせて）
celltype_order <- c(
  "Granule", "Oligodendrocytes", "Bergmann_glia", "Astrocytes",
  "Interneurons", "Endothelial", "Purkinje", "OPC",
  "Microglia", "Fibroblasts"
)

# 存在しない細胞タイプを除外
celltype_order <- celltype_order[celltype_order %in% unique(prop_data$CellType)]

# 因子として設定
prop_data$CellType <- factor(prop_data$CellType, levels = celltype_order)
prop_data$Group <- factor(prop_data$Group, levels = c("DTA", "CTL"))

cat("✓ 構成比計算完了\n")

# =============================================================
# 3. 色の定義 ---------------------------------------------------
# 画像に合わせた色設定
celltype_colors <- c(
  "Granule" = "#E31A1C",           # 赤
  "Oligodendrocytes" = "#FF7F00",   # オレンジ
  "Bergmann_glia" = "#B2DF8A",     # 黄緑
  "Astrocytes" = "#33A02C",         # 緑
  "Interneurons" = "#1F78B4",       # ダークグリーン
  "Endothelial" = "#A6CEE3",        # ティール
  "Purkinje" = "#6A3D9A",           # ライトブルー
  "OPC" = "#CAB2D6",                # 紫
  "Microglia" = "#FB9A99",          # ピンク
  "Fibroblasts" = "#E7298A"         # ダークピンク
)

# =============================================================
# 4. Panel C: 積み上げ棒グラフ（構成比） ------------------------
cat("\nPanel C: 積み上げ棒グラフを作成中...\n")

p_proportion <- ggplot(prop_data, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = celltype_colors, name = "Cell Type") +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    x = "",
    y = "Proportion",
    title = "C"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# 拡大図用のデータ（0.75以上のみ）
prop_data_zoom <- prop_data |>
  filter(Proportion >= 0.75 | (Group == "DTA" & Proportion < 0.75) | (Group == "CTL" & Proportion < 0.75))

# 拡大図（inset）
p_proportion_zoom <- ggplot(prop_data, aes(x = Group, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  scale_y_continuous(
    labels = scales::percent_format(),
    limits = c(0.75, 1.0),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(x = "", y = "") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    plot.background = element_rect(color = "black", linewidth = 1.5, fill = "white")
  )

cat("✓ Panel C作成完了\n")

# =============================================================
# 5. Panel D: Log2(FC) Ridge Plot ------------------------------
cat("\nPanel D: Log2(FC) Ridge Plotを作成中...\n")

# Log2(FC)を計算
# FC = DTAの割合 / CTLの割合
# Log2(FC) = log2(DTAの割合 / CTLの割合)

# Log2(FC)を計算するためのデータ準備
log2fc_data <- prop_data |>
  select(Group, CellType, Count, Proportion) |>
  pivot_wider(names_from = Group, values_from = c(Count, Proportion)) |>
  mutate(
    Log2FC = log2(Proportion_DTA / Proportion_CTL),
    DTA_total = sum(prop_data$Count[prop_data$Group == "DTA"]),
    CTL_total = sum(prop_data$Count[prop_data$Group == "CTL"])
  )

# ブートストラップサンプリングで分布を生成
set.seed(123)  # 再現性のため
log2fc_samples_list <- list()

for (i in 1:nrow(log2fc_data)) {
  row <- log2fc_data[i, ]
  
  # 各グループで1000回ブートストラップサンプリング
  n_bootstrap <- 1000
  dta_samples <- rbinom(n_bootstrap, row$DTA_total, row$Proportion_DTA) / row$DTA_total
  ctl_samples <- rbinom(n_bootstrap, row$CTL_total, row$Proportion_CTL) / row$CTL_total
  
  # ゼロ除算を避ける
  fc_samples <- dta_samples / (ctl_samples + 1e-10)
  log2fc_samples <- log2(fc_samples)
  
  log2fc_samples_list[[i]] <- data.frame(
    CellType = row$CellType,
    Log2FC = log2fc_samples
  )
}

# サンプルデータを結合
log2fc_expanded <- bind_rows(log2fc_samples_list)

# log2fc_expandedは上で作成済み

# 細胞タイプの順序を逆順に（画像では下から上へ）
log2fc_expanded$CellType <- factor(
  log2fc_expanded$CellType,
  levels = rev(celltype_order)
)

# Ridge plotを作成
p_ridge <- ggplot(log2fc_expanded, aes(x = Log2FC, y = CellType, fill = CellType)) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 1.2,
    rel_min_height = 0.01
  ) +
  scale_fill_manual(values = celltype_colors, guide = "none") +
  scale_x_continuous(
    limits = c(-2, 2),
    breaks = seq(-2, 2, by = 0.5),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  annotate(
    "text",
    x = -1.5, y = length(celltype_order) + 0.5,
    label = "Decreased in DTA",
    size = 5,
    fontface = "bold"
  ) +
  annotate(
    "segment",
    x = -1.5, xend = -1.8, y = length(celltype_order) + 0.5, yend = length(celltype_order) + 0.5,
    arrow = arrow(length = unit(0.3, "cm")),
    linewidth = 1
  ) +
  annotate(
    "text",
    x = 1.5, y = length(celltype_order) + 0.5,
    label = "Increased in DTA",
    size = 5,
    fontface = "bold"
  ) +
  annotate(
    "segment",
    x = 1.5, xend = 1.8, y = length(celltype_order) + 0.5, yend = length(celltype_order) + 0.5,
    arrow = arrow(length = unit(0.3, "cm")),
    linewidth = 1
  ) +
  labs(
    x = "Log2(FC)",
    y = "",
    title = "D"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

cat("✓ Panel D作成完了\n")

# =============================================================
# 6. 図を結合して保存 -------------------------------------------
cat("\n図を結合して保存中...\n")

# 図を結合
p_combined <- p_proportion + p_ridge +
  plot_layout(ncol = 2, widths = c(1, 1.5))

# 保存（明示的なパスを指定）
output_file1 <- file.path(output_dir, "細胞構成比比較.png")
output_file2 <- file.path(output_dir, "細胞構成比_積み上げ棒グラフ.png")
output_file3 <- file.path(output_dir, "細胞構成比_Log2FC_RidgePlot.png")

ggsave(
  filename = output_file1,
  plot = p_combined,
  width = 16,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("✓ 図を保存しました:", output_file1, "\n")

# 個別にも保存
ggsave(
  filename = output_file2,
  plot = p_proportion,
  width = 8,
  height = 8,
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = output_file3,
  plot = p_ridge,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("✓ 個別図も保存しました:\n")
cat("  -", output_file2, "\n")
cat("  -", output_file3, "\n")

# =============================================================
# 7. 統計情報を出力 --------------------------------------------
cat("\n=== 細胞タイプ構成比の統計 ===\n")
cat("\n【DTA群】\n")
for (i in 1:nrow(prop_data[prop_data$Group == "DTA", ])) {
  row <- prop_data[prop_data$Group == "DTA", ][i, ]
  cat(sprintf("  %-18s: %6d cells (%5.2f%%)\n",
              row$CellType, row$Count, row$Proportion * 100))
}

cat("\n【CTL群】\n")
for (i in 1:nrow(prop_data[prop_data$Group == "CTL", ])) {
  row <- prop_data[prop_data$Group == "CTL", ][i, ]
  cat(sprintf("  %-18s: %6d cells (%5.2f%%)\n",
              row$CellType, row$Count, row$Proportion * 100))
}

cat("\n【Log2(FC)】\n")
log2fc_summary <- log2fc_data |>
  select(CellType, Log2FC) |>
  arrange(desc(Log2FC))
for (i in 1:nrow(log2fc_summary)) {
  cat(sprintf("  %-18s: %6.3f\n",
              log2fc_summary$CellType[i], log2fc_summary$Log2FC[i]))
}

cat("\n解析完了！\n")

