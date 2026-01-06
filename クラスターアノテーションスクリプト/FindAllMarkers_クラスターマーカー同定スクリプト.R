# =============================================================
# FindAllMarkersによるクラスターマーカー遺伝子同定スクリプト
# =============================================================

# 0. ライブラリ -------------------------------------------------
library(Seurat)
library(dplyr)
library(openxlsx)
library(future)

# 並列処理の設定（処理を高速化）
plan("multicore", workers = 4)  # 4コアを使用（環境に応じて調整）
options(future.globals.maxSize = 8000 * 1024^2)  # 8GB

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
  stop("エラー: seu_dta_processed.rdsが見つかりません。\n")
}

if (!file.exists(ctl_file)) {
  stop("エラー: seu_ctl_processed.rdsが見つかりません。\n")
}

cat("✓ データファイルを検出:\n")
cat("  DTA:", dta_file, "\n")
cat("  CTL:", ctl_file, "\n")

seu_dta <- readRDS(dta_file)
seu_ctl <- readRDS(ctl_file)

cat("✓ DTAデータ読み込み完了:", ncol(seu_dta), "cells\n")
cat("✓ CTLデータ読み込み完了:", ncol(seu_ctl), "cells\n")

# seurat_clustersが存在するか確認
if (!"seurat_clusters" %in% colnames(seu_dta@meta.data)) {
  stop("エラー: DTAのseurat_clustersがメタデータに存在しません。\n")
}

if (!"seurat_clusters" %in% colnames(seu_ctl@meta.data)) {
  stop("エラー: CTLのseurat_clustersがメタデータに存在しません。\n")
}

# =============================================================
# 2. FindAllMarkersの実行 --------------------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("FindAllMarkersによるマーカー遺伝子同定\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# パラメータ設定（論文の記述に合わせる）
logfc_threshold <- 0.25  # 最小log2 fold change
min_pct <- 0.1           # 最小発現細胞割合
only_pos <- FALSE        # 正のマーカーのみかどうか

# 高速化のための追加パラメータ
max_cells_per_ident <- 200  # 各クラスターから最大200細胞を使用（計算量削減）

# DTAのマーカー遺伝子同定
cat("\n--- DTA群のマーカー遺伝子同定 ---\n")
Idents(seu_dta) <- "seurat_clusters"

cat("FindAllMarkers実行中（logfc.threshold =", logfc_threshold, 
    ", min.pct =", min_pct, ")...\n")
cat("  並列処理を使用して高速化中...\n")
cat("  これは時間がかかる場合があります（数分〜数十分）...\n")

# 高速化のため、各クラスターからランダムサンプリング（オプション）
# 実際のデータでは全細胞を使用する方が正確ですが、時間短縮のため
markers_dta <- FindAllMarkers(
  object = seu_dta,
  only.pos = only_pos,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold,
  test.use = "wilcox",  # Wilcoxon rank sum test
  max.cells.per.ident = max_cells_per_ident,  # 計算量削減
  verbose = TRUE
)

cat("✓ DTA: マーカー遺伝子同定完了\n")
cat("  クラスター数:", length(unique(markers_dta$cluster)), "\n")
cat("  マーカー遺伝子総数:", nrow(markers_dta), "\n")

# CTLのマーカー遺伝子同定
cat("\n--- CTL群のマーカー遺伝子同定 ---\n")
Idents(seu_ctl) <- "seurat_clusters"

cat("FindAllMarkers実行中（logfc.threshold =", logfc_threshold, 
    ", min.pct =", min_pct, ")...\n")
cat("  並列処理を使用して高速化中...\n")
cat("  これは時間がかかる場合があります（数分〜数十分）...\n")

markers_ctl <- FindAllMarkers(
  object = seu_ctl,
  only.pos = only_pos,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold,
  test.use = "wilcox",
  max.cells.per.ident = max_cells_per_ident,  # 計算量削減
  verbose = TRUE
)

cat("✓ CTL: マーカー遺伝子同定完了\n")
cat("  クラスター数:", length(unique(markers_ctl$cluster)), "\n")
cat("  マーカー遺伝子総数:", nrow(markers_ctl), "\n")

# =============================================================
# 3. 結果のサマリー表示 ----------------------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("結果サマリー\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# DTAのクラスターごとのマーカー数
cat("DTA群 - クラスターごとのマーカー遺伝子数:\n")
dta_summary <- markers_dta %>%
  group_by(cluster) %>%
  summarise(
    n_markers = n(),
    n_pos_markers = sum(avg_log2FC > 0),
    n_neg_markers = sum(avg_log2FC < 0),
    top_gene = gene[which.max(avg_log2FC)],
    max_log2FC = max(avg_log2FC),
    .groups = 'drop'
  ) %>%
  arrange(as.numeric(cluster))
print(dta_summary)

cat("\nCTL群 - クラスターごとのマーカー遺伝子数:\n")
ctl_summary <- markers_ctl %>%
  group_by(cluster) %>%
  summarise(
    n_markers = n(),
    n_pos_markers = sum(avg_log2FC > 0),
    n_neg_markers = sum(avg_log2FC < 0),
    top_gene = gene[which.max(avg_log2FC)],
    max_log2FC = max(avg_log2FC),
    .groups = 'drop'
  ) %>%
  arrange(as.numeric(cluster))
print(ctl_summary)

# =============================================================
# 4. 文献ベースのマーカー遺伝子との比較 ------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("文献ベースのマーカー遺伝子との比較\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 文献ベースのマーカー遺伝子リスト
literature_markers <- list(
  Granule = c("Rims1", "Grm4", "Gabra6"),
  Purkinje = c("Itpr1", "Calb1", "Car8"),
  Bergmann_glia = c("Tubb2b", "Aqp4"),
  Interneurons = c("Gad1", "Pvalb"),
  Astrocytes = c("Ttn", "Aqp1", "Gfap"),
  Oligodendrocytes = c("Plp1", "Mbp"),
  OPC = c("Pdgfra", "Olig1"),
  Microglia = c("Cd74", "Csf1r"),
  Endothelial = c("Cldn5", "Vwf"),
  Fibroblasts = c("Dcn", "Apod")
)

all_literature_markers <- unlist(literature_markers)
names(all_literature_markers) <- NULL

cat("文献ベースのマーカー遺伝子総数:", length(all_literature_markers), "\n")
cat("  遺伝子リスト:", paste(all_literature_markers, collapse = ", "), "\n\n")

# DTAで同定されたマーカーとの重複確認
dta_marker_genes <- unique(markers_dta$gene)
dta_overlap <- intersect(dta_marker_genes, all_literature_markers)
cat("DTA群:\n")
cat("  同定されたマーカー遺伝子数:", length(dta_marker_genes), "\n")
cat("  文献マーカーとの重複数:", length(dta_overlap), "\n")
cat("  重複率:", round(length(dta_overlap) / length(all_literature_markers) * 100, 1), "%\n")
if (length(dta_overlap) > 0) {
  cat("  重複遺伝子:", paste(dta_overlap, collapse = ", "), "\n")
} else {
  cat("  重複遺伝子: なし\n")
}

# CTLで同定されたマーカーとの重複確認
ctl_marker_genes <- unique(markers_ctl$gene)
ctl_overlap <- intersect(ctl_marker_genes, all_literature_markers)
cat("\nCTL群:\n")
cat("  同定されたマーカー遺伝子数:", length(ctl_marker_genes), "\n")
cat("  文献マーカーとの重複数:", length(ctl_overlap), "\n")
cat("  重複率:", round(length(ctl_overlap) / length(all_literature_markers) * 100, 1), "%\n")
if (length(ctl_overlap) > 0) {
  cat("  重複遺伝子:", paste(ctl_overlap, collapse = ", "), "\n")
} else {
  cat("  重複遺伝子: なし\n")
}

# =============================================================
# 5. 結果の保存 ------------------------------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("結果を保存中...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

output_dir <- dirname(dta_file)
if (output_dir == ".") {
  output_dir <- getwd()
}

# CSVファイルとして保存
write.csv(markers_dta, 
          file.path(output_dir, "FindAllMarkers_DTA_結果.csv"),
          row.names = FALSE)
cat("✓ DTAマーカー遺伝子を保存: FindAllMarkers_DTA_結果.csv\n")

write.csv(markers_ctl, 
          file.path(output_dir, "FindAllMarkers_CTL_結果.csv"),
          row.names = FALSE)
cat("✓ CTLマーカー遺伝子を保存: FindAllMarkers_CTL_結果.csv\n")

# Excelファイルとして保存（クラスターごとにシート分け）
wb <- createWorkbook()

# DTAのクラスターごとにシート作成
dta_clusters <- sort(unique(markers_dta$cluster))
for (clust in dta_clusters) {
  clust_data <- markers_dta %>%
    filter(cluster == clust) %>%
    arrange(desc(avg_log2FC))
  
  sheet_name <- paste0("DTA_Cluster_", clust)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, clust_data)
}

# CTLのクラスターごとにシート作成
ctl_clusters <- sort(unique(markers_ctl$cluster))
for (clust in ctl_clusters) {
  clust_data <- markers_ctl %>%
    filter(cluster == clust) %>%
    arrange(desc(avg_log2FC))
  
  sheet_name <- paste0("CTL_Cluster_", clust)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, clust_data)
}

# サマリーシート
addWorksheet(wb, "Summary_DTA")
writeData(wb, "Summary_DTA", dta_summary)

addWorksheet(wb, "Summary_CTL")
writeData(wb, "Summary_CTL", ctl_summary)

# 文献マーカーとの比較シート
comparison_df <- data.frame(
  CellType = names(literature_markers),
  Literature_Markers = sapply(literature_markers, function(x) paste(x, collapse = ", ")),
  stringsAsFactors = FALSE
)
addWorksheet(wb, "Literature_Markers")
writeData(wb, "Literature_Markers", comparison_df)

saveWorkbook(wb, 
             file.path(output_dir, "FindAllMarkers_結果_クラスター別.xlsx"),
             overwrite = TRUE)
cat("✓ Excelファイルを保存: FindAllMarkers_結果_クラスター別.xlsx\n")

# =============================================================
# 6. 上位マーカー遺伝子の表示 ----------------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("各クラスターの上位5マーカー遺伝子（log2FC順）\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

cat("DTA群:\n")
for (clust in dta_clusters) {
  top_markers <- markers_dta %>%
    filter(cluster == clust) %>%
    arrange(desc(avg_log2FC)) %>%
    head(5)
  
  cat("\nCluster", clust, ":\n")
  for (i in 1:nrow(top_markers)) {
    cat(sprintf("  %d. %s (log2FC=%.3f, p_val=%.2e)\n",
                i, top_markers$gene[i], 
                top_markers$avg_log2FC[i],
                top_markers$p_val[i]))
  }
}

cat("\nCTL群:\n")
for (clust in ctl_clusters) {
  top_markers <- markers_ctl %>%
    filter(cluster == clust) %>%
    arrange(desc(avg_log2FC)) %>%
    head(5)
  
  cat("\nCluster", clust, ":\n")
  for (i in 1:nrow(top_markers)) {
    cat(sprintf("  %d. %s (log2FC=%.3f, p_val=%.2e)\n",
                i, top_markers$gene[i], 
                top_markers$avg_log2FC[i],
                top_markers$p_val[i]))
  }
}

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("✓ 解析完了！\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

