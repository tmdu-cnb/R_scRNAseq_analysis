# =============================================================
# 細胞タイプマーカー遺伝子ヒートマップ作成スクリプト（クラスター版）
# 横軸をCTLのクラスターに変更したバージョン
# =============================================================

# 0. ライブラリ -------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

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

# 出力ディレクトリを設定
output_dir <- dirname(dta_file)
if (output_dir == ".") {
  output_dir <- getwd()
}
cat("✓ 出力ディレクトリ:", output_dir, "\n")

seu_ctl <- readRDS(ctl_file)
cat("✓ CTLデータ読み込み完了:", ncol(seu_ctl), "cells\n")

# seurat_clustersが存在するか確認
if (!"seurat_clusters" %in% colnames(seu_ctl@meta.data)) {
  stop("エラー: seurat_clustersがメタデータに存在しません。\n",
       "      UMAPメイン作成スクリプト_データ読み込み.Rを先に実行してください。")
}

# クラスター数を確認
# seurat_clustersは文字列として保存されている可能性があるため、文字列として扱う
# 数値順にソートする
clusters_raw <- unique(as.character(seu_ctl$seurat_clusters))
clusters <- sort(as.numeric(clusters_raw))
clusters <- as.character(clusters)  # 文字列に戻す
cat("✓ CTLクラスター数:", length(clusters), "\n")
cat("  クラスター番号:", paste(clusters, collapse = ", "), "\n\n")

# =============================================================
# 2. 指定されたマーカー遺伝子を使用 --------------------------
cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("指定されたマーカー遺伝子を使用します...\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# 細胞タイプの順序を定義
celltype_order <- c(
  "Granule", "Oligodendrocytes", "Bergmann_glia", "Astrocytes",
  "Interneurons", "Endothelial", "Purkinje", "OPC",
  "Microglia", "Fibroblasts"
)

# 指定されたマーカー遺伝子（論文準拠）
top_markers_list <- list(
  "Granule" = c("Etv1", "Rims1", "Scn2a", "Grm4"),
  "Oligodendrocytes" = c("Plp1", "Mbp", "Mobp", "Cnp"),
  "Bergmann_glia" = c("Slc1a3", "Ndrg2", "Aqp4", "Tubb2b"),
  "Astrocytes" = c("Ttn", "Gfap", "Aldh1l1", "S100b"),
  "Interneurons" = c("Kit", "Gad1", "Gad2", "Pvalb"),
  "Endothelial" = c("Cldn5", "Pecam1", "Vwf", "Flt1"),
  "Purkinje" = c("Itpr1", "Calb1", "Pcp4", "Car8"),
  "OPC" = c("Ptprz1", "Vcan", "Tnr", "Pdgfra"),
  "Microglia" = c("Cd74", "C1qb", "C3", "Csf1r"),
  "Fibroblasts" = c("Dcn", "Apod", "Tagln", "Fbln1")
)

# データセット内の実際の遺伝子名を確認
available_genes <- rownames(GetAssayData(seu_ctl, layer = "data"))

# 遺伝子名のマッピング（大文字小文字を考慮）
map_gene_name <- function(gene_name, available_genes) {
  # 完全一致
  if (gene_name %in% available_genes) {
    return(gene_name)
  }
  # 大文字小文字を無視して検索
  matched <- available_genes[tolower(available_genes) == tolower(gene_name)]
  if (length(matched) > 0) {
    return(matched[1])
  }
  # 見つからない場合は元の名前を返す（警告は後で出す）
  return(gene_name)
}

# 各マーカー遺伝子を実際の遺伝子名にマッピング
top_markers_list_mapped <- list()
for (celltype in names(top_markers_list)) {
  mapped_genes <- sapply(top_markers_list[[celltype]], function(g) {
    map_gene_name(g, available_genes)
  })
  # 利用可能な遺伝子のみを保持
  mapped_genes <- mapped_genes[mapped_genes %in% available_genes]
  top_markers_list_mapped[[celltype]] <- mapped_genes
  
  if (length(mapped_genes) < length(top_markers_list[[celltype]])) {
    missing <- setdiff(top_markers_list[[celltype]], mapped_genes)
    cat(sprintf("警告: %sの以下のマーカー遺伝子が見つかりません: %s\n", 
                celltype, paste(missing, collapse = ", ")))
  }
  cat(sprintf("  %-18s: %d markers (指定: %d)\n", 
              celltype, length(mapped_genes), length(top_markers_list[[celltype]])))
}

# 選択されたマーカー遺伝子のリストを結合
top_markers <- unlist(top_markers_list_mapped)
cat("\n✓ 合計", length(top_markers), "個のマーカー遺伝子を使用します\n\n")

# =============================================================
# 3. ヒートマップ用のデータを準備（クラスター単位） ------------
cat("ヒートマップ用のデータを準備中（クラスター単位）...\n")

# Identsをクラスターに設定
Idents(seu_ctl) <- "seurat_clusters"

# 正規化データを取得
data_norm <- GetAssayData(seu_ctl, layer = "data")

# 選択されたマーカー遺伝子のみを含むデータを抽出
available_markers <- top_markers[top_markers %in% rownames(data_norm)]
if (length(available_markers) < length(top_markers)) {
  cat("警告: 一部のマーカー遺伝子がデータに存在しません\n")
  cat("  選択された:", length(top_markers), "個\n")
  cat("  利用可能:", length(available_markers), "個\n")
}

# 各クラスターごとに平均発現量を計算
heatmap_data <- matrix(0, nrow = length(available_markers), ncol = length(clusters))
rownames(heatmap_data) <- available_markers
colnames(heatmap_data) <- paste0("Cluster_", clusters)

for (cluster_id in clusters) {
  cluster_name <- paste0("Cluster_", cluster_id)
  # seurat_clustersは文字列として扱う
  cells_in_cluster <- WhichCells(seu_ctl, idents = cluster_id)
  if (length(cells_in_cluster) > 0) {
    # そのクラスターの平均発現量を計算
    mean_expr <- rowMeans(data_norm[available_markers, cells_in_cluster, drop = FALSE])
    heatmap_data[, cluster_name] <- mean_expr
    cat(sprintf("  %s: %d cells\n", cluster_name, length(cells_in_cluster)))
  }
}

# z-scoreスケーリング（中央揃えおよびスケール表示）
# 各行（遺伝子）について、平均を引いて標準偏差で割る
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# マーカー遺伝子を細胞タイプごとにグループ化して並び替え
gene_order <- c()
for (celltype in celltype_order) {
  if (celltype %in% names(top_markers_list_mapped)) {
    markers <- top_markers_list_mapped[[celltype]]
    markers_available <- markers[markers %in% available_markers]
    gene_order <- c(gene_order, markers_available)
  }
}

# マーカー遺伝子の順序を反映
heatmap_data_scaled <- heatmap_data_scaled[gene_order, ]

cat("✓ ヒートマップ用データの準備完了\n")
cat("  データサイズ:", nrow(heatmap_data_scaled), "genes x", ncol(heatmap_data_scaled), "clusters\n\n")

# =============================================================
# 4. クラスターの細胞タイプ構成を確認 --------------------------
cat("各クラスターの細胞タイプ構成を確認中...\n")

# celltype_10classが存在するか確認
if ("celltype_10class" %in% colnames(seu_ctl@meta.data)) {
  cat("\n各クラスターの主要な細胞タイプ:\n")
  for (cluster_id in clusters) {
    cluster_cells <- as.character(seu_ctl$seurat_clusters) == cluster_id
    if (sum(cluster_cells) > 0) {
      celltype_counts <- table(seu_ctl$celltype_10class[cluster_cells])
      celltype_prop <- prop.table(celltype_counts) * 100
      top_celltype <- names(celltype_counts)[which.max(celltype_counts)]
      top_prop <- celltype_prop[top_celltype]
      cat(sprintf("  Cluster %s: %s (%.1f%%, n=%d)\n", 
                  cluster_id, top_celltype, top_prop, sum(cluster_cells)))
    }
  }
} else {
  cat("警告: celltype_10classがメタデータに存在しません。\n")
  cat("      クラスターの細胞タイプ情報は表示されません。\n")
}

cat("\n")

# =============================================================
# 5. ヒートマップを作成（pheatmapを使用） ---------------------
cat("ヒートマップを作成中...\n")

# 色のパレット（青-白-赤、0で真っ白）
color_palette <- colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(100)

# 列のアノテーション（クラスター番号）
col_ann <- data.frame(
  Cluster = factor(colnames(heatmap_data_scaled), levels = colnames(heatmap_data_scaled))
)
rownames(col_ann) <- colnames(heatmap_data_scaled)

# クラスターごとの色を定義（デフォルトはパステルカラー）
n_clusters <- length(clusters)
cluster_colors <- rainbow(n_clusters, alpha = 0.7)
names(cluster_colors) <- colnames(heatmap_data_scaled)

ann_colors <- list(Cluster = cluster_colors)

# ヒートマップを作成
png(
  filename = file.path(output_dir, "細胞タイプマーカー遺伝子ヒートマップ_クラスター版.png"),
  width = 14,
  height = 16,
  units = "in",
  res = 300
)

# カラースケールの範囲を設定（0で真っ白になるように）
value_range <- range(heatmap_data_scaled, na.rm = TRUE)
max_abs_value <- max(abs(value_range))

# 対称的な範囲でカラースケールを設定
breaks <- seq(-max_abs_value, max_abs_value, length.out = 101)

pheatmap(
  heatmap_data_scaled,
  color = color_palette,
  breaks = breaks,  # 0で真っ白になるように対称的なbreaksを設定
  cluster_rows = FALSE,  # 遺伝子の順序は固定
  cluster_cols = FALSE,  # クラスターの順序は固定（数値順）
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize = 10,
  fontsize_row = 8,
  fontsize_col = 12,
  annotation_col = col_ann,
  annotation_colors = ann_colors,
  gaps_row = cumsum(sapply(celltype_order, function(x) {
    if (x %in% names(top_markers_list_mapped)) {
      sum(top_markers_list_mapped[[x]] %in% available_markers)
    } else {
      0
    }
  }))[-length(celltype_order)],  # 各細胞タイプのマーカー遺伝子の間に区切り線
  gaps_col = NULL,
  main = "Marker Gene Expression by Cluster (CTL)",
  border_color = NA,
  scale = "none"  # 既にスケーリング済み
)

dev.off()

cat("✓ ヒートマップを保存しました\n")
cat("  -", file.path(output_dir, "細胞タイプマーカー遺伝子ヒートマップ_クラスター版.png"), "\n\n")

# =============================================================
# 6. クラスターごとの細胞タイプ構成をCSVに保存 ------------------
if ("celltype_10class" %in% colnames(seu_ctl@meta.data)) {
  cat("クラスターごとの細胞タイプ構成をCSVに保存中...\n")
  
  cluster_celltype_summary <- data.frame()
  for (cluster_id in clusters) {
    cluster_cells <- as.character(seu_ctl$seurat_clusters) == cluster_id
    if (sum(cluster_cells) > 0) {
      celltype_counts <- table(seu_ctl$celltype_10class[cluster_cells])
      celltype_prop <- prop.table(celltype_counts) * 100
      
      for (celltype in names(celltype_counts)) {
        cluster_celltype_summary <- rbind(
          cluster_celltype_summary,
          data.frame(
            Cluster = cluster_id,
            CellType = celltype,
            Count = as.numeric(celltype_counts[celltype]),
            Proportion = as.numeric(celltype_prop[celltype]),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  write.csv(
    cluster_celltype_summary,
    file = file.path(output_dir, "クラスター別細胞タイプ構成.csv"),
    row.names = FALSE
  )
  
  cat("✓ クラスター別細胞タイプ構成を保存しました\n")
  cat("  -", file.path(output_dir, "クラスター別細胞タイプ構成.csv"), "\n\n")
}

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("解析完了！\n")
cat("※ ヒートマップは、横軸がCTLのクラスター、縦軸がマーカー遺伝子です。\n")
cat("   発現量はz-scoreスケーリング（中央揃えおよびスケール表示）されています。\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

