# =============================================================
# 細胞タイプマーカー遺伝子ヒートマップ作成スクリプト
# 各細胞タイプごとの上位5つのマーカー遺伝子をヒートマップで表示
# 論文のPanel Bに類似した形式
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

# celltype_10classが存在しない場合は自動分類
if (!"celltype_10class" %in% colnames(seu_ctl@meta.data)) {
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
  
  if (!"celltype_10class" %in% colnames(seu_ctl@meta.data)) {
    cat("  CTLデータを分類中...\n")
    result_ctl <- classify_cerebellar_cells(seu_ctl)
    seu_ctl$celltype_10class <- result_ctl$celltype
    cat("  ✓ CTL分類完了\n")
  }
  
  cat("✓ 自動分類完了\n\n")
}

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
# 大文字小文字を確認して実際の遺伝子名にマッピング
top_markers_list <- list(
  "Granule" = c("Etv1", "Rims1", "Scn2a", "Grm4"),  # Camk4を削除
  "Oligodendrocytes" = c("Plp1", "Mbp", "Mobp", "Cnp"),  # Trfを削除
  "Bergmann_glia" = c("Slc1a3", "Ndrg2", "Aqp4", "Tubb2b"),  # Apoeを削除
  "Astrocytes" = c("Ttn", "Gfap", "Aldh1l1", "S100b"),  # Aqp1→ALDH1L1, Gja1→S100β(S100b)
  "Interneurons" = c("Kit", "Gad1", "Gad2", "Pvalb"),  # Homer3を削除
  "Endothelial" = c("Cldn5", "Pecam1", "Vwf", "Flt1"),  # B2m→CD31(Pecam1)
  "Purkinje" = c("Itpr1", "Calb1", "Pcp4", "Car8"),  # Aldocを削除
  "OPC" = c("Ptprz1", "Vcan", "Tnr", "Pdgfra"),  # Olig1を削除して4個に
  "Microglia" = c("Cd74", "C1qb", "C3", "Csf1r"),  # Aif1を削除
  "Fibroblasts" = c("Dcn", "Apod", "Tagln", "Fbln1")  # C7を削除
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
# 4. ヒートマップ用のデータを準備 -----------------------------
cat("ヒートマップ用のデータを準備中...\n")

# Identsを設定（重要！）
Idents(seu_ctl) <- "celltype_10class"

# 正規化データを取得
data_norm <- GetAssayData(seu_ctl, layer = "data")

# 選択されたマーカー遺伝子のみを含むデータを抽出
available_markers <- top_markers[top_markers %in% rownames(data_norm)]
if (length(available_markers) < length(top_markers)) {
  cat("警告: 一部のマーカー遺伝子がデータに存在しません\n")
  cat("  選択された:", length(top_markers), "個\n")
  cat("  利用可能:", length(available_markers), "個\n")
}

# 各細胞タイプごとに平均発現量を計算
heatmap_data <- matrix(0, nrow = length(available_markers), ncol = length(celltype_order))
rownames(heatmap_data) <- available_markers
colnames(heatmap_data) <- celltype_order

for (celltype in celltype_order) {
  if (celltype %in% seu_ctl$celltype_10class) {
    cells_in_cluster <- WhichCells(seu_ctl, idents = celltype)
    if (length(cells_in_cluster) > 0) {
      # その細胞タイプの平均発現量を計算
      mean_expr <- rowMeans(data_norm[available_markers, cells_in_cluster, drop = FALSE])
      heatmap_data[, celltype] <- mean_expr
    }
  }
}

# z-scoreスケーリング（中央揃えおよびスケール表示）
# 各行（遺伝子）について、平均を引いて標準偏差で割る
heatmap_data_scaled <- t(scale(t(heatmap_data)))

# マーカー遺伝子を細胞タイプごとにグループ化して並び替え
# 各細胞タイプのマーカー遺伝子をグループ化
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
cat("  データサイズ:", nrow(heatmap_data_scaled), "genes x", ncol(heatmap_data_scaled), "cell types\n\n")

# =============================================================
# 5. ヒートマップを作成（pheatmapを使用） ---------------------
cat("ヒートマップを作成中...\n")

# 色のパレット（青-白-赤、0で真っ白）
# カラースケールは0で真っ白になるように設定
# 対称的なカラースケール（例：-2から2までの範囲で、0が白）
color_palette <- colorRampPalette(c("#2166AC", "#FFFFFF", "#B2182B"))(100)

# 列のアノテーション（細胞タイプの色分け）
col_ann <- data.frame(
  CellType = factor(colnames(heatmap_data_scaled), levels = celltype_order)
)
rownames(col_ann) <- colnames(heatmap_data_scaled)

# 細胞タイプごとの色を定義
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

ann_colors <- list(CellType = celltype_colors[celltype_order])

# ヒートマップを作成
png(
  filename = file.path(output_dir, "細胞タイプマーカー遺伝子ヒートマップ.png"),
  width = 12,
  height = 16,
  units = "in",
  res = 300
)

# カラースケールの範囲を設定（0で真っ白になるように）
# z-scoreスケーリング後の値の範囲を取得
value_range <- range(heatmap_data_scaled, na.rm = TRUE)
max_abs_value <- max(abs(value_range))

# 対称的な範囲でカラースケールを設定（例：-2.5から2.5、0が白）
breaks <- seq(-max_abs_value, max_abs_value, length.out = 101)

pheatmap(
  heatmap_data_scaled,
  color = color_palette,
  breaks = breaks,  # 0で真っ白になるように対称的なbreaksを設定
  cluster_rows = FALSE,  # 遺伝子の順序は固定
  cluster_cols = FALSE,  # 細胞タイプの順序は固定
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
  main = "B Expression",
  border_color = NA,
  scale = "none"  # 既にスケーリング済み
)

dev.off()

cat("✓ ヒートマップを保存しました\n")
cat("  -", file.path(output_dir, "細胞タイプマーカー遺伝子ヒートマップ.png"), "\n\n")

# =============================================================
# 6. マーカー遺伝子のリストをCSVに保存 -------------------------
cat("マーカー遺伝子のリストをCSVに保存中...\n")

# 指定されたマーカー遺伝子のリストを保存
markers_summary <- data.frame()
for (celltype in celltype_order) {
  if (celltype %in% names(top_markers_list_mapped)) {
    markers <- top_markers_list_mapped[[celltype]]
    if (length(markers) > 0) {
      markers_df <- data.frame(
        gene = markers,
        cluster = celltype,
        order = 1:length(markers),
        stringsAsFactors = FALSE
      )
      markers_summary <- rbind(markers_summary, markers_df)
    }
  }
}

write.csv(
  markers_summary,
  file = file.path(output_dir, "細胞タイプマーカー遺伝子リスト.csv"),
  row.names = FALSE
)

cat("✓ マーカー遺伝子リストを保存しました\n")
cat("  -", file.path(output_dir, "細胞タイプマーカー遺伝子リスト.csv"), "\n\n")

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("解析完了！\n")
cat("※ ヒートマップは、各細胞タイプごとの上位5つのマーカー遺伝子を表示しています。\n")
cat("   発現量はz-scoreスケーリング（中央揃えおよびスケール表示）されています。\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

