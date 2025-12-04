# =============================================================
# 0. ライブラリ -------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)
library(SingleR)
library(celldex)
library(SingleCellExperiment)

# =============================================================
# 1. データ読み込み ---------------------------------------------
mat_dir <- "/Users/godakyosuke/Desktop/goda sc RNA seq/yellow-aggr_outs/count/filtered_feature_bc_matrix"
seu <- Read10X(data.dir = mat_dir) |> CreateSeuratObject(project = "yellow_aggr")

# =============================================================
# 2. サンプル情報付与 -------------------------------------------
seu$lib_idx <- sub(".*-(\\d+)$", "\\1", colnames(seu))
lib_ids <- c(
  "7928-DTA_dta_1", "7935-Kir_control_1", "7936-Kir_control_2",
  "7937-Kir_control_3", "7938-Kir_kir_1", "7939-Kir_kir_2",
  "7940-Kir_kir_3", "8112-DD3", "8693-DC1", "8694-DC2",
  "8695-DC3", "8696-DD2"
)
seu$sample_id <- lib_ids[as.integer(seu$lib_idx)]

# --- グルーピング ---
seu$condition <- case_when(
  grepl("DTA|DD2|DD3", seu$sample_id) ~ "DTA",
  grepl("DC1|DC2|DC3", seu$sample_id) ~ "CTL",
  TRUE ~ "OTHER"
)

# --- Kir群除外 ---
seu <- subset(seu, subset = condition %in% c("DTA", "CTL"))
cat("✓ 対象サンプル数:", length(unique(seu$sample_id)), "（DTAとCTLのみ使用）\n")

# =============================================================
# 3. QCフィルタリング -------------------------------------------
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
seu <- subset(seu,
              subset = nFeature_RNA > 200 &
                nFeature_RNA < 5000 &
                nCount_RNA > 1000 &
                nCount_RNA < 30000 &
                percent.mt < 15)
cat("✓ QCフィルタリング後:", ncol(seu), "cells\n")

# =============================================================
# 4. 群ごとのUMAP作成 ------------------------------------------
seu_list <- SplitObject(seu, split.by = "condition")

# ---- DTA ----
seu_dta <- seu_list[["DTA"]] |> NormalizeData() |> FindVariableFeatures() |> ScaleData() |> RunPCA()
seu_dta <- FindNeighbors(seu_dta, dims = 1:10)
seu_dta <- FindClusters(seu_dta, resolution = 0.5)
seu_dta <- RunUMAP(seu_dta, dims = 1:10)

# ---- CTL ----
seu_ctl <- seu_list[["CTL"]] |> NormalizeData() |> FindVariableFeatures() |> ScaleData() |> RunPCA()
seu_ctl <- FindNeighbors(seu_ctl, dims = 1:10)
seu_ctl <- FindClusters(seu_ctl, resolution = 0.5)
seu_ctl <- RunUMAP(seu_ctl, dims = 1:10)

# =============================================================
# 5. SingleRによる自動アノテーション ----------------------------
cat("\n✓ SingleRによる細胞型アノテーション開始...\n")

# 参照データセット取得（マウス）
ref_data <- celldex::MouseRNAseqData()

# SingleCellExperimentオブジェクトに変換
sce_dta <- as.SingleCellExperiment(seu_dta)
sce_ctl <- as.SingleCellExperiment(seu_ctl)

# SingleRでアノテーション（各細胞レベル）
# label.main: 大まかな分類
pred_dta_main <- SingleR(test = sce_dta, ref = ref_data, labels = ref_data$label.main)
pred_ctl_main <- SingleR(test = sce_ctl, ref = ref_data, labels = ref_data$label.main)

# label.fine: 詳細な分類（Neuronをさらに細分化）
pred_dta_fine <- SingleR(test = sce_dta, ref = ref_data, labels = ref_data$label.fine)
pred_ctl_fine <- SingleR(test = sce_ctl, ref = ref_data, labels = ref_data$label.fine)

# Seuratオブジェクトに結果を追加
seu_dta$celltype_main <- pred_dta_main$labels
seu_dta$celltype_fine <- pred_dta_fine$labels
seu_dta$celltype_score <- pred_dta_main$scores |> apply(1, max)
seu_ctl$celltype_main <- pred_ctl_main$labels
seu_ctl$celltype_fine <- pred_ctl_fine$labels
seu_ctl$celltype_score <- pred_ctl_main$scores |> apply(1, max)

cat("✓ DTA: 推定された細胞型種類数（main） =", length(unique(pred_dta_main$labels)), "\n")
cat("✓ CTL: 推定された細胞型種類数（main） =", length(unique(pred_ctl_main$labels)), "\n")
cat("✓ DTA: 推定された細胞型種類数（fine） =", length(unique(pred_dta_fine$labels)), "\n")
cat("✓ CTL: 推定された細胞型種類数（fine） =", length(unique(pred_ctl_fine$labels)), "\n")

cat("\nDTA 細胞型の内訳（main）:\n")
print(table(seu_dta$celltype_main))
cat("\nCTL 細胞型の内訳（main）:\n")
print(table(seu_ctl$celltype_main))

cat("\n=== Neuronクラスターの詳細分類（fine）===\n")
cat("DTA Neurons（fine分類）:\n")
print(table(seu_dta$celltype_fine[seu_dta$celltype_main == "Neurons"]))
cat("\nCTL Neurons（fine分類）:\n")
print(table(seu_ctl$celltype_fine[seu_ctl$celltype_main == "Neurons"]))

# =============================================================
# 6. Purkinje候補細胞のタグ付け --------------------------------
cnt_dta <- GetAssayData(seu_dta, layer = "counts")
cnt_ctl <- GetAssayData(seu_ctl, layer = "counts")

seu_dta$Purkinje_tag <- ifelse(cnt_dta["Calb1", ] > 30 & cnt_dta["Car8", ] > 30, "Purkinje", "Other")
seu_ctl$Purkinje_tag <- ifelse(cnt_ctl["Calb1", ] > 30 & cnt_ctl["Car8", ] > 30, "Purkinje", "Other")

cat("DTA Purkinje:", sum(seu_dta$Purkinje_tag == "Purkinje"), "\n")
cat("CTL Purkinje:", sum(seu_ctl$Purkinje_tag == "Purkinje"), "\n")

# =============================================================
# 7. 論文ベースの小脳細胞10タイプ分類 -------------------------
cat("\n✓ 論文ベースのマーカー遺伝子による10タイプ分類開始...\n")
cat("参考: Nature Communications 2023 - Cerebellar Cell Atlas\n")

# マーカー遺伝子定義（論文準拠）
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
result_dta <- classify_cerebellar_cells(seu_dta)
result_ctl <- classify_cerebellar_cells(seu_ctl)

seu_dta$celltype_10class <- result_dta$celltype
seu_dta$celltype_10class_score <- apply(result_dta$scores, 1, max)
seu_ctl$celltype_10class <- result_ctl$celltype
seu_ctl$celltype_10class_score <- apply(result_ctl$scores, 1, max)

cat("\n=== 10種類の細胞タイプ分類結果 ===\n")
cat("\nDTA 細胞タイプ内訳:\n")
print(table(seu_dta$celltype_10class))
cat("\nCTL 細胞タイプ内訳:\n")
print(table(seu_ctl$celltype_10class))

# 各細胞タイプの割合
cat("\nDTA 細胞タイプ割合:\n")
print(round(prop.table(table(seu_dta$celltype_10class)) * 100, 2))
cat("\nCTL 細胞タイプ割合:\n")
print(round(prop.table(table(seu_ctl$celltype_10class)) * 100, 2))

# =============================================================
# 8. UMAP可視化: 多角的なアノテーション表示 --------------------

# ---- DTA ----
umap_dta_df <- Embeddings(seu_dta, "umap") |> as.data.frame()
colnames(umap_dta_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_dta_df$Purkinje_tag <- seu_dta$Purkinje_tag

p_dta <- DimPlot(seu_dta, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.4) +
  geom_point(
    data = subset(umap_dta_df, Purkinje_tag == "Purkinje"),
    aes(x = UMAP_1, y = UMAP_2),
    color = "red3", size = 1.2, inherit.aes = FALSE
  ) +
  ggtitle("DTA: Clusters + Purkinje Highlighted (Calb1>30 & Car8>30)") +
  theme_classic()

# ---- CTL ----
umap_ctl_df <- Embeddings(seu_ctl, "umap") |> as.data.frame()
colnames(umap_ctl_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_ctl_df$Purkinje_tag <- seu_ctl$Purkinje_tag

p_ctl <- DimPlot(seu_ctl, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.4) +
  geom_point(
    data = subset(umap_ctl_df, Purkinje_tag == "Purkinje"),
    aes(x = UMAP_1, y = UMAP_2),
    color = "red3", size = 1.2, inherit.aes = FALSE
  ) +
  ggtitle("CTL: Clusters + Purkinje Highlighted (Calb1>30 & Car8>30)") +
  theme_classic()

# ---- SingleRアノテーション結果のUMAP（main） ----
p_dta_main <- DimPlot(seu_dta, reduction = "umap", group.by = "celltype_main", 
                       pt.size = 0.4, label = TRUE, repel = TRUE) +
  ggtitle("DTA: SingleR (Main Category)") +
  theme_classic() +
  theme(legend.position = "right")

p_ctl_main <- DimPlot(seu_ctl, reduction = "umap", group.by = "celltype_main", 
                       pt.size = 0.4, label = TRUE, repel = TRUE) +
  ggtitle("CTL: SingleR (Main Category)") +
  theme_classic() +
  theme(legend.position = "right")

# ---- SingleRアノテーション結果のUMAP（fine） ----
p_dta_fine <- DimPlot(seu_dta, reduction = "umap", group.by = "celltype_fine", 
                       pt.size = 0.4, label = TRUE, repel = TRUE, label.size = 3) +
  ggtitle("DTA: SingleR (Fine Category)") +
  theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 8))

p_ctl_fine <- DimPlot(seu_ctl, reduction = "umap", group.by = "celltype_fine", 
                       pt.size = 0.4, label = TRUE, repel = TRUE, label.size = 3) +
  ggtitle("CTL: SingleR (Fine Category)") +
  theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 8))

# ---- 論文ベース10種類の細胞タイプUMAP ----
# カラーパレット設定（見やすい色）
cell_colors <- c(
  "Granule" = "#E41A1C",
  "Purkinje" = "#377EB8",
  "Bergmann_glia" = "#4DAF4A",
  "Interneurons" = "#984EA3",
  "Astrocytes" = "#FF7F00",
  "Oligodendrocytes" = "#FFFF33",
  "OPC" = "#A65628",
  "Microglia" = "#F781BF",
  "Endothelial" = "#999999",
  "Fibroblasts" = "#66C2A5",
  "Unclassified" = "#CCCCCC"
)

p_dta_10class <- DimPlot(seu_dta, reduction = "umap", group.by = "celltype_10class", 
                          pt.size = 0.5, label = TRUE, repel = TRUE, label.size = 3.5) +
  scale_color_manual(values = cell_colors) +
  ggtitle("DTA: 10 Major Cell Types (Paper-based)") +
  theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 9))

p_ctl_10class <- DimPlot(seu_ctl, reduction = "umap", group.by = "celltype_10class", 
                          pt.size = 0.5, label = TRUE, repel = TRUE, label.size = 3.5) +
  scale_color_manual(values = cell_colors) +
  ggtitle("CTL: 10 Major Cell Types (Paper-based)") +
  theme_classic() +
  theme(legend.position = "right", legend.text = element_text(size = 9))

# ---- 並列表示 ----
cat("\n=== 1. クラスタとPurkinje細胞のハイライト ===\n")
print(p_dta + p_ctl)

cat("\n=== 2. SingleRによる細胞型アノテーション（Main） ===\n")
print(p_dta_main + p_ctl_main)

cat("\n=== 3. SingleRによる細胞型アノテーション（Fine - Neuron詳細） ===\n")
print(p_dta_fine + p_ctl_fine)

cat("\n=== 4. 論文ベース10種類の細胞タイプ分類 ===\n")
print(p_dta_10class + p_ctl_10class)

# ---- アノテーションスコアの分布確認 ----
cat("\n=== SingleR アノテーションスコアの分布 ===\n")
p_score_dta <- ggplot(seu_dta@meta.data, aes(x = celltype_score)) +
  geom_histogram(bins = 50, fill = "skyblue", color = "black") +
  ggtitle("DTA: SingleR Score Distribution") +
  xlab("Annotation Score") +
  theme_classic()

p_score_ctl <- ggplot(seu_ctl@meta.data, aes(x = celltype_score)) +
  geom_histogram(bins = 50, fill = "salmon", color = "black") +
  ggtitle("CTL: SingleR Score Distribution") +
  xlab("Annotation Score") +
  theme_classic()

print(p_score_dta + p_score_ctl)

# =============================================================
# 9. マーカー遺伝子の発現確認（FeaturePlot） ------------------
cat("\n✓ 主要マーカー遺伝子の発現パターン確認...\n")

# 各細胞タイプの代表マーカー
key_markers <- c("Calb1", "Itpr1",    # Purkinje
                 "Rims1", "Gabra6",   # Granule
                 "Aqp4", "Gfap",      # Bergmann/Astrocytes
                 "Pvalb", "Gad1",     # Interneurons
                 "Mbp", "Plp1",       # Oligodendrocytes
                 "Pdgfra",            # OPC
                 "Csf1r")             # Microglia

# 利用可能なマーカーのみを選択
available_markers <- key_markers[key_markers %in% rownames(seu_dta)]

if (length(available_markers) > 0) {
  cat("発現確認マーカー:", paste(available_markers, collapse = ", "), "\n")
  
  # 主要な4つのマーカーをプロット
  top_markers <- head(available_markers, 4)
  
  p_markers_dta <- FeaturePlot(seu_dta, features = top_markers, 
                                ncol = 2, pt.size = 0.3) &
    theme_classic() &
    theme(legend.position = "right")
  
  p_markers_ctl <- FeaturePlot(seu_ctl, features = top_markers, 
                                ncol = 2, pt.size = 0.3) &
    theme_classic() &
    theme(legend.position = "right")
  
  cat("\n=== 5. 主要マーカー遺伝子の発現パターン ===\n")
  cat("DTA:\n")
  print(p_markers_dta)
  cat("\nCTL:\n")
  print(p_markers_ctl)
}

# =============================================================
# 10. ニューロンタイプのサブクラスタリング（オプション） ------
cat("\n✓ ニューロン系細胞（Granule, Purkinje, Interneurons）の詳細解析...\n")

# ニューロン系のみ抽出
neuron_types <- c("Granule", "Purkinje", "Interneurons")
seu_dta_neurons <- subset(seu_dta, subset = celltype_10class %in% neuron_types)
seu_ctl_neurons <- subset(seu_ctl, subset = celltype_10class %in% neuron_types)

if (ncol(seu_dta_neurons) > 50) {  # 十分な細胞数がある場合のみ
  # 再正規化・再クラスタリング
  seu_dta_neurons <- seu_dta_neurons |> 
    NormalizeData() |> 
    FindVariableFeatures() |> 
    ScaleData() |> 
    RunPCA(verbose = FALSE)
  
  seu_dta_neurons <- FindNeighbors(seu_dta_neurons, dims = 1:15, verbose = FALSE)
  seu_dta_neurons <- FindClusters(seu_dta_neurons, resolution = 0.3, verbose = FALSE)
  seu_dta_neurons <- RunUMAP(seu_dta_neurons, dims = 1:15, verbose = FALSE)
  
  cat("✓ DTA Neurons: サブクラスター数 =", length(unique(seu_dta_neurons$seurat_clusters)), "\n")
  cat("   - Granule:", sum(seu_dta_neurons$celltype_10class == "Granule"), "cells\n")
  cat("   - Purkinje:", sum(seu_dta_neurons$celltype_10class == "Purkinje"), "cells\n")
  cat("   - Interneurons:", sum(seu_dta_neurons$celltype_10class == "Interneurons"), "cells\n")
  
  # サブクラスターのUMAP
  p_dta_neuron_subcluster <- DimPlot(seu_dta_neurons, reduction = "umap", 
                                      group.by = "seurat_clusters", 
                                      pt.size = 1, label = TRUE) +
    ggtitle("DTA: Neuron Sub-clustering") +
    theme_classic()
  
  # ニューロンタイプ別の色分け
  p_dta_neuron_type <- DimPlot(seu_dta_neurons, reduction = "umap", 
                                group.by = "celltype_10class", 
                                pt.size = 1) +
    scale_color_manual(values = cell_colors) +
    ggtitle("DTA: Neuron Types (Granule/Purkinje/Interneurons)") +
    theme_classic()
  
  cat("\n=== 6. Neuronのみの再クラスタリング（DTA） ===\n")
  print(p_dta_neuron_subcluster + p_dta_neuron_type)
}

if (ncol(seu_ctl_neurons) > 50) {
  seu_ctl_neurons <- seu_ctl_neurons |> 
    NormalizeData() |> 
    FindVariableFeatures() |> 
    ScaleData() |> 
    RunPCA(verbose = FALSE)
  
  seu_ctl_neurons <- FindNeighbors(seu_ctl_neurons, dims = 1:15, verbose = FALSE)
  seu_ctl_neurons <- FindClusters(seu_ctl_neurons, resolution = 0.3, verbose = FALSE)
  seu_ctl_neurons <- RunUMAP(seu_ctl_neurons, dims = 1:15, verbose = FALSE)
  
  cat("✓ CTL Neurons: サブクラスター数 =", length(unique(seu_ctl_neurons$seurat_clusters)), "\n")
  cat("   - Granule:", sum(seu_ctl_neurons$celltype_10class == "Granule"), "cells\n")
  cat("   - Purkinje:", sum(seu_ctl_neurons$celltype_10class == "Purkinje"), "cells\n")
  cat("   - Interneurons:", sum(seu_ctl_neurons$celltype_10class == "Interneurons"), "cells\n")
  
  p_ctl_neuron_subcluster <- DimPlot(seu_ctl_neurons, reduction = "umap", 
                                      group.by = "seurat_clusters", 
                                      pt.size = 1, label = TRUE) +
    ggtitle("CTL: Neuron Sub-clustering") +
    theme_classic()
  
  p_ctl_neuron_type <- DimPlot(seu_ctl_neurons, reduction = "umap", 
                                group.by = "celltype_10class", 
                                pt.size = 1) +
    scale_color_manual(values = cell_colors) +
    ggtitle("CTL: Neuron Types (Granule/Purkinje/Interneurons)") +
    theme_classic()
  
  cat("\n=== 7. Neuronのみの再クラスタリング（CTL） ===\n")
  print(p_ctl_neuron_subcluster + p_ctl_neuron_type)
}

# =============================================================
# 11. 細胞タイプ別統計サマリー -----------------------------------
cat("\n" , paste(rep("=", 60), collapse = ""), "\n")
cat("✓ 全ての解析が完了しました！\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("=== 最終結果サマリー ===\n\n")

cat("【10種類の細胞タイプ分類】\n")
cat("論文ベースのマーカー遺伝子による分類:\n")
cat("  1. Granule cells (RIMS1, GRM4)\n")
cat("  2. Purkinje cells (ITPR1, CALB1)\n")
cat("  3. Bergmann glia (TUBB2B, AQP4)\n")
cat("  4. Interneurons (GAD1, PVALB)\n")
cat("  5. Astrocytes (TTN, AQP1)\n")
cat("  6. Oligodendrocytes (PLP1, MBP)\n")
cat("  7. OPC (PDGFRA, OLIG1)\n")
cat("  8. Microglia (CD74, CSF1R)\n")
cat("  9. Endothelial (CLDN5, VWF)\n")
cat(" 10. Fibroblasts (DCN, APOD)\n\n")

cat("【DTA群】\n")
cat("総細胞数:", ncol(seu_dta), "\n")
dta_tab <- table(seu_dta$celltype_10class)
dta_prop <- prop.table(dta_tab) * 100
for (i in 1:length(dta_tab)) {
  cat(sprintf("  %-18s: %5d cells (%5.2f%%)\n", 
              names(dta_tab)[i], dta_tab[i], dta_prop[i]))
}

cat("\n【CTL群】\n")
cat("総細胞数:", ncol(seu_ctl), "\n")
ctl_tab <- table(seu_ctl$celltype_10class)
ctl_prop <- prop.table(ctl_tab) * 100
for (i in 1:length(ctl_tab)) {
  cat(sprintf("  %-18s: %5d cells (%5.2f%%)\n", 
              names(ctl_tab)[i], ctl_tab[i], ctl_prop[i]))
}

cat("\n" , paste(rep("=", 60), collapse = ""), "\n")
cat("保存されているメタデータ列:\n")
cat("  - celltype_main: SingleR大分類\n")
cat("  - celltype_fine: SingleR詳細分類\n")
cat("  - celltype_10class: 論文ベース10種類分類 ★推奨★\n")
cat("  - celltype_10class_score: 分類スコア\n")
cat("  - Purkinje_tag: Purkinje細胞タグ\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("主要オブジェクト:\n")
cat("  - seu_dta: DTAグループのSeuratオブジェクト\n")
cat("  - seu_ctl: CTLグループのSeuratオブジェクト\n")
if (exists("seu_dta_neurons")) {
  cat("  - seu_dta_neurons: DTAニューロン細胞のみ（再クラスタリング済み）\n")
}
if (exists("seu_ctl_neurons")) {
  cat("  - seu_ctl_neurons: CTLニューロン細胞のみ（再クラスタリング済み）\n")
}

cat("\n次のステップ:\n")
cat("  - DimPlot(seu_dta, group.by = 'celltype_10class') で確認\n")
cat("  - FeaturePlot(seu_dta, features = 'Calb1') でマーカー確認\n")
cat("  - saveRDS(seu_dta, 'seu_dta_annotated.rds') で保存\n")

# =============================================================
# 12. 図の保存 --------------------------------------------------
cat("\n✓ 図をPNGファイルとして保存中...\n")

# 保存先ディレクトリの作成
output_dir <- "/Users/godakyosuke/Desktop/goda sc RNA seq/UMAP_outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 1. クラスター + Purkinjeハイライト
ggsave(file.path(output_dir, "01_Clusters_Purkinje.png"), 
       p_dta + p_ctl, width = 14, height = 6, dpi = 300)

# 2. SingleR Main分類
ggsave(file.path(output_dir, "02_SingleR_Main.png"), 
       p_dta_main + p_ctl_main, width = 14, height = 6, dpi = 300)

# 3. SingleR Fine分類
ggsave(file.path(output_dir, "03_SingleR_Fine.png"), 
       p_dta_fine + p_ctl_fine, width = 16, height = 6, dpi = 300)

# 4. 10種類細胞タイプ（論文ベース）★メイン★
ggsave(file.path(output_dir, "04_10CellTypes_PaperBased.png"), 
       p_dta_10class + p_ctl_10class, width = 16, height = 6, dpi = 300)

# 5. アノテーションスコア分布
ggsave(file.path(output_dir, "05_Annotation_Scores.png"), 
       p_score_dta + p_score_ctl, width = 12, height = 5, dpi = 300)

# 6. マーカー遺伝子発現パターン（もし作成されていれば）
if (exists("p_markers_dta")) {
  ggsave(file.path(output_dir, "06_Markers_DTA.png"), 
         p_markers_dta, width = 10, height = 8, dpi = 300)
}

if (exists("p_markers_ctl")) {
  ggsave(file.path(output_dir, "07_Markers_CTL.png"), 
         p_markers_ctl, width = 10, height = 8, dpi = 300)
}

# 7. Neuronサブクラスタリング（もし作成されていれば）
if (exists("p_dta_neuron_subcluster")) {
  ggsave(file.path(output_dir, "08_Neurons_Subcluster_DTA.png"), 
         p_dta_neuron_subcluster + p_dta_neuron_type, 
         width = 14, height = 6, dpi = 300)
}

if (exists("p_ctl_neuron_subcluster")) {
  ggsave(file.path(output_dir, "09_Neurons_Subcluster_CTL.png"), 
         p_ctl_neuron_subcluster + p_ctl_neuron_type, 
         width = 14, height = 6, dpi = 300)
}

# PDFでも保存（高品質・全図まとめ）
pdf(file.path(output_dir, "All_Figures.pdf"), width = 16, height = 6)
print(p_dta + p_ctl)
print(p_dta_main + p_ctl_main)
print(p_dta_fine + p_ctl_fine)
print(p_dta_10class + p_ctl_10class)
print(p_score_dta + p_score_ctl)
if (exists("p_markers_dta")) print(p_markers_dta)
if (exists("p_markers_ctl")) print(p_markers_ctl)
if (exists("p_dta_neuron_subcluster")) {
  print(p_dta_neuron_subcluster + p_dta_neuron_type)
}
if (exists("p_ctl_neuron_subcluster")) {
  print(p_ctl_neuron_subcluster + p_ctl_neuron_type)
}
dev.off()

cat("✓ 図の保存完了:", output_dir, "\n")
cat("  保存ファイル:\n")
cat("    - 01_Clusters_Purkinje.png\n")
cat("    - 02_SingleR_Main.png\n")
cat("    - 03_SingleR_Fine.png\n")
cat("    - 04_10CellTypes_PaperBased.png ★メイン★\n")
cat("    - 05_Annotation_Scores.png\n")
cat("    - All_Figures.pdf（全図まとめ）\n")

cat("\n解析完了！\n")
