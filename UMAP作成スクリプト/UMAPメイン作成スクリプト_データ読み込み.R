# =============================================================
# UMAPメイン作成スクリプト
# メインクラスタリングとUMAP図作成まで
# =============================================================

# 0. ライブラリ -------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)
library(ggforce)  # グループを囲む破線用

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
cat("\n✓ DTA群の処理開始...\n")
seu_dta <- seu_list[["DTA"]] |> 
  NormalizeData() |> 
  FindVariableFeatures() |> 
  ScaleData() |> 
  RunPCA()
seu_dta <- FindNeighbors(seu_dta, dims = 1:10)
seu_dta <- FindClusters(seu_dta, resolution = 0.5)
seu_dta <- RunUMAP(seu_dta, dims = 1:10)
cat("✓ DTA: クラスター数 =", length(unique(seu_dta$seurat_clusters)), "\n")

# ---- CTL ----
cat("\n✓ CTL群の処理開始...\n")
seu_ctl <- seu_list[["CTL"]] |> 
  NormalizeData() |> 
  FindVariableFeatures() |> 
  ScaleData() |> 
  RunPCA()
seu_ctl <- FindNeighbors(seu_ctl, dims = 1:10)
seu_ctl <- FindClusters(seu_ctl, resolution = 0.5)
seu_ctl <- RunUMAP(seu_ctl, dims = 1:10)
cat("✓ CTL: クラスター数 =", length(unique(seu_ctl$seurat_clusters)), "\n")

# =============================================================
# 4.5. 処理済みデータの保存 ------------------------------------
cat("\n✓ 処理済みデータを保存中...\n")
saveRDS(seu_dta, file = "/Users/godakyosuke/Desktop/goda sc RNA seq/seu_dta_processed.rds")
saveRDS(seu_ctl, file = "/Users/godakyosuke/Desktop/goda sc RNA seq/seu_ctl_processed.rds")
cat("✓ 保存完了: seu_dta_processed.rds, seu_ctl_processed.rds\n")
cat("  次回からは「UMAPメイン作成スクリプト_図作成.R」を使用すると高速です\n")

# =============================================================
# 5. UMAP可視化: クラスタリング結果 ----------------------------

# グループ化情報（オプション：必要に応じて編集）
# 例：cluster_groups_dta <- list(
#   "Group1" = c("0", "1", "2"),
#   "Group2" = c("3", "4", "5")
# )
# デフォルトではグループ化なし
cluster_groups_dta <- NULL
cluster_groups_ctl <- NULL

# カスタムテーマ設定（画像のスタイルに合わせる）
# 枠線の太さを統一（1.5に統一）
border_width <- 1.5

custom_theme <- theme_classic() +
  theme(
    # 軸ラベルとタイトルのフォントサイズを大きく
    axis.title = element_text(size = 18, face = "bold"),
    # 縦軸横軸の数値のフォント太さを統一（太字に）
    axis.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    # クラスターラベルのフォントサイズを大きく、太字
    legend.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 18, face = "bold"),
    # 軸の線を太く（枠線の太さと統一）
    axis.line = element_line(linewidth = border_width),
    # 背景を白に
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # 図を四角く囲む（枠線の太さと統一）
    panel.border = element_rect(color = "black", fill = NA, linewidth = border_width),
    plot.margin = margin(20, 20, 20, 20, "pt")
  )

# ---- DTA: クラスタリング結果のUMAP ----
# UMAP座標を取得
umap_dta <- Embeddings(seu_dta, "umap") |> as.data.frame()
colnames(umap_dta) <- c("UMAP_1", "UMAP_2")
umap_dta$cluster <- seu_dta$seurat_clusters

# スケール統一のため、範囲を計算（後でCTLと統合）
dta_x_range <- range(umap_dta$UMAP_1)
dta_y_range <- range(umap_dta$UMAP_2)

# クラスターの中心座標を計算（ラベル位置用）
cluster_centers_dta <- umap_dta |>
  group_by(cluster) |>
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

p_dta_clusters <- ggplot(umap_dta, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.5, alpha = 0.7)

# グループを囲む太い赤い破線（オプション）
if (!is.null(cluster_groups_dta)) {
  for (group_name in names(cluster_groups_dta)) {
    group_clusters <- cluster_groups_dta[[group_name]]
    group_data <- umap_dta |> filter(cluster %in% group_clusters)
    if (nrow(group_data) > 2) {
      p_dta_clusters <- p_dta_clusters +
        geom_mark_hull(
          data = group_data,
          aes(x = UMAP_1, y = UMAP_2),
          color = "red",
          fill = NA,
          linetype = "dashed",
          linewidth = 2,
          expand = unit(2, "mm"),
          label = group_name,
          label.fontsize = 14,
          label.fontface = "bold",
          inherit.aes = FALSE
        )
    }
  }
}

# ---- CTL: クラスタリング結果のUMAP ----
umap_ctl <- Embeddings(seu_ctl, "umap") |> as.data.frame()
colnames(umap_ctl) <- c("UMAP_1", "UMAP_2")
umap_ctl$cluster <- seu_ctl$seurat_clusters

# スケール統一のため、範囲を計算
ctl_x_range <- range(umap_ctl$UMAP_1)
ctl_y_range <- range(umap_ctl$UMAP_2)

# DTAとCTLの共通範囲を計算（スケール統一）
common_x_range <- range(c(dta_x_range, ctl_x_range))
common_y_range <- range(c(dta_y_range, ctl_y_range))

# 少し余白を追加
x_margin <- diff(common_x_range) * 0.05
y_margin <- diff(common_y_range) * 0.05
common_x_range <- common_x_range + c(-x_margin, x_margin)
common_y_range <- common_y_range + c(-y_margin, y_margin)

p_dta_clusters <- p_dta_clusters +
  scale_color_discrete(name = "Cluster") +
  ggtitle("DTA: UMAP Clusters") +
  custom_theme +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  # スケールを統一
  coord_cartesian(xlim = common_x_range, ylim = common_y_range)

cluster_centers_ctl <- umap_ctl |>
  group_by(cluster) |>
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  )

p_ctl_clusters <- ggplot(umap_ctl, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.5, alpha = 0.7)

# グループを囲む太い赤い破線（オプション）
if (!is.null(cluster_groups_ctl)) {
  for (group_name in names(cluster_groups_ctl)) {
    group_clusters <- cluster_groups_ctl[[group_name]]
    group_data <- umap_ctl |> filter(cluster %in% group_clusters)
    if (nrow(group_data) > 2) {
      p_ctl_clusters <- p_ctl_clusters +
        geom_mark_hull(
          data = group_data,
          aes(x = UMAP_1, y = UMAP_2),
          color = "red",
          fill = NA,
          linetype = "dashed",
          linewidth = 2,
          expand = unit(2, "mm"),
          label = group_name,
          label.fontsize = 14,
          label.fontface = "bold",
          inherit.aes = FALSE
        )
    }
  }
}

p_ctl_clusters <- p_ctl_clusters +
  scale_color_discrete(name = "Cluster") +
  ggtitle("CTL: UMAP Clusters") +
  custom_theme +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  # スケールを統一（DTAと同じ範囲）
  coord_cartesian(xlim = common_x_range, ylim = common_y_range)

# ---- 並列表示 ----
cat("\n=== UMAP クラスタリング結果 ===\n")
print(p_dta_clusters + p_ctl_clusters)

# =============================================================
# 6. 図の保存 --------------------------------------------------
cat("\n✓ 図をPNGファイルとして保存中...\n")

# 保存先ディレクトリの作成
output_dir <- "/Users/godakyosuke/Desktop/goda sc RNA seq/UMAP_outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# クラスタリング結果のUMAP図を保存
ggsave(file.path(output_dir, "UMAP_Clusters_Main.png"), 
       p_dta_clusters + p_ctl_clusters, 
       width = 16, 
       height = 7, 
       dpi = 300,
       bg = "white")

cat("✓ 図の保存完了:", output_dir, "\n")
cat("  保存ファイル: UMAP_Clusters_Main.png\n")

cat("\n解析完了！\n")

