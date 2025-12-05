# =============================================================
# UMAPメイン作成スクリプト_図作成版
# 処理済みデータから図を作成（高速）
# =============================================================

# 0. ライブラリ -------------------------------------------------
library(Seurat)
library(dplyr)
library(ggplot2)
library(tibble)
library(patchwork)
library(ggforce)  # グループを囲む破線用

# =============================================================
# 1. 処理済みデータの読み込み ------------------------------------
cat("✓ 処理済みデータを読み込み中...\n")
seu_dta <- readRDS("/Users/godakyosuke/Desktop/goda sc RNA seq/seu_dta_processed.rds")
seu_ctl <- readRDS("/Users/godakyosuke/Desktop/goda sc RNA seq/seu_ctl_processed.rds")
cat("✓ DTA: クラスター数 =", length(unique(seu_dta$seurat_clusters)), "\n")
cat("✓ CTL: クラスター数 =", length(unique(seu_ctl$seurat_clusters)), "\n")

# =============================================================
# 2. UMAP可視化: クラスタリング結果 ----------------------------

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
    # 凡例の背景と枠線を削除
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    # 軸の線を無効化（panel.borderで統一するため）
    axis.line = element_blank(),
    # 背景を白に
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    # 図を四角く囲む（全辺で枠線の太さを統一）- これが唯一の枠線
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
  ) |>
  # クラスター番号を1から始まるように変換（0→1, 1→2, ...）
  # factor型の場合も正しく変換するため、文字列経由で変換
  mutate(cluster_label = as.integer(as.character(cluster)) + 1)

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

# CTLのデータ範囲を計算
ctl_x_range <- range(umap_ctl$UMAP_1)
ctl_y_range <- range(umap_ctl$UMAP_2)

# CTLのデータの中心を計算
ctl_x_center <- mean(ctl_x_range)
ctl_y_center <- mean(ctl_y_range)

# CTLのデータ範囲の大きさを計算
ctl_x_size <- diff(ctl_x_range)
ctl_y_size <- diff(ctl_y_range)
ctl_max_size <- max(ctl_x_size, ctl_y_size)

# CTL用の軸範囲と目盛りを設定（データを中央に配置）
# データの中心を0の周りに配置するように範囲を調整
ctl_x_range_set <- ctl_x_center + c(-ctl_max_size/2 - 2, ctl_max_size/2 + 2)  # 少し余白を追加
ctl_y_range_set <- ctl_y_center + c(-ctl_max_size/2 - 2, ctl_max_size/2 + 2)

# キリのいい数字に調整（範囲を丸める）
ctl_x_range_set <- c(floor(min(ctl_x_range_set)/5)*5, ceiling(max(ctl_x_range_set)/5)*5)
ctl_y_range_set <- c(floor(min(ctl_y_range_set)/5)*5, ceiling(max(ctl_y_range_set)/5)*5)

# 目盛りを設定
ctl_x_breaks <- seq(min(ctl_x_range_set), max(ctl_x_range_set), by = 5)
ctl_y_breaks <- seq(min(ctl_y_range_set), max(ctl_y_range_set), by = 5)

# DTA用の軸範囲と目盛りを設定（個別設定可能）
dta_x_range_set <- c(-10, 20)  # DTAのX軸範囲
dta_y_range_set <- c(-15, 15)  # DTAのY軸範囲
dta_x_breaks <- seq(-10, 20, by = 5)  # DTAのX軸目盛り
dta_y_breaks <- seq(-15, 15, by = 5)  # DTAのY軸目盛り

p_dta_clusters <- p_dta_clusters +
  # クラスター番号をクラスターの上に表示
  geom_text(
    data = cluster_centers_dta,
    aes(x = UMAP_1, y = UMAP_2, label = cluster_label),
    color = "black",
    size = 5,  # タイトルより小さめ（タイトルは20、これは約14pt相当）
    fontface = "bold",  # タイトルと同じ太字
    inherit.aes = FALSE
  ) +
  scale_color_discrete(
    name = "Cluster",
    labels = function(x) as.character(as.integer(as.character(x)) + 1)  # 凡例も1から始まるように
  ) +
  # DTA用の縦横軸の目盛りを設定
  scale_x_continuous(breaks = dta_x_breaks) +
  scale_y_continuous(breaks = dta_y_breaks) +
  ggtitle("DTA: UMAP Clusters") +
  custom_theme +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  # DTA用の範囲で正方形に固定
  coord_fixed(ratio = 1, xlim = dta_x_range_set, ylim = dta_y_range_set)

cluster_centers_ctl <- umap_ctl |>
  group_by(cluster) |>
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2),
    .groups = "drop"
  ) |>
  # クラスター番号を1から始まるように変換（0→1, 1→2, ...）
  # factor型の場合も正しく変換するため、文字列経由で変換
  mutate(cluster_label = as.integer(as.character(cluster)) + 1)

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
  # クラスター番号をクラスターの上に表示
  geom_text(
    data = cluster_centers_ctl,
    aes(x = UMAP_1, y = UMAP_2, label = cluster_label),
    color = "black",
    size = 5,  # タイトルより小さめ（タイトルは20、これは約14pt相当）
    fontface = "bold",  # タイトルと同じ太字
    inherit.aes = FALSE
  ) +
  scale_color_discrete(
    name = "Cluster",
    labels = function(x) as.character(as.integer(as.character(x)) + 1)  # 凡例も1から始まるように
  ) +
  # CTL用の縦横軸の目盛りを設定
  scale_x_continuous(breaks = ctl_x_breaks) +
  scale_y_continuous(breaks = ctl_y_breaks) +
  ggtitle("CTL: UMAP Clusters") +
  custom_theme +
  guides(color = guide_legend(override.aes = list(size = 5, alpha = 1))) +
  # CTL用の範囲で正方形に固定
  coord_fixed(ratio = 1, xlim = ctl_x_range_set, ylim = ctl_y_range_set)

# ---- 個別表示 ----
cat("\n=== UMAP クラスタリング結果 ===\n")
cat("DTA:\n")
print(p_dta_clusters)
cat("\nCTL:\n")
print(p_ctl_clusters)

# =============================================================
# 3. 図の保存（1枚ずつ） --------------------------------------
cat("\n✓ 図をPNGファイルとして保存中...\n")

# 保存先ディレクトリの作成
output_dir <- "/Users/godakyosuke/Desktop/goda sc RNA seq/UMAP_outputs"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# DTAのUMAP図を個別に保存
ggsave(file.path(output_dir, "UMAP_Clusters_DTA.png"), 
       p_dta_clusters, 
       width = 8, 
       height = 8, 
       dpi = 300,
       bg = "white")

# CTLのUMAP図を個別に保存
ggsave(file.path(output_dir, "UMAP_Clusters_CTL.png"), 
       p_ctl_clusters, 
       width = 8, 
       height = 8, 
       dpi = 300,
       bg = "white")

cat("✓ 図の保存完了:", output_dir, "\n")
cat("  保存ファイル:\n")
cat("    - UMAP_Clusters_DTA.png\n")
cat("    - UMAP_Clusters_CTL.png\n")

cat("\n解析完了！\n")

