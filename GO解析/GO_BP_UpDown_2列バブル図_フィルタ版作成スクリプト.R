#!/usr/bin/env Rscript
# =============================================================
# GO(BP) Up/Down 比較バブル図（横軸2列）- フィルタ版
# - 横軸: Up / Down
# - 縦軸: GO term（BP）
# - バブルサイズ: GeneRatio
# - 色: p.adjust
# - フィルタ: Count >= 3
# - "positive regulation of muscle hypertrophy" までを表示
# =============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(forcats)
})

out_dir <- "thesis_figures_GO_BP"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ego_down_rds <- "GO_BP_EnrichmentMap_outputs/ego_down_BP.rds"
ego_up_rds   <- "GO_BP_EnrichmentMap_outputs/ego_up_BP.rds"

if (!file.exists(ego_down_rds) || !file.exists(ego_up_rds)) {
  stop("ego_*_BP.rds が見つかりません。先に GO_BP_EnrichmentMap作成スクリプト.R を実行してください。")
}

ego_down <- readRDS(ego_down_rds)
ego_up   <- readRDS(ego_up_rds)

df_down <- as.data.frame(ego_down) %>% mutate(Direction = "Down")
df_up   <- as.data.frame(ego_up)   %>% mutate(Direction = "Up")

df <- bind_rows(df_down, df_up) %>%
  # Count>=3に限定
  filter(Count >= 3) %>%
  # GeneRatioを数値化（"5/61" -> 0.0819）
  mutate(
    GeneRatio_num = {
      parts <- str_split(GeneRatio, "/", simplify = TRUE)
      as.numeric(parts[, 1]) / as.numeric(parts[, 2])
    }
  )

if (nrow(df) == 0) stop("Count>=3 のtermがありません（フィルタ条件を確認してください）。")

# 縦軸termの並び：最小p.adjustで並べる（Up/Downどちらかで強いものを上に）
term_order <- df %>%
  group_by(Description) %>%
  summarise(best_p = min(p.adjust, na.rm = TRUE), best_ratio = max(GeneRatio_num, na.rm = TRUE), .groups = "drop") %>%
  arrange(best_p, desc(best_ratio)) %>%
  pull(Description)

# "positive regulation of muscle hypertrophy" の位置を探す
target_term <- "positive regulation of muscle hypertrophy"
if (target_term %in% term_order) {
  # 該当termの位置（インデックス）を取得
  target_idx <- which(term_order == target_term)
  # その位置まで（含む）のtermのみを保持
  term_order_filtered <- term_order[1:target_idx]
  message("✓ '", target_term, "' までを表示（", length(term_order_filtered), " terms）")
} else {
  # 見つからない場合は警告して全て表示
  warning("'", target_term, "' が見つかりません。全てのtermを表示します。")
  term_order_filtered <- term_order
}

# フィルタリングされたtermのみを保持
df <- df %>%
  filter(Description %in% term_order_filtered) %>%
  mutate(
    Direction = factor(Direction, levels = c("Down", "Up")),
    Description = factor(Description, levels = rev(term_order_filtered))
  )

# 出力データも保存
write_csv(df %>% select(Direction, ID, Description, GeneRatio, GeneRatio_num, Count, p.adjust, qvalue, geneID),
          file.path(out_dir, "Figure_GO_BP_UpDown_2col_Count3_filtered_source.csv"))

p <- ggplot(df, aes(x = Direction, y = Description)) +
  geom_point(aes(size = GeneRatio_num, color = p.adjust), alpha = 0.95) +
  scale_color_viridis_c(option = "C", direction = -1, name = "p.adjust", limits = c(0, 0.05)) +
  scale_size_continuous(name = "GeneRatio", range = c(4, 14)) +
  labs(
    title = "GO BP enrichment (Up vs Down)",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 11, face = "bold"),
    plot.margin = margin(10, 10, 10, 18)
  )

# 高さを動的に調整
n_terms <- length(term_order_filtered)
plot_height <- max(6, min(10, n_terms * 0.35 + 2))

ggsave(file.path(out_dir, "Figure_GO_BP_UpDown_2col_Count3_filtered.png"), p, width = 9, height = plot_height, dpi = 300)
ggsave(file.path(out_dir, "Figure_GO_BP_UpDown_2col_Count3_filtered.pdf"), p, width = 9, height = plot_height)

message("✓ 出力完了: ", out_dir)

