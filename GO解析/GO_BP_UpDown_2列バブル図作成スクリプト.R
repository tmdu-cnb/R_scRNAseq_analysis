#!/usr/bin/env Rscript
# =============================================================
# GO(BP) Up/Down 比較バブル図（横軸2列）
# - 横軸: Up / Down
# - 縦軸: GO term（BP）
# - バブルサイズ: GeneRatio
# - 色: p.adjust
# - フィルタ: Count >= 3
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

df <- df %>%
  mutate(
    Direction = factor(Direction, levels = c("Down", "Up")),
    Description = factor(Description, levels = rev(term_order))
  )

# 出力データも保存
write_csv(df %>% select(Direction, ID, Description, GeneRatio, GeneRatio_num, Count, p.adjust, qvalue, geneID),
          file.path(out_dir, "Figure_GO_BP_UpDown_2col_Count3_source.csv"))

p <- ggplot(df, aes(x = Direction, y = Description)) +
  geom_point(aes(size = GeneRatio_num, color = p.adjust), alpha = 0.95) +
  scale_color_viridis_c(option = "C", direction = -1, name = "p.adjust") +
  scale_size_continuous(name = "GeneRatio", range = c(2.5, 10)) +
  labs(
    title = "GO BP enrichment (Up vs Down)",
    subtitle = "Filtered to Count ≥ 3; bubble size = GeneRatio",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 11, face = "bold")
  )

ggsave(file.path(out_dir, "Figure_GO_BP_UpDown_2col_Count3.png"), p, width = 9, height = 10, dpi = 300)
ggsave(file.path(out_dir, "Figure_GO_BP_UpDown_2col_Count3.pdf"), p, width = 9, height = 10)

message("✓ 出力完了: ", out_dir)


