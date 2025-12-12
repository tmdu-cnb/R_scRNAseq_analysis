# KEGG_GO_遺伝子リストを表形式の画像として出力するスクリプト

library(dplyr)
library(gridExtra)
library(grid)
library(ggplot2)

# データの読み込み
data <- read.csv("KEGG_GO_遺伝子リスト.csv", stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# データの前処理
# 遺伝子リストが長い場合は省略表示
data$Gene_Symbol_Short <- ifelse(
  nchar(data$Gene_Symbol) > 50,
  paste0(substr(data$Gene_Symbol, 1, 47), "..."),
  data$Gene_Symbol
)

# カテゴリーごとに分割
go_down <- data %>% filter(Category == "GO_Down")
go_up <- data %>% filter(Category == "GO_Up")
kegg_down <- data %>% filter(Category == "KEGG_Down")
kegg_up <- data %>% filter(Category == "KEGG_Up")

# 表を作成する関数
create_table <- function(df, title) {
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  # 表示用のデータフレームを作成
  display_df <- df %>%
    select(ID, Description, Gene_Symbol_Short, Count, p.adjust) %>%
    mutate(
      p.adjust = sprintf("%.4f", p.adjust),
      Description = ifelse(nchar(Description) > 40, 
                          paste0(substr(Description, 1, 37), "..."), 
                          Description)
    )
  
  # 列名を日本語に変更
  colnames(display_df) <- c("ID", "説明", "遺伝子", "遺伝子数", "p.adjust")
  
  # テーブルグロブを作成
  table_grob <- tableGrob(
    display_df,
    rows = NULL,
    theme = ttheme_default(
      base_size = 8,
      padding = unit(c(4, 4), "mm"),
      colhead = list(
        bg_params = list(fill = "#4472C4", alpha = 0.7),
        fg_params = list(col = "white", fontface = "bold")
      ),
      core = list(
        fg_params = list(col = "black"),
        bg_params = list(fill = c("#F2F2F2", "white"))
      )
    )
  )
  
  # タイトルを追加
  title_grob <- textGrob(
    title,
    gp = gpar(fontsize = 14, fontface = "bold"),
    y = unit(0.5, "npc"),
    vjust = 0.5
  )
  
  # タイトルとテーブルを結合
  combined <- gTree(
    children = gList(
      title_grob,
      table_grob
    ),
    vp = viewport(
      layout = grid.layout(
        nrow = 2,
        ncol = 1,
        heights = unit(c(0.5, 1), c("cm", "null"))
      )
    )
  )
  
  return(combined)
}

# 各カテゴリーの表を作成
cat("表画像を作成中...\n")

# 統合表を作成（4つの表を1つの画像に）
create_combined_table <- function() {
  tables <- list()
  
  if (nrow(go_down) > 0) {
    display_df_down <- go_down %>%
      select(ID, Description, Gene_Symbol_Short, Count, p.adjust) %>%
      mutate(
        p.adjust = sprintf("%.4f", p.adjust),
        Description = ifelse(nchar(Description) > 35, 
                            paste0(substr(Description, 1, 32), "..."), 
                            Description),
        Category = "GO_Down"
      )
    colnames(display_df_down) <- c("ID", "説明", "遺伝子", "遺伝子数", "p.adjust", "Category")
    tables[["GO_Down"]] <- display_df_down
  }
  
  if (nrow(go_up) > 0) {
    display_df_up <- go_up %>%
      select(ID, Description, Gene_Symbol_Short, Count, p.adjust) %>%
      mutate(
        p.adjust = sprintf("%.4f", p.adjust),
        Description = ifelse(nchar(Description) > 35, 
                            paste0(substr(Description, 1, 32), "..."), 
                            Description),
        Category = "GO_Up"
      )
    colnames(display_df_up) <- c("ID", "説明", "遺伝子", "遺伝子数", "p.adjust", "Category")
    tables[["GO_Up"]] <- display_df_up
  }
  
  if (nrow(kegg_down) > 0) {
    display_df_kegg_down <- kegg_down %>%
      select(ID, Description, Gene_Symbol_Short, Count, p.adjust) %>%
      mutate(
        p.adjust = sprintf("%.4f", p.adjust),
        Description = ifelse(nchar(Description) > 35, 
                            paste0(substr(Description, 1, 32), "..."), 
                            Description),
        Category = "KEGG_Down"
      )
    colnames(display_df_kegg_down) <- c("ID", "説明", "遺伝子", "遺伝子数", "p.adjust", "Category")
    tables[["KEGG_Down"]] <- display_df_kegg_down
  }
  
  if (nrow(kegg_up) > 0) {
    display_df_kegg_up <- kegg_up %>%
      select(ID, Description, Gene_Symbol_Short, Count, p.adjust) %>%
      mutate(
        p.adjust = sprintf("%.4f", p.adjust),
        Description = ifelse(nchar(Description) > 35, 
                            paste0(substr(Description, 1, 32), "..."), 
                            Description),
        Category = "KEGG_Up"
      )
    colnames(display_df_kegg_up) <- c("ID", "説明", "遺伝子", "遺伝子数", "p.adjust", "Category")
    tables[["KEGG_Up"]] <- display_df_kegg_up
  }
  
  return(tables)
}

# 統合表を作成
tables <- create_combined_table()

# ページサイズを設定（A4横向き）
png("KEGG_GO_遺伝子リスト_表.png", 
    width = 16, height = 11, units = "in", res = 300)

grid.newpage()

# レイアウトを設定
pushViewport(viewport(layout = grid.layout(
  nrow = 2, ncol = 2,
  heights = unit(c(1, 1), "null"),
  widths = unit(c(1, 1), "null"),
  respect = TRUE
)))

# GO解析（ダウンレギュレーション）
if (!is.null(tables[["GO_Down"]])) {
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  title_grob1 <- textGrob("GO解析（ダウンレギュレーション）", 
                          gp = gpar(fontsize = 12, fontface = "bold", col = "#4472C4"),
                          y = unit(0.98, "npc"), vjust = 1)
  grid.draw(title_grob1)
  
  table_grob1 <- tableGrob(
    tables[["GO_Down"]][, 1:5],
    rows = NULL,
    theme = ttheme_default(
      base_size = 7,
      padding = unit(c(2, 3), "mm"),
      colhead = list(
        bg_params = list(fill = "#4472C4", alpha = 0.8),
        fg_params = list(col = "white", fontface = "bold")
      ),
      core = list(
        fg_params = list(col = "black"),
        bg_params = list(fill = c("#E8F0F8", "white"))
      )
    )
  )
  grid.draw(table_grob1)
  popViewport()
}

# GO解析（アップレギュレーション）
if (!is.null(tables[["GO_Up"]])) {
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
  title_grob2 <- textGrob("GO解析（アップレギュレーション）", 
                          gp = gpar(fontsize = 12, fontface = "bold", col = "#4472C4"),
                          y = unit(0.98, "npc"), vjust = 1)
  grid.draw(title_grob2)
  
  table_grob2 <- tableGrob(
    tables[["GO_Up"]][, 1:5],
    rows = NULL,
    theme = ttheme_default(
      base_size = 7,
      padding = unit(c(2, 3), "mm"),
      colhead = list(
        bg_params = list(fill = "#4472C4", alpha = 0.8),
        fg_params = list(col = "white", fontface = "bold")
      ),
      core = list(
        fg_params = list(col = "black"),
        bg_params = list(fill = c("#E8F0F8", "white"))
      )
    )
  )
  grid.draw(table_grob2)
  popViewport()
}

# KEGG解析（ダウンレギュレーション）
if (!is.null(tables[["KEGG_Down"]])) {
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  title_grob3 <- textGrob("KEGG解析（ダウンレギュレーション）", 
                          gp = gpar(fontsize = 12, fontface = "bold", col = "#70AD47"),
                          y = unit(0.98, "npc"), vjust = 1)
  grid.draw(title_grob3)
  
  table_grob3 <- tableGrob(
    tables[["KEGG_Down"]][, 1:5],
    rows = NULL,
    theme = ttheme_default(
      base_size = 7,
      padding = unit(c(2, 3), "mm"),
      colhead = list(
        bg_params = list(fill = "#70AD47", alpha = 0.8),
        fg_params = list(col = "white", fontface = "bold")
      ),
      core = list(
        fg_params = list(col = "black"),
        bg_params = list(fill = c("#E8F5E9", "white"))
      )
    )
  )
  grid.draw(table_grob3)
  popViewport()
}

# KEGG解析（アップレギュレーション）
if (!is.null(tables[["KEGG_Up"]])) {
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
  title_grob4 <- textGrob("KEGG解析（アップレギュレーション）", 
                          gp = gpar(fontsize = 12, fontface = "bold", col = "#70AD47"),
                          y = unit(0.98, "npc"), vjust = 1)
  grid.draw(title_grob4)
  
  table_grob4 <- tableGrob(
    tables[["KEGG_Up"]][, 1:5],
    rows = NULL,
    theme = ttheme_default(
      base_size = 7,
      padding = unit(c(2, 3), "mm"),
      colhead = list(
        bg_params = list(fill = "#70AD47", alpha = 0.8),
        fg_params = list(col = "white", fontface = "bold")
      ),
      core = list(
        fg_params = list(col = "black"),
        bg_params = list(fill = c("#E8F5E9", "white"))
      )
    )
  )
  grid.draw(table_grob4)
  popViewport()
}

popViewport()
dev.off()
cat("✓ KEGG_GO_遺伝子リスト_表.pngを保存しました\n")

# より詳細な表（遺伝子リスト全体を含む）を別ファイルとして作成
# 各カテゴリーごとに個別の画像を作成
create_detailed_table <- function(df, title, filename) {
  if (nrow(df) == 0) {
    return(NULL)
  }
  
  # 表示用のデータフレームを作成
  display_df <- df %>%
    select(ID, Description, Gene_Symbol, Count, p.adjust) %>%
    mutate(
      p.adjust = sprintf("%.4f", p.adjust)
    )
  
  # 列名を日本語に変更
  colnames(display_df) <- c("ID", "説明", "遺伝子リスト", "遺伝子数", "p.adjust")
  
  # カテゴリーに応じた色を決定
  if (grepl("GO", title)) {
    header_color <- "#4472C4"
    bg_color <- c("#E8F0F8", "white")
    title_color <- "#4472C4"
  } else {
    header_color <- "#70AD47"
    bg_color <- c("#E8F5E9", "white")
    title_color <- "#70AD47"
  }
  
  # テーブルグロブを作成（フォントサイズを小さく）
  table_grob <- tableGrob(
    display_df,
    rows = NULL,
    theme = ttheme_default(
      base_size = 7,
      padding = unit(c(3, 3), "mm"),
      colhead = list(
        bg_params = list(fill = header_color, alpha = 0.8),
        fg_params = list(col = "white", fontface = "bold")
      ),
      core = list(
        fg_params = list(col = "black"),
        bg_params = list(fill = bg_color)
      )
    )
  )
  
  # PNGファイルとして保存
  png(filename, width = 14, height = max(8, nrow(df) * 0.5), units = "in", res = 300)
  grid.newpage()
  grid.draw(table_grob)
  title_grob <- textGrob(
    title,
    gp = gpar(fontsize = 16, fontface = "bold", col = title_color),
    y = unit(0.98, "npc"),
    vjust = 1
  )
  grid.draw(title_grob)
  dev.off()
  
  cat(sprintf("✓ %sを保存しました\n", filename))
}

# 詳細な表を個別に作成
if (nrow(go_down) > 0) {
  create_detailed_table(go_down, "GO解析（ダウンレギュレーション）", 
                       "GO解析_ダウンレギュレーション_詳細表.png")
}

if (nrow(go_up) > 0) {
  create_detailed_table(go_up, "GO解析（アップレギュレーション）", 
                       "GO解析_アップレギュレーション_詳細表.png")
}

if (nrow(kegg_down) > 0) {
  create_detailed_table(kegg_down, "KEGG解析（ダウンレギュレーション）", 
                       "KEGG解析_ダウンレギュレーション_詳細表.png")
}

if (nrow(kegg_up) > 0) {
  create_detailed_table(kegg_up, "KEGG解析（アップレギュレーション）", 
                       "KEGG解析_アップレギュレーション_詳細表.png")
}

cat("\n完了！\n")

