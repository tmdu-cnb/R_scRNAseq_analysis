# KEGG解析とGO解析の結果から、各項目と対応する遺伝子リストを作成するスクリプト

library(dplyr)
library(tidyr)

# データの読み込み
go_down <- read.csv("GO解析_ダウンレギュレーション_BP.csv", stringsAsFactors = FALSE)
go_up <- read.csv("GO解析_アップレギュレーション_BP.csv", stringsAsFactors = FALSE)
kegg_down <- read.csv("KEGG解析_ダウンレギュレーション.csv", stringsAsFactors = FALSE)
kegg_up <- read.csv("KEGG解析_アップレギュレーション.csv", stringsAsFactors = FALSE)

# 遺伝子IDをリストに変換する関数
parse_gene_ids <- function(gene_id_string) {
  if (is.na(gene_id_string) || gene_id_string == "") {
    return(character(0))
  }
  # スラッシュで区切られた遺伝子IDを分割
  genes <- strsplit(gene_id_string, "/")[[1]]
  return(genes)
}

# ENTREZIDを遺伝子シンボルに変換する関数（org.Mm.eg.dbを使用）
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  cat("警告: org.Mm.eg.dbパッケージがインストールされていません。ENTREZIDのまま出力します。\n")
  convert_to_symbol <- function(entrez_ids) {
    return(entrez_ids)
  }
} else {
  library(org.Mm.eg.db)
  convert_to_symbol <- function(entrez_ids) {
    if (length(entrez_ids) == 0) {
      return(character(0))
    }
    # ENTREZIDを遺伝子シンボルに変換
    symbols <- tryCatch({
      mapIds(org.Mm.eg.db, 
             keys = as.character(entrez_ids),
             column = "SYMBOL",
             keytype = "ENTREZID",
             multiVals = "first")
    }, error = function(e) {
      return(entrez_ids)
    })
    # NAを元のENTREZIDに置き換え
    symbols[is.na(symbols)] <- entrez_ids[is.na(symbols)]
    return(symbols)
  }
}

# GO解析（ダウンレギュレーション）の処理（上位10項目のみ）
cat("=== GO解析（ダウンレギュレーション）の遺伝子リスト作成中（上位10項目）... ===\n")
# p.adjustでソートして上位10項目を取得
go_down_sorted <- go_down %>% arrange(p.adjust) %>% head(10)
go_down_list <- lapply(1:nrow(go_down_sorted), function(i) {
  term_id <- go_down_sorted$ID[i]
  description <- go_down_sorted$Description[i]
  gene_ids <- parse_gene_ids(go_down_sorted$geneID[i])
  gene_symbols <- convert_to_symbol(gene_ids)
  
  data.frame(
    Category = "GO_Down",
    ID = term_id,
    Description = description,
    Gene_ENTREZID = paste(gene_ids, collapse = "; "),
    Gene_Symbol = paste(gene_symbols, collapse = "; "),
    Count = length(gene_ids),
    p.adjust = go_down_sorted$p.adjust[i],
    stringsAsFactors = FALSE
  )
})
go_down_df <- bind_rows(go_down_list)

# GO解析（アップレギュレーション）の処理（上位10項目のみ）
cat("=== GO解析（アップレギュレーション）の遺伝子リスト作成中（上位10項目）... ===\n")
# p.adjustでソートして上位10項目を取得
go_up_sorted <- go_up %>% arrange(p.adjust) %>% head(10)
go_up_list <- lapply(1:nrow(go_up_sorted), function(i) {
  term_id <- go_up_sorted$ID[i]
  description <- go_up_sorted$Description[i]
  gene_ids <- parse_gene_ids(go_up_sorted$geneID[i])
  gene_symbols <- convert_to_symbol(gene_ids)
  
  data.frame(
    Category = "GO_Up",
    ID = term_id,
    Description = description,
    Gene_ENTREZID = paste(gene_ids, collapse = "; "),
    Gene_Symbol = paste(gene_symbols, collapse = "; "),
    Count = length(gene_ids),
    p.adjust = go_up_sorted$p.adjust[i],
    stringsAsFactors = FALSE
  )
})
go_up_df <- bind_rows(go_up_list)

# KEGG解析（ダウンレギュレーション）の処理
cat("=== KEGG解析（ダウンレギュレーション）の遺伝子リスト作成中... ===\n")
kegg_down_list <- lapply(1:nrow(kegg_down), function(i) {
  pathway_id <- kegg_down$ID[i]
  description <- kegg_down$Description[i]
  gene_ids <- parse_gene_ids(kegg_down$geneID[i])
  gene_symbols <- convert_to_symbol(gene_ids)
  
  data.frame(
    Category = "KEGG_Down",
    ID = pathway_id,
    Description = description,
    Gene_ENTREZID = paste(gene_ids, collapse = "; "),
    Gene_Symbol = paste(gene_symbols, collapse = "; "),
    Count = length(gene_ids),
    p.adjust = kegg_down$p.adjust[i],
    stringsAsFactors = FALSE
  )
})
kegg_down_df <- bind_rows(kegg_down_list)

# KEGG解析（アップレギュレーション）の処理
cat("=== KEGG解析（アップレギュレーション）の遺伝子リスト作成中... ===\n")
kegg_up_list <- lapply(1:nrow(kegg_up), function(i) {
  pathway_id <- kegg_up$ID[i]
  description <- kegg_up$Description[i]
  gene_ids <- parse_gene_ids(kegg_up$geneID[i])
  gene_symbols <- convert_to_symbol(gene_ids)
  
  data.frame(
    Category = "KEGG_Up",
    ID = pathway_id,
    Description = description,
    Gene_ENTREZID = paste(gene_ids, collapse = "; "),
    Gene_Symbol = paste(gene_symbols, collapse = "; "),
    Count = length(gene_ids),
    p.adjust = kegg_up$p.adjust[i],
    stringsAsFactors = FALSE
  )
})
kegg_up_df <- bind_rows(kegg_up_list)

# 全ての結果を統合
all_results <- bind_rows(
  go_down_df,
  go_up_df,
  kegg_down_df,
  kegg_up_df
)

# p.adjustでソート
all_results <- all_results %>%
  arrange(Category, p.adjust)

# CSVファイルとして保存
write.csv(all_results, "KEGG_GO_遺伝子リスト.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("✓ KEGG_GO_遺伝子リスト.csvを保存しました\n")

# カテゴリー別に分割して保存
write.csv(go_down_df, "GO解析_ダウンレギュレーション_遺伝子リスト.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("✓ GO解析_ダウンレギュレーション_遺伝子リスト.csvを保存しました\n")

write.csv(go_up_df, "GO解析_アップレギュレーション_遺伝子リスト.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("✓ GO解析_アップレギュレーション_遺伝子リスト.csvを保存しました\n")

write.csv(kegg_down_df, "KEGG解析_ダウンレギュレーション_遺伝子リスト.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("✓ KEGG解析_ダウンレギュレーション_遺伝子リスト.csvを保存しました\n")

write.csv(kegg_up_df, "KEGG解析_アップレギュレーション_遺伝子リスト.csv", row.names = FALSE, fileEncoding = "UTF-8")
cat("✓ KEGG解析_アップレギュレーション_遺伝子リスト.csvを保存しました\n")

# サマリーを表示
cat("\n=== サマリー ===\n")
cat(sprintf("GO解析（ダウンレギュレーション）: %d個のGO term（図に表示されている上位10項目）\n", nrow(go_down_df)))
cat(sprintf("GO解析（アップレギュレーション）: %d個のGO term（図に表示されている上位10項目）\n", nrow(go_up_df)))
cat(sprintf("KEGG解析（ダウンレギュレーション）: %d個のパスウェイ\n", nrow(kegg_down_df)))
cat(sprintf("KEGG解析（アップレギュレーション）: %d個のパスウェイ\n", nrow(kegg_up_df)))
cat(sprintf("合計: %d個の項目\n", nrow(all_results)))

