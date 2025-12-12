# =============================================================
# DEGデータのFDRフィルタリングスクリプト
# Pval(FDR)が0.1以下の行のみを抽出
# =============================================================

# ライブラリ
library(dplyr)

# 入力ファイル
input_file <- "scRNAseqからidepでDEG.csv"
output_file <- "scRNAseqからidepでDEG_FDR0.1以下.csv"

cat("=============================================================\n")
cat("DEGデータのFDRフィルタリング開始\n")
cat("=============================================================\n\n")

# データ読み込み
cat("データを読み込み中...\n")
data <- read.csv(input_file, stringsAsFactors = FALSE)
cat("✓ 読み込み完了: ", nrow(data), "行, ", ncol(data), "列\n\n")

# CTL-DTA_adjPvalカラムの確認（Rはハイフンをドットに変換するため）
col_name <- "CTL.DTA_adjPval"
if (!col_name %in% colnames(data)) {
  cat("利用可能なカラム名（最初の10個）:\n")
  print(head(colnames(data), 10))
  stop("エラー: CTL.DTA_adjPvalカラムが見つかりません。")
}

# カラムにアクセスするためのインデックスを使用
adjpval_col <- which(colnames(data) == col_name)

# フィルタリング前の統計
cat("フィルタリング前の統計:\n")
cat("  総行数: ", nrow(data), "\n")
na_count <- sum(is.na(data[[adjpval_col]]))
cat("  NAの数: ", na_count, "\n")
below_01 <- sum(!is.na(data[[adjpval_col]]) & data[[adjpval_col]] <= 0.1)
above_01 <- sum(!is.na(data[[adjpval_col]]) & data[[adjpval_col]] > 0.1)
cat("  0.1以下の数: ", below_01, "\n")
cat("  0.1より大きい数: ", above_01, "\n\n")

# FDR <= 0.1の行のみを抽出（NAは除外）
data_filtered <- data[!is.na(data[[adjpval_col]]) & data[[adjpval_col]] <= 0.1, ]

cat("フィルタリング後の統計:\n")
cat("  残存行数: ", nrow(data_filtered), "\n")
cat("  除外行数: ", nrow(data) - nrow(data_filtered), "\n")
cat("  残存率: ", round(nrow(data_filtered) / nrow(data) * 100, 2), "%\n\n")

# 結果を保存
cat("結果を保存中...\n")
write.csv(data_filtered, file = output_file, row.names = FALSE)
cat("✓ 保存完了: ", output_file, "\n\n")

cat("=============================================================\n")
cat("フィルタリング完了！\n")
cat("=============================================================\n")

