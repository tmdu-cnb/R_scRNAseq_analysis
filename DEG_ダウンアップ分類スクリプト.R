# =============================================================
# DEGデータのダウン/アップレギュレーション分類スクリプト
# FDR ≤ 0.1のデータを、ダウンレギュレーションとアップレギュレーションに分類
# =============================================================

# ライブラリ
library(dplyr)

# 入力ファイル
input_file <- "scRNAseqからidepでDEG_FDR0.1以下.csv"
output_file_down <- "scRNAseqからidepでDEG_FDR0.1以下_ダウンレギュレーション.csv"
output_file_up <- "scRNAseqからidepでDEG_FDR0.1以下_アップレギュレーション.csv"

cat("=============================================================\n")
cat("DEGデータのダウン/アップレギュレーション分類開始\n")
cat("=============================================================\n\n")

# データ読み込み
cat("データを読み込み中...\n")
data <- read.csv(input_file, stringsAsFactors = FALSE)
cat("✓ 読み込み完了: ", nrow(data), "行, ", ncol(data), "列\n\n")

# Log2FCカラムの確認
log2fc_col <- "CTL.DTA_log2FC"
if (!log2fc_col %in% colnames(data)) {
  cat("利用可能なカラム名（最初の10個）:\n")
  print(head(colnames(data), 10))
  stop("エラー: CTL.DTA_log2FCカラムが見つかりません。")
}

log2fc_idx <- which(colnames(data) == log2fc_col)

# Log2FCの統計
cat("Log2FCの統計:\n")
cat("  NAの数: ", sum(is.na(data[[log2fc_idx]])), "\n")
cat("  Log2FC > 0 (CTL > DTA, DTAでダウン): ", sum(!is.na(data[[log2fc_idx]]) & data[[log2fc_idx]] > 0), "\n")
cat("  Log2FC < 0 (CTL < DTA, DTAでアップ): ", sum(!is.na(data[[log2fc_idx]]) & data[[log2fc_idx]] < 0), "\n")
cat("  Log2FC = 0: ", sum(!is.na(data[[log2fc_idx]]) & data[[log2fc_idx]] == 0), "\n\n")

# ダウンレギュレーション（DTAで発現が減少 = CTL > DTA = Log2FC > 0）
data_down <- data[!is.na(data[[log2fc_idx]]) & data[[log2fc_idx]] > 0, ]

# アップレギュレーション（DTAで発現が増加 = CTL < DTA = Log2FC < 0）
data_up <- data[!is.na(data[[log2fc_idx]]) & data[[log2fc_idx]] < 0, ]

cat("分類結果:\n")
cat("  ダウンレギュレーション（DTAで減少）: ", nrow(data_down), "遺伝子\n")
cat("  アップレギュレーション（DTAで増加）: ", nrow(data_up), "遺伝子\n")
cat("  合計: ", nrow(data_down) + nrow(data_up), "遺伝子\n\n")

# 結果を保存
cat("結果を保存中...\n")
write.csv(data_down, file = output_file_down, row.names = FALSE)
cat("✓ ダウンレギュレーション: ", output_file_down, "\n")

write.csv(data_up, file = output_file_up, row.names = FALSE)
cat("✓ アップレギュレーション: ", output_file_up, "\n\n")

# サマリー表示
if (nrow(data_down) > 0) {
  cat("ダウンレギュレーション遺伝子（上位10個）:\n")
  data_down_sorted <- data_down[order(data_down[[log2fc_idx]], decreasing = TRUE), ]
  print(data_down_sorted[1:min(10, nrow(data_down_sorted)), c("symbol", log2fc_col, "CTL.DTA_adjPval")])
  cat("\n")
}

if (nrow(data_up) > 0) {
  cat("アップレギュレーション遺伝子（上位10個）:\n")
  data_up_sorted <- data_up[order(data_up[[log2fc_idx]]), ]
  print(data_up_sorted[1:min(10, nrow(data_up_sorted)), c("symbol", log2fc_col, "CTL.DTA_adjPval")])
  cat("\n")
}

cat("=============================================================\n")
cat("分類完了！\n")
cat("=============================================================\n")

