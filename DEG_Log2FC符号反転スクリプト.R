# =============================================================
# DEGデータのLog2FC符号反転スクリプト
# CTL.DTA_log2FC → DTA.CTL_log2FC（符号を反転）
# =============================================================

# ライブラリ
library(dplyr)

# 入力ファイル
input_file_down <- "scRNAseqからidepでDEG_FDR0.1以下_ダウンレギュレーション.csv"
input_file_up <- "scRNAseqからidepでDEG_FDR0.1以下_アップレギュレーション.csv"
input_file_all <- "scRNAseqからidepでDEG_FDR0.1以下.csv"

# 出力ファイル
output_file_down <- "scRNAseqからidepでDEG_FDR0.1以下_ダウンレギュレーション.csv"
output_file_up <- "scRNAseqからidepでDEG_FDR0.1以下_アップレギュレーション.csv"
output_file_all <- "scRNAseqからidepでDEG_FDR0.1以下.csv"

cat("=============================================================\n")
cat("Log2FC符号反転スクリプト\n")
cat("CTL.DTA_log2FC → DTA.CTL_log2FC（符号を反転）\n")
cat("=============================================================\n\n")

# 関数：Log2FCを反転
invert_log2fc <- function(data) {
  if ("CTL.DTA_log2FC" %in% colnames(data)) {
    # 符号を反転
    data$DTA.CTL_log2FC <- -data$CTL.DTA_log2FC
    
    # カラム名を変更（オプション：元のカラムを残すか削除するか）
    # ここでは両方残す
    cat("  Log2FCを反転しました\n")
    cat("  CTL.DTA_log2FC → DTA.CTL_log2FC（符号反転）\n")
  } else {
    cat("  警告: CTL.DTA_log2FCカラムが見つかりません\n")
  }
  return(data)
}

# 1. 全データ
if (file.exists(input_file_all)) {
  cat("1. 全データを処理中...\n")
  data_all <- read.csv(input_file_all, stringsAsFactors = FALSE)
  cat("  読み込み完了: ", nrow(data_all), "行\n")
  
  data_all <- invert_log2fc(data_all)
  
  write.csv(data_all, file = output_file_all, row.names = FALSE)
  cat("  ✓ 保存完了: ", output_file_all, "\n\n")
}

# 2. ダウンレギュレーション
if (file.exists(input_file_down)) {
  cat("2. ダウンレギュレーションデータを処理中...\n")
  data_down <- read.csv(input_file_down, stringsAsFactors = FALSE)
  cat("  読み込み完了: ", nrow(data_down), "行\n")
  
  cat("  反転前のLog2FC範囲: ", 
      min(data_down$CTL.DTA_log2FC, na.rm = TRUE), " ～ ", 
      max(data_down$CTL.DTA_log2FC, na.rm = TRUE), "\n")
  
  data_down <- invert_log2fc(data_down)
  
  cat("  反転後のLog2FC範囲: ", 
      min(data_down$DTA.CTL_log2FC, na.rm = TRUE), " ～ ", 
      max(data_down$DTA.CTL_log2FC, na.rm = TRUE), "\n")
  
  write.csv(data_down, file = output_file_down, row.names = FALSE)
  cat("  ✓ 保存完了: ", output_file_down, "\n\n")
}

# 3. アップレギュレーション
if (file.exists(input_file_up)) {
  cat("3. アップレギュレーションデータを処理中...\n")
  data_up <- read.csv(input_file_up, stringsAsFactors = FALSE)
  cat("  読み込み完了: ", nrow(data_up), "行\n")
  
  cat("  反転前のLog2FC範囲: ", 
      min(data_up$CTL.DTA_log2FC, na.rm = TRUE), " ～ ", 
      max(data_up$CTL.DTA_log2FC, na.rm = TRUE), "\n")
  
  data_up <- invert_log2fc(data_up)
  
  cat("  反転後のLog2FC範囲: ", 
      min(data_up$DTA.CTL_log2FC, na.rm = TRUE), " ～ ", 
      max(data_up$DTA.CTL_log2FC, na.rm = TRUE), "\n")
  
  write.csv(data_up, file = output_file_up, row.names = FALSE)
  cat("  ✓ 保存完了: ", output_file_up, "\n\n")
}

cat("=============================================================\n")
cat("符号反転完了！\n")
cat("\n新しい解釈:\n")
cat("  DTA.CTL_log2FC > 0（プラス） → DTA > CTL → DTAでアップレギュレーション\n")
cat("  DTA.CTL_log2FC < 0（マイナス） → DTA < CTL → DTAでダウンレギュレーション\n")
cat("=============================================================\n")

