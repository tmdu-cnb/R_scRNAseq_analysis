# =============================================================
# DEG遺伝子の機能解析スクリプト
# GO解析とパスウェイ解析による遺伝子分類
# =============================================================

# ライブラリ
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)  # マウス用
library(enrichplot)
library(ggplot2)

# データ読み込み
cat("=============================================================\n")
cat("DEG遺伝子の機能解析開始\n")
cat("=============================================================\n\n")

data_down <- read.csv("scRNAseqからidepでDEG_FDR0.1以下_ダウンレギュレーション.csv", stringsAsFactors = FALSE)
data_up <- read.csv("scRNAseqからidepでDEG_FDR0.1以下_アップレギュレーション.csv", stringsAsFactors = FALSE)

cat("ダウンレギュレーション遺伝子数:", nrow(data_down), "\n")
cat("アップレギュレーション遺伝子数:", nrow(data_up), "\n\n")

# 遺伝子シンボルを取得
genes_down <- data_down$symbol
genes_up <- data_up$symbol

cat("遺伝子シンボル例（ダウン）:", paste(head(genes_down, 5), collapse = ", "), "\n")
cat("遺伝子シンボル例（アップ）:", paste(head(genes_up, 5), collapse = ", "), "\n\n")

# 遺伝子シンボルをEntrez IDに変換
cat("遺伝子ID変換中...\n")
gene_mapping_down <- bitr(genes_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
gene_mapping_up <- bitr(genes_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

cat("✓ ダウン: ", nrow(gene_mapping_down), "/", length(genes_down), "遺伝子が変換されました\n")
cat("✓ アップ: ", nrow(gene_mapping_up), "/", length(genes_up), "遺伝子が変換されました\n\n")

# GO解析（Biological Process）
cat("GO解析（Biological Process）を実行中...\n")

ego_down <- enrichGO(
  gene = gene_mapping_down$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_up <- enrichGO(
  gene = gene_mapping_up$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

cat("✓ GO解析完了\n")
cat("  ダウン: ", nrow(ego_down), "個のGO term\n")
cat("  アップ: ", nrow(ego_up), "個のGO term\n\n")

# 結果を保存
if (nrow(ego_down) > 0) {
  write.csv(as.data.frame(ego_down), 
            file = "GO解析_ダウンレギュレーション_BP.csv", 
            row.names = FALSE)
  cat("✓ ダウンレギュレーションGO解析結果を保存\n")
}

if (nrow(ego_up) > 0) {
  write.csv(as.data.frame(ego_up), 
            file = "GO解析_アップレギュレーション_BP.csv", 
            row.names = FALSE)
  cat("✓ アップレギュレーションGO解析結果を保存\n")
}

# KEGG解析
cat("\nKEGG解析を実行中...\n")

# ダウンレギュレーション側は閾値を緩めて実行（p.adjust < 0.1）
kegg_down <- enrichKEGG(
  gene = gene_mapping_down$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.1,  # 閾値を緩める
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

# p.adjust < 0.1でフィルタリング
kegg_down_filtered <- kegg_down[kegg_down$p.adjust < 0.1, ]

# アップレギュレーション側も閾値を緩めて実行（p.adjust < 0.1）
kegg_up <- enrichKEGG(
  gene = gene_mapping_up$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.1,  # 閾値を緩める
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2
)

# p.adjust < 0.1でフィルタリング
kegg_up_filtered <- kegg_up[kegg_up$p.adjust < 0.1, ]

cat("✓ KEGG解析完了\n")
cat("  ダウン（p.adjust < 0.1）: ", nrow(kegg_down_filtered), "個のパスウェイ\n")
cat("  アップ（p.adjust < 0.1）: ", nrow(kegg_up_filtered), "個のパスウェイ\n\n")

# 結果を保存
if (nrow(kegg_down_filtered) > 0) {
  write.csv(as.data.frame(kegg_down_filtered), 
            file = "KEGG解析_ダウンレギュレーション.csv", 
            row.names = FALSE)
  cat("✓ ダウンレギュレーションKEGG解析結果を保存（p.adjust < 0.1）\n")
}

if (nrow(kegg_up_filtered) > 0) {
  write.csv(as.data.frame(kegg_up_filtered), 
            file = "KEGG解析_アップレギュレーション.csv", 
            row.names = FALSE)
  cat("✓ アップレギュレーションKEGG解析結果を保存（p.adjust < 0.1）\n")
}

# 可視化
cat("\n可視化を作成中...\n")

# GO解析結果の可視化（上位10個、GeneRatioを横軸に）
if (nrow(ego_down) > 0) {
  png("GO解析_ダウンレギュレーション_BP.png", width = 12, height = 8, units = "in", res = 300)
  p1 <- barplot(ego_down, showCategory = 10, x = "GeneRatio", title = "GO Biological Process - Downregulated") +
    scale_fill_gradient(low = "#CC0000", high = "#FFE5E5", 
                       name = "p.adjust",
                       guide = guide_colorbar(reverse = TRUE))
  print(p1)
  dev.off()
  cat("✓ ダウンレギュレーションGO解析図を保存\n")
}

if (nrow(ego_up) > 0) {
  png("GO解析_アップレギュレーション_BP.png", width = 12, height = 8, units = "in", res = 300)
  p2 <- barplot(ego_up, showCategory = 10, x = "GeneRatio", title = "GO Biological Process - Upregulated") +
    scale_fill_gradient(low = "#CC0000", high = "#FFE5E5", 
                       name = "p.adjust",
                       guide = guide_colorbar(reverse = TRUE))
  print(p2)
  dev.off()
  cat("✓ アップレギュレーションGO解析図を保存\n")
}

# KEGG解析結果の可視化（GeneRatioを横軸に、青系のカラースケール）
if (nrow(kegg_down_filtered) > 0) {
  png("KEGG解析_ダウンレギュレーション.png", width = 12, height = 8, units = "in", res = 300)
  # enrichResultオブジェクトを使用してbarplotを作成
  # kegg_downからp.adjust < 0.1の項目のみを表示
  p3_base <- barplot(kegg_down, showCategory = min(10, nrow(kegg_down_filtered)), 
                     x = "GeneRatio", title = "KEGG Pathway - Downregulated (p.adjust < 0.1)")
  # 青系のスケールを追加（既存のスケールを上書き、濃薄を逆転）
  # lowとhighを逆にして、小さいp値（より有意）を濃い青に
  p3 <- p3_base + scale_fill_gradient(low = "#0066CC", high = "#E5F0FF", 
                                      name = "p.adjust",
                                      guide = guide_colorbar(reverse = TRUE))
  print(p3)
  dev.off()
  cat("✓ ダウンレギュレーションKEGG解析図を保存（p.adjust < 0.1）\n")
}

if (nrow(kegg_up_filtered) > 0) {
  png("KEGG解析_アップレギュレーション.png", width = 12, height = 8, units = "in", res = 300)
  # enrichResultオブジェクトを使用してbarplotを作成
  # kegg_upからp.adjust < 0.1の項目のみを表示
  p4_base <- barplot(kegg_up, showCategory = min(10, nrow(kegg_up_filtered)), 
                     x = "GeneRatio", title = "KEGG Pathway - Upregulated (p.adjust < 0.1)")
  # 青系のスケールを追加（既存のスケールを上書き、濃薄を逆転）
  p4 <- p4_base + scale_fill_gradient(low = "#0066CC", high = "#E5F0FF", 
                                      name = "p.adjust",
                                      guide = guide_colorbar(reverse = TRUE))
  print(p4)
  dev.off()
  cat("✓ アップレギュレーションKEGG解析図を保存（p.adjust < 0.1）\n")
}

cat("\n=============================================================\n")
cat("機能解析完了！\n")
cat("=============================================================\n")

# 結果のサマリーを表示
cat("\n【ダウンレギュレーション - 上位GO term】\n")
if (nrow(ego_down) > 0) {
  print(head(as.data.frame(ego_down)[, c("Description", "GeneRatio", "p.adjust")], 10))
}

cat("\n【アップレギュレーション - 上位GO term】\n")
if (nrow(ego_up) > 0) {
  print(head(as.data.frame(ego_up)[, c("Description", "GeneRatio", "p.adjust")], 10))
}

cat("\n【ダウンレギュレーション - KEGGパスウェイ】\n")
if (nrow(kegg_down) > 0) {
  print(head(as.data.frame(kegg_down)[, c("Description", "GeneRatio", "p.adjust")], 10))
}

cat("\n【アップレギュレーション - KEGGパスウェイ（p.adjust < 0.1）】\n")
if (nrow(kegg_up_filtered) > 0) {
  print(head(as.data.frame(kegg_up_filtered)[, c("Description", "GeneRatio", "p.adjust")], 10))
}

