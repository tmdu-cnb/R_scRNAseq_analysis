# =============================================================
# 3. Purkinje細胞を上書きで赤点・枠付き表示 ----------------------

# UMAP座標を抽出して列名を揃える
umap_dta <- Embeddings(seu_dta, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")

umap_ctl <- Embeddings(seu_ctl, "umap") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")

# UMAP座標列が1,2列目にあることを保証
if (!all(c("UMAP_1", "UMAP_2") %in% colnames(umap_dta))) {
  colnames(umap_dta)[2:3] <- c("UMAP_1", "UMAP_2")
}
if (!all(c("UMAP_1", "UMAP_2") %in% colnames(umap_ctl))) {
  colnames(umap_ctl)[2:3] <- c("UMAP_1", "UMAP_2")
}

# Purkinje細胞の座標だけ抽出
purk_dta <- umap_dta[seu_dta$Purkinje_tag == "Purkinje", c("UMAP_1", "UMAP_2")]
purk_ctl <- umap_ctl[seu_ctl$Purkinje_tag == "Purkinje", c("UMAP_1", "UMAP_2")]

# Purkinje上書き
p_dta_highlight <- p_dta_base +
  geom_point(data = purk_dta, aes(x = UMAP_1, y = UMAP_2),
             color = "red3", size = 1.8, shape = 21, stroke = 0.7, fill = NA) +
  ggtitle("DTA: UMAP with highlighted Purkinje (Calb1>10 & Car8>10)")

p_ctl_highlight <- p_ctl_base +
  geom_point(data = purk_ctl, aes(x = UMAP_1, y = UMAP_2),
             color = "red3", size = 1.8, shape = 21, stroke = 0.7, fill = NA) +
  ggtitle("CTL: UMAP with highlighted Purkinje (Calb1>10 & Car8>10)")

# =============================================================
# 4. 並列で表示 ------------------------------------------------
library(patchwork)
p_dta_highlight + p_ctl_highlight
