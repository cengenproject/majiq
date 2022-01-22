# test UMAP and heatmaps on gene AND splicing


# Init ----
library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(281)

#~ Retrieve gene for individual sample ----
gene_expr <- read_delim("data/genes/bsn9_featureCoounts_ab2828_112821.txt",
           comment = "#") |>
  select(-c("Chr","Start","End","Strand","Length")) |>
  pivot_longer(-Geneid,
               names_to = "replicate",
               values_to = "count") |>
  mutate(replicate = str_match(replicate, "^\\.\\./bsn9_bams/([A-Z0-9ef]+r[0-9]+t?[12]?)\\.bam$")[,2],
         sample = str_match(replicate, "^([A-Z0-9ef]+r[0-9]+)t?[12]?$")[,2]) |>
  group_by(Geneid, sample) |>
  summarize(count = sum(count),
            .groups = "drop") |>
  pivot_wider(names_from = "sample",
              values_from = "count") |>
  column_to_rownames("Geneid")

# saveRDS(gene_expr, "data/genes/bsn9_featureCoounts_ab2828_112821.tibble.rds")
gene_expr <- readRDS("data/genes/bsn9_featureCoounts_ab2828_112821.tibble.rds")

m_gene_expr <- gene_expr |>
  select(-c("RICr133","PVMr122")) |>
  as.matrix()





# select variable features

all_mads <- matrixStats::rowMads(m_gene_expr)
hist(log1p(all_mads), breaks = 150); abline(v=4, col="red3")

table(log1p(all_mads) >= 4)
table(all_mads >= exp(4)-1)
gene_expr_variable <- gene_expr[all_mads >= exp(4)-1,]
dim(gene_expr_variable)

pheatmap::pheatmap(log1p(gene_expr_variable))
plot(hclust(dist(log1p(t(gene_expr_variable)))))

#

#


library(Seurat)
seu <- CreateSeuratObject(m_gene_expr)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)


top10 <- head(VariableFeatures(seu), 100)
i2s(top10, gids)
VariableFeaturePlot(seu)




seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, npcs=80, verbose = FALSE)

DimPlot(seu, reduction = "pca")
ElbowPlot(seu, ndims = 80)



seu <- FindNeighbors(seu, dims = 1:40)
seu <- FindClusters(seu, resolution = 1.5)
seu <- RunUMAP(seu, dims = 1:40,n.neighbors = 3)
DimPlot(seu, reduction = "umap")

seu[["neuron"]] <- str_match(colnames(seu), "^([A-Z0-9ef]+)r[0-9]{1,3}$")[,2]

DimPlot(seu, group.by = "neuron",label = TRUE,repel = TRUE,pt.size = 2) + NoLegend()



tibble(seurat_clusters = seu$seurat_clusters,
       neuron = seu$neuron) |>
  ggplot() +
  theme_bw() +
  geom_jitter(aes(x=neuron,y=seurat_clusters), alpha=.3, width = .2, height = .2) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ))




#~ Retrieve signif PSI ----
dpsi <- readRDS("intermediates/211130_dpsi.rds")

variable_lsvs <- dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  pull(lsv_id) |>
  unique()

lsv_belonging_to_cassette <- readRDS("intermediates/220118_dmod.rds") |>
  filter(cassette > 0, complex == FALSE) |>
  pull(lsv_id) |>
  map(~ str_split(.x, ";")) |>
  unlist() |>
  unique()

lsvs_to_analyze <- intersect(variable_lsvs, lsv_belonging_to_cassette)


# Retrieve psi for neurons quantified
get_psi <- function(lsv_id, neuron, gene){
  tryCatch(rhdf5::h5read(paste0("data/2021-11-30_outs/psi/",neuron,".psi.voila"),
                         paste0("/lsvs/",gene,"/",lsv_id,"/means"))[1],
           error = \(e) NA_real_)
}

gene_expr <- read.delim("data/genes/aggr_ave_integrant_GeTMM_011822_v2.tsv")
neurs_integrated <- colnames(gene_expr)

lsv_psi <- expand_grid(lsv_id = lsvs_to_analyze,
                       neuron = neurs_integrated) |>
  mutate(gene = str_split_fixed(lsv_id, ":", 2)[,1])

lsv_psi$psi <- pmap_dbl(lsv_psi, get_psi)





# Find "interesting" LSVs: highest MAD; keep single best LSV per gene
xx <- lsv_psi |>
  group_by(lsv_id) |>
  summarize(mad = mad(psi, na.rm = TRUE),
            sd = sd(psi, na.rm = TRUE),
            nb = sum(!is.na(psi))) |>
  mutate(gene_id = str_split_fixed(lsv_id,":",2)[,1]) |>
  group_by(gene_id) |>
  slice(which.max(mad)) |>
  ungroup() |>
  filter(nb >= 20,   #avoid those that are measured only in few neurons
         mad > .15) 
  

interesting_lsvs <- xx$lsv_id

m_interesting_lsv <- lsv_psi |>
  filter(lsv_id %in% interesting_lsvs) |>
  select(-gene) |>
  pivot_wider(names_from = "neuron",
              values_from = psi) |>
  column_to_rownames("lsv_id") |>
  as.matrix()

# remove neurons with too many NA
m_interesting_lsv <- m_interesting_lsv[,colSums(is.na(m_interesting_lsv)) < 15]

rownames(m_interesting_lsv) <- rownames(m_interesting_lsv) |>
  (\(.x) str_split_fixed(.x, ":", 2)[,1])() |>
  i2s(gids) |>
  make.unique(sep = "_")


hc_cols <- hclust(as.dist(1-cor(m_interesting_lsv, use = "pairwise.complete.obs")))
hc_rows <- hclust(as.dist(1-cor(t(m_interesting_lsv), use = "pairwise.complete.obs")))



pheatmap::pheatmap(m_interesting_lsv,
                   # scale = 'row',
                   cluster_rows = hc_rows,
                   cluster_cols = hc_cols,
                   cutree_rows = 2,
                   cutree_cols = 2,
                   # filename = "presentations/heatmap_ds_de/top_spl_cassettes.png",
                   width = 8,
                   height = 6
                   )


genes_interesting_lsv2 <- as.matrix(gene_expr)[s2i(hc_rows$labels, gids),
                                              hc_cols$labels]
rownames(genes_interesting_lsv2) <- i2s(rownames(genes_interesting_lsv2), gids)

pheatmap::pheatmap(log10(genes_interesting_lsv2),
                   cluster_rows = hc_rows,
                   cluster_cols = hc_cols,
                   cutree_rows = 2,
                   cutree_cols = 2,
                   # filename = "presentations/heatmap_ds_de/gene_for_top_spl_cassettes.pdf",
                   width = 8,
                   height = 6
                   )

# recluster by gene
ghc_cols <- hclust(as.dist(1-cor(genes_interesting_lsv, use = "pairwise.complete.obs")))
ghc_rows <- hclust(as.dist(1-cor(t(genes_interesting_lsv), use = "pairwise.complete.obs")))

pheatmap::pheatmap(genes_interesting_lsv,
                   scale = "row",
                   cluster_rows = ghc_rows,
                   cluster_cols = ghc_cols,
                   cutree_rows = 2,
                   cutree_cols = 2)










# HClust ----

#~ select variable genes ----
m_gene_expr <- as.matrix(gene_expr)
all_mads <- matrixStats::rowMads(m_gene_expr)
hist(log1p(all_mads), breaks = 150); abline(v=3, col="red3")

table(log1p(all_mads) >= 3)
table(all_mads >= exp(3)-1)
gene_expr_variable <- gene_expr[all_mads >= exp(3)-1,]
dim(gene_expr_variable)

plot(hclust(dist(log1p(t(gene_expr_variable)))))


#~ select variable LSVs ----

variable_lsvs <- dpsi |>
  filter(p20 >.95 & p05 < .001) |>
  pull(lsv_id) |>
  unique()


lsv_psi_var <- expand_grid(lsv_id = variable_lsvs,
                       neuron = neurs_integrated) |>
  mutate(gene = str_split_fixed(lsv_id, ":", 2)[,1])

lsv_psi_var$psi <- pmap_dbl(lsv_psi_var, get_psi)


m_variable_lsv <- lsv_psi_var |>
  select(-gene) |>
  pivot_wider(names_from = "neuron",
              values_from = psi) |>
  column_to_rownames("lsv_id") |>
  as.matrix()


all_mads_lsv <- matrixStats::rowMads(m_variable_lsv, na.rm = TRUE)
hist(log1p(all_mads_lsv), breaks = 150); abline(v=.1, col="red3")

table(log1p(all_mads_lsv) >= .1)
table(all_mads_lsv >= exp(.1)-1)
m_variable_lsv_variable <- m_variable_lsv[which(all_mads_lsv >= exp(.1)-1),]
dim(m_variable_lsv_variable)

plot(hclust(dist(log1p(t(m_variable_lsv_variable)))))


weights <- set_names(neurs_integrated,
                     exp(seq_along(neurs_integrated))) |>
  sort() |>
  names() |>
  as.numeric()

hm_callback <- function(hc, ...){
  as.hclust(reorder(as.dendrogram(hc), wts = weights))
}

hc_gene <- hclust(dist(log1p(t(gene_expr_variable))))
hc_lsv <- hclust(dist(log1p(t(m_variable_lsv_variable))))

plot(hm_callback(hc_lsv), main = "spl")
plot(hm_callback(hc_gene), main = "gene")

pdf("presentations/hclust_ds_de/splicing.pdf",width = 7, height=4)
  plot(hm_callback(hc_lsv), main = "spl", hang = -1)
dev.off()
pdf("presentations/hclust_ds_de/genes.pdf",width = 7, height=4)
  plot(hm_callback(hc_gene), main = "gene", hang = -1)
dev.off()






dendextend::entanglement(dendextend::untangle(as.dendrogram(hc_gene), as.dendrogram(hc_lsv), method = "step2side"))
dendextend::tanglegram(dendextend::untangle(as.dendrogram(hc_gene),
                                            as.dendrogram(hc_lsv),
                                            method = "step2side"),
                       columns_width = c(3,3,3),
                       main_left = "Gene expression",
                       main_right = "Splicing",
                       margin_outer = 2,
                       axes = FALSE,
                       lab.cex = 1.1,
                       highlight_branches_lwd = FALSE)


pdf("presentations/hclust_ds_de/tanglegram.pdf",width = 7, height=6)
dendextend::tanglegram(dendextend::untangle(as.dendrogram(hc_gene),
                                            as.dendrogram(hc_lsv),
                                            method = "step2side"),
                       columns_width = c(3,3,3),
                       main_left = "Gene expression",
                       main_right = "Splicing",
                       margin_outer = 2,
                       axes = FALSE,
                       lab.cex = 1.1,
                       highlight_branches_lwd = FALSE)
dev.off()
png("presentations/hclust_ds_de/tanglegram.png",width = 900, height=650)
dendextend::tanglegram(dendextend::untangle(as.dendrogram(hc_gene),
                                            as.dendrogram(hc_lsv),
                                            method = "step2side"),
                       columns_width = c(3,3,3),
                       main_left = "Gene expression",
                       main_right = "Splicing",
                       margin_outer = 2,
                       axes = FALSE,
                       lab.cex = 1.1,
                       highlight_branches_lwd = FALSE)
dev.off()






