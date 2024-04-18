# Load deltas


# Initializations ----
library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(289)


export_dir <- "presentations/2024-03_rerun/"

neuron_properties <- read_csv("data/neuron_properties.csv")





#~ Load ----
#~ deltapsi ----
data_dir <- "data/2024-03-04_outs/deltapsi/"

files_dpsi <- list.files(data_dir,
                         pattern = "\\.tsv$",
                         full.names = FALSE)

dpsidta <- map_dfr(files_dpsi,
                  ~read_tsv(file.path(data_dir, .x),
                            col_names = c("gene_id", "lsv_id", "lsv_type",
                                          "dpsi", "p20", "p05", "psiA", "psiB",
                                          "nb_sj", "nb_exons",
                                          "sj_coords", "ir_coords"),
                            col_types = cols(
                              gene_id = col_character(),
                              lsv_id = col_character(),
                              lsv_type = col_character(),
                              dpsi = col_character(),
                              p20 = col_character(),
                              p05 = col_character(),
                              psiA = col_character(),
                              psiB = col_character(),
                              nb_sj = col_double(),
                              nb_exons = col_double(),
                              sj_coords = col_character(),
                              ir_coords = col_character()
                            ),
                            skip = 1,
                            na = "na",
                            progress = FALSE) |>
                    add_column(neurA = str_split(.x,"[-\\.]")[[1]][1],
                               neurB = str_split(.x,"[-\\.]")[[1]][2]),
                  .progress = TRUE)

# neurs_here <- neuron_properties |>
#   filter(include == "yes") |>
#   pull(Neuron_type)
# 
# dpsidta <- dpsidta |>
#   filter(neurA %in% neurs_here,
#          neurB %in% neurs_here)

# separate data in E(PSI) and Std(PSI) fields
# We get one row per junction
dpsi <- dpsidta |>
  separate_rows(dpsi, p20, p05, psiA, psiB,
                sep = ";") |>
  group_by(neurA, neurB, lsv_id) |>
  mutate(junction_id = row_number()) |>
  ungroup()  |>
  mutate(across(c(dpsi,p20, p05, psiA, psiB), as.double),
         junction_id = factor(junction_id))

# qs::qsave(dpsi, "intermediates/240304_dpsi.qs")





#~ Load PSI ----
data_dir <- "data/2024-03-04_outs/psi/"

files_psi <- list.files(data_dir,
                         pattern = "\\.tsv$",
                         full.names = FALSE)

# xx <- read_tsv(file.path(data_dir, files_psi[[2]]),na = "na")

psidta <- map_dfr(files_psi,
                   ~read_tsv(file.path(data_dir, .x),
                             col_names = c("gene_id", "lsv_id", "lsv_type",
                                           "mean_psi_per_lsv_junction", "stdev_psi_per_lsv_junction",
                                           "nb_sj", "nb_exons",
                                           "sj_coords", "ir_coords"),
                             col_types = cols(
                               gene_id = col_character(),
                               lsv_id = col_character(),
                               lsv_type = col_character(),
                               mean_psi_per_lsv_junction = col_character(),
                               stdev_psi_per_lsv_junction = col_character(),
                               nb_sj = col_double(),
                               nb_exons = col_double(),
                               sj_coords = col_character(),
                               ir_coords = col_character()
                             ),
                             skip = 1,
                             na = "na",
                             progress = FALSE) |>
                     add_column(neuron = str_split_1(.x, "\\.")[[1]]),
                  .progress = TRUE)

psi <- psidta |>
  separate_rows(mean_psi_per_lsv_junction, stdev_psi_per_lsv_junction,
                sep = ";") |>
  group_by(neuron, lsv_id) |>
  mutate(junction_id = row_number()) |>
  ungroup()  |>
  mutate(across(c(mean_psi_per_lsv_junction,stdev_psi_per_lsv_junction), as.double),
         junction_id = factor(junction_id))
# qs::qsave(psi, "intermediates/240304_psi.qs")





# Load precomputed ----

psi <- qs::qread("intermediates/240304_psi.qs")
dpsi <- qs::qread("intermediates/240304_dpsi.qs")
all_neurs_sequenced <- unique(c(dpsi$neurA, dpsi$neurB))

stopifnot(identical(sort(all_neurs_sequenced), sort(unique(psi$neuron))))


#~ Alec genes expressed ----
# Use Alec's integrated GeTMMs to determine what genes are expressed in each neuron type
# Update 2022: use more recent data, already binary
gene_expression <- 1L * (cengenDataSC::cengen_sc_3_bulk > 0)
# gene_expression <- read.delim("data/2024-03-05_alec_integration/bsn12_subtracted_integrated_binarized_expression_withVDDD_FDR0.05_030424.tsv")

neurs_with_known_expr <- colnames(gene_expression)

# stopifnot(all.equal(sort(neurs_with_known_expr),
#                     sort(neurons_here) |> setdiff(c("ADF", "M4", "AWC", "DD", "VD"))))
# 
# 
# table(gene_expression_bk)
# genes_both <- intersect(rownames(gene_expression_bk), rownames(gene_expression))
# neurs_both <- intersect(colnames(gene_expression_bk), colnames(gene_expression))
# 
# table(sc_threshold = gene_expression_bk[genes_both, neurs_both],
#       integr = as.matrix(gene_expression)[genes_both, neurs_both])


#~ filter dpsi based on expression ----
gene_expr_tib <- gene_expression |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  as_tibble() |>
  pivot_longer(-gene_id,
               names_to = "neuron_id",
               values_to = "expressed")

# find nb of coexpressed genes which are DS
dpsi_filt <- dpsi |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurA = "neuron_id")) |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurB = "neuron_id")) |>
  filter(expressed.x > 0, expressed.y > 0) |>
  select(-expressed.x, -expressed.y)

psi_filt <- psi |>
  left_join(gene_expr_tib,
            by = c("gene_id", neuron = "neuron_id")) |>
  filter(expressed > 0) |>
  select(-expressed)



# Alec's integrations

degs <- read.delim("data/2021-11-30_alec_integration/Total_integrated_DEGS_pairwise_113021.tsv") |>
  as_tibble() |>
  rowwise() |>
  mutate(pair = c_across(starts_with("cell_")) |> sort() |> paste0(collapse = "-"))

neurs_with_deg_data <- unique(union(degs$cell_A, degs$cell_B))


genes_degs <- readRDS("data/genes/integrated_significant_genes_pairwise_bsn9_012722.rds") |>
  map(as_tibble, rownames = "gene_id") |>
  imap_dfr(~add_column(.x, pair = .y)) |>
  mutate(pair = str_remove_all(pair, "[()]"))




# >>>>>  Analysis <<<<< ----


# DS genes  ----


ds_genes <- dpsi |>
       filter(p20 >.50 & p05 < .05) |>
       pull(gene_id) |> unique()


dsf_genes <- dpsi_filt |>
  filter(p20 >.50 & p05 < .05) |>
  pull(gene_id) |> unique()

neur_genes <- wormDatasets::genes_by_pattern


table(ds_genes %in% neur_genes$present_in_neurons)

ds_genes |>
  intersect(neur_genes$present_in_neurons) |>
  head() |>
  i2s(gids)
dsf_genes |>
  intersect(neur_genes$present_in_neurons) |>
  head() |>
  i2s(gids)

ds_genes |>
  intersect(neur_genes$nonneuronal) |>
  head() |>
  i2s(gids)

dsf_genes |>
  intersect(neur_genes$nonneuronal) |>
  head() |>
  i2s(gids)




# Nb of genes with DS ----

# nb genes tested
length(unique(dpsi$gene_id))
length(unique(dpsi_filt$gene_id))

# total DS
dpsi_filt |>
  filter(p20 >.50 & p05 < .05) |>
  pull(gene_id) |> unique() |>
  length()

#~ By neur classes ----
nb_signif_genes_by_test <- dpsi_filt |>
  filter(p20 >.50 & p05 < .05) |>
  select(gene_id, neurA, neurB) |>
  distinct() |>
  group_by(neurA, neurB) |>
  summarize(nb_DS_genes = n(),
             .groups = 'drop')

nb_signif_genes_by_test |>
  rename(neurA = neurB,
         neurB = neurA) |>
  bind_rows(nb_signif_genes_by_test) |>
  ggplot() +
  theme_classic() +
  geom_tile(aes(x = neurA, y=neurB, fill = nb_DS_genes))





# Compare literature ----

bib_all <- read_tsv("data/biblio/bib_ds_genes.tsv")$gene_id


# Restrict to genes expressed in more than 2 neurons in our dataset

nb_neurs_where_gene_expr <- rowSums(gene_expression[,all_neurs_sequenced])
genes_in_our_neur_sample <- names(nb_neurs_where_gene_expr)[nb_neurs_where_gene_expr > 2]

plot(eulerr::euler(list(literature = bib_all,
                        `this work` = genes_in_our_neur_sample)),
     quantities = TRUE,
     main = "Number of detectable genes")


bib_all <- intersect(bib_all, genes_in_our_neur_sample)



# ds_genes <- dpsi |>
#   filter(p20 >.50 & p05 < .05) |>
#   pull(gene_id) |>
#   unique()

length(ds_genes)

(litt_plot <- plot(eulerr::euler(list(literature = bib_all,
                                      `this work` = ds_genes)),
                   quantities = TRUE,
                   main = "Number of Differentially Spliced genes"))

# ggsave("compare_literature.pdf", path = export_dir, plot = litt_plot,
#        width = 12, height = 9, units = "cm")


# Gene type ----

#~ GO with background ----

writeLines(genes_in_our_neur_sample, "intermediates/240304_background_genes.txt")
writeClipboard(ds_genes[ds_genes %in% genes_in_our_neur_sample])
# -> use Wormbase enrichment analysis



#~ Gene families ----

bind_rows(wormDatasets::gene_families |>
            select(gene_id, gene_name, family),
          wormDatasets::list_rbp_tamburino2013 |>
              select(gene_id, gene_name) |>
              mutate(family = "RBP")) |>
  filter(gene_id %in% genes_in_our_neur_sample) |>
  mutate(family = recode(
    family,
    DEG_ENaC = "DEG/ENaC channel",
    downstream_GPCR = "GPCR signaling",
    GPCR = "GPCR signaling",
    neuropept_metabo = "neuropeptide signaling",
    nt_degradation = "neuropeptide signaling",
    nt_synthesis = "neuropeptide signaling",
    nt_transporter = "neuropeptide signaling",
    potassium_channel = "potassium channel",
    ribosome = "ribosome subunit",
    RBP = "RNA-binding protein",
    synaptic_ves = "synaptic vesicle",
    TF = "transcription factor",
    trp_channel = "TRP channel"
  )) |>
  mutate(is_ds = gene_id %in% ds_genes) |>
  summarize(nb_ds = sum(is_ds),
            nb_tot = n(),
            prop_ds = nb_ds / nb_tot,
            label = paste0(nb_ds, "/", nb_tot),
            .by = "family") |>
  ggplot() +
  theme_classic() +
  geom_col(aes(x = family, y = prop_ds)) +
  geom_text(aes(x = family, y = .9,
                label = label),
            angle = 90)+
  xlab(NULL) + ylab("Proportion DS genes") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1), limits = c(0,1))

# ggsave("ds_by_family.pdf", path = export_dir,
#        width = 15, height = 15, units = "cm")























# Heatmap ----


annot_df <- neuron_properties |>
  select(-include) |>
  column_to_rownames("Neuron_type")


gene_expr_tib <- gene_expression |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  as_tibble() |>
  pivot_longer(-gene_id,
               names_to = "neuron",
               values_to = "expressed")

# find nb of coexpressed genes which are DS
coexpr_ds <- dpsi |>
  filter(neurA %in% neurs_with_known_expr,
         neurB %in% neurs_with_known_expr) |>
  mutate(ds = (p20 >.50 & p05 < .05)) |>
  summarize(has_ds = any(ds),
            .by = c(neurA, neurB, gene_id)) |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurA = "neuron")) |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurB = "neuron")) |>
  mutate(coexpressed = expressed.x & expressed.y) |>
  filter(! is.na(coexpressed)) |>   # a number of annotated pseudogenes have splicing quantified but not expression
  select(neurA, neurB, gene_id, has_ds, coexpressed) |>
  group_by(neurA, neurB) |>
  summarize(nb_ds = sum(has_ds),
            nb_coexpr = sum(coexpressed),
            nb_both = sum(has_ds & coexpressed),
            nb_total = n(),
            .groups = "drop")

# make matrix of proportion coexpr genes that are DS
hm_coexpr <- coexpr_ds |>
  mutate(prop_ds = 100*nb_both/nb_coexpr) |>
  select(neurA,neurB,prop_ds)


ds_mat <- hm_coexpr |>
  rename(neurA = neurB,
         neurB = neurA) |>
  bind_rows(hm_coexpr) |>
  pivot_wider(names_from = neurB,
              values_from = prop_ds) |>
  arrange(neurA) |>
  (\(.x) select(.x, order(colnames(.x))))() |>
  column_to_rownames("neurA") |> 
  as.matrix()

# Show neurons clustered by similarity of prop of DS
hc1 <- hclust(dist(ds_mat, method = "euclidean"), method = "complete")
plot(hc1)
pheatmap::pheatmap(ds_mat,
                   color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name =
                                                                       "Blues"))(100),
                   scale = "none",
                   breaks = (0:100)/5,
                   cluster_rows = hc1, cluster_cols = hc1,
                   cutree_rows = 4,
                   cutree_cols = 4,
                   main = "Proportion of coexpressed genes DS",
                   # filename = file.path(export_dir, "heatmap_ds.png"),
                   width = 8,
                   height = 7
)


# Show neurons clustered by prop of DS
dist_mat_prop <- ds_mat/max(ds_mat, na.rm = TRUE)
heatmap(dist_mat_prop, Rowv = NA, Colv = NA)
hc <- hclust(as.dist(dist_mat_prop), method = "ward.D")
# plot(hc)

pheatmap::pheatmap(ds_mat,
                   color = colorRampPalette(RColorBrewer::brewer.pal(n = 9, name =
                                                                       "Blues"))(151),
                   scale = "none",
                   breaks = (0:150)/5,
                   cluster_rows = hc,
                   cluster_cols = hc,
                   # cutree_rows = 4,
                   # cutree_cols = 4,
                   # main = "Proportion of DS",
                   legend_breaks = c(0,5,10,15,20, 25, 30),
                   legend_labels = c("0%","5%","10%","15%","20%","25%","30%"),
                   # annotation_row = annot_df,
                   annotation_col = annot_df,
                   border_color = NA,
                   # filename = file.path(export_dir, "heatmap_ds.pdf"),
                   width = 10.1,
                   height = 7.7
)



# Show neurons clustered by nb of DS
nb_ds_mat <-  coexpr_ds |>
  select(neurA = neurB, neurB = neurA, nb_both) |>
  bind_rows(coexpr_ds |> select(neurA, neurB, nb_both)) |>
  pivot_wider(names_from = neurB,
              values_from = nb_both) |>
  arrange(neurA) |>
  (\(.x) select(.x, order(colnames(.x))))() |>
  column_to_rownames("neurA") |>
  as.matrix()


dist_mat_nb <- nb_ds_mat/max(nb_ds_mat, na.rm = TRUE)
heatmap(dist_mat_nb, Rowv = NA, Colv = NA)
plot(hclust(as.dist(dist_mat_nb), method = "complete"))





# Show neurons clustered by similarity of profile

ds_mat |> 
  dist(method = "canberra") |>
  hclust(method = "complete") |>
  plot()

nb_ds_mat |> 
  dist(method = "canberra") |>
  hclust(method = "complete") |>
  plot()


# Pairs of similar neurons ----
# fig S2C
source("R/utils.R")


plot_hc_clustered(hc, h = .2, cex = .7,
                  xlab = "Neuron type", sub = NA, main = "Neuron similarity (by proportion of genes dAS)")
abline(h = .2, lty = 'dashed', col = 'grey80')

# pdf(file.path(export_dir, "dendrogram_similarity.pdf"),
#     width = 10, height = 9)
# 
# plot_hc_clustered(hc, h = .2, cex = .8,
#                   xlab = "Neuron type", sub = NA, main = "Neuron similarity (by proportion of genes dAS)")
# abline(h = .2, lty = 'dashed', col = 'grey80')
# dev.off()








# MDS ----
dist_mat_prop[1:3,1:3]
diag(dist_mat_prop) <- 0
mds_ds <- cmdscale(dist_mat_prop, eig = TRUE)

ggplot(tibble(eig = mds_ds$eig[seq_len(100)], k = seq(along = eig)),
       aes(x = k, y = eig)) + theme_minimal() +
  scale_x_discrete("k", limits = as.factor(seq_len(100))) + 
  geom_bar(stat = "identity", width = 0.5, fill = "#ffd700", col = "#0057b7")

mds_ds_df <- mds_ds$points |> as.data.frame() |> rownames_to_column("neuron")

ggplot(mds_ds_df,
       aes(x = V1, y = V2, label = neuron)) +
  geom_point() +
  ggrepel::geom_text_repel(col = "#0057b7") + coord_fixed() 



#~ Cluster by PCA ----

mat_psi <- psi |>
  mutate(lsv_sj_id = paste0(lsv_id, "_", junction_id)) |>
  select(neuron, lsv_sj_id, mean_psi_per_lsv_junction) |>
  pivot_wider(names_from = "lsv_sj_id",
              values_from = "mean_psi_per_lsv_junction") |>
  column_to_rownames("neuron") |>
  as.matrix()
mat_psi2 <- mat_psi[, colSums(is.na(mat_psi))/nrow(mat_psi) < .5]
mat_psi_imp <- impute::impute.knn(mat_psi2)$data

pca <- prcomp(mat_psi_imp)
plot(pca$x)
library(ggfortify)
autoplot(pca, label = TRUE, shape = FALSE)+ theme_classic()

km <- kmeans(mat_psi_imp, 5)

autoplot(km, data = mat_psi_imp, label = TRUE, shape = FALSE) + theme_classic()

d <- dist(mat_psi_imp)
hc_psi <- hclust(d)

plot(hc_psi)



# Cluster analysis...
# need to redefine clusters first

table(cutree(hc1, k=4))
high_ds <- hc2$labels[cutree(hc2, k=2) == 1]
low_ds <- hc2$labels[cutree(hc2, k=2) == 2]

high_ds_mat <- ds_mat[high_ds, high_ds]
high_ds_mat[lower.tri(high_ds_mat)] <- NA
xx <- as.numeric(high_ds_mat)
xx <- xx[!is.na(xx)]
range(xx)
median(xx)


low_ds_mat <- ds_mat[low_ds, low_ds]
low_ds_mat[lower.tri(low_ds_mat)] <- NA
xx <- as.numeric(low_ds_mat)
xx <- xx[!is.na(xx)]
range(xx)
median(xx)









# # Event broadness ----
# Probably wrong approach: we would prefer to work on PSI, not dPSI
# as it's more interesting to say "this spliceform is un/common" than to
# look at a number of pairs (could be 1 single neuron that differs)
# 
# lsv_dpsi <- dpsi |>
#   mutate(is_ds_jct = (p20 >.50 & p05 < .05)) |>
#   group_by(lsv_id, neurA, neurB) |>
#   summarize(is_ds_lsv = any(is_ds_jct),
#             .groups = 'drop')
# 
# 
# lsv_ds <- lsv_dpsi |>
#   summarize(nb_times_ds = sum(is_ds_lsv),
#             n = length(is_ds_lsv),
#             .by = "lsv_id") |>
#   mutate(prop_ds = nb_times_ds/n)
# 
# hist(lsv_ds$n, breaks = 100)
# hist(lsv_ds$prop_ds, breaks = 100)
# hist(lsv_ds$nb_times_ds, breaks = 100)
# 
# lsv_ds |>
#   mutate(category = case_when(
#     n <= 2 ~ "unmeasured",
#     nb_times_ds == 0 ~ "not_ds",
#     nb_times_ds <= 3 ~ "rare",
#     .default = "common"
#   ))
# 
# 
# dpsi |>
#   filter(lsv_id == "WBGene00000006:s:2573999-2574267") |>
#   View()
# 
# 
# 
# dpsi |>
#   mutate(is_ds = (p20 >.50 & p05 < .05)) |>
#   group_by(neurA, neurB) |>
#   summarize(nb_times_ds = sum(is_ds),
#             n = length(is_ds))




# Dominant spliceform ----


# by lsv

dominant_jct_pattern <- function(jct_id, threshold = .7){
  
  tab <- sort(table(jct_id), decreasing = TRUE)
  
  stopifnot(sum(tab) == length(jct_id))
  prop_tab <- tab/sum(tab)
  
  if(prop_tab[[1]] >= threshold){
    return("single dominant form")
    
  } else if(sum(prop_tab[1:2]) >= threshold){
    return("two main forms")
    
  } else if(sum(prop_tab[1:3]) >= threshold){
    return("three main forms")
    
  } else{
    return("complex")
  }
  
}




dominant_jct <- psi |>
  filter(mean_psi_per_lsv_junction > .3) |>
  summarize(jct_pattern = dominant_jct_pattern(junction_id) |>
              factor(levels = c("single dominant form", "two main forms", "three main forms", "complex")),
            gene_id = gene_id[[1]],
            .by = c(lsv_id))





dominant_jct$jct_pattern |> table()

gg_dominant_jct <- dominant_jct |>
  mutate(jct_pattern = if_else(jct_pattern == "three main forms", "complex", jct_pattern) |>
           factor(levels = c("single dominant form", "two main forms", "complex"))) |>
  mutate(jct_pattern = recode_factor(
    jct_pattern,
    `single dominant form` = "single\ndominant\nform",
    `two main forms` = "two\nmain\nforms"
  )) |>
  ggplot() +
  theme_classic() +
  geom_bar(aes(x = jct_pattern)) +
  # scale_y_log10() +
  ylab("Number of LSVs") +
  xlab("Spliceform expression pattern")

# ggsave("dominant_spliceform_by_lsv.pdf", path = export_dir,
#        width = 9, height = 12, units = "cm")


# by gene

dominant_gene_pattern <- function(jct_pat){
  
  jct_pat <- as.integer(jct_pat)
  
  nb_tot <- length(jct_pat)
  nb_single <- sum(jct_pat == 1L)
  nb_two <- sum(jct_pat == 2L)
  nb_three <- sum(jct_pat == 3L)
  nb_more <- sum(jct_pat > 3L)
  
  if(nb_single == nb_tot){
    return("single dominant form")
    
  } else if(nb_two == 1 && (nb_two + nb_single) == nb_tot){
    return("two main forms")
  } else if(nb_three == 1 && (nb_two + nb_single + nb_three) == nb_tot){
    return("three main forms")
  } else{
    return("complex")
  }
}


dominant_jct_gene <- dominant_jct |>
  summarize(gene_pattern = dominant_gene_pattern(jct_pattern) |>
              factor(levels = c("single dominant form", "two main forms", "three main forms", "complex")),
            .by = gene_id)



dominant_jct_gene$gene_pattern |> table()

gg_dominant_jct_gene <- dominant_jct_gene |>
  mutate(gene_pattern = if_else(gene_pattern == "three main forms", "complex", gene_pattern) |>
           factor(levels = c("single dominant form", "two main forms", "complex"))) |>
  mutate(gene_pattern = recode_factor(
    gene_pattern,
    `single dominant form` = "single\ndominant\nform",
    `two main forms` = "two\nmain\nforms"
  )) |>
  ggplot() +
  theme_classic() +
  geom_bar(aes(x = gene_pattern)) +
  # scale_y_log10() +
  ylab("Number of genes") +
  xlab("Spliceform expression pattern")

# ggsave("dominant_spliceform_by_gene.pdf", path = export_dir,
#        width = 8, height = 12, units = "cm")

patchwork::wrap_plots(gg_dominant_jct, gg_dominant_jct_gene)

ggsave("dominant_spliceforms.pdf", path = export_dir,
       width = 16, height = 10, units = "cm")








# DE vs DS ----

#~ neuron level ----



dsgs <- nb_signif_genes_by_test |>
  filter(neurA %in% neurs_with_deg_data,
         neurB %in% neurs_with_deg_data) |>
  rowwise() |>
  mutate(pair = c_across(starts_with("neur")) |> sort() |> paste0(collapse = "-"))


de_ds <- full_join(degs, dsgs, by = "pair")

ggplot(de_ds, aes(x = total_integrated_DEGs, y = nb_DS_genes)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm") +
  ggrepel::geom_text_repel(aes(label = pair),
                           max.overlaps = 5) +
  scale_x_log10() + scale_y_log10() +
  xlab("Number of DE genes (log)") +
  ylab("Number of DAS genes (log)")
# ggsave("DS_vs_DE.pdf", path = export_dir,
#        width = 15, height = 13, units = "cm")

mod <- lm(log(nb_DS_genes)~log(total_integrated_DEGs), data = de_ds)
mod <- lm(nb_DS_genes~total_integrated_DEGs, data = de_ds)

plot(fitted.values(mod), residuals(mod))
summary(mod)
qqnorm(residuals(mod)); qqline(residuals(mod))




#~ Top DEGs vs top DSGs ----


top_degs <- genes_degs |>
  filter(p.hmp < .05) |>
  group_by(pair) |>
  slice_max(order_by = abs(Average_logFC),
            n = 100) |>
  summarize(top_DEGs = list(as.character(gene_id)))


# check on app:
# top_degs$top_DEGs[[150]] |> clipr::write_clip()
# top_degs$pair[[150]]

top_dsgs <- dpsi |>
  filter(neurA %in% neurs_with_deg_data,
         neurB %in% neurs_with_deg_data) |>
  filter(p20 >.50 & p05 < .05) |>
  summarize(max_dpsi = max(dpsi),
            .by = c(neurA, neurB, gene_id)) |>
  group_by(neurA, neurB) |>
  slice_max(order_by = max_dpsi,
            n = 100) |>
  summarize(top_DSGs = list(as.character(gene_id)),
            .groups = 'drop') |>
  rowwise() |>
  mutate(pair = c_across(starts_with("neur")) |>
           sort() |>
           paste0(collapse = "-")) |>
  select(-neurA, -neurB)

sapply(top_degs$top_DEGs, length) |> hist()
sapply(top_dsgs$top_DSGs, length) |> hist()
table(lengths(top_dsgs$top_DSGs) == 100)

tops <- left_join(top_degs, top_dsgs,
          by = "pair") |>
  filter(lengths(top_DEGs) == 100,
         lengths(top_DSGs) == 100) |>
  rowwise() |>
  mutate(overlap = length(intersect(top_DEGs, top_DSGs)),
         nb_degs = length(top_DEGs),
         nb_dsgs = length(top_DSGs))

table(tops$overlap)


empty_eulerr_degs_dsgs <- eulerr::euler(list(`Top 100 DE genes` = 1:5,
                   `Top 100 DS genes` = 5:9)) |>
  plot()

# ggsave("empty_eulerr_degs_dsgs.pdf", plot = empty_eulerr_degs_dsgs,
#        path = export_dir,
#        width = 9, height = 6, units = "cm")

ggplot(tops) +
  theme_classic() +
  geom_bar(aes(x = overlap)) +
  scale_x_continuous(breaks = 0:7) +
  xlab("Number of overlapping genes") +
  ylab("Number of neuron pairs")

# ggsave("hist_degs_dsgs.pdf", path = export_dir,
#        width = 9, height = 6, units = "cm")



#~ Only RBPs ----


rbps <- wormDatasets::list_rbp_tamburino2013 |>
  filter(gene_id %in% rownames(cengenDataSC::cengen_proportion_bulk)) |>
  pull(gene_id)

genes_degs <- readRDS("data/genes/integrated_significant_genes_pairwise_bsn9_012722.rds") |>
  map(as_tibble, rownames = "gene_id") %>%
  map2_dfr(., names(.), ~add_column(.x, pair = .y))

nb_rbps_deg <- genes_degs |>
  filter(gene_id %in% rbps) |>
  group_by(pair) |>
  summarize(nb_degs = n()) |>
  mutate(pair = str_remove_all(pair, "[()]"))

stopifnot(sum(nb_rbps_deg$pair %in% dsgs$pair) == sum(dsgs$pair %in% nb_rbps_deg$pair))
stopifnot(sum(nb_rbps_deg$pair %in% dsgs$pair) == nrow(nb_rbps_deg))
stopifnot(sum(nb_rbps_deg$pair %in% dsgs$pair) == nrow(dsgs))

de_ds_rbps <- full_join(nb_rbps_deg, dsgs, by = "pair")


ggplot(de_ds_rbps, aes(x = nb_degs, y = nb_DS_genes)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm") +
  # ggrepel::geom_text_repel(aes(label = pair)) +
  scale_x_log10() + scale_y_log10() +
  xlab("Number of DE RBPs (log)") +
  ylab("Number of DAS genes (log)")
# ggsave("DS_vs_DE_rbps.pdf", path = export_dir,
#        width = 4, height = 4, units = "in")

mod <- lm(log(nb_DS_genes)~log(nb_degs), data = de_ds_rbps)

plot(fitted.values(mod), residuals(mod))
summary(mod)
qqnorm(residuals(mod));qqline(residuals(mod))





# Subsample neurons ----

# How many neurons do we need to sequence to detect that many DS genes?
# Subsamples



#~ single replicate ----
nb_ds_genes_sub <- function(n){
  neurset <- sample(all_neurs_sequenced, n)
  
  dpsi |>
    filter(neurA %in% neurset,
           neurB %in% neurset) |>
    filter(p20 >.50 & p05 < .05) |>
    pull(gene_id) |> unique() |>
    length()
}

# ==> use version with replicates below <==

# sub_nb_genes_ds <- tibble(nb_neurs = 1:length(all_neurs_sequenced),
#                           nb_ds_genes = map_int(nb_neurs, nb_ds_genes_sub))
# 
# 
# 
# # Fit Michaelis-Menten model
# mod_micment <- nls(nb_ds_genes ~ Gmax * nb_neurs/(b + nb_neurs),
#                    data = sub_nb_genes_ds,
#                    start = list(Gmax = 1000, b = 10))
# summary(mod_micment)
# plot(fitted.values(mod_micment), residuals(mod_micment))
# qqnorm(residuals(mod_micment))
# qqline(residuals(mod_micment))
# 
# 
# sub_nb_genes_ds <- cbind(sub_nb_genes_ds, fit=predict(mod_micment))
# 
# ggplot(sub_nb_genes_ds) +
#   theme_classic() +
#   geom_line(aes(x=nb_neurs, y=fit), color = "red3", linetype = "dashed") +
#   geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
#   geom_hline(aes(yintercept = coef(mod_micment)[["Gmax"]]), color = "gray50", linetype = "dotted") +
#   xlab("Number of neurons sampled") +
#   ylab("Number of genes detected as DS")
# # ggsave("nb_genes_ds.pdf", path = export_dir,
# #        width = 5, height = 5/1.375)





#~ 10 replicates per subsample ----

#~~ DAS ----
sub_nb_genes_ds_reps <- expand_grid(nb_neurs = 1:length(all_neurs_sequenced),
                                    rep = 1:10) |>
  mutate(nb_ds_genes = map_int(nb_neurs, nb_ds_genes_sub,
                               .progress = TRUE))



#~~ expr ----

nb_expr_genes_sub <- function(n){
  neurset <- sample(all_neurs_sequenced, n)
  
  gene_expression[,neurset] |> matrixStats::rowAnys() |> sum()
}

sub_nb_genes_expr_reps <- expand_grid(nb_neurs = 2:length(all_neurs_sequenced),
                                      rep = 1:10) |>
  mutate(nb_expr_genes = map_int(nb_neurs, nb_expr_genes_sub,
                                 .progress = TRUE))


ggplot(sub_nb_genes_ds_reps) +
  theme_classic() +
  geom_boxplot(aes(x = nb_neurs, y = nb_ds_genes, group = nb_neurs))+
  # geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  xlab("Number of neurons") +
  ylab("Number of genes detected as DAS")

ggplot(sub_nb_genes_expr_reps) +
  theme_classic() +
  geom_boxplot(aes(x = nb_neurs, y = nb_expr_genes, group = nb_neurs))+
  # geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  xlab("Number of neurons") +
  ylab("Number of genes detected as DEG")






#~ Fit on log-linear ----

#~~ DAS ----
mod_lm <- lm(nb_ds_genes ~ log(nb_neurs),
             data = sub_nb_genes_ds_reps)
plot(fitted.values(mod_lm), residuals(mod_lm))
summary(mod_lm)



sub_nb_genes_ds_reps$fit_lm <- predict(mod_lm)

ggplot(sub_nb_genes_ds_reps) +
  theme_classic() +
  geom_boxplot(aes(x = nb_neurs, y = nb_ds_genes, group = nb_neurs)) +
  geom_line(aes(x=nb_neurs, y=fit_lm), color = "red3", linetype = "dashed") +
  # geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  xlab("Number of neurons sampled") +
  ylab("Number of genes detected as DAS") +
  scale_x_continuous(limits = c(0,NA)) +
  scale_y_continuous(limits = c(0,NA), labels = scales::label_comma())
# ggsave("nb_genes_das.pdf", path = export_dir,
#        width = 16, height = 8, units = "cm")

predict(mod_lm, newdata = data.frame(nb_neurs = 119))



#~~ expr ----

gmod_lm <- lm(nb_expr_genes ~ log(nb_neurs),
              data = sub_nb_genes_expr_reps)
plot(fitted.values(gmod_lm), residuals(gmod_lm))
summary(gmod_lm)
predict(gmod_lm, newdata = data.frame(nb_neurs = 119))
sub_nb_genes_expr_reps$gmod_lm <- predict(gmod_lm, newdata = data.frame(nb_neurs=sub_nb_genes_expr_reps$nb_neurs))

ggplot(sub_nb_genes_expr_reps) +
  theme_classic() +
  geom_boxplot(aes(x = nb_neurs, y = nb_expr_genes, group = nb_neurs)) +
  geom_line(aes(x=nb_neurs, y=gmod_lm), color = "red3", linetype = "dashed") +
  # geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  # geom_hline(aes(yintercept = coef(mod_micment)[["Gmax"]]), color = "gray50", linetype = "dotted") +
  xlab("Number of neurons sampled") +
  ylab("Number of genes detected as expressed")



# Prop
predict(mod_lm, newdata = data.frame(nb_neurs = 119))/predict(gmod_lm, newdata = data.frame(nb_neurs = 119))




#~ micment ----

#~~ DAS ----


dmod_micment <- nls(nb_ds_genes ~ Gmax * nb_neurs/(b + nb_neurs),
                   data = sub_nb_genes_ds_reps,
                   start = list(Gmax = 3000, b = 10))
summary(dmod_micment)
plot(fitted.values(dmod_micment), residuals(dmod_micment))
predict(dmod_micment, newdata = data.frame(nb_neurs = 119))


sub_nb_genes_ds_reps <- cbind(sub_nb_genes_ds_reps, fit=predict(dmod_micment))

ggplot(sub_nb_genes_ds_reps) +
  theme_classic() +
  geom_boxplot(aes(x = nb_neurs, y = nb_ds_genes, group = nb_neurs)) +
  geom_line(aes(x=nb_neurs, y=fit), color = "red3", linetype = "dashed") +
  # geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  geom_hline(aes(yintercept = coef(dmod_micment)[["Gmax"]]), color = "gray50", linetype = "dotted") +
  xlab("Number of neurons sampled") +
  ylab("Number of genes detected as DS")


#~~ expr ----



gmod_micment <- nls(nb_expr_genes ~ Gmax * nb_neurs/(b + nb_neurs),
                    data = sub_nb_genes_expr_reps,
                    start = list(Gmax = 8000, b = 10))
summary(gmod_micment)
plot(fitted.values(gmod_micment), residuals(gmod_micment))
predict(gmod_micment, newdata = data.frame(nb_neurs = 119))

# Prop
predict(dmod_micment, newdata = data.frame(nb_neurs = 119))/predict(gmod_micment, newdata = data.frame(nb_neurs = 119))






# the rest below isn't used as of 2024-03

# ~~~~~~~~~~~~~~ ----


#~ By proportion of expressed genes ----
# Say a gene is expressed if threshold 2 from sc detects it in at least 2 of the sampled classes
# note: this is somewhat equivalent to the heatmap (but with subsampling), not very interesting

prop_ds_genes_sub <- function(n){
  neurset <- sample(all_neurs_sequenced, n)
  
  nb_neurs_where_gene_expr <- rowSums(gene_expr_int[,neurset])
  nb_genes_expr_in_neurset <- sum(nb_neurs_where_gene_expr > 1)
  
  
  nb_genes_ds_in_neurset <- dpsi |>
    filter(neurA %in% neurset,
           neurB %in% neurset) |>
    filter(p20 >.50 & p05 < .05) |>
    pull(gene_id) |> unique() |>
    length()
  
  nb_genes_ds_in_neurset/nb_genes_expr_in_neurset
}


gene_coexpr_by_neur_pair <- dpsi |>
  mutate(ds = (p20 >.50 & p05 < .05)) |>
  select(gene_id, ds, neurA, neurB) |>
  summarize(has_ds = any(ds),
            .by = c(neurA, neurB, gene_id)) |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurA = "neuron")) |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurB = "neuron")) |>
  mutate(coexpressed = expressed.x & expressed.y) |>
  filter(! is.na(coexpressed)) |>   # a number of annotated pseudogenes have splicing quantified but not expression
  select(neurA, neurB, gene_id, has_ds, coexpressed)


prop_ds_genes_sub <- function(n){
  neurset <- sample(all_neurs_sequenced, n)
  
  gene_coexpr_by_neur_pair |>
    filter(neurA %in% neurset,
           neurB %in% neurset,
           coexpressed) |>
    group_by(neurA,neurB) |>
    summarize(prop_ds_in_pair = mean(has_ds),
              .groups = "drop") |>
    summarize(mean_prop_ds = mean(prop_ds_in_pair)) |>
    pull(mean_prop_ds)
}

sub_prop_reps <- expand_grid(nb_neurs = 3:length(all_neurs_sequenced),
                        rep = 1:10) |>
  mutate(prop_ds_genes = map_dbl(nb_neurs, prop_ds_genes_sub))

plot(sub_prop_reps$nb_neurs, sub_prop_reps$prop_ds_genes)
fit <- lm(prop_ds_genes ~ nb_neurs, data = sub_prop_reps)
summary(fit)

plot(fitted.values(fit), residuals(fit))
qqnorm(residuals(fit)); qqline(residuals(fit))



ggplot(sub_prop_reps, aes(x=nb_neurs, y=prop_ds_genes)) +
  theme_classic() +
  # geom_line(aes(x=nb_neurs, y=fit), color = "red3", linetype = "dashed") +
  # geom_smooth(method = "lm", se = FALSE) +
  geom_point() +
  geom_abline(intercept = coef(fit)[[1]], slope = coef(fit)[[2]],
              linetype = "dashed", color = "red3") +
  xlab("Number of neurons sampled") +
  ylab("Mean proportion detected as DS") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1))

# ggsave("subsample_proportion_of_ds.pdf", path = export_dir,
#        width = 16, height = 8, units = "cm")




#~ proportion of DS not subsampling ----
# alternatively, subsampling may not be useful here, asking across neuron pairs




gene_coexpr_by_neur_pair <- dpsi |>
  mutate(ds = (p20 >.50 & p05 < .05)) |>
  select(gene_id, ds, neurA, neurB) |>
  summarize(has_ds = any(ds),
            .by = c(neurA, neurB, gene_id)) |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurA = "neuron")) |>
  left_join(gene_expr_tib,
            by = c("gene_id", neurB = "neuron")) |>
  mutate(coexpressed = expressed.x & expressed.y) |>
  filter(! is.na(coexpressed)) |>   # a number of annotated pseudogenes have splicing quantified but not expression
  select(neurA, neurB, gene_id, has_ds, coexpressed)


ds_among_coexpr_by_pair <- gene_coexpr_by_neur_pair |>
  filter(coexpressed) |>
  summarize(prop_ds_in_pair = mean(has_ds),
            nb_ds_in_pair = sum(has_ds),
            nb_coexpr_genes = n(),
            .by = c(neurA,neurB))

ds_among_coexpr_by_pair |>
  ggplot() + theme_classic() +
  geom_point(aes(x = nb_coexpr_genes, y = prop_ds_in_pair))
#> no obvious bias

ds_among_coexpr_by_pair |>
  ggplot() +
  theme_classic() +
  geom_histogram(aes(x = prop_ds_in_pair),
                 bins = 60, color = 'white') +
  geom_vline(xintercept = mean(ds_among_coexpr_by_pair$prop_ds_in_pair),
              linetype = "dashed", color = "red3") +
  xlab("Proportion DS among genes coexpressed in pair") +
  ylab("Number of pairs") +
  scale_x_continuous(breaks = seq(5, 25, by = 5)/100,
                     labels = scales::label_percent(accuracy = 1))


ds_among_coexpr_by_pair |> filter(neurA %in% c("OLL","RME"),
                                  neurB %in% c("OLL","RME"))

                     
coexpr_ds








#~ older version, using sc ----

# expr_sc <- cengenDataSC::cengen_sc_3_bulk > 0


prop_ds_genes_sub <- function(n){
  neurset <- sample(all_neurs_sequenced, n)
  
  nb_neurs_where_gene_expr <- rowSums(gene_expr_int[,neurset])
  nb_genes_expr_in_neurset <- sum(nb_neurs_where_gene_expr > 1)
  
  
  nb_genes_ds_in_neurset <- dpsi |>
    filter(neurA %in% neurset,
           neurB %in% neurset) |>
    filter(p20 >.50 & p05 < .05) |>
    pull(gene_id) |> unique() |>
    length()
  
  nb_genes_ds_in_neurset/nb_genes_expr_in_neurset
}


sub_prop_genes_ds <- tibble(nb_neurs = 3:length(all_neurs_sequenced),
                            prop_ds_genes = map_dbl(nb_neurs, prop_ds_genes_sub))


# Fit Michaelis-Menten model
mod_micment <- nls(prop_ds_genes ~ Gmax * nb_neurs/(b + nb_neurs),
                   data = sub_prop_genes_ds,
                   start = list(Gmax = 3000, b = 10))

sub_prop_genes_ds <- cbind(sub_prop_genes_ds, fit=predict(mod_micment))

ggplot(sub_prop_genes_ds) +
  theme_classic() +
  geom_line(aes(x=nb_neurs, y=fit), color = "red3", linetype = "dashed") +
  geom_point(aes(x=nb_neurs, y=prop_ds_genes)) +
  geom_hline(aes(yintercept = coef(mod_micment)[["Gmax"]]), color = "gray50", linetype = "dotted") +
  scale_y_continuous(limits = c(0,.3), labels = \(xx) scales::percent(xx, accuracy = 1)) +
  xlab("Number of neurons sampled") +
  ylab("Proportion of genes detected as DS")
# ggsave("prop_genes_ds.pdf", path = export_dir,
#        width = 5, height = 5/1.375)



# event DS broadness ----

tau <- rowSums(1 - as.matrix(gene_expr_int))/(ncol(gene_expr_int) - 1)
tau[tau > 1] <- NA


tibble(gene = genes_in_our_neur_sample,
       tau = tau[genes_in_our_neur_sample]) |>
  mutate(is_ds = gene %in% ds_genes) |>
  ggplot(aes(x = is_ds, y = tau, fill = is_ds)) +
  theme_classic() +
  geom_boxplot() +
  ggbeeswarm::geom_quasirandom()



neurs_integrated_noD <- neurs_integrated |> setdiff(c("VD","DD"))

degs <- read.delim("data/2021-11-30_alec_integration/Total_integrated_DEGS_pairwise_113021.tsv") |>
  as_tibble() |>
  rowwise() |>
  mutate(pair = c_across(starts_with("cell_")) |> sort() |> paste0(collapse = "-"))

dsgs <- nb_signif_genes_by_test |>
  ungroup() |>
  filter(neurA %in% neurs_integrated_noD,
         neurB %in% neurs_integrated_noD) |>
  rowwise() |>
  mutate(pair = c_across(starts_with("neur")) |> sort() |> paste0(collapse = "-"))



library(edgeR)
x <- DGEList(read.delim("data/2021-11-30_alec_integration/average_integration_GeTMM_113021.tsv"))
x$samples$group <- str_split_i(colnames(x), "r", 1) |> as.factor()

cpm <- cpm(x, log = TRUE)

pca_genes <- prcomp(t(cpm))

autoplot(pca_genes, label = TRUE, shape = FALSE) + theme_classic()

keep.exprs <- filterByExpr(x, min.total.count = 30, min.prop = .9, min.count = 20)
table(keep.exprs)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

plotMDS(cpm)

Glimma::glMDSPlot(x, labels = x$samples$group, groups = x$samples$group)

design <- model.matrix(~0+group, data = x$samples)

v <- voom(x, design, plot = TRUE)

vfit <- lmFit(v, design)
plotSA(vfit)

# pairwise tests
all_groups <- colnames(design)

expand_grid(neurA = all_groups,
            neurB = all_groups)
contr.matrix <- combn(all_groups, 2, simplify = FALSE) |>
  map_chr(paste0, collapse = "-")|>
  makeContrasts(contrasts = _, levels = all_groups)

vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")




# one neur vs all others
all_groups <- colnames(design)

all_contrasts <- map_chr(all_groups,
        ~paste0( .x, " - (", paste(all_groups |> setdiff(.x), collapse = "+"), ")/(",length(all_groups)-1,")"))
  
contr.matrix <- makeContrasts(contrasts = all_contrasts, levels = all_groups)

vfit <- contrasts.fit(vfit, contrasts = contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

decideTests(efit) |> summary() |> dim()

tfit <- treat(vfit, lfc = 1)
dt <- decideTests(tfit)
summary(dt)
colnames(dt) <- str_match(colnames(dt), "^group([[:alnum:]]+) - ")[,2]




#~ Check degs vs dsgs ----


genes_de <- abs(as.matrix(decideTests(efit))) |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  as_tibble() |>
  pivot_longer(cols = -gene_id,
               names_to = "test",
               values_to = "deg") |>
  filter(deg == 1) |>
  select(-deg) |>
  distinct() |>
  mutate(test = str_remove_all(test, "group"))

genes_ds <- dpsi |>
  filter(p20 >.50 & p05 < .05 & abs(dpsi) > .5) |>
  mutate(test = paste0(neurA, "-", neurB)) |>
  select(gene_id, test) |>
  distinct()

#~~ per neuron pair ----
# reproducing above


#' Ex: order_test("RIA-ADL")
#' returns it alphabetically sorted: ADL-RIA
order_test <- function(neur_pair_str){
  neurons <- str_split_1(neur_pair_str, "-")
  paste0(sort(neurons), collapse = "-")
}

degs2 <- genes_de |>
  summarize(nb_degs = n(),
            .by = "test") |>
  rowwise() |>
  mutate(test = order_test(test)) |>
  ungroup()

dsgs2 <- genes_ds |>
  summarize(nb_dsgs = n(),
            .by = "test") |>
  rowwise() |>
  mutate(test = order_test(test)) |>
  ungroup()




de_ds <- full_join(degs2, dsgs2, by = "test")

ggplot(de_ds, aes(x = nb_degs, y = nb_dsgs)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm") +
  ggrepel::geom_text_repel(aes(label = test)) +
  # scale_x_log10() +
  # scale_y_log10() +
  xlab("Number of DE genes (log)") +
  ylab("Number of DAS genes (log)")



#~~ gene-level ----

degs3 <- genes_de |>
  summarize(nb_pairs_de = n(),
            .by = gene_id)

dsgs3 <- genes_ds |>
  summarize(nb_pairs_ds = n(),
            .by = gene_id)

pair_de_ds <- full_join(degs3, dsgs3, by = "gene_id")

ggplot(pair_de_ds, aes(x = nb_pairs_de, y = nb_pairs_ds)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm") +
  # ggrepel::geom_text_repel(aes(label = gene_id)) +
  # scale_x_log10() +
  # scale_y_log10() +
  xlab("Number of pairs where DE") +
  ylab("Number of pairs where DS")


pair_de_ds$nb_pairs_de |> hist(breaks = 50)
pair_de_ds$nb_pairs_ds |> hist(breaks = 50)

pair_de_ds |>
  filter(!is.na(nb_pairs_ds) & !is.na(nb_pairs_de)) |>
  ggplot(aes(x = nb_pairs_de, y = nb_pairs_ds)) +
  theme_classic() +
  geom_point() +
  geom_smooth(method = "lm") +
  # ggrepel::geom_text_repel(aes(label = gene_id)) +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Number of pairs where DE") +
  ylab("Number of pairs where DS")

pair_de_ds |>
  filter(nb_pairs_de < 300,
         nb_pairs_ds > 10)

#~ Dominant event per sample ----
psi




# ~~~~~~~~~~~~~~~~~~~~~~~~ ----


# Long exons with DS ----

# Consider only genes witrh DS
genes_ds <- dpsi |>
  filter(neurA %in% all_neurs_sequenced,
         neurB %in% all_neurs_sequenced) |>
  filter(p20 >.50 & p05 < .05) |>
  pull(gene_id) |> unique()


# of those, how many have an exon with > 500 bp

exons <- wb_load_exon_coords(281)

big_ex_ds <- exons |>
  filter(gene_id %in% genes_ds) |>
  mutate(exon_length = end-start+1) |>
  filter(exon_length > 1000) |>
  mutate(gene_name = i2s(gene_id, gids))




#~ old Heatmap normd by nb_genes_coexpressedping genes ----


hm_mat <- nb_signif_genes_by_test |>
  rename(neurA = neurB,
         neurB = neurA) |>
  bind_rows(nb_signif_genes_by_test) |>
  pivot_wider(names_from = neurB,
              values_from = nb_DS_genes) |>
  column_to_rownames("neurA") |> 
  as.matrix()

hm_mat <- hm_mat[sort(rownames(hm_mat)), sort(colnames(hm_mat))]

heatmap_annot <- neuron_properties |>
  column_to_rownames("Neuron_type")

hc <- hclust(dist(hm_mat, method = "canberra"), method = "complete")

nb_genes_coexpressed <- t(cengenDataSC::cengen_sc_3_bulk[,colnames(hm_mat)] >0) %*% (cengenDataSC::cengen_sc_3_bulk[,colnames(hm_mat)]>0)
all.equal(colnames(hm_mat), colnames(nb_genes_coexpressed))
all.equal(rownames(hm_mat), rownames(nb_genes_coexpressed))



# To set the order (roughly)
weights <- set_names(c("OLL","IL1","VD","DD","PHA","AWC","AVH","DA","PVC","I5","AVM","VC","AIN","OLQ",
                       "RIM","RMD","ASER","ASI","AFD","BAG","ASEL","ASG","AWA","NSM","CAN","ADL",
                       "ASK","IL2","RIA","RIC","AVK","SMD","AVA","AVE","AVG","AWB","PVD","PVM","RIS",
                       "VB","AIY"),
                     exp(1:41)) |>
  sort() |>
  names() |>
  as.numeric()


hm_callback <- function(hc, ...){
  as.hclust(reorder(as.dendrogram(hc), wts = weights))
}

hc2 <- hm_callback(hc)

pheatmap::pheatmap((hm_mat/nb_genes_coexpressed),
                   scale = "none",
                   cluster_rows = hc2,
                   cluster_cols = hc2,
                   cutree_rows = 2,
                   cutree_cols = 2,
                   main = "Normalized by number of nb_genes_coexpressedping genes",
                   # filename = file.path(export_dir, "heatmap_ds_bynb_genes_coexpressed.pdf"),
                   # width = 10,
                   # height = 9
)

table(cutree(hc2, k=2))




#~ other Heatmaps ----



pheatmap::pheatmap(hm_mat, scale = "none",
                   clustering_distance_rows = "canberra",
                   clustering_distance_cols = "canberra",
                   clustering_method = "complete",
                   annotation_row = heatmap_annot,
                   annotation_col = heatmap_annot,
                   # filename = file.path(export_dir, "211111_heatmap_ds.pdf"),
                   # width = 10,
                   # height = 9
                   )






pheatmap::pheatmap((hm_mat/nb_genes_coexpressed),
                   scale = "none",
                   clustering_distance_rows = "canberra",
                   clustering_distance_cols = "canberra",
                   clustering_method = "complete",
                   clustering_callback = hm_callback,
                   main = "Normalized by number of nb_genes_coexpressedping genes",
                   # filename = file.path(export_dir, "heatmap_ds_bynb_genes_coexpressed_no_annot.pdf"),
                   # width = 10,
                   # height = 9
)


#~ Heatmap normalized by neurons ----

hc <- hclust(dist(hm_mat, method = "canberra"), method = "complete")
hm_mat_by_row <- hm_mat/colSums(hm_mat, na.rm = TRUE)
hm_mat_by_col <- t(hm_mat/colSums(hm_mat, na.rm = TRUE))

pheatmap::pheatmap(hm_mat_by_row, scale = "none",
                   cluster_cols = hc,
                   cluster_rows = hc,
                   annotation_row = heatmap_annot,
                   annotation_col = heatmap_annot,
                   main = "Normalized by row",
                   # filename = file.path(export_dir, "211111_heatmap_ds_byrow.pdf"),
                   # width = 10,
                   # height = 9
)
# 
# pheatmap::pheatmap(hm_mat_by_col, scale = "none",
#                    cluster_cols = hc,
#                    cluster_rows = hc,
#                    annotation_row = heatmap_annot,
#                    annotation_col = heatmap_annot)
# 
# pheatmap::pheatmap(t(scale(t(hm_mat), center = FALSE)), scale = "none",
#                    cluster_cols = hc,
#                    cluster_rows = hc,
#                    annotation_row = heatmap_annot,
#                    annotation_col = heatmap_annot)
# 
# pheatmap::pheatmap(hm_mat, scale = "row",
#                    cluster_cols = hc,
#                    cluster_rows = hc,
#                    annotation_row = heatmap_annot,
#                    annotation_col = heatmap_annot)


#~ Normalization comparison ----

#~~ bargraph ----

tot_genes_per_neur <- nb_signif_genes_by_test |>
  pivot_longer(cols = c("neurA","neurB")) |>
  group_by(value) |>
  summarize(tot_neur = sum(nb_DS_genes),
            .groups = "drop")

nb_signif_genes_by_test |>
  ungroup() |>
  left_join(tot_genes_per_neur,
            by = c(neurA = "value")) |>
  left_join(tot_genes_per_neur,
            by = c(neurB = "value"),
            suffix = c("A","B")) |>
  mutate(pair = paste0(neurA, "-", neurB),
         nb_by_A = nb_DS_genes/tot_neurA,
         nb_by_B = nb_DS_genes/tot_neurB) |>
  select(pair, nb_by_A, nb_by_B) |>
  pivot_longer(-pair, names_to = "norm_by", values_to = "norm_nb_DS_genes") |>
  ggplot() +
  geom_col(aes(x = pair, y = norm_nb_DS_genes, fill = norm_by))

#~~ scatterplot ----

nb_signif_genes_by_test |>
  ungroup() |>
  left_join(tot_genes_per_neur,
            by = c(neurA = "value")) |>
  left_join(tot_genes_per_neur,
            by = c(neurB = "value"),
            suffix = c("A","B")) |>
  mutate(pair = paste0(neurA, "-", neurB),
         nb_by_A = nb_DS_genes/tot_neurA,
         nb_by_B = nb_DS_genes/tot_neurB) |>
  select(pair, nb_by_A, nb_by_B) |>
  ggplot(aes(x = nb_by_A, y = nb_by_B)) +
  geom_point() +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Number of DS genes for neuron pair normalized by first neuron") +
  ylab("Number of DS genes for neuron pair normalized by second neuron") +
  ggrepel::geom_text_repel(aes(label = pair))

# ggsave("nb_ds_genes_norm_by_both_neurs.pdf", path = export_dir)
# ggsave("nb_ds_genes_norm_by_both_neurs_labels.pdf", path = export_dir)

# How many tests this DS gene appears in?
nb_tests_per_ds_gene <- dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  select(gene_id, neurA, neurB) |>
  distinct() |>
  group_by(gene_id) |>
  summarize(nb_tests = n())

hist(nb_tests_per_ds_gene$nb_tests, breaks = 100)

nb_tests_per_ds_gene |>
    filter(nb_tests > 600) |>
  mutate(gene_name = i2s(gene_id, gids))

dpsi |>
  filter(p20 >.50 & p05 < .05,
         neurA == "AIY" | neurB == "AIY") |>
  select(gene_id, neurA, neurB) |>
  distinct() |>
  group_by(gene_id) |>
  summarize(nb_tests = n()) |>
  pull(nb_tests) |>
  hist(breaks = 100)


xx <- dpsi |>
  filter(neurA == "AIY" | neurB == "AIY") |>
  mutate(has_ds = p20 > .5 & p05 < .05) |>
  group_by(neurA, neurB, gene_id) |>
  summarize(nb_ds = sum(has_ds),
            .groups = "drop") |>
  mutate(neur = if_else(neurA == "AIY", neurB, neurA)) |>
  select(-neurA, -neurB) |>
  group_by(gene_id) |>
  summarize(prop_tests_ds = mean(nb_ds > 0))

hist(xx$prop_tests_ds, breaks = 50)



xx <- dpsi |>
  mutate(has_ds = p20 > .5 & p05 < .05) |>
  group_by(neurA, neurB, gene_id) |>
  summarize(nb_ds = sum(has_ds),
            .groups = "drop") |>
  group_by(gene_id) |>
  summarize(prop_tests_ds = sum(nb_ds > 0)/nrow(select(dpsi, neurA, neurB) |> distinct()))

hist(xx$prop_tests_ds, breaks = 50)

xx |>
  filter(prop_tests_ds > .8) |>
  mutate(gene_name = i2s(gene_id, gids))


# Nb of LSVs with DS ----


dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  pull(lsv_id) |> unique() |>
  length()
#> 4611 / 10,232

#~ By neur classes ----
signif_lsvs <- dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  select(lsv_id, neurA, neurB) |>
  distinct() |>
  group_by(neurA, neurB) |>
  summarize(nb_ds_lsvs = n(),
            .groups = "drop")

signif_lsvs |>
  rename(neurA = neurB,
         neurB = neurA) |>
  bind_rows(signif_lsvs) |>
  ggplot() +
  theme_classic() +
  geom_tile(aes(x = neurA, y=neurB, fill = nb_ds_lsvs))



#~ Heatmap ----
hm_mat_lsvs <- signif_lsvs |>
  rename(neurA = neurB,
         neurB = neurA) |>
  bind_rows(signif_lsvs) |>
  pivot_wider(names_from = neurB,
              values_from = nb_ds_lsvs) |>
  column_to_rownames("neurA") |> 
  as.matrix()

hm_mat_lsvs <- hm_mat_lsvs[sort(rownames(hm_mat_lsvs)), sort(colnames(hm_mat_lsvs))]


pheatmap::pheatmap(hm_mat_lsvs, scale = "none",
                   clustering_distance_rows = "canberra",
                   clustering_distance_cols = "canberra",
                   clustering_method = "complete",
                   annotation_row = heatmap_annot,
                   annotation_col = heatmap_annot,
                   # filename = file.path(export_dir, "211111_heatmap_lsvs.pdf"),
                   # width = 10,
                   # height = 9
)



#~ Heatmap normd by nb_genes_coexpressedping genes ----


nb_genes_coexpressed <- t(cengenDataSC::cengen_sc_3_bulk[,colnames(hm_mat_lsvs)]) %*% cengenDataSC::cengen_sc_3_bulk[,colnames(hm_mat_lsvs)]
all.equal(colnames(hm_mat_lsvs), colnames(nb_genes_coexpressed))
all.equal(rownames(hm_mat_lsvs), rownames(nb_genes_coexpressed))



# To set the order (roughly)
weights <- set_names(c("OLL","OLQ","PHA","AWC","IL1","M4","VD","DD","AVH","DA","PVC","I5","AVM","VC","AIN","RIM","RMD","ASER","ASI","AFD","BAG","ADF","ASK","IL2","RIA","RIC","AVK","SMD","ASEL","ASG","AWA","NSM","CAN","ADL","AVA","AVE","AVG","AWB","PVD","PVM","RIS","VB","AIY"),
                     exp(1:43)) |>
  sort() |>
  names() |>
  as.numeric()


hm_callback <- function(hc, ...){
  as.hclust(reorder(as.dendrogram(hc), wts = weights))
}


pheatmap::pheatmap((hm_mat_lsvs/nb_genes_coexpressed),
                   scale = "none",
                   clustering_distance_rows = "canberra",
                   clustering_distance_cols = "canberra",
                   clustering_method = "complete",
                   clustering_callback = hm_callback,
                   annotation_row = heatmap_annot,
                   annotation_col = heatmap_annot,
                   main = "Normalized by number of nb_genes_coexpressedping LSVs",
                   # filename = file.path(export_dir, "heatmap_lsvs_bynb_genes_coexpressed.pdf"),
                   # width = 10,
                   # height = 9
)




#~ Heatmap normalized by neurons ----

hc <- hclust(dist(hm_mat_lsvs, method = "canberra"), method = "complete")
hm_mat_by_row <- hm_mat_lsvs/colSums(hm_mat_lsvs, na.rm = TRUE)
hm_mat_by_col <- t(hm_mat_lsvs/colSums(hm_mat_lsvs, na.rm = TRUE))

pheatmap::pheatmap(hm_mat_by_row, scale = "none",
                   cluster_cols = hc,
                   cluster_rows = hc,
                   annotation_row = heatmap_annot,
                   annotation_col = heatmap_annot,
                   main = "Normalized by row",
                   # filename = file.path(export_dir, "211111_heatmap_lsvs_byrow.pdf"),
                   # width = 10,
                   # height = 9
)
# 
# pheatmap::pheatmap(hm_mat_by_col, scale = "none",
#                    cluster_cols = hc,
#                    cluster_rows = hc,
#                    annotation_row = heatmap_annot,
#                    annotation_col = heatmap_annot)
# 
# pheatmap::pheatmap(t(scale(t(hm_mat), center = FALSE)), scale = "none",
#                    cluster_cols = hc,
#                    cluster_rows = hc,
#                    annotation_row = heatmap_annot,
#                    annotation_col = heatmap_annot)
# 
# pheatmap::pheatmap(hm_mat, scale = "row",
#                    cluster_cols = hc,
#                    cluster_rows = hc,
#                    annotation_row = heatmap_annot,
#                    annotation_col = heatmap_annot)



#~ Normalization comparison ----
#~~ bargraph ----

tot_lsvs_per_neur <- signif_lsvs |>
  pivot_longer(cols = c("neurA","neurB")) |>
  group_by(value) |>
  summarize(tot_neur = sum(nb_ds_lsvs),
            .groups = "drop")

signif_lsvs |>
  ungroup() |>
  left_join(tot_lsvs_per_neur,
            by = c(neurA = "value")) |>
  left_join(tot_lsvs_per_neur,
            by = c(neurB = "value"),
            suffix = c("A","B")) |>
  mutate(pair = paste0(neurA, "-", neurB),
         nb_by_A = nb_ds_lsvs/tot_neurA,
         nb_by_B = nb_ds_lsvs/tot_neurB) |>
  select(pair, nb_by_A, nb_by_B) |>
  pivot_longer(-pair, names_to = "norm_by", values_to = "norm_nb_DS_lsvs") |>
  ggplot() +
  geom_col(aes(x = pair, y = norm_nb_DS_lsvs, fill = norm_by))


#~~ scatterplot ----

signif_lsvs |>
  ungroup() |>
  left_join(tot_lsvs_per_neur,
            by = c(neurA = "value")) |>
  left_join(tot_lsvs_per_neur,
            by = c(neurB = "value"),
            suffix = c("A","B")) |>
  mutate(pair = paste0(neurA, "-", neurB),
         nb_by_A = nb_ds_lsvs/tot_neurA,
         nb_by_B = nb_ds_lsvs/tot_neurB) |>
  select(pair, nb_by_A, nb_by_B) |>
  ggplot(aes(x = nb_by_A, y = nb_by_B)) +
  geom_point() +
  scale_x_continuous(labels = scales::percent) +
  scale_y_continuous(labels = scales::percent) +
  xlab("Number of DS LSVs for neuron pair normalized by first neuron") +
  ylab("Number of DS LSVs for neuron pair normalized by second neuron") +
  ggrepel::geom_text_repel(aes(label = pair))

# ggsave("nb_ds_lsvs_norm_by_both_neurs.pdf", path = export_dir)
# ggsave("nb_ds_lsvs_norm_by_both_neurs_labels.pdf", path = export_dir)


# How many tests this DS gene appears in?
nb_tests_per_ds_gene <- dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  select(gene_id, neurA, neurB) |>
  distinct() |>
  group_by(gene_id) |>
  summarize(nb_tests = n())

hist(nb_tests_per_ds_gene$nb_tests, breaks = 100)

nb_tests_per_ds_gene |>
  filter(nb_tests > 600) |>
  mutate(gene_name = i2s(gene_id, gids))

dpsi |>
  filter(p20 >.50 & p05 < .05,
         neurA == "AIY" | neurB == "AIY") |>
  select(gene_id, neurA, neurB) |>
  distinct() |>
  group_by(gene_id) |>
  summarize(nb_tests = n()) |>
  pull(nb_tests) |>
  hist(breaks = 100)


xx <- dpsi |>
  filter(neurA == "AIY" | neurB == "AIY") |>
  mutate(has_ds = p20 > .5 & p05 < .05) |>
  group_by(neurA, neurB, gene_id) |>
  summarize(nb_ds = sum(has_ds),
            .groups = "drop") |>
  mutate(neur = if_else(neurA == "AIY", neurB, neurA)) |>
  select(-neurA, -neurB) |>
  group_by(gene_id) |>
  summarize(prop_tests_ds = mean(nb_ds > 0))

hist(xx$prop_tests_ds, breaks = 50)



xx <- dpsi |>
  mutate(has_ds = p20 > .5 & p05 < .05) |>
  group_by(neurA, neurB, gene_id) |>
  summarize(nb_ds = sum(has_ds),
            .groups = "drop") |>
  group_by(gene_id) |>
  summarize(prop_tests_ds = sum(nb_ds > 0)/nrow(select(dpsi, neurA, neurB) |> distinct()))

hist(xx$prop_tests_ds, breaks = 50)

xx |>
  filter(prop_tests_ds > .8) |>
  mutate(gene_name = i2s(gene_id, gids))





# ~~~~~~~~~~~  ----














