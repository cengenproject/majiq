# Load deltas


# Initializations ----
library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(281)


export_dir <- "presentations/2021-11_for_preprint_december/"

neuron_properties <- read_csv("data/neuron_properties.csv")


#~ Load ----
data_dir <- "data/2021-11-10_outs/deltapsi/"

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
                            na = "na") |>
                    add_column(neurA = str_split(.x,"[-\\.]")[[1]][1],
                               neurB = str_split(.x,"[-\\.]")[[1]][2]))

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

# saveRDS(dpsi, "intermediates/211110_dpsi.rds")







dpsi <- readRDS("intermediates/211110_dpsi.rds")




# Nb of genes with DS ----

dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  pull(gene_id) |> unique() |>
  length()

#~ By neur classes ----
signif_genes <- dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  select(gene_id, neurA, neurB) |>
  distinct() |>
  group_by(neurA, neurB) |>
  summarize(nb_DS_genes = n())

signif_genes |>
  rename(neurA = neurB,
         neurB = neurA) |>
  bind_rows(signif_genes) |>
  ggplot() +
  theme_classic() +
  geom_tile(aes(x = neurA, y=neurB, fill = nb_DS_genes))



#~ Heatmap ----
hm_mat <- signif_genes |>
  rename(neurA = neurB,
         neurB = neurA) |>
  bind_rows(signif_genes) |>
  pivot_wider(names_from = neurB,
              values_from = nb_DS_genes) |>
  column_to_rownames("neurA") |> 
  as.matrix()

hm_mat <- hm_mat[sort(rownames(hm_mat)), sort(colnames(hm_mat))]

heatmap_annot <- neuron_properties |>
  column_to_rownames("Neuron_type")


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

#~ Heatmap normd by overlapping genes ----

hc <- hclust(dist(hm_mat, method = "canberra"), method = "complete")

overlap <- t(cengenDataSC::cengen_sc_3_bulk[,colnames(hm_mat)] >0) %*% (cengenDataSC::cengen_sc_3_bulk[,colnames(hm_mat)]>0)
all.equal(colnames(hm_mat), colnames(overlap))
all.equal(rownames(hm_mat), rownames(overlap))



# To set the order (roughly)
weights <- set_names(c("OLL","OLQ","PHA","AWC","IL1","M4","VD","DD","AVH","DA","PVC","I5","AVM","VC","AIN","RIM","RMD","ASER","ASI","AFD","BAG","ADF","ASEL","ASG","AWA","NSM","CAN","ADL","ASK","IL2","RIA","RIC","AVK","SMD","AVA","AVE","AVG","AWB","PVD","PVM","RIS","VB","AIY"),
                     exp(1:43)) |>
  sort() |>
  names() |>
  as.numeric()


hm_callback <- function(hc, ...){
  as.hclust(reorder(as.dendrogram(hc), wts = weights))
}


pheatmap::pheatmap((hm_mat/overlap),
                   scale = "none",
                   clustering_distance_rows = "canberra",
                   clustering_distance_cols = "canberra",
                   clustering_method = "complete",
                   clustering_callback = hm_callback,
                   annotation_row = heatmap_annot,
                   annotation_col = heatmap_annot,
                   main = "Normalized by number of overlapping genes",
                   # filename = file.path(export_dir, "heatmap_ds_byoverlap.pdf"),
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

tot_genes_per_neur <- signif_genes |>
  pivot_longer(cols = c("neurA","neurB")) |>
  group_by(value) |>
  summarize(tot_neur = sum(nb_DS_genes),
            .groups = "drop")

signif_genes |>
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

signif_genes |>
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



#~ Heatmap normd by overlapping genes ----


overlap <- t(cengenDataSC::cengen_sc_3_bulk[,colnames(hm_mat_lsvs)]) %*% cengenDataSC::cengen_sc_3_bulk[,colnames(hm_mat_lsvs)]
all.equal(colnames(hm_mat_lsvs), colnames(overlap))
all.equal(rownames(hm_mat_lsvs), rownames(overlap))



# To set the order (roughly)
weights <- set_names(c("OLL","OLQ","PHA","AWC","IL1","M4","VD","DD","AVH","DA","PVC","I5","AVM","VC","AIN","RIM","RMD","ASER","ASI","AFD","BAG","ADF","ASK","IL2","RIA","RIC","AVK","SMD","ASEL","ASG","AWA","NSM","CAN","ADL","AVA","AVE","AVG","AWB","PVD","PVM","RIS","VB","AIY"),
                     exp(1:43)) |>
  sort() |>
  names() |>
  as.numeric()


hm_callback <- function(hc, ...){
  as.hclust(reorder(as.dendrogram(hc), wts = weights))
}


pheatmap::pheatmap((hm_mat_lsvs/overlap),
                   scale = "none",
                   clustering_distance_rows = "canberra",
                   clustering_distance_cols = "canberra",
                   clustering_method = "complete",
                   clustering_callback = hm_callback,
                   annotation_row = heatmap_annot,
                   annotation_col = heatmap_annot,
                   main = "Normalized by number of overlapping LSVs",
                   # filename = file.path(export_dir, "heatmap_lsvs_byoverlap.pdf"),
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


# Subsample neurons ----

# How many neurons do we need to sequence to detect that many DS genes?
# Subsamples

all_neurs <- unique(c(dpsi$neurA, dpsi$neurB))


nb_ds_genes_sub <- function(n){
  neurset <- sample(all_neurs, n)
  
  dpsi |>
    filter(neurA %in% neurset,
           neurB %in% neurset) |>
    filter(p20 >.50 & p05 < .05) |>
    pull(gene_id) |> unique() |>
    length()
}

sub_nb_genes_ds <- tibble(nb_neurs = 1:length(all_neurs),
                          nb_ds_genes = map_int(nb_neurs, nb_ds_genes_sub))

ggplot(sub_nb_genes_ds) +
  theme_classic() +
  geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  xlab("Number of neurons") +
  ylab("Number of genes detected as DS")

# Fit Michaelis-Menten model
mod_micment <- nls(nb_ds_genes ~ Gmax * nb_neurs/(b + nb_neurs),
    data = sub_nb_genes_ds,
    start = list(Gmax = 3000, b = 10))
summary(mod_micment)

sub_nb_genes_ds <- cbind(sub_nb_genes_ds, fit=predict(mod_micment))

ggplot(sub_nb_genes_ds) +
  theme_classic() +
  geom_line(aes(x=nb_neurs, y=fit), color = "red3", linetype = "dashed") +
  geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  geom_hline(aes(yintercept = coef(mod_micment)[["Gmax"]]), color = "gray50", linetype = "dotted") +
  xlab("Number of neurons sampled") +
  ylab("Number of genes detected as DS")
# ggsave("nb_genes_ds.pdf", path = export_dir,
#        width = 5, height = 5/1.375)


#~ Fit on log-linear ----
ggplot(sub_nb_genes_ds) +
  theme_classic() +
  geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  xlab("Number of neurons") +
  ylab("Number of genes detected as DS") +
  scale_x_log10()

mod_lm <- lm(nb_ds_genes ~ log(nb_neurs),
                   data = sub_nb_genes_ds)
summary(mod_lm)

sub_nb_genes_ds<- cbind(sub_nb_genes_ds, fit_lm=predict(mod_lm))

ggplot(sub_nb_genes_ds) +
  theme_classic() +
  geom_line(aes(x=nb_neurs, y=fit_lm), color = "red3", linetype = "dashed") +
  geom_point(aes(x=nb_neurs, y=nb_ds_genes)) +
  xlab("Number of neurons") +
  ylab("Number of genes detected as DS") +
  scale_x_log10()

predict(mod_lm, newdata = data.frame(nb_neurs = 120))
predict(mod_micment, newdata = data.frame(nb_neurs = 120))


#~ By proportion of expressed genes ----
# Say a gene is expressed if threshold 2 from sc detects it in at least 2 of the sampled classes

expr_sc <- cengenDataSC::cengen_sc_3_bulk > 0


prop_ds_genes_sub <- function(n){
  neurset <- sample(all_neurs, n)
  
  nb_neurs_where_gene_expr <- rowSums(expr_sc[,neurset])
  nb_genes_expr_in_neurset <- sum(nb_neurs_where_gene_expr > 1)
  
  
  nb_genes_ds_in_neurset <- dpsi |>
    filter(neurA %in% neurset,
           neurB %in% neurset) |>
    filter(p20 >.50 & p05 < .05) |>
    pull(gene_id) |> unique() |>
    length()
  
  nb_genes_ds_in_neurset/nb_genes_expr_in_neurset
}


sub_prop_genes_ds <- tibble(nb_neurs = 3:length(all_neurs),
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




# Compare litterature ----

bib_by_sf <- readRDS("../../../bulk/psi_methods/data/biblio/bib_by_SF.rds")


bib_all <- unique(unlist(unlist(bib_by_sf)))
length(bib_all)
i2s(head(bib_all), gids)

# Restrict to genes expressed in at least 2 neurons in our dataset
expr_sc <- cengenDataSC::cengen_sc_3_bulk > 0

nb_neurs_where_gene_expr <- rowSums(expr_sc[,all_neurs])
genes_in_our_neur_sample <- names(nb_neurs_where_gene_expr)[nb_neurs_where_gene_expr>2]

table(bib_all %in% genes_in_our_neur_sample)
bib_all <- intersect(bib_all, genes_in_our_neur_sample)


all_signif_genes <- dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  pull(gene_id) |>
  unique()

length(all_signif_genes)

(litt_plot <- plot(eulerr::euler(list(litterature = bib_all,
                        `this work` = all_signif_genes)),
     quantities = TRUE))

# ggsave("compare_litterature2.pdf", path = export_dir, plot = litt_plot,
#        width = 4, height = 3, units = "in")



# Classify events ----
dpsi |>
  select(gene_id, lsv_id, lsv_type, nb_sj, nb_exons, sj_coords, ir_coords) |>
  filter(gene_id == s2i("unc-36", gids)) |>
  distinct() |> View()








