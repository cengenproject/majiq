# Load quantifications from MAJIQ and use them for general descriptions


## Initializations ----
library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(289)



# Load and format data ----

# gene_patterns <- readRDS("data/cengen_sc/211028_genes_categorized_by_pattern.rds")
gene_patterns <- wormDatasets::genes_by_pattern


path_psi <- "data/2024-03-04_outs/psi/"
files_psi <- list.files(path_psi, pattern = "*.tsv")

# Load file, one row per LSV
mjqdta <- map_dfr(files_psi,
                  ~read_tsv(file.path(path_psi, .x),
                            col_types = cols(
                              gene_id = col_character(),
                              lsv_id = col_character(),
                              lsv_type = col_character(),
                              mean_psi_per_lsv_junction = col_character(),
                              stdev_psi_per_lsv_junction = col_character(),
                              num_junctions = col_integer(),
                              num_exons = col_integer(),
                              junctions_coords = col_character(),
                              ir_coords = col_character()
                            ),
                            na = "na") |>
                    add_column(neuron = str_split(.x,"\\.")[[1]][1]))

# separate data in E(PSI) and Std(PSI) fields
# We get one row per junction
mjq <- mjqdta |>
  separate_rows(mean_psi_per_lsv_junction,
                stdev_psi_per_lsv_junction,
                sep = ";") |>
  group_by(neuron, lsv_id) |>
  mutate(junction_id = row_number()) |>
  ungroup() |>
  dplyr::select(neuron,
                gene_id,
                lsv_id,
                junction_id,
                psi = mean_psi_per_lsv_junction,
                sd = stdev_psi_per_lsv_junction,
                jct_coords = junctions_coords,
                ir_coords) |>
  mutate(psi = as.double(psi),
         sd = as.double(sd),
         junction_id = factor(junction_id))







# Global LSV characterizations ----

nb_jctions_per_lsv <- mjq |>
  group_by(gene_id, neuron, lsv_id) |>
  summarize(nb_jctions = n()) |> 
  group_by(gene_id, lsv_id) |>
  summarize(max_nb_jctions = max(nb_jctions))

ggplot(nb_jctions_per_lsv) +
  theme_classic() +
  geom_histogram(aes(x=max_nb_jctions),bins=60, color = "white") +
  scale_x_continuous(trans = "log1p", breaks = c(0,1:10,12,15,20,30,40,50,100)) +
  scale_y_continuous(trans = "log1p", breaks = c(0,1,10,100,1000,10000)) +
  xlab("Number of junctions per LSV")

# A look at the LSVs with a lot of junctions
nb_jctions_per_lsv |>
  filter(max_nb_jctions >= 10) |>
  group_by(gene_id) |>
  summarize(nb_jctions = paste0(max_nb_jctions, collapse=", "),
            max_max = max(max_nb_jctions)) |>
  arrange(desc(max_max)) |>
  mutate(gene_name = i2s(gene_id, gids)) |>
  select(gene_name,
         `Number of junctions per LSV` = nb_jctions) |>
  head()
  # clipr::write_clip()



nb_lsv_per_gene <- mjq |>
  group_by(neuron, gene_id)|>
  summarize(nb_lsv = n()) |>
  group_by(gene_id) |>
  summarize(max_nb_lsv = max(nb_lsv))

ggplot(nb_lsv_per_gene) +
  theme_classic() +
  geom_histogram(aes(x=max_nb_lsv),bins=60, color = "white") +
  scale_x_continuous(trans = "log1p", breaks = c(0,1:10,12,15,20,30,40,50,100,200)) +
  scale_y_continuous(trans = "log1p", breaks = c(0,1,10,100,1000,2000)) +
  xlab("Number of LSV per gene")

# Look at genes with many LSVs
nb_lsv_per_gene |>
  filter(max_nb_lsv >= 30) |>
  arrange(desc(max_nb_lsv)) |>
  mutate(gene_name = i2s(gene_id, gids)) |>
  head()
  # clipr::write_clip()


## Call differential splicing ----


# Examples to develop method

mjq |>
  filter(gene_id == s2i("unc-40", gids),
         lsv_id == "WBGene00006776:s:5682284-5682373") |>
  ggplot(aes(x=neuron,y=psi, color=junction_id)) +
  theme_classic() +
  geom_point() +
  geom_errorbar(aes(ymin = psi - sd , ymax = psi + sd)) +
  ylab("E(PSI) +/- SD") + xlab(NULL) +
  # ylab("E(PSI)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mjq |>
  filter(gene_id == s2i("unc-40", gids)) |>
  group_by(lsv_id, junction_id) |>
  mutate(is_DS = max(psi-sd) > min(psi+sd) |  min(psi+sd) < max(psi-sd)) |> 
  pull(is_DS)


mjq |>
  filter(lsv_id == "WBGene00000065:s:11074415-11074549") |>
  group_by(lsv_id, junction_id) |>
  mutate(is_DS = max(psi-sd) > .25 & min(psi+sd) < .25 |
                 min(psi+sd) < .75 & max(psi-sd) > .75) |> 
  select(1,3,4,5,6,9)




mjq |>
  filter(lsv_id == "WBGene00000065:s:11074415-11074549") |>
  ggplot(aes(x=neuron,y=psi, color=junction_id)) +
  theme_classic() +
  geom_point() +
  geom_errorbar(aes(ymin = psi - sd , ymax = psi + sd)) +
  ylab("E(PSI) +/- SD") + xlab(NULL) +
  # ylab("E(PSI)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




mjq |>
  filter(gene_id == s2i("unc-40", gids)) |> 
  pull(lsv_id) |> unique()



all_ds |>
  filter(lsv_id == "WBGene00000065:s:11074415-11074549")


# pmk-2 all over the place
mjq |>
  filter(lsv_id == "WBGene00004056:s:8147927-8148131") |>
  mutate(m_min_sd = psi - 2*sd,
         m_plus_sd = psi + 2*sd) |>
  group_by(lsv_id, junction_id) |>
  mutate(is_ds = min(m_plus_sd)<max(m_min_sd)) |>
  ggplot(aes(x=neuron,y=psi, color=junction_id)) +
  theme_classic() +
  geom_point() +
  geom_errorbar(aes(ymin = psi - sd , ymax = psi + sd)) +
  ylab("E(PSI) +/- SD") + xlab(NULL) +
  scale_y_continuous(limits = c(0,1)) +
  # ylab("E(PSI)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# cah-6 none
mjq |>
  filter(lsv_id == "WBGene00000284:s:3630246-3630476") |>
  mutate(m_min_sd = psi - 2*sd,
         m_plus_sd = psi + 2*sd) |>
  group_by(lsv_id, junction_id) |>
  mutate(is_ds = max(m_min_sd)-min(m_plus_sd)>.1) |>
  ggplot(aes(x=neuron,y=psi, color=junction_id)) +
  theme_classic() +
  geom_point() +
  geom_errorbar(aes(ymin = psi - 2*sd , ymax = psi + 2*sd)) +
  ylab("E(PSI) +/- SD") + xlab(NULL) +
  scale_y_continuous(limits = c(0,1)) +
  # ylab("E(PSI)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




# Compute DS ----

# all_ds1 <- mjq |>
#   group_by(lsv_id, junction_id) |>
#   mutate(is_DS = max(psi-sd) > .25 & min(psi+sd) < .25 |
#            min(psi+sd) < .75 & max(psi-sd) > .75) |>
#   group_by(gene_id, lsv_id, junction_id) |>
#   summarize(is_DS = any(is_DS),
#             .groups = "drop") |>
#   group_by(gene_id, lsv_id) |>
#   summarize(nb_DS = sum(is_DS),
#             .groups = "drop") |>
#   group_by(gene_id) |>
#   summarize(has_DS = any(nb_DS > 0),
#             .groups = "drop")


all_ds <- mjq |>
  mutate(m_min_sd = psi - 2*sd,
         m_plus_sd = psi + 2*sd) |>
  group_by(gene_id, lsv_id, junction_id) |>
  summarize(is_ds = max(m_min_sd)-min(m_plus_sd)>.1,
            .groups = "drop") |>
  group_by(gene_id, lsv_id) |>
  summarize(nb_ds = sum(is_ds),
            .groups = "drop") |>
  group_by(gene_id) |>
  summarize(has_ds = any(nb_ds > 0),
            .groups = "drop") |>
  mutate(gene_name = i2s(gene_id, gids))



table(all_ds$has_ds)





# Neuronal enrichment ----


table(has_ds = all_ds$has_ds,
      is_neuronal = all_ds$gene_id %in% gene_patterns$present_in_neurons)

table(gene_patterns$nondetected %in% all_ds$gene_id)
table(gene_patterns$nonneuronal %in% all_ds$gene_id)
table(gene_patterns$present_in_neurons %in% all_ds$gene_id)


gene_patterns$nondetected[gene_patterns$nondetected %in% all_ds$gene_id] |>
  i2s(gids)
#> essentially unstudied

tibble(gene_type = names(gene_patterns)) |>
  mutate(prop_splicing = map_dbl(gene_type,
                                 ~ mean(all_ds$has_ds[all_ds$gene_id %in% gene_patterns[[.x]]])),
         nb_genes = map_int(gene_type,
                            ~ length(gene_patterns[[.x]])),
         nb_diff_spliced = map_dbl(gene_type,
                              ~ sum(all_ds$has_ds[all_ds$gene_id %in% gene_patterns[[.x]]])),
         nb_non_diff_spliced = map_dbl(gene_type,
                                 ~ sum(! all_ds$has_ds[all_ds$gene_id %in% gene_patterns[[.x]]]))) |> 
  mutate(nb_unquantified = nb_genes - (nb_diff_spliced + nb_non_diff_spliced)) |>
  select(-prop_splicing,
         sc_pattern = gene_type) |> 
  pivot_longer(-c(sc_pattern , nb_genes),
               names_to = "majiq",
               names_prefix = "nb_",
               values_to = "nb") |>
  mutate(majiq = factor(majiq, levels = c("diff_spliced","non_diff_spliced","unquantified")),
         sc_pattern = factor(sc_pattern, levels = c("nondetected","nonneuronal",
                                                    "present_in_neurons","broadly_neuronal",
                                                    "panneuronal","neuron_specific","ubiquitous"))) |> 
  ggplot() +
  theme_classic() +
  geom_col(aes(x = sc_pattern, y = nb, fill = majiq), position = "fill") +
  geom_text(aes(x= sc_pattern, y = .1, label = nb_genes)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Expression pattern in scRNA-Seq") +
  ylab(NULL) +
  scale_y_continuous(labels = scales::percent)



# Build a new list slightly different
fig_gene_patterns <- list(in_every_neuron = character(),
                          specific_to_neurons = character(),
                          other_neuronal = character(),
                          panneuronal = character(),
                          ubiquitous = character())

fig_gene_patterns$in_every_neuron <- gene_patterns$panneuronal |>
                            union(gene_patterns$broadly_neuronal) |>
                            union(gene_patterns$ubiquitous)

fig_gene_patterns$specific_to_neurons = gene_patterns$neuron_specific |>
  union(gene_patterns$panneuronal) |>
  union(gene_patterns$broadly_neuronal)

fig_gene_patterns$other_neuronal <- gene_patterns$present_in_neurons |>
  setdiff(fig_gene_patterns$in_all_neurons) |>
  setdiff(fig_gene_patterns$specific_to_neurons)

fig_gene_patterns$panneuronal <- gene_patterns$panneuronal

fig_gene_patterns$ubiquitous <- gene_patterns$ubiquitous



tibble(gene_type = names(fig_gene_patterns)) |>
  mutate(prop_splicing = map_dbl(gene_type,
                                 ~ mean(all_ds$has_ds[all_ds$gene_id %in% fig_gene_patterns[[.x]]])),
         nb_genes = map_int(gene_type,
                            ~ length(fig_gene_patterns[[.x]])),
         nb_diff_spliced = map_dbl(gene_type,
                                   ~ sum(all_ds$has_ds[all_ds$gene_id %in% fig_gene_patterns[[.x]]])),
         nb_non_diff_spliced = map_dbl(gene_type,
                                       ~ sum(! all_ds$has_ds[all_ds$gene_id %in% fig_gene_patterns[[.x]]]))) |> 
  mutate(nb_unquantified = nb_genes - (nb_diff_spliced + nb_non_diff_spliced)) |>
  select(-prop_splicing,
         sc_pattern = gene_type) |> 
  pivot_longer(-c(sc_pattern , nb_genes),
               names_to = "majiq",
               names_prefix = "nb_",
               values_to = "nb") |>
  mutate(majiq = factor(majiq, levels = c("diff_spliced","non_diff_spliced","unquantified")),
         sc_pattern = fct_inorder(sc_pattern)) |> 
  ggplot() +
  theme_classic() +
  geom_col(aes(x = sc_pattern, y = nb, fill = majiq), position = "fill") +
  geom_text(aes(x= sc_pattern, y = .1, label = nb_genes)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab("Expression pattern in scRNA-Seq") +
  ylab(NULL) +
  scale_y_continuous(labels = scales::percent)


tibble(sc_pattern = names(fig_gene_patterns)) |>
  mutate(nb_genes = map_int(sc_pattern,
                            ~ length(fig_gene_patterns[[.x]])),
         nb_diff_spliced = map_dbl(sc_pattern,
                                   ~ sum(all_ds$has_ds[all_ds$gene_id %in% fig_gene_patterns[[.x]]]))) |>
  mutate(sc_pattern = fct_inorder(sc_pattern),
         prop = nb_diff_spliced/nb_genes,
         sep = sqrt(prop*(1-prop)/nb_genes)) |>
  ggplot(aes(x = sc_pattern, y = prop)) +
  theme_classic() +
  geom_col() +
  geom_errorbar(aes(ymin = prop-1.96*sep, ymax = prop+1.96*sep),
                width = .2) +
  scale_y_continuous(limits = c(0,1), labels = scales::percent) +
  xlab("Expression pattern in scRNA-Seq") +
  ylab("Proportion of genes AS (+/- 95% CI)") +
  xlab(NULL)



#~ Genes expressed in all current neurons ----

# we see that genes expressed in all neurons are more likely to be DS. But that could be because
# we are not measuring every neuron: perhaps 40% of all genes are DS, but 15% of the events happen
# in unsequenced neurons and not detected here.

cur_neurs <- unique(mjq$neuron)

prop_neurs_with_expr <- rowMeans(cengenDataSC::cengen_sc_1_bulk[,cur_neurs])
fig_gene_patterns$in_sequenced_neurons <- names(prop_neurs_with_expr)[prop_neurs_with_expr > .9] |>
  setdiff(fig_gene_patterns$in_every_neuron)



#~ Examine some subcategories ----


# Genes that are panneurs in sc but not quantified by majiq

setdiff(c(gene_patterns$broadly_neuronal, gene_patterns$panneuronal),
        all_ds$gene_id) |> 
  # clipr::write_clip()
  i2s(gids)


# Genes that are non-neuronal (sc) and unquantified (majiq)

setdiff(gene_patterns$nonneuronal,
        all_ds$gene_id) |> #length()
  sample(20) |> i2s(gids)
  # clipr::write_clip()
#> mostly just noise; a few are high (possibly bc conta or miscategorized in sc), but no DS



# Genes that are non-neuronal (sc) and DS (majiq)
all_ds |> 
  filter(gene_id %in% gene_patterns$nonneuronal,    #(\(.x) table(.x$has_DS))()
         has_DS) |>
  pull(gene_id) |> #length()
  sample(20) |> 
  clipr::write_clip()

#> most have high level, and the splicing looks true. Could be conta or miscategorization



# Genes that are panneur and non-diff spliced (though quantified)
all_ds |>
  filter(gene_id %in% gene_patterns$panneuronal,
         !has_DS) |>
  pull(gene_id) |> #length()
  sample(20) |> 
  clipr::write_clip()



# Ubiquitous genes

# ubiq and unquantif
setdiff(gene_patterns$ubiquitous,
        all_ds$gene_id) |> #length()
  sample(20) |> 
  clipr::write_clip()
#> seems to be highly expressed but no AS to quantify (e.g. rpl genes)


# Multi-exon genes
all_exons <- wb_load_exon_coords(277)
nb_ex <- all_exons |>
  group_by(gene_id) |>
  summarize(n = n(),
            .groups = "drop")
hist(nb_ex$n, breaks = 150)
table(nb_ex$n)

plot(eulerr::euler(list(ds = all_ds$gene_id[all_ds$has_ds],
                        multi_ex = nb_ex$gene_id[nb_ex$n > 1])))
table(nb_ex$gene_id %in% all_ds$gene_id)
table(all_ds$gene_id %in% nb_ex$gene_id)




# Gene DS in a single neuron type ----


ds_in_single_neur_type <- mjq |>
  group_by(gene_id, lsv_id, junction_id) |>
  summarize(is_ds = (sum(psi + sd < .25) >= 10 & sum(psi - sd > .25) == 1) |
              (sum(psi - sd > .75) >= 10 & sum(psi + sd < .75) == 1),
            .groups = "drop") |>
  group_by(gene_id, lsv_id) |>
  summarize(nb_ds = sum(is_ds),
            .groups = "drop") |>
  group_by(gene_id) |>
  summarize(has_ds = any(nb_ds > 0),
            .groups = "drop") |>
  mutate(gene_name = i2s(gene_id, gids))

ds_in_single_neur_type2 <- mjq |>
  group_by(gene_id, lsv_id, junction_id) |>
  summarize(is_ds = (sum(psi + sd < .35) >= 10 & sum(psi - sd > .25) == 1) |
              (sum(psi - sd > .65) >= 10 & sum(psi + sd < .75) == 1),
            .groups = "drop") |>
  group_by(gene_id, lsv_id) |>
  summarize(nb_ds = sum(is_ds),
            .groups = "drop") |>
  group_by(gene_id) |>
  summarize(has_ds = any(nb_ds > 0),
            .groups = "drop") |>
  mutate(gene_name = i2s(gene_id, gids))




ds_in_single_neur_type_lsv <- mjq |>
  group_by(gene_id, lsv_id, junction_id) |>
  summarize(is_ds = (sum(psi + sd < .25) >= 10 & sum(psi - sd > .25) == 1) |
              (sum(psi - sd > .75) >= 10 & sum(psi + sd < .75) == 1),
            .groups = "drop") |>
  group_by(gene_id, lsv_id) |>
  summarize(nb_ds = sum(is_ds),
            .groups = "drop") |>
  mutate(gene_name = i2s(gene_id, gids))


ds_in_single_neur_type_lsv2 <- mjq |>
  group_by(gene_id, lsv_id, junction_id) |>
  summarize(is_ds = (sum(psi + sd < .35) >= 10 & sum(psi - sd > .25) == 1) |
              (sum(psi - sd > .65) >= 10 & sum(psi + sd < .75) == 1),
            .groups = "drop") |>
  group_by(gene_id, lsv_id) |>
  summarize(nb_ds = sum(is_ds),
            .groups = "drop") |>
  mutate(gene_name = i2s(gene_id, gids))


ds_in_single_neur_type2[ds_in_single_neur_type2$has_ds & !ds_in_single_neur_type$has_ds,]


ds_in_single_neur_type_lsv2[ds_in_single_neur_type_lsv2$gene_name == "nhr-128",]

ds_in_single_neur_type3 <- mjq |>
  group_by(gene_id, lsv_id, junction_id) |>
  summarize(highest = max(psi - sd),
            highest_var_inter = sort(psi, decreasing = TRUE)[2] + 2*sd(psi),
            lowest = min(psi + sd),
            lowest_var_inter  = sort(psi, decreasing = FALSE)[2] - 2*sd(psi),
            nb_neurs = n(),
            .groups = "drop") |>
  mutate(is_ds = nb_neurs > 5 & (highest - highest_var_inter > 0.05 |
                                   lowest - lowest_var_inter < -0.05)) |>
  group_by(gene_id, lsv_id) |>
  summarize(nb_ds = sum(is_ds),
            .groups = "drop") |>
  group_by(gene_id) |>
  summarize(has_ds = any(nb_ds > 0),
            .groups = "drop") |>
  mutate(gene_name = i2s(gene_id, gids))


ds_in_single_neur_type_lsv3 <- mjq |>
  group_by(gene_id, lsv_id, junction_id) |>
  summarize(highest = max(psi - sd),
            highest_var_inter = sort(psi, decreasing = TRUE)[2] + 2*sd(psi),
            lowest = min(psi + sd),
            lowest_var_inter  = sort(psi, decreasing = FALSE)[2] - 2*sd(psi),
            nb_neurs = n(),
            .groups = "drop") |>
  mutate(is_ds = nb_neurs > 5 & (highest - highest_var_inter > 0.05 |
                                   lowest - lowest_var_inter < -0.05)) |>
  group_by(gene_id, lsv_id) |>
  summarize(nb_ds = sum(is_ds),
            .groups = "drop") |>
  mutate(gene_name = i2s(gene_id, gids))


my_gene <- sample(ds_in_single_neur_type3$gene_id[which(ds_in_single_neur_type3$has_ds)], 1)
ds_in_single_neur_type_lsv3 |>
  filter(gene_id == my_gene)




mjq |>
  filter(lsv_id == "WBGene00003166:t:5582959-5583062") |>
  group_by(gene_id, lsv_id, junction_id) |>
  summarize(highest = max(psi - sd),
            highest_var_inter = sort(psi, decreasing = TRUE)[2] + 2*sd(psi),
            lowest = min(psi + sd),
            lowest_var_inter  = sort(psi, decreasing = FALSE)[2] - 2*sd(psi),
            .groups = "drop") |>
  mutate(is_ds = highest > highest_var_inter |
           lowest < lowest_var_inter) 


# AS event type ----

# 




# Keep only LSVs with DS
# all_ds_lsv <- mjq |>
#   group_by(lsv_id, junction_id) |>
#   mutate(is_DS = max(psi-sd) > .25 & min(psi+sd) < .25 |
#            min(psi+sd) < .75 & max(psi-sd) > .75) |>
#   group_by(gene_id, lsv_id, junction_id) |>
#   summarize(is_DS = any(is_DS),
#             .groups = "drop") |>
#   group_by(gene_id, lsv_id) |>
#   summarize(nb_DS = sum(is_DS),
#             .groups = "drop")


all_ds_lsv <- mjq |>
  mutate(m_min_sd = psi - 2*sd,
         m_plus_sd = psi + 2*sd) |>
  group_by(gene_id, lsv_id, junction_id) |>
  summarize(is_ds = max(m_min_sd)-min(m_plus_sd)>.1,
            .groups = "drop") |>
  group_by(gene_id, lsv_id) |>
  summarize(nb_ds = sum(is_ds),
            .groups = "drop") |>
  mutate(gene_name = i2s(gene_id, gids))

# Compare two ways to call DS
# xx <- all_ds_lsv2[all_ds_lsv2$nb_ds > 0 & all_ds_lsv$nb_DS == 0,]
# # my_lsv <- sample(xx$lsv_id, 1)
# 
# # xx <- all_ds_lsv2[all_ds_lsv$nb_DS > 0 & all_ds_lsv2$nb_ds == 0,]
# # my_lsv <- sample(xx$lsv_id, 1)


table(all_ds_lsv$nb_DS > 0)
table(all_ds_lsv2$nb_ds > 0)
table(all_ds_lsv2$nb_ds > 0 & all_ds_lsv$nb_DS == 0)
table(all_ds_lsv$nb_DS > 0 & all_ds_lsv2$nb_ds == 0)




mjqdta2 <- mjqdta |>
  filter(lsv_id %in% all_ds_lsv$lsv_id[all_ds_lsv$nb_ds > 0])


# local type
mjqdta2 |>
  select(A5SS,A3SS,ES) |>
  eulerr::euler()



# Find the first and last exon of every tx (note, first/last in genomic coordinates, so inverted for - strand)
txdb <- wb_load_TxDb(277)
all_exons <- exonsBy(txdb, "tx", use.names = TRUE)

all_first_start <- lapply(start(all_exons), \(.x) .x[[1]])
all_last_start <- lapply(start(all_exons), \(.x) .x[[length(.x)]])

all_first_end <- lapply(end(all_exons), \(.x) .x[[1]])
all_last_end <- lapply(end(all_exons), \(.x) .x[[length(.x)]])
rm(all_exons)
rm(txdb)

mjqdta2 |>
  mutate(is_first = str_detect(all_first_end[startsWith(names(all_first_end),
                                                        wb_id2seq(xx$`Gene ID`, gids))],
                               xx$`Junctions coords`))

xx


mjqdta2 |>
  slice_sample(n=1) |> 
  as.data.frame() ->xx
xx
writeClipboard(xx$`Gene ID`)
all_first_end[startsWith(names(all_first_end),
                         wb_id2seq(xx$`Gene ID`, gids))]
str_detect(xx$`Junctions coords`,
           paste0(all_first_end[startsWith(names(all_first_end),
                                           wb_id2seq(xx$`Gene ID`, gids))],
                  collapse = "|"))


all_first_end[]




# Looking at some events ----
mjqdta2 |>
  slice_sample(n=1) |> 
  as.data.frame() ->xx


table(mjqdta$A3SS, mjqdta$A5SS)

mjqdta |>
  # filter(A3SS & A5SS) |>
  slice_sample(n=1) |> 
  as.data.frame()



format_sj_coords <- function(s){
  patterns <- str_match_all(s, "([\\d]+)-([\\d]+)")[[1]]
  all_starts <- format(as.integer(patterns[,2]), big.mark = ",")
  all_ends <- format(as.integer(patterns[,3]), big.mark = ",")
  paste(all_starts, all_ends, sep = "-", collapse = "; ")
}


xx <- mjqdta2 |>
  filter(A3SS, A5SS) |>
  slice_sample(n=1) |>
  mutate(`Junctions coords` = format_sj_coords(`Junctions coords`)) |>
  as.data.frame()



xx <- mjqdta2 |>
  select(`LSV ID`, A5SS, A3SS, ES) |>
  mutate(across(where(is_logical), as.integer))

UpSetR::upset(xx)

xx <- mjqdta2 |>
  select(A5SS, A3SS, ES)

xx2 <- xx[apply(xx,1,any),]

eulerr::euler(list(#all = unique(mjqdta2$`LSV ID`),
                   `alternative 3' SS` = unique(mjqdta2$`LSV ID`[mjqdta2$A3SS]),
                   `alternative 5' SS` = unique(mjqdta2$`LSV ID`[mjqdta2$A5SS]),
                   `exon skipping` = unique(mjqdta2$`LSV ID`[mjqdta2$ES]))) |>
  plot()


xx <- list(all = unique(mjqdta2$`LSV ID`),
           a3 = unique(mjqdta2$`LSV ID`[mjqdta2$A3SS]),
           a5 = unique(mjqdta2$`LSV ID`[mjqdta2$A5SS]),
           es = unique(mjqdta2$`LSV ID`[mjqdta2$ES]))
length(unique(unlist(xx)))
length(unique(xx$es))
table(xx$es %in% xx$all)
which(! xx$es %in% xx$all)


# AS event type ----

# based on Calarco's paper:
# cassette-type, alternative 3'/5' splice site selection, alternative 
# start/terminal exons, mutually exclusive exons and intron retention
# 
# Cassette-type exons will have two LSVs 
# with two junctions each, involving the same three exons; one junction unique to each LSV and one 
# shared. Alternative 3'/5' splice site selection will have one LSV comprised of two junctions and
# two  exons.  Similarly,  alternative  start  or  terminal  exons  will  have  two  junctions,  both  with  the 
# same donor site on a common exon, and two more exons that have no upstream or downstream 
# junction connectivity, respectively. Mutually exclusive exons have two LSVs with two junctions 
# and  three  exons  each.  The  two  exons  that  are  shared  by  both  LSVs  are  the  mutually  exclusive 
# exons,  which  have  no  connecting  junction.  Intron  retention  events  were  explicitly  detected  by 
# Majiq and filtered for sufficient read coverage and reproducibility. Complex events were classified 
# as LSVs with more than three junctions connecting more than four exons. Lastly, for the splicing 
# events  of  all  classes,  it  was  verified  that  no  junctions  external  to  the  relevant  LSVs  overlapped 
# within the boundaries of each categorized splicing event.

