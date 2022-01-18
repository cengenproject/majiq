library(tidyverse)
library(wbData)

# Norris et al 2014 ----
# Adam D. Norris, Shangbang Gao, Megan L. Norris, Debashish Ray, Arun K. Ramani,
# Andrew G. Fraser, Quaid Morris, Timothy R. Hughes, Mei Zhen, John A. Calarco
# A Pair of RNA-Binding Proteins Controls Networks of Splicing Events Contributing to
# Specialization of Neural Cell Types
# 2014
# https://www.cell.com/molecular-cell/fulltext/S1097-2765(14)00398-0

norris_2014 <- readRDS("../../../bulk/psi_methods/data/biblio/norris_2014_a_pair/targets.rds") |>
  unlist() |> unique() |>
  enframe(value = "gene_id") |>
  filter(startsWith(gene_id, "WBGene")) |>
  mutate(name = "Norris et al., 2014")



# Tan and Fraser 2017 ----
# The combinatorial control of alternative splicing in C. elegans
# https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007033

tan_fraser_2017 <- readRDS("../../../bulk/psi_methods/data/biblio/tan_fraser_2016_targets.rds") |>
  unlist() |> unique() |>
  enframe(value = "gene_id") |>
  mutate(name = "Tan and Fraser, 2017")



# Norris, Gracida and Calarco 2017 ----
# https://elifesciences.org/articles/28129

gids_280 <- wb_load_gene_ids(280) # WBcel235

norris_2017_table_s4  <- readxl::read_excel("../../../bulk/psi_methods/data/biblio/norris_2017_elife-28129-supp4-v4.xlsx") %>%
  select(gene_name = `gene name`, ends_with("Pvalue")) %>%
  mutate(gene_id = s2i(gene_name, gids_280, TRUE),
         gene_id = recode(gene_name,
                          F32A5.2 = "WBGene00017968",
                          C54G10.4 = "WBGene00008320",
                          C53B4.4 = "WBGene00008274",
                          F28F8.5 = "WBGene00009223",
                          Y55B1BR.5 = "WBGene00044401",
                          R17.3 = "WBGene00011269",
                          W09D12.1 = "WBGene00012364",
                          C18C4.5 = "WBGene00015968",
                          R11A8.7 = "WBGene00011240",
                          F28C6.4 = "WBGene00009204",
                          Y39G10AR.7 = "WBGene00021465",
                          F53C3.13 = "WBGene00018756",
                          Y40B1A.3 = "WBGene00012734",
                          C32C4.1 = "WBGene00007862",
                          F37D6.2 = "WBGene00009508",
                          M01A8.2 = "WBGene00010796",
                          Y52B11A.2 = "WBGene00013122",
                          F21A10.2 = "WBGene00008999",
                          F55A4.8 = "WBGene00018860",
                          F33A8.6 = "WBGene00009354",
                          .default = gene_id))


bib_norris_2017 <- list(`exc-7` = unique(norris_2017_table_s4$gene_id[norris_2017_table_s4$`wild type versus exc-7 Pvalue` < 0.05]),
                        `mbl-1` = unique(norris_2017_table_s4$gene_id[norris_2017_table_s4$`wild type versus mbl-1 Pvalue` < 0.05]))

norris_gracida_calarco_2017 <- bib_norris_2017 |>
  unlist() |> unique() |>
  enframe(value = "gene_id") |>
  mutate(name = "Norris, Gracida and Calarco 2017")


# Barberan-Soler et al. 2011 ----
# https://doi.org/10.1093/nar/gkq767
# Note, here ratio is log2(junction/skipped)

barberan_table_s1 <- readxl::read_excel("../../../bulk/psi_methods/data/biblio/barberan-soler_2011_SuppTable1.xls",
                                        sheet = 1,
                                        skip = 1) %>%
  select(-notes)

bib_barberan_2011 <- select(barberan_table_s1,-"probe", -gene_id) %>%
  map(~ set_names(.x,barberan_table_s1$gene_id)) %>%
  map(~ names(which(abs(.x) > 1))) %>%
  set_names(tolower(names(.))) %>%
  # set_names(str_remove(tolower(names(.)), "-")) %>%
  rlist::list.filter(length(.) > 1) %>%
  rlist::list.remove("hrpf-1/sym-2") %>%
  set_names(., recode(names(.),
                      `hrp-1` = "hrpa-1"))

barberan_soler_2011 <- bib_barberan_2011 |>
  unlist() |> unique() |>
  enframe(value = "gene_id") |>
  mutate(name = "Barberan-Soler et al., 2011")




# Others ----

# Galvin (2011): spk-1 likely regulates splicing of ced-4 (indirectly, as a kinase targeting SF)

# Kuroyanagi (2006): egl-15 regulated by asd-1 and fox-1 redundantly

# Calixto (2010): mec-8 regulates mec-2 splicing (probably unc-52 and possibly fbn-1 too)

# Ohno (2012): sup-12 and asd-2 regulate unc-60 in muscles (and probably neurons)
# Ohno (2008): asd-2 also regulates let-2 (in muscles)

# Kabat (2009) and Heintz (2016), cited by Kotugama (2019): hrp-2 (=hrpr-1) regulates unc-52, lin-10 and ret-1

# Wang (2010): grld-1 may regulate grl-1 splicing, see paper https://www.nature.com/articles/nn.2667

others <- tibble(gene_id = c("ced-4","egl-15", "unc-60", "let-2",
                             "unc-52","lin-10","ret-1","grl-1"),
                 name = c("Galvin et al., 2011","Kuroyanagu et al., 2006", "Ohno et al., 2012", "Ohno et al., 2008",
                          "Kabat et al., 2009;Kotugama et al., 2019","Kabat et al., 2009;Kotugama et al., 2019","Kabat et al., 2009;Kotugama et al., 2019",
                          "Wang et al., 2010")) |>
  mutate(gene_id = s2i(gene_id, gids_280, warn_missing = TRUE))


bind_rows(norris_2014, tan_fraser_2017, norris_gracida_calarco_2017, barberan_soler_2011,others) |>
  rename(References = name) |>
  group_by(gene_id) |>
  summarize(References = paste(References, collapse=";")) |>
  mutate(gene_symbol = i2s(gene_id, gids_280, warn_missing = TRUE)) |> View()





