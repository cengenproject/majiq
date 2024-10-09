library(tidyverse)
library(wbData)

gids <- wb_load_gene_ids(289)

export_dir <- "presentations/2024-03_rerun/"


# Load deltaPSI
dpsi <- qs::qread("intermediates/240304_dpsi.qs")

signif_lsv <- dpsi |>
  filter(p20 >.50 & p05 < .05) |>
  pull(lsv_id) |> unique()


# Load modulizer ----
modules <- read_tsv("data/2024-03-22_voila_modulizer/heatmap.tsv", comment = "#")

table(cplx = modules$complex, composite = str_detect(modules$module_event_combination, "\\|"))

modules |>
  filter(complex,
         !str_detect(modules$module_event_combination, "\\|")) |> View()


modules |>
  filter(!complex,
         str_detect(modules$module_event_combination, "\\|")) |> View()


head(modules$module_event_combination)


tot_nb_modules <- nrow(modules)

# filter
modules |>
  select(module_id, lsv_id, complex, module_event_combination) |>
  mutate(is_ds = lsv_id %in% signif_lsv) |>
  mutate(category = if_else(complex,
                        "complex",
                        str_split_i(module_event_combination, "\\^", 1))) |>
  mutate(category_lumped = case_match(
    category,
    "afe" ~ "Alt. first exon",
    "ale" ~ "Alt. last exon",
    "alt3" ~ "Alt. 3' ss",
    "alt5" ~ "Alt. 5' ss",
    "alt3_5" ~ "Alt. 3' and 5' ss",
    "cassette" ~ "Cassette exon",
    "complex" ~ "Complex",
    "ir" ~ "Intron retention",
    "mxe" ~ "Multi-exon",
    "putative_afe" ~ "Alt. first exon",
    "putative_ale" ~ "Alt. last exon",
    "putative_alt3" ~ "Alt. 3' ss",
    "putative_alt5" ~ "Alt. 5' ss",
    "tandem_cassette" ~ "Multi-exon"
    )) |>
  filter(! category_lumped == "Alt. 3' and 5' ss") |>
  mutate(label = paste0(sum(is_ds), "/", n()),
         .by = category_lumped) |>
  # summarize(nb_modules_tot = n(),
  #           nb_modules_ds = sum(is_ds),
  #           nb_modules_nonds = nb_modules_tot - nb_modules_ds,
  #           .by = category_lumped) |>
  # mutate(percent = round(100*nb_modules_tot/tot_nb_modules),
  #        percent = paste0(percent, "%")) |>
  ggplot(aes(x = category_lumped,fill = is_ds)) +
  theme_classic() +
  geom_bar() +
  geom_text(aes(label = label, y = 1000), nudge_y = 25) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  ylab("Number of modules") + xlab(NULL)
  

my_format <- \(x) format(x, big.mark = ",")
justify_left <- \(s){
  len <- nchar(s)
  filling <- sapply(max(len) - len,
         \(.x) paste0(rep(" ", .x), collapse = ""))
  paste0(s, filling)
}


modules |>
  select(module_id, lsv_id, complex, module_event_combination) |>
  mutate(is_ds = lsv_id %in% signif_lsv) |>
  mutate(category = if_else(complex,
                            "complex",
                            str_split_i(module_event_combination, "\\^", 1))) |>
  mutate(category_lumped = case_match(
    category,
    "afe" ~ "Alt. first exon",
    "ale" ~ "Alt. last exon",
    "alt3" ~ "Alt. 3' ss",
    "alt5" ~ "Alt. 5' ss",
    "alt3_5" ~ "Alt. 3' and 5' ss",
    "cassette" ~ "Cassette exon",
    "complex" ~ "Complex",
    "ir" ~ "Intron retention",
    "mxe" ~ "Multi-exon",
    "putative_afe" ~ "Alt. first exon",
    "putative_ale" ~ "Alt. last exon",
    "putative_alt3" ~ "Alt. 3' ss",
    "putative_alt5" ~ "Alt. 5' ss",
    "tandem_cassette" ~ "Multi-exon"
  )) |>
  filter(! category_lumped == "Alt. 3' and 5' ss") |>
  summarize(n_tot = n(),
            nb_modules_ds = sum(is_ds),
            nb_modules_nonds = n_tot - nb_modules_ds,
            .by = category_lumped) |>
  mutate(lab = paste0(my_format(nb_modules_ds), "/", my_format(n_tot)),
         .by = category_lumped) |>
  mutate(lab = justify_left(lab)) |>
  pivot_longer(-c(category_lumped, lab, n_tot),
               names_to = "ds",
               values_to = "n",
               names_prefix = "nb_modules_") |>
  mutate(ds = factor(ds, levels = c("nonds","ds"))) |>
  ggplot(aes(x = category_lumped, y = n, fill = ds)) +
  theme_classic() +
  geom_col() +
  geom_text(aes(label = lab, y = n_tot), nudge_y = 250, angle = 90) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_manual(values = c("grey30", "darkred")) +
  ylab("Number of modules") + xlab(NULL) +
  theme(legend.position = 'none') +
  scale_y_continuous(limits = c(0, 2500))


# ggsave("DS_module_category2.pdf", path = export_dir,
#        width = 5, height = 6, units = "in")





# Compare AF with suppa's definition

suppa_psi <- read.delim("../suppa_events/data/240301b_psiPerEvent.psi") |>
  rownames_to_column("event_id") |>
  as_tibble() |>
  separate_wider_regex(event_id,
                       patterns = c(gene_id = "^WBGene[0-9]{8}", ";",
                                    event_type = "[SEA53MXRIFL]{2}", "\\:",
                                    event_coordinates = "[IXV]+\\:[0-9:\\-]+:[+-]$"),
                       cols_remove = FALSE)



majiq_afe <- modules |>
  mutate(category = if_else(complex,
                            "complex",
                            str_split_i(module_event_combination, "\\^", 1)))  |>
  filter(category == "afe") |>
  select(gene_id, gene_name, lsv_id)



suppa_afe <- suppa_psi |>
  filter(event_type == "AF") |>
  select(gene_id, event_coordinates)


majiq_afe |>
  summarize(n_ev = n(),
            .by = "gene_id") |>
  arrange(desc(n_ev))

suppa_afe |>
  summarize(n_ev_suppa = n(),
            .by = "gene_id") |>
  full_join(majiq_afe |>
              summarize(n_ev_majiq = n(),
                        .by = "gene_id"),
            by = "gene_id") |>
  arrange(desc(n_ev_suppa)) |> View()



majiq_afe |>
  filter(gene_id == "WBGene00020294")
suppa_afe |>
  filter(gene_id == "WBGene00006792")

modules |>
  filter(gene_id == "WBGene00020294") |> View()





# Check cassette
suppa_se <- suppa_psi |>
  filter(event_type == "SE") |>
  select(gene_id, event_coordinates)



majiq_se <- modules |>
  mutate(category = if_else(complex,
                            "complex",
                            str_split_i(module_event_combination, "\\^", 1)))  |>
  filter(category == "cassette") |>
  select(gene_id, gene_name, lsv_id)


suppa_se |>
  summarize(n_ev_suppa = n(),
            .by = "gene_id") |>
  full_join(majiq_se |>
              summarize(n_ev_majiq = n(),
                        .by = "gene_id"),
            by = "gene_id") |>
  mutate(n_ev_suppa = replace_na(n_ev_suppa, 0),
         n_ev_majiq = replace_na(n_ev_majiq, 0)) |>
  arrange(desc(n_ev_suppa)) |>
  ggplot() +
  geom_jitter(aes(x = n_ev_majiq,
                 y = n_ev_suppa),
             alpha = .2,
             width = .1,
             height = .1)



















