
mjq |>
  filter(gene_id == s2i("vab-10", gids)) |> 
  pull(lsv_id) |> unique()

mjq |>
  filter(gene_id == s2i("vab-10", gids)) |>
  filter(lsv_id == unique(lsv_id)[1]) |>
  ggplot(aes(x=neuron,y=psi, color=junction_id)) +
  theme_classic() +
  geom_point() +
  geom_errorbar(aes(ymin = psi - sd , ymax = psi + sd)) +
  ylab("Junction usage (%)") + xlab(NULL) +
  # ylab("E(PSI) +/- SD") + xlab(NULL) +
  # ylab("E(PSI)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(limits=c(0,1), labels = scales::percent)




mjq |>
  filter(gene_id == s2i("unc-40", gids)) |>
  filter(lsv_id == unique(lsv_id)[2],
         junction_id == 1) |>
  ggplot(aes(x=neuron,y=psi)) +
  theme_classic() +
  geom_point() +
  geom_errorbar(aes(ymin = psi - sd , ymax = psi + sd), width = .5) +
  ylab("PSI +/- SD") + xlab(NULL) +
  # ylab("E(PSI)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(limits=c(NA,1), labels = scales::percent)




# gene_id in variable ----
mjq |>
  filter(gene_id == my_gene) |> 
  pull(jct_coords) |> unique()


mjq |>
  filter(gene_id == my_gene) |> 
  pull(lsv_id) |> unique()



mjq |>
  filter(gene_id == my_gene) |>
  filter(lsv_id == unique(lsv_id)[1]) |>
  ggplot(aes(x=neuron,y=psi, color=junction_id)) +
  theme_classic() +
  geom_point() +
  geom_errorbar(aes(ymin = psi - sd , ymax = psi + sd)) +
  scale_y_continuous(limits = c(0,1)) +
  ylab("E(PSI) +/- SD") + xlab(NULL) +
  # ylab("E(PSI)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

writeClipboard(my_gene)

mjq |>
  filter(gene_id == my_gene) |>
  mutate(m_min_sd = psi - 2*sd,
         m_plus_sd = psi + 2*sd) |>
  group_by(gene_id, lsv_id, junction_id) |>
  summarize(is_ds = max(m_min_sd)-min(m_plus_sd),
            .groups = "drop") 



#~ by lsv ----


mjq |>
  filter(lsv_id == my_lsv) |>
  ggplot(aes(x=neuron,y=psi, color=junction_id)) +
  theme_classic() +
  geom_point() +
  geom_errorbar(aes(ymin = psi - sd , ymax = psi + sd)) +
  scale_y_continuous(limits = c(0,1)) +
  ylab("E(PSI) +/- SD") + xlab(NULL) +
  # ylab("E(PSI)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_hline(aes(yintercept = .25),linetype="dashed",color="grey") +
  geom_hline(aes(yintercept = .75),linetype="dashed",color="grey")

writeClipboard(strsplit(my_lsv, ":")[[1]][[1]])

mjq |>
  filter(lsv_id == my_lsv) |>
  pull(jct_coords) |> unique()


mjq |>
  filter(lsv_id == my_lsv) |>
  mutate(m_min_sd = psi - 2*sd,
         m_plus_sd = psi + 2*sd) |>
  group_by(gene_id, lsv_id, junction_id) |>
  summarize(is_ds = max(m_min_sd)-min(m_plus_sd),
            .groups = "drop") 





mjq |>
  filter(gene_id == s2i("unc-40", gids),
         jct_coords == "5681447-5682284;5681447-5681876") |>
  ggplot(aes(x=neuron,y=psi, color=junction_id)) +
  theme_classic() +
  geom_point() +
  geom_errorbar(aes(ymin = psi - sd , ymax = psi + sd)) +
  ylab("E(PSI) +/- SD") + xlab(NULL) +
  # ylab("E(PSI)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# same, more streamlined

mjq |>
  select(gene_id, lsv_id, jct_coords, ir_coords) |> 
  filter(gene_id == s2i("skn-1", gids)) |> 
  distinct() |>
  View()
  
lsv_ids <- mjq |>
  filter(gene_id == s2i("skn-1", gids)) |> 
  pull(lsv_id) |> unique()


mjq |>
  filter(gene_id == s2i("skn-1", gids),
         lsv_id == lsv_ids[[4]]) |>
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







all_ds |>
  filter(lsv_id == "WBGene00000065:s:11074415-11074549")




# Nice fig for unc-36


mjq |>
  filter(gene_id == s2i("unc-36", gids),
         neuron %in% c("ADL","AVH")) |>
  mutate(lsv_name = if_else(lsv_id == "WBGene00006772:t:8197487-8197828", "downstream","upstream"),
         junction_name = case_when(lsv_name == "upstream" & junction_id == "1" ~ "inclusion",
                                   lsv_name == "upstream" & junction_id == "2" ~ "exclusion",
                                   lsv_name == "downstream" & junction_id == "1" ~ "exclusion",
                                   lsv_name == "downstream" & junction_id == "2" ~ "inclusion"),
         junction_name = factor(junction_name, levels = c("inclusion","exclusion")),
         lsv_name = factor(lsv_name, levels = c("upstream","downstream"))) |>
  ggplot(aes(x=neuron,y=psi, color=neuron)) +
  theme_bw() +
  facet_grid(junction_name ~ lsv_name) +
  geom_point() +
  geom_errorbar(aes(ymin = psi - sd , ymax = psi + sd), width = .2) +
  ylab("Percent Splied In (+/- SD)") + xlab(NULL) +
  # ylab("E(PSI)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




# Nice fig for skn-1
mjq |>
  select(gene_id, lsv_id, jct_coords, ir_coords) |> 
  filter(gene_id == s2i("skn-1", gids)) |> 
  distinct() |>
  View()

lsv_ids <- mjq |>
  filter(gene_id == s2i("skn-1", gids)) |> 
  pull(lsv_id) |> unique()


mjq |>
  filter(gene_id == s2i("skn-1", gids),
         lsv_id == lsv_ids[[4]],
         junction_id %in% c(1,3)) |>
  ggplot(aes(x=neuron,y=psi, color=junction_id)) +
  theme_classic() +
  geom_point() +
  geom_errorbar(aes(ymin = psi - sd , ymax = psi + sd), width = .5) +
  ylab("Percent Spliced-In (+/- SD)") + xlab(NULL) +
  # ylab("E(PSI)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



# Plot sc expression ----

xx <- cengenDataSC::cengen_sc_3_bulk[s2i("unc-36", gids),] * 
  cengenDataSC::cengen_TPM_bulk[s2i("unc-36", gids),colnames(cengenDataSC::cengen_sc_2_bulk)]

enframe(xx) |>
  mutate(value = if_else(value == 0, NA_real_, value),) |>
  filter(!is.na(value),
         name %in% all_neurs) |>
  mutate(neuron = fct_inorder(name),
         neuron = factor(neuron, levels = rev(levels(neuron)))) |>
  ggplot() +
  theme_minimal() +
  geom_tile(aes(x = 1,y = neuron, fill = log10(value))) +
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_fill_gradient2(name = "TPM (log)",
                      mid = "white",
                      high = "red4",
                      midpoint = .1,
                      limits = c(0,3))

# ggsave("heatmap_unc36.pdf", path = export_dir,
#        width = 2, height = 10, units = "in")  





# Plot bulk integrated expression ----
gene_expr <- read.delim("data/2021-11-30_alec_integration/aggr_ave_integrant_GeTMM_113021.tsv")
# chosen to match a FDR of 0.105 (threshold 3 of sc paper) -> 256
# chosen to match FDR = 14 (same as threshold 2 in sc props) is 73
threshold <- 73
gene_expr_bin <- gene_expr > threshold

xx <- (gene_expr[s2i("unc-36", gids),] * gene_expr_bin[s2i("unc-36", gids),]) |>
  (\(.x) set_names(as.numeric(.x), names(.x)))()


enframe(xx) |>
  mutate(value = if_else(value == 0, NA_real_, value),) |>
  filter(!is.na(value),
         name %in% all_neurs) |>
  mutate(neuron = fct_inorder(name),
         neuron = factor(neuron, levels = rev(levels(neuron)))) |>
  ggplot() +
  theme_minimal() +
  geom_tile(aes(x = 1,y = neuron, fill = log10(value))) +
  xlab(NULL) + ylab(NULL) +
  theme(axis.text.x = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_fill_gradient2(name = "GeTMM (log)",
                       mid = "white",
                       high = "red4",
                       midpoint = .1,
                       limits = c(0,3))

ggsave("heatmap_unc36.pdf", path = export_dir,
       width = 2, height = 10, units = "in")  



# Make scale for the IGV representation

tibble(CPM = 0:7) |>
  ggplot() +
  geom_tile(aes(x=1,y=CPM, fill = log1p(CPM))) +
  scale_fill_gradient2(name = "CPM",
                       low = hsv(0,0,1),
                       mid = hsv(240/360,81/100,69/100),
                       high = hsv(240/360,100/100,69/100),
                       midpoint = 1,
                       limits = c(0,1))
# ggsave("headmap_scale.pdf", path = export_dir,
#        width = 2, height = 10, units = "in")  
