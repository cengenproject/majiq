
norris_2014 <- readRDS("data/biblio/norris_2014_a_pair/targets.rds") |>
  unlist() |> unique() |>
  enframe(value = "gene_id") |>
  mutate(name = "Norris et al., 2014")
tan_fraser_2017 <- readRDS("data/biblio/tan_fraser_2016_targets.rds") |>
  unlist() |> unique() |>
  enframe(value = "gene_id") |>
  mutate(name = "Tan and Fraser, 2017")
norris_gracida_calarco_2017 <- bib_norris_2017 |>
  unlist() |> unique() |>
  enframe(value = "gene_id") |>
  mutate(name = "Norris, Gracida and Calarco 2017")
barberan_soler_2011 <- bib_barberan_2011 |>
  unlist() |> unique() |>
  enframe(value = "gene_id") |>
  mutate(name = "Barberan-Soler et al., 2011")

others <- tibble(gene_id = c("ced-4"),
                 name = c("Galvin et al., 2011"))
