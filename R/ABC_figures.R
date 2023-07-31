library(slimr)
library(readr)
library(dplyr)
library(fastcluster)
library(NbClust)
library(stringr)
library(ggplot2)
library(GGally)
library(purrr)

abc_sims <- read_rds("data/abc_param_sims_full.rds")
abc_params <- read_rds("data/abc_param_sims_full_params.rds")

## remove weird outlier
outlier <- which(abc_params$mut_rate > 0.0002)
abc_sims <- abc_sims[-outlier]
abc_params <- abc_params[-outlier, ]

ggpairs(abc_params)

abc_params_scaled <- scale(abc_params)

abc_clust <- hclust(dist(abc_params_scaled), "ward.D2")
split_clust <- cutree(abc_clust, h = 15) 

num_splits <- max(split_clust)

abc_clust_k <- kmeans(abc_params_scaled, centers = 2)

abc_clust_nest <- NbClust(abc_params_scaled, method = "kmeans")

abc_clust_nest2 <- NbClust(abc_params_scaled, method = "ward.D2")

fst_res <- abc_sims[[1]]

extract_fst <- function(fst_res) {
  fst_dat <- fst_res$output_data %>%
        dplyr::filter(name == "all_fsts") %>%
        dplyr::select(generation, data) %>%
        dplyr::mutate(data = str_split(data, " ")) %>%
        tidyr::unnest_longer(data, values_to = "fst", indices_to = "index") %>%
        dplyr::mutate(fst = as.numeric(fst), `subpop\npair` = c("BL-BR", "BL-TR", "BR-TR")[index]) 
  fst_dat
}

fst_dat <- extract_fst(abc_sims[[5]])

ggplot(fst_dat %>%
         dplyr::filter(generation > 250), aes(generation, fst)) +
  geom_path(aes(colour = `subpop\npair`)) +
  scale_y_sqrt()

extract_sample_fst <- function(fst_res, rep) {
  fst_dat <- fst_res$output_data %>%
    filter(name != "all_fsts") %>%
    group_by(name) %>%
    mutate(fst = as.numeric(data),
           year = rep(c("2006", "2007", "2008"), 100)[rank(generation)]) %>%
    group_by(year) %>%
    summarise(mean_fst = mean(fst, na.rm = TRUE)) %>%
    mutate(rep = rep)
}

all_fsts <- imap_dfr(abc_sims, ~ extract_sample_fst(.x, .y))

all_fsts <- all_fsts %>%
  filter(mean_fst < 0.2)

fst_df <- read_rds("data/fst_df.rds") %>%
  group_by(year) %>%
  summarise(mean_fst = mean(fst)) %>%
  mutate(rep = 999999)

ggplot(all_fsts) +
  geom_line(aes(year, mean_fst, group = rep), alpha = 0.1,
            position = position_jitter(height = 0.01, width = 0)) +
  geom_line(aes(year, mean_fst), data = fst_df, colour = "red", group = 1, inherit.aes = FALSE) +
  geom_point(aes(year, mean_fst), data = fst_df, colour = "red", size = 3, inherit.aes = FALSE) +
  ylab("Mean Fst") +
  xlab("") +
  theme_minimal()
