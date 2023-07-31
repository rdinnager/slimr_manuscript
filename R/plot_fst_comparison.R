plot_fst_comparison <- function(sim_result, observed_fst, tech_reps = 6L) {
  slim_genlights <- slim_extract_genlight(sim_result, name = "pop_sample",
                                       by = "generation")
  
  names(sim_result) <- paste0("rep_", seq_along(sim_result))
  
  slim_subpops <- imap_dfr(sim_result, ~ .x$output_data %>%
                            filter(name == "subpops") %>%
                            filter(!duplicated(generation)) %>% ## deals with an intermittent bug
                            slim_results_to_data() %>%
                            unnest_longer(col = data, values_to = "subpop",
                                          indices_to = "genome_num") %>%
                            mutate(individual_num = c(rep(1:sum(samp_sizes[ ,
                                                                            1]), each = 2),
                                                      rep(1:sum(samp_sizes[ , 2]), each = 2),
                                                      rep(1:sum(samp_sizes[ , 3]), each = 2)) %>%
                                     rep(tech_reps),
                                   rep = .y) %>%
                            group_by(generation, individual_num, rep) %>%
                            summarise(subpop = subpop[1],
                                      .groups = "drop"))
  
  ## reorder by individual label
  slim_genlights$genlight <- map(slim_genlights$genlight,
                                 ~.x[order(as.numeric(.x$ind.names)), ])
  ## add subpop data
  slim_genlights$genlight <- map2(slim_genlights$genlight, 
                                  slim_subpops %>%
                                    group_split(rep, generation),
                                  ~{pop(.x) <- .y$subpop; .x})
  
  
  fst_by_year_sim <- map(slim_genlights$genlight,
                     ~pairFST(.x) %>%
                       as.dist() %>%
                       as.matrix())
  years <- rep(c("2006", "2007", "2008"), tech_reps * 6)
  
  
  fst_sim_df <- imap_dfr(fst_by_year_sim,
                     ~combn(c("Subpopulation<p1>", 
                              "Subpopulation<p2>", 
                              "Subpopulation<p3>"), 2) %>%
                       t() %>%
                       as.data.frame() %>%
                       rename(pop1 = V1, pop2 = V2) %>%
                       mutate(fst = .x[cbind(pop1, pop2)],
                              year = years[.y],
                              rep = slim_genlights$rep[.y],
                              generation = slim_genlights$generation[.y],
                              pop_combo = paste(pop1, pop2, sep = " to "))) %>%
    group_by(rep, generation, year) %>%
    summarise(fst = mean(fst), .groups = "drop") %>%
    mutate(iter = rep(1:(tech_reps*6), each = 3))
  
  all_summ <- fst_sim_df %>%
    group_by(year) %>%
    summarise(fst_mean = mean(fst), fst_sd = sd(fst),
              fst = fst_mean,
              .groups = "drop")
  
  fst_obs_df <- observed_fst %>%
    group_by(year) %>%
    summarise(fst_mean = mean(fst),
              fst = fst_mean,
              iter = 99999,
              .groups = "drop")
  
  fst_plot_df1 <- bind_rows(fst_sim_df %>% mutate(type = "simulated"),
                           fst_obs_df %>% mutate(type = "observed")) 
  
  fst_plot_df2 <- bind_rows(all_summ %>% mutate(type = "simulated"),
                           fst_obs_df %>% mutate(type = "observed")) 
  
  p1 <- ggplot(fst_plot_df1, aes(year, fst)) +
    geom_errorbar(aes(ymin = fst_mean - fst_sd,
                      ymax = fst_mean + fst_sd),
                  width = 0.1,
                  data = all_summ) +
    geom_point(aes(colour = type, alpha = type, size = type)) +
    geom_path(aes(colour = type, alpha = type, group = iter)) +
    geom_point(aes(year, fst),  size = 2,
               data = all_summ, inherit.aes = FALSE) +
    geom_path(aes(year, fst), group = 1, data = all_summ, 
              inherit.aes = FALSE) +
    scale_alpha_manual(values = c(1, 0.5)) +
    scale_size_manual(values = c(2, 1)) +
    theme_minimal()
    
  p1
  
}