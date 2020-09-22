library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(mapview)
library(purrr)
library(conflicted)
library(future)
library(furrr)
library(lubridate)
library(zoo)
library(ggplot2)
library(patchwork)
library(sf)
library(adegenet)

library(ggsfc)

library(slimr)

conflict_prefer("filter", "dplyr")

gen <- readr::read_rds("data/herm.rdata")
abund <- readr::read_csv("data/mammal_captures.csv")

gen_meta <- gen@other$ind.metrics

sites <- tribble(~pop, ~SiteName,
                 "CS", "Carlo",
                 "FRN", "Field River North",
                 "FRS", "Field River South",
                 "KSE", "Kunnamuka Swamp East",
                 "MC", "Main Camp",
                 "SS", "South Site",
                 "WS", "Way Site"
)

gen_meta <- gen_meta %>%
  mutate(pop = case_when(pop == "FR" ~ substr(as.character(SiteGrid), 1, 3),
                         pop == "KS" ~ "KSE",
                         TRUE ~ as.character(pop))) %>%
  left_join(sites)

coords <- gen_meta %>%
  group_by(pop) %>%
  summarise(lon = mean(lon),
            lat = mean(lat))

sites1 <- sf::st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
#mapview(coords, xcol = "lon", ycol = "lat", label = coords$pop)
mapview(sites1, label = coords$pop)

gen_meta <- gen_meta %>%
  mutate(three_pop = case_when(pop %in% c("MC", "SS", "WS") ~ "BR",
                               pop %in% c("FRN", "FRS") ~ "BL",
                               pop %in% c("KSE", "CS") ~ "TR")) %>%
  group_by(three_pop) %>%
  mutate(lon_3 = mean(lon),
         lat_3 = mean(lat))

coords_3 <- gen_meta %>%
  group_by(three_pop) %>%
  summarise(lon = mean(lon_3),
            lat = mean(lat_3))

sites_3 <- sf::st_as_sf(coords_3, coords = c("lon", "lat"), crs = 4326)
#mapview(coords_3, xcol = "lon", ycol = "lat", label = coords_3$three_pop)
mapview(sites_3, label = coords_3$three_pop)

pop(gen) <- gen_meta$three_pop
gen@other$ind.metrics <- gen_meta

gen <- gen[gen@other$ind.metrics$pop != "WS"]

abund_by_site <- abund %>%
  left_join(sites) %>%
  mutate(three_pop = case_when(pop %in% c("MC", "SS", "WS") ~ "BR",
                               pop %in% c("FRN", "FRS") ~ "BL",
                               pop %in% c("KSE", "CS") ~ "TR")) %>%
  drop_na(three_pop) %>%
  separate(MonthYear, c("Month", "year"), "\\.") %>%
  mutate(Month = sapply(Month, function(x) which(month.abb == x))) %>%
  mutate(Year_month = date(as.yearmon(paste(Year, Month, sep = "-")))) %>%
  group_by(Year_month, three_pop) %>%
  summarise(abund = mean(`Pseudomys hermannsburgensis`)) %>%
  arrange(Year_month) %>%
  group_by(Year_month) %>%
  mutate(abund_mean = mean(abund))

ggplot(abund_by_site, aes(Year_month, abund)) +
  geom_path(aes(colour = three_pop)) +
  scale_x_date(date_labels = "%b - %Y",
               date_breaks = "6 months") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))


######### simulation ###########
## params:
##   mut_rate: mutation rate of simulated population
##   genome_size: size of simulated genome
##   selection_strength: mean strength of selection specified as a standard deviation of selection coeffs 
##   migration_rates: rate of migration amongst three pops when abundance is high, must be between 0 and 1
##   abund_threshold: abundance (before scaling) above which migration between populations is "turned" on,
##     makes sense this should be between 3 (5 better) min and 12 max
##   recomb_rate: recombination rate, a nuisance parameter
##   popsize_scaling: multiply observed abundances by this value to get total subpop size
##  

pop_abunds <- abund_by_site %>%
  filter(Year_month < "2009-10-01" & Year_month > "1995-03-01") %>%
  select(-abund_mean) %>%
  pivot_wider(names_from = Year_month, values_from = abund) %>%
  as.matrix()

pop_abunds <- rbind(approx(pop_abunds[1, ], n = 50)$y,
                    approx(pop_abunds[2, ], n = 50)$y,
                    approx(pop_abunds[3, ], n = 50)$y)

## add small constant to any zeroes because setting a population size to zero in SLiM destroys that population,
## requiring you to recreate it, which is annoying. Plus we don't think any population actually goes exinct ever in
## real life
pop_abunds[pop_abunds == 0] <- 0.02

plot(pop_abunds[2, ], type = "l", col = "blue")
lines(pop_abunds[1, ], col = "red")
lines(pop_abunds[3, ], col = "green")

max_generations <- 1000

## first get data only from 2008 (panmictic low rainfall, year after big rainfall)
## fill-in na values with a multinomial draw with prob based based on non-missing data (cats are 0, 1, 2)
gen_2008 <- gen[gen@other$ind.metrics$year == 2008, ]
snp_mat <- as.matrix(gen_2008)
#mostly_nas <- apply(snp_mat, 2, function(x) sum(is.finite(x)) < 2)
nas <- apply(snp_mat, 2, function(x) !any(!is.finite(x)))
snp_mat[ , !nas] <- apply(snp_mat[ , !nas], 2, function(x) {
    tab <- table(x[is.finite(x)]);
    x[!is.finite(x)] <- as.integer(sample(names(tab), sum(!is.finite(x)), replace = TRUE, prob = tab / sum(tab)));
    x
  })

dim(snp_mat)

genome_size <- 500000

sexes <- as.character(gen_2008@other$ind.metrics$Sex)
sexes[sexes == ""] <- sample(c("M", "F"), sum(sexes == ""), replace = TRUE)
sexes[] <- c("M", "F")

init_pop <- slimr::slim_make_pop_input(snp_mat, sim_gen = 1,
                                       ind_pops = gen_2008@other$ind.metrics$three_pop,
                                       ind_sex = sexes,
                                       mut_pos = sample.int(genome_size, nLoc(gen_2008)),
                                       mut_prev = apply(as.matrix(snp_mat), 2, sum))

wsl_file <- slimr:::convert_to_wsl_path(init_pop)
#wsl_file <- init_pop

init_popsize <- matrix(table(gen_2008@other$ind.metrics$three_pop), nrow = 1)

pop_sim <- slim_script(
  
  slim_function("o<Subpopulation>$ subpop1", "o<Subpopulation>$ subpop2",
                name = "calcFST",
                return_type = "f$", body = {
                  ## Calculate the FST between two subpopulations
                  p1_p = sim.mutationFrequencies(subpop1);
                  p2_p = sim.mutationFrequencies(subpop2);
                  mean_p = (p1_p + p2_p) / 2.0;
                  H_t = 2.0 * mean_p * (1.0 - mean_p);
                  H_s = p1_p * (1.0 - p1_p) + p2_p * (1.0 - p2_p);
                  fst = 1.0 - H_s/H_t;
                  fst = fst[isFinite(fst)]; ## exclude muts where mean_p is 0.0 or 1.0
                  return(mean(fst));
                }),
  
  slim_block(initialize(), {
    
    initializeMutationRate(slimr_template("mut_rate", 1e-6));
    initializeMutationType("m1", 0.5, "n", 0, slimr_template("selection_strength", 0.1));
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, 0, slimr_template("genome_size", !!genome_size) - 1);
    initializeRecombinationRate(slimr_template("recomb_rate", 1e-8));
    initializeSex("A");
    defineConstant("abund", slimr_inline(pop_abunds));
    
  }),
  slim_block(1, {
    
    init_pop = slimr_inline(init_popsize)
    
    ## set populations to initial size
    sim.addSubpop("p1", asInteger(init_pop[0, 0]));
    sim.addSubpop("p2", asInteger(init_pop[0, 1]));
    sim.addSubpop("p3", asInteger(init_pop[0, 2]));
    
  }),
  
  slim_block(1, late(), {
    sim.readFromPopulationFile(!!wsl_file);
    ## migration on or off flags for pops 1-3 (using tag)
    p1.tag = 0;
    p2.tag = 0;
    p3.tag = 0;
  }),
  
  slim_block(1, !!max_generations, late(), {
    
    ## update generation number
    gen = sim.generation %% 50
    if(gen == 0) {
      gen = 50
    }
    
    ## set population size to observed levels
    p1.setSubpopulationSize(asInteger(ceil(abund[0, gen - 1] * slimr_template("popsize_scaling", 1))));
    p2.setSubpopulationSize(asInteger(ceil(abund[1, gen - 1] * ..popsize_scaling..)));
    p3.setSubpopulationSize(asInteger(ceil(abund[2, gen - 1] * ..popsize_scaling..)));
    
    ## increase migration when above abundance threshold
    if(p1.tag == 0 & abund[0, gen - 1] > slimr_template("abund_threshold", 5)) {
      p2.setMigrationRates(p1, slimr_template("migration_rate", 0))
      p3.setMigrationRates(p1, ..migration_rate..)
      p1.tag = 1;
    } 
    if(p1.tag == 1 & abund[0, gen - 1] <= ..abund_threshold..) {
      p2.setMigrationRates(p1, 0)
      p3.setMigrationRates(p1, 0)
      p1.tag = 0;
    }
    
    if(p2.tag == 0 & abund[1, gen - 1] > ..abund_threshold..) {
      p1.setMigrationRates(p2, ..migration_rate..)
      p3.setMigrationRates(p2, ..migration_rate..)
      p2.tag = 1;
    } 
    if(p2.tag == 1 & abund[1, gen - 1] <= ..abund_threshold..) {
      p1.setMigrationRates(p2, 0)
      p3.setMigrationRates(p2, 0)
      p2.tag = 0;
    }    
    
    if(p3.tag == 0 & abund[2, gen - 1] > ..abund_threshold..) {
      p1.setMigrationRates(p3, ..migration_rate..)
      p2.setMigrationRates(p3, ..migration_rate..)
      p3.tag = 1;
    } 
    if(p3.tag == 1 & abund[2, gen - 1] <= ..abund_threshold..) {
      p1.setMigrationRates(p3, 0)
      p2.setMigrationRates(p3, 0)
      p3.tag = 0;
    }
    
    muts = sim.mutations[sim.mutationCounts(NULL) > 1]
    
   
    ## output data to visualise
    slimr_output(c(calcFST(p1, p2), calcFST(p1, p3), calcFST(p2, p3)), "fsts");
    slimr_output(sim.subpopulations.individualCount, "vis_data_N");
    slimr_output(sapply(sim.subpopulations, "sim.mutationFrequencies(applyValue, muts);"), "vis_data_freq");
    slimr_output(muts.mutationType, "vis_data_types");
    slimr_output(muts.position, "vis_data_pos");
    slimr_output(muts.selectionCoeff, "vis_data_sel")
    
  })
  
)

pop_sim

scaling <- 100
test_default <- slimr_script_render(pop_sim, template = list(popsize_scaling = scaling))
test_default

# test_run <- slim_run(test_default, capture_output = "|")
# out_dat <- test_run$output_data

extract_vis_data <- function(out_dat) {
  
  gens <- table(out_dat$generation)
  
  if(any(gens == 6)) {
  
    curr_gen <- max(as.numeric(names(gens))[gens == 6])
    
    out_dat <- out_dat %>%
      dplyr::filter(generation == curr_gen)
    
    N_dat <- out_dat %>%
      dplyr::filter(name == "vis_data_N") %>%
      dplyr::select(data) %>%
      dplyr::mutate(data = stringr::str_split(data, " ")) %>%
      tidyr::unnest_longer(data, values_to = "N", indices_to = "pop") %>%
      dplyr::mutate(N = as.numeric(N), generation = curr_gen)
    
    n_pop <- max(N_dat$pop)
    
    type_dat <- out_dat %>%
      dplyr::filter(name == "vis_data_types") %>%
      dplyr::select(data) %>%
      dplyr::mutate(data = stringr::str_split(data, " ")) %>%
      tidyr::unnest_longer(data, values_to = "type", indices_to = "index") %>%
      dplyr::mutate(type = stringr::str_match(type, "MutationType<(.*?)>")[, 2])
    
    pos_dat <- out_dat %>%
      dplyr::filter(name == "vis_data_pos") %>%
      dplyr::select(data) %>%
      dplyr::mutate(data = stringr::str_split(data, " ")) %>%
      tidyr::unnest_longer(data, values_to = "pos", indices_to = "index") %>%
      dplyr::mutate(pos = as.integer(pos) + 1L)
    
    sel_dat <- out_dat %>%
      dplyr::filter(name == "vis_data_sel") %>%
      dplyr::select(data) %>%
      dplyr::mutate(data = stringr::str_split(data, " ")) %>%
      tidyr::unnest_longer(data, values_to = "sel", indices_to = "index")
    
    n_mut <- max(pos_dat$index)
    
    freq_dat <- out_dat %>%
      dplyr::filter(name == "vis_data_freq") %>%
      dplyr::select(data) %>%
      dplyr::mutate(data = stringr::str_split(data, " ")) %>%
      tidyr::unnest_longer(data, values_to = "freq", indices_to = "index") %>%
      dplyr::mutate(pop = rep(seq_len(n_pop), each = n_mut),
                    index = index - (pop - 1) * n_mut,
                    freq = as.numeric(freq)) %>%
      dplyr::left_join(type_dat, by = "index") %>%
      dplyr::left_join(pos_dat, by = "index") %>%
      dplyr::left_join(sel_dat, by = "index")
    
    fst_dat <- out_dat %>%
      dplyr::filter(name == "fsts") %>%
      dplyr::select(data) %>%
      dplyr::mutate(data = stringr::str_split(data, " ")) %>%
      tidyr::unnest_longer(data, values_to = "fst", indices_to = "index") %>%
      dplyr::mutate(fst = as.numeric(fst)) %>%
      dplyr::summarise(fst = mean(fst), generation = curr_gen)
  
    list(fst = fst_dat, N = N_dat, freq = freq_dat)
  } else {
    return(NULL)
  }
  
}

make_vis_function <-  function(ylims = c(N = NA, fst = NA)) {
  if(is.na(ylims["N"])) {
    max_pop <- max(pop_abunds * scaling) + 1 * scaling
  } else {
    max_pop <- ylims["N"] 
  }
  if(is.na(ylims["fst"])) {
    max_fst <- 1
  } else {
    max_fst <- ylims["fst"]
  }
  fst_dat <- dplyr::tibble(fst = 0, generation = 0)
  N_dat <- dplyr::tibble(N = 0, pop = c(1L, 2L, 3L), generation = 0)
  N_base <- N_dat %>%
    ggplot(., aes(generation, N)) +
    scale_x_continuous(limits = c(0, max_generations)) +
    scale_y_continuous(limits = c(0, max_pop)) +
    scale_colour_brewer(name = "population", type = "qual") +
    theme_minimal() +
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  fst_base <- fst_dat %>%
    ggplot(., aes(generation, fst)) +
    scale_x_continuous(limits = c(0, max_generations)) +
    scale_y_continuous(limits = c(0, max_fst)) +
    scale_colour_brewer(name = "population", type = "qual") +
    theme_minimal()
  function(data) {
    grDevices::dev.hold()
    all_dat <- extract_vis_data(data)
    #print(is.null(all_dat))
    if(is.null(all_dat)) {
      grDevices::dev.flush()
      return(invisible(NULL))
    }
    fst_dat <<- dplyr::bind_rows(fst_dat, all_dat$fst)
    N_dat <<- dplyr::bind_rows(N_dat, all_dat$N)
    p1 <- {N_base +
        geom_path(aes(colour = factor(pop)), data = N_dat)} +
    {fst_base +
        geom_path(data = fst_dat)} +
      plot_layout(ncol = 1)
    
    #print("is it here?")
    plot(p1)
    
    invisible(grDevices::dev.flush())
    
  }
}

vis_fun <- make_vis_function()

# grDevices::X11()
# grDevices::dev.hold()
# vis_fun(out_dat)
# grDevices::dev.flush()

grDevices::X11(width = 14, height = 8)

sim_res <- slim_run(test_default, capture_output = "|", 
                    callbacks = list(vis_fun))




test_2 <- slimr_script_render(pop_sim, template = list(popsize_scaling = scaling,
                                                             migration_rate = 0.5))

vis_fun <- make_vis_function(ylims = c(fst = 0.2, N = NULL))

grDevices::X11(width = 14, height = 8)

sim_res <- slim_run(test_2, capture_output = "|", 
                    callbacks = list(vis_fun))



test_3 <- slimr_script_render(pop_sim, template = list(popsize_scaling = scaling,
                                                       migration_rate = 0.5,
                                                       abund_threshold = 0,
                                                       selection_strength = 0.000000001))

vis_fun <- make_vis_function(ylims = c(fst = 0.12, N = NULL))

grDevices::X11(width = 14, height = 8)

sim_res <- slim_run(test_3, capture_output = "|", 
                    callbacks = list(vis_fun))


test_4 <- slimr_script_render(pop_sim, template = list(popsize_scaling = scaling,
                                                       migration_rate = 0.5,
                                                       abund_threshold = 0,
                                                       selection_strength = 0.2))

vis_fun <- make_vis_function(ylims = c(fst = 0.1, N = NULL))

grDevices::X11(width = 14, height = 8)

sim_res <- slim_run(test_4, capture_output = "|", 
                    callbacks = list(vis_fun))




test_5 <- slimr_script_render(pop_sim, template = list(popsize_scaling = scaling,
                                                       migration_rate = 0.25))

vis_fun <- make_vis_function(ylims = c(fst = 0.4, N = NULL))

#grDevices::X11(width = 14, height = 8)

sim_res <- slim_run(test_5, capture_output = "|", 
                    callbacks = list(vis_fun),
                    new_grdev = TRUE)



out_dat <- sim_res$output_data

data <- extract_vis_data(sim_res$output_data)

extract_vis_data <- function(out_dat) {
  
  gens <- table(out_dat$generation)
  
  if(any(gens == 6)) {
    
    curr_gen <- max(as.numeric(names(gens))[gens == 6])
    
    out_dat <- out_dat %>%
      dplyr::filter(generation == curr_gen)
    
    N_dat <- out_dat %>%
      dplyr::filter(name == "vis_data_N") %>%
      dplyr::select(data) %>%
      dplyr::mutate(data = stringr::str_split(data, " ")) %>%
      tidyr::unnest_longer(data, values_to = "N", indices_to = "pop") %>%
      dplyr::mutate(N = as.numeric(N), generation = curr_gen)
    
    n_pop <- max(N_dat$pop)
    
    type_dat <- out_dat %>%
      dplyr::filter(name == "vis_data_types") %>%
      dplyr::select(data) %>%
      dplyr::mutate(data = stringr::str_split(data, " ")) %>%
      tidyr::unnest_longer(data, values_to = "type", indices_to = "index") %>%
      dplyr::mutate(type = stringr::str_match(type, "MutationType<(.*?)>")[, 2])
    
    pos_dat <- out_dat %>%
      dplyr::filter(name == "vis_data_pos") %>%
      dplyr::select(data) %>%
      dplyr::mutate(data = stringr::str_split(data, " ")) %>%
      tidyr::unnest_longer(data, values_to = "pos", indices_to = "index") %>%
      dplyr::mutate(pos = as.integer(pos) + 1L)
    
    sel_dat <- out_dat %>%
      dplyr::filter(name == "vis_data_sel") %>%
      dplyr::select(data) %>%
      dplyr::mutate(data = stringr::str_split(data, " ")) %>%
      tidyr::unnest_longer(data, values_to = "sel", indices_to = "index") %>%
      dplyr::mutate(sel = as.numeric(sel))
    
    n_mut <- max(pos_dat$index)
    
    freq_dat <- out_dat %>%
      dplyr::filter(name == "vis_data_freq") %>%
      dplyr::select(data) %>%
      dplyr::mutate(data = stringr::str_split(data, " ")) %>%
      tidyr::unnest_longer(data, values_to = "freq", indices_to = "index") %>%
      dplyr::mutate(pop = rep(seq_len(n_pop), each = n_mut),
                    index = index - (pop - 1) * n_mut,
                    freq = as.numeric(freq)) %>%
      dplyr::left_join(type_dat, by = "index") %>%
      dplyr::left_join(pos_dat, by = "index") %>%
      dplyr::left_join(sel_dat, by = "index")
    
    fst_dat <- out_dat %>%
      dplyr::filter(name == "fsts") %>%
      dplyr::select(data) %>%
      dplyr::mutate(data = stringr::str_split(data, " ")) %>%
      tidyr::unnest_longer(data, values_to = "fst", indices_to = "index") %>%
      dplyr::mutate(fst = as.numeric(fst)) %>%
      dplyr::summarise(fst = mean(fst), generation = curr_gen)
    
    list(fst = fst_dat, N = N_dat, freq = freq_dat)
  } else {
    return(NULL)
  }
  
}

make_vis_function <-  function(ylims = c(N = NA, fst = NA)) {
  if(is.na(ylims["N"])) {
    max_pop <- max(pop_abunds * scaling) + 1 * scaling
  } else {
    max_pop <- ylims["N"] 
  }
  if(is.na(ylims["fst"])) {
    max_fst <- 1
  } else {
    max_fst <- ylims["fst"]
  }
  fst_dat <- dplyr::tibble(fst = 0, generation = 0)
  N_dat <- dplyr::tibble(N = 0, pop = c(1L, 2L, 3L), generation = 0)
  N_base <- N_dat %>%
    ggplot(., aes(generation, N)) +
    scale_x_continuous(limits = c(0, max_generations)) +
    scale_y_continuous(limits = c(0, max_pop)) +
    scale_colour_brewer(name = "population", type = "qual") +
    theme_minimal() +
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  fst_base <- fst_dat %>%
    ggplot(., aes(generation, fst)) +
    scale_x_continuous(limits = c(0, max_generations)) +
    scale_y_continuous(limits = c(0, max_fst)) +
    scale_colour_brewer(name = "population", type = "qual") +
    theme_minimal()
  function(data) {
    #grDevices::dev.hold()
    all_dat <- extract_vis_data(data)
    #print(is.null(all_dat))
    if(is.null(all_dat)) {
      #grDevices::dev.flush()
      return(invisible(NULL))
    }
    fst_dat <<- dplyr::bind_rows(fst_dat, all_dat$fst)
    N_dat <<- dplyr::bind_rows(N_dat, all_dat$N)
    p1 <- {N_base +
        geom_path(aes(colour = factor(pop)), data = N_dat)} +
      {fst_base +
          geom_path(data = fst_dat)} +
      plot_layout(ncol = 1)
    
    #print("is it here?")
    plot(p1)
    
    #invisible(grDevices::dev.flush())
    return(invisible(NULL))
    
  }
}

vis_fun <- make_vis_function()





test_fst <- out_dat %>%
  dplyr::filter(name == "fsts")

dat <- extract_vis_data(test_run$output_data)

limits <- c(1, genome_size)
vis_flowsnake <- function(dat, group, pos, colour, limits) {
  dat_split <- dat %>%
    dplyr::group_by(pop) %>%
    dplyr::group_split()
  
  flowsnake <- ggsfc::sfcurve("flowsnake", 4)
  
  plot(as.list(flowsnake), xlim = c(-1, 3), ylim = c(-1, 3), col = "green", asp = 1, type = "l")
  flowsnake2 <- flowsnake
  flowsnake2$x <- flowsnake2$x + 1.66666
  flowsnake2$y <- flowsnake2$y + 0.52
  lines(flowsnake2, col = "red")
  
  flowsnake3 <- flowsnake
  flowsnake3$x <- flowsnake3$x + 0.4
  flowsnake3$y <- flowsnake3$y + 1.79
  lines(flowsnake3, col = "blue")
  
  ang <- asin(sqrt(3) / (5*sqrt(7)))
  ang <- -atan2(flowsnake$y[1] * sqrt(3), flowsnake$x[1])
  flowsnake_rot <- as.matrix(flowsnake) %*% matrix(c(cos(ang), sin(ang), -sin(ang), cos(ang)), nrow = 2)
  
  plot(flowsnake_rot, type = "l", ylim = c(-1, 3), xlim = c(-1, 3), col = "green", asp = 1)
  flowsnake_up <- flowsnake_rot
  flowsnake_up[ , 2] <- flowsnake_up[ , 2] + 1.8
  lines(flowsnake_up, col = "red")
  
  flowsnake_right <- flowsnake_rot
  flowsnake_right[ , 2] <- flowsnake_right[ , 2] + 0.93
  flowsnake_right[ , 1] <- flowsnake_right[ , 1] + 1.48
  lines(flowsnake_right, col = "blue")
  
  flowsnakes_base <- purrr::map(1:3)
  
  plot(flowsnake, col = rainbow(nrow(flowsnake)))
}


########### write a callback to visualize results ###############
data <- test_run$output_data[100, ]
calc_fsts <- function(data) {
  tt <- slim_outputFull_extract(data, "genomes")
  gen <- slim_output_genlight(data, "pop_sample")
}



library(gapminder)
library(ggplot2)
makeplot <- function(){
  datalist <- split(gapminder, gapminder$year)
  lapply(datalist, function(data){
    p <- ggplot(data, aes(gdpPercap, lifeExp, size = pop, color = continent)) +
      scale_size("population", limits = range(gapminder$pop)) + geom_point() + ylim(20, 90) +
      scale_x_log10(limits = range(gapminder$gdpPercap)) + ggtitle(data$year) + theme_classic()
    print(p)
    grDevices::dev.copy(which = disp)
    grDevices::dev.flush()
    grDevices::dev.set(grDevices::dev.prev())
  })
}

# Play 1 plot per sec, and use an interpolation filter to convert into 10 fps
video_file <- file.path(tempdir(), 'output.mp4')

grDevices::X11(width = 14, height = 8)
disp <- grDevices::dev.cur()

av_capture_graphics(makeplot(), video_file, 1280, 720, res = 144, vfilter = 'framerate=fps=10')
av::av_media_info(video_file)

utils::browseURL(video_file)



################ reprex for av issue ############
library(av)

make_vid <- function(file = tempfile(fileext = ".mp4")) {
  plot_fun <- function(num) plot(cbind(runif(num), runif(num)))
  args <- list(expr = rlang::expr(plot_fun(100)),
               output = file)
  do.call(av::av_capture_graphics, args)
}

make_vid()