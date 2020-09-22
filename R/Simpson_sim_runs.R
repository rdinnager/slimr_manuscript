library(readr)
#library(dartR)
library(adegenet)
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

gen_meta <- gen_meta %>%
  mutate(three_pop = case_when(pop %in% c("MC", "SS", "WS") ~ "BR",
                               pop %in% c("FRN", "FRS") ~ "BL",
                               pop %in% c("KSE", "CS") ~ "TR")) %>%
  group_by(three_pop) %>%
  mutate(lon_3 = mean(lon),
         lat_3 = mean(lat))

pop(gen) <- gen_meta$three_pop
gen@other$ind.metrics <- gen_meta

gen <- gen[gen@other$ind.metrics$pop != "WS"]

#fst_overall <- gl.fst.pop(gen)

# plan(multisession(workers = 3))
# fst_by_year <- future_map(unique(gen_meta$year),
#                    ~gl.fst.pop(gen[gen@other$ind.metrics$year == .x]))


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

sample_times <- c("2006" = 40, "2007" = 45, "2008" = 49)

plot(pop_abunds[1, ], type = "l", col = "red")
lines(pop_abunds[2, ], col = "blue")
lines(pop_abunds[3, ], col = "green")
abline(v = sample_times)

max_generations <- 1000

num_cycles <- max_generations / 50

years_of_sim <- 14 * num_cycles

## setup generations to sample
## we will just sample the three years corresponding to the data, but do it during the last three cycles, as "technical" replicates
sample_these <- purrr::map(c(17:19),
                           ~50*.x + sample_times) %>%
  purrr::reduce(union)

plot(rep(pop_abunds[1, ], 20), type = "l", col = "red", xlim = c(800, 1000))
lines(rep(pop_abunds[2, ], 20), col = "blue")
lines(rep(pop_abunds[3, ], 20), col = "green")
abline(v = sample_these, col = "yellow")

## first get data only from 2008 (panmictic low rainfall, year after big rainfall)
## fill-in na values with a bernoulli draw with prob based based on non-missing data
gen_2008 <- gen[gen@other$ind.metrics$year == 2008, ]
snp_mat <- as.matrix(gen_2008)
#mostly_nas <- apply(snp_mat, 2, function(x) sum(is.finite(x)) < 2)
nas <- apply(snp_mat, 2, function(x) !any(!is.finite(x)))
snp_mat[ , !nas] <- apply(snp_mat[ , !nas], 2, function(x) {
  tab <- table(x[is.finite(x)]);
  x[!is.finite(x)] <- as.integer(sample(names(tab), sum(!is.finite(x)), replace = TRUE, prob = tab / sum(tab)));
  x
})

genome_size <- 500000

sexes <- as.character(gen_2008@other$ind.metrics$Sex)
sexes[sexes == ""] <- sample(c("M", "F"), sum(sexes == ""), replace = TRUE)

## sex ratio of sample is skewed, so replace actually sexes with perfectly even alternation
## otherwise simulation crashes because it is unable to sample enough females (need to look into this issue more)
sexes[] <- c("M", "F")
init_pop <- slimr::slim_make_pop_input(snp_mat, sim_gen = 1,
                                       ind_pops = gen_2008@other$ind.metrics$three_pop,
                                       ind_sex = sexes,
                                       mut_pos = sample.int(genome_size, nLoc(gen_2008)),
                                       mut_prev = apply(as.matrix(snp_mat), 2, sum))

#wsl_file <- slimr:::convert_to_wsl_path(init_pop)
wsl_file <- init_pop

init_popsize <- c(table(gen_2008@other$ind.metrics$three_pop))

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
                  fst = fst[!isNAN(fst)]; ## exclude muts where mean_p is 0.0 or 1.0
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
    defineConstant("sample_these", slimr_inline(sample_these));
    
  }),
  slim_block(1, {
    
    init_pop = slimr_inline(init_popsize)
    
    ## set populations to initial size
    sim.addSubpop("p1", asInteger(init_pop[0]));
    sim.addSubpop("p2", asInteger(init_pop[1]));
    sim.addSubpop("p3", asInteger(init_pop[2]));
    
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
    
    if(any(match(sample_these, sim.generation) >= 0)) {
      slimr_output(sim.outputFull(), "pop_sample", do_every = 1, send_to = "file",
                  file_name = "..file_name..", format = "csv");
    }
    
  }),
  
  slim_block(!!max_generations, {
    sim.simulationFinished();
  })
  
)

pop_sim

test_default <- slimr_script_render(pop_sim, template = list(popsize_scaling = 100, file_name = "/mnt/f/Projects/slimr_manuscript/test.csv.gz"))
test_default

test_run <- slim_run(test_default, capture_output = "file")

#################### setup parameters ##################

n_sims <- 10000
params_df <- dplyr::tibble(mut_rate = runif(n_sims, 1e-7, 1e-4),
                           selection_strength = runif(n_sims, 0.01, 0.2),
                           genome_size = genome_size,
                           recomb_rate = runif(n_sims, 0.5 * 1e-8, 2 * 1e-8),
                           popsize_scaling = runif(n_sims, 10, 1000),
                           abund_threshold = runif(n_sims, 0, 12),
                           migration_rate = runif(n_sims, 0, 0.5)) %>%
  dplyr::mutate(file_name = paste0("/mnt/f/Projects/slimr_manuscript/data/sims/sim_", 1:dplyr::n(), "csv.gz"))

readr::write_csv(params_df, "data/sims_params_df.csv")

many_sims <- slimr_script_render(pop_sim, template = params_df)

readr::write_rds(many_sims, "data/sims_rendered_scripts.rds", compress = "gz")

future::plan(future::multisession(workers = 6))
sim_res <- slim_run(many_sims, capture_output = "|", parallel = TRUE, progress = TRUE)
