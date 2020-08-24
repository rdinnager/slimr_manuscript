library(slimr)
library(adegenet)
library(dartR)
library(igraph)
library(ggplot2)
library(dplyr)
library(purrr)
library(future)
library(furrr)

set.seed(123456)


disp_mat <- matrix(c(1, 1, 0.78, 
                     1, 2, 0.10, 
                     1, 3, 0.12, 
                     2, 1, 0.01, 
                     2, 2, 0.96, 
                     2, 3, 0.03, 
                     3, 1, 0.33, 
                     3, 2, 0.17, 
                     3, 3, 0.5),
                   ncol = 3, byrow = TRUE)

slim_script(
  slim_block(initialize(), {
    setSeed(123456);
    initializeMutationRate(2e-6);
    initializeMutationType("m1", 0.5, "f", 0.0);
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, 0, 99999);
    initializeRecombinationRate(1e-8);
  }),
  slim_block(1, {
    for (i in 1:3) {
      sim.addSubpop(i, 100);
    }
    subpops = sim.subpopulations;
    
    disp_mat = slimr_inline(disp_mat)
    
    for (line in seqLen(nrow(disp_mat))) {
      i = asInteger(disp_mat[line, 0]);
      j = asInteger(disp_mat[line, 1]);
      m = asFloat(disp_mat[line, 2]);
      if (i != j) {
        p_i = subpops[subpops.id == i];
        p_j = subpops[subpops.id == j];
        p_j.setMigrationRates(p_i, m);
      }
    }
  }),
  slim_block(100000, late(), {
    slimr_output(sim.outputFull(), "final_output");
  })
) -> script_1

script_1



########## run script ###############

results <- slim_run(script_1, progress = FALSE)

########### convert to genlight ###########

inds <- slim_output_genlight(results, "final_output")
plot(inds)

########### calculate Fst ##############

p_fst <- dartR::gl.fst.pop(inds, nboots = 1)
p_fst

########### compare with average migration rates ###########

overall_mig_rates <- disp_mat %>%
  as.data.frame() %>%
  igraph::graph_from_data_frame() %>%
  igraph::as.undirected(mode = "collapse", edge.attr.comb = "mean") %>%
  igraph::as_adjacency_matrix(type = "lower", attr = "V3", sparse = FALSE)

diag(overall_mig_rates) <- NA
overall_mig_rates[upper.tri(overall_mig_rates)] <- NA

p_fst
overall_mig_rates

plot(as.vector(overall_mig_rates), as.vector(p_fst))

########## make simpler script to examine Fst further ###########

## put modified slim script here: 
## change num gens to 1000

slim_script(
  slim_block(initialize(), {
    initializeMutationRate(2e-6);
    initializeMutationType("m1", 0.5, "f", 0.0);
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, 0, 99999);
    initializeRecombinationRate(1e-8);
  }),
  slim_block(1, {
    sim.addSubpop("p1", slimr_template("pop_size", default = 100));
    sim.addSubpop("p2", ..pop_size..);
    p1.setMigrationRates(p2, slimr_template("m", default = 0.1));
    p2.setMigrationRates(p1, ..m..);
  }),
  slim_block(10000, late(), {
    slimr_output(sim.outputFull(), "final_output");
  })
) -> script_2

script_2

########## render template ###########
## edit below to render
slimr_script_render(script_2, template = list(pop_size = 100,
                                              m = 0.1))

## you can also use a data.frame
n_iters <- 100
param_df <- data.frame(pop_size = sample(c(5:25, 75:100), n_iters, replace = TRUE),
                       m = runif(n_iters, 0, 0.707)^2)

param_df

plot(param_df$pop_size, param_df$m)

script_list <- slimr_script_render(script_2,
                                   template = param_df)


script_list

########### run all scripts in parallel ##############

plan(multisession)

all_results <- slim_run(script_list, parallel = TRUE, progress = TRUE)


########## check for errors ######

errors <- map_lgl(all_results,
                  ~length(.x$error) > 0)
sum(errors)


######### convert all to genlight ##########

all_gl <- future_map(all_results,
                     ~suppressMessages(try(slim_output_genlight(.x, "final_output"))),
                     .progress = TRUE)

######### any conversion errors #############

errored <- which(sapply(all_gl, function(x) class(x) == "try-error"))
errored

map(errored,
    ~all_results[[.x]]$output_data$data)

if(length(errored) > 0) {
  all_gl <- all_gl[-errored]
}

######### calculate all Fsts ##############

all_fsts <- future_map_dbl(all_gl,
                           ~dartR::gl.fst.pop(.x, nboots = 1)[2, 1],
                           .progress = TRUE)


######## visualize results ##############
## set negative Fsts to zero
all_fsts[all_fsts < 0] <- 0.0

hist(all_fsts, breaks = 100)

## add results to parameter data
if(length(errored) > 0) {
  param_df <- param_df[-errored, ]
}

fst_dat <- param_df %>%
  mutate(fst = all_fsts,
         pop_size_st = (pop_size - mean(pop_size)) / sd(pop_size),
         m_st = (m - mean(m)) / sd(m))

## run model
mod <- lm(fst ~ pop_size_st*m_st, data = fst_dat)
summary(mod)


## plot it!
## split data into 2 population size bins and plot each with separate fitted lines
fst_dat <- fst_dat %>%
  mutate(ps_bin = ifelse(pop_size <= 25, "Pop. Size 5 to 25", "Pop. Size 75 to 100")) %>%
  mutate(ps_bin = factor(ps_bin, c("Pop. Size 5 to 25", "Pop. Size 75 to 100"), ordered = TRUE))

ggplot(fst_dat, aes(m, fst)) +
  geom_point(aes(colour = ps_bin)) +
  geom_smooth(aes(colour = ps_bin, fill = ps_bin)) +
  scale_x_sqrt() +
  theme_minimal()



