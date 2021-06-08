library(slimr)
library(adegenet)

slim_script(
  slim_block(initialize(),
             {
               initializeMutationRate(1e-7);
               initializeMutationType("m1", 0.5, "f", 0.0);
               initializeGenomicElementType("g1", m1, 1.0);
               initializeGenomicElement(g1, 0, 99999);
               initializeRecombinationRate(1e-8);
             }),
  slim_block(1,
             {
               #slimr_track_progress();
               sim%.%.SS$addSubpop(subpopID = "p1", size = 500)
             }),
  slim_block(10000,
             {
               slimr_output_full();
               sim.simulationFinished();
             })
  ) -> script_1
suppressWarnings({

  
options(warn=-1)  
  
  
    
run_1 <- slim_run(script_1)
data_1 <- slim_output_genlight(run_1)

plot(data_1)




})

slim_script(
  slim_block(initialize(),
             {
               initializeMutationRate(1e-7);
               initializeMutationType("m1", 0.5, "f", 0.0);
               initializeGenomicElementType("g1", m1, 1.0);
               initializeGenomicElement(g1, 0, 99999);
               initializeRecombinationRate(1e-8);
             }),
  slim_block(1,
             {
               #slimr_track_progress();
               sim%.%.SS$addSubpop("p1", 500)
             }),
  slim_block(10000,
             {
               slimr_output(p1, "p1");
               slimr_output(sim.subpopulations.id, "ids")
               slimr_output_full();
               sim.simulationFinished();
             })
) -> script_2

test <- slim_run(script_2, capture_output = "|")


## slimr_script to generate simulation of three
## populations with migration matrix from R
library(slimr)

disp_mat <- matrix(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 
                     1, 2, 3, 1, 2, 3, 1, 2, 3, 
                     0.78, 0.1, 0.12, 0.01, 0.96, 
                     0.03, 0.33, 0.17, 0.50), 
                   ncol = 3)

slim_script(
  ## minimal initialize block
  slimr_block_init_minimal( 
    ## template mut rate, genome size, and recomb
    mut = slimr_template("mut_rate", 1e-7),
    gen = slimr_template("genome_size", 99999),
    recomb = slimr_template("recomb_rate", 1e-8)
  ),
  
  ## setup pops and migration rates in first gen
  slim_block(1, {
    
    for (i in 1:3) {
      sim.addSubpop(i, 100);
    }
    subpops = sim.subpopulations;
    ## pull in migration rate matrix from R
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
    ## output full sim output at gen 100000
    slimr_output(sim.outputFull(), 
                 "final_output");
  })
  
) -> script_1

script_1

slimr_script_render(script_1, data.frame(mut_rate = c(1e-6, 1e-8),
                                         genome_size = c(1e5, 1e6)))





pop_sim <- slim_script(
  
  slim_block(initialize(), {
    
    initializeMutationRate(slimr_template("mut_rate", 1e-6));
    initializeMutationType("m1", 0.5, "n", 0, slimr_template("selection_strength", 0.1));
    initializeGenomicElementType("g1", m1, 1.0);
    initializeGenomicElement(g1, 0, slimr_template("genome_size", 50000) - 1);
    initializeRecombinationRate(slimr_template("recomb_rate", 1e-8));
    initializeSex("A");
    defineConstant("abund", slimr_inline(pop_abunds, delay = TRUE));
    defineConstant("sample_these", slimr_inline(sample_these, delay = TRUE));
    
  }),
  slim_block(1, {
    
    init_pop = slimr_inline(init_popsize, delay = TRUE)
    
    ## set populations to initial size
    sim.addSubpop("p1", asInteger(init_pop[0]));
    sim.addSubpop("p2", asInteger(init_pop[1]));
    sim.addSubpop("p3", asInteger(init_pop[2]));
    
  }),
  
  slim_block(1, late(), {
    ## get starting population from a file which we will fill-in later
    sim.readFromPopulationFile(slimr_inline(starting_pop, delay = TRUE));
    ## migration on or off flags for pops 1-3 (using tag)
    p1.tag = 0;
    p2.tag = 0;
    p3.tag = 0;
  }),
  
  slim_block(1, 1000, late(), {
    
    ## update generation number
    gen = sim.generation %% 50
    if(gen == 0) {
      gen = 50
    }
    
    ## set population size to observed levels
    p1.setSubpopulationSize(asInteger(ceil(abund[0, gen - 1] * slimr_template("popsize_scaling", 100))));
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
      ind_sample = sample(sim.subpopulations.individuals, 50)
      slimr_output(ind_sample.genomes.output(), "pop_sample", do_every = 1);
    }
    
  }),
  
  slim_block(1000, late(), {
    sim.simulationFinished()
  })
  
)

result <- slim_run(slimr_script_render(pop_sim))


