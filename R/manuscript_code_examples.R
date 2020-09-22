library(slimr)

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
               slimr_track_progress();
               sim%.%.SS$addSubpop(subpopID = "P1", 500)
             }),
  slim_block(10000,
             {
               slimr_output_full();
               sim.simulationFinished();
             })
  ) -> script_1
  

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
