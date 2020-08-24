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

results <- slim_run(script_1, progress = FALSE)

inds <- slim_output_genlight(results, "final_output")
plot(inds)




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
) -> script_3

script_3
