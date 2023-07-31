library(slimr)

slim_script(
  slim_block(initialize(),
             {
               ## set the overall mutation rate
               initializeMutationRate(1e-7); 
               ## m1 mutation type: neutral
               initializeMutationType("m1", 0.5, "f", 0.0);
               ## g1 genomic element type: uses m1 for all mutations
               initializeGenomicElementType("g1", m1, 1.0);
               ## uniform chromosome of length 100 kb
               initializeGenomicElement(g1, 0, 99999);
               ## uniform recombination along the chromosome
               initializeRecombinationRate(1e-8);
             }),
  slim_block_progress(10),
  slim_block(1,
             {
               sim%.%SS$addSubpop("p1", 500);
             }),
  slim_block(500000,
             {
               sim%.%SS$simulationFinished();
             })
) -> script_1

out <- slim_run(script_1)


slim_script(
  slim_block_init_minimal(),
  slim_block(1,
             {
               sim%.%SS$addSubpop("p1", 500)
             }),
  slim_block(10000,
             {
               sim%.%SS$simulationFinished();
             })
) -> script_1

script_1




slim_load_globals()

slim_script(
  slim_block_init_minimal(),
  slim_block(1,
             {
               sim$addSubpop("p1", 500)
             }),
  slim_block(10000,
             {
               sim$simulationFinished();
             })
) -> script_1

script_1

slim_script(
  slim_block_init_minimal(),
  slim_block(1,
             {
               sim$addSubpop("p1", 500)
             }),
  slim_block(10000,
             {
               sim$simulationFinished();
             })
) -> script_1







slim_block(1, { sim%.%Sp$addSubpop("p1", 500); })

slim_block(1, { sim$addSubpop("p1", 500); })












   slim_block(1, { sim%.%Sp$ad })




   
   
   
   


   
   
   
   
   
   
   slim_load_globals()
   slim_block(1, { sim$addSubpop() })




   
   
   
   
   
   
   
   
   

