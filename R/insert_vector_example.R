library(slimr)

vec <- replicate(3, runif(10), simplify = FALSE)

test_script <- slim_script(
  
  slim_block(initialize(), {
    initializeMutationRate(1e-07); 
    initializeMutationType("m1", 0.5, "f", 0); 
    initializeGenomicElementType("g1", m1, 1); 
    initializeGenomicElement(g1, 0, 1e+05 - 1); 
    initializeRecombinationRate(1e-08)
    
    defineConstant("test_vector", slimr_inline(vec[[slimr_template("index", 0)]], delay = TRUE));
    
  }),
  
  slim_block(1, {
    cat(test_vector)
  })
  
)

test_script

test_render <- slimr_script_render(test_script,
                                   template = data.frame(index = 1:3))

test_render

slim_run(test_render, simple_run = TRUE)