FirstGeneration = function(maxDepthIni, populationSize){
  method = "grow"
  
  population = vector('list', populationSize) 
  
  for (i in 1:populationSize){
    
    method = ifelse(method == "grow", "full", "grow")
    actualDepth = ifelse(method == "grow", sample(0:maxDepthIni, 1), maxDepthIni)
    population[[i]] = CreateIndividual(actualDepth, method, restricted = FALSE)
    
  }
  
  population
  
}