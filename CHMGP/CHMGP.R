source('CHMGP_GlobalVar.R')
source('CHMGP_FirstGeneration.R')
source('CHMGP_Individual.R')
source('CHMGP_ComputationFunctions.R')
source('CHMGP_ComputeIndividual.R')
source('CHMGP_VariationOperators.R')
source('CHMGP_SuperflexFunctions.R')
#source('CHMGP_DEoptim.R')
source('CHMGP_DDS.R')
source('New_Marrmot_functions.R')

#system('R CMD SHLIB CHMGP_SuperflexFunctions.f90', ignore.stdout = FALSE)
#dyn.load('CHMGP_SuperflexFunctions.so')

CHMGP = function(DependentVariable, IndependentVariables, DataSet, Weights = NULL,
                 FunctionSet = c('+', '-', '*', '/', 'sqrt','^'), ConstantRange = c(-1, 1), StrictFunctionsArgumentsRange = c(lower = 0, upper = 1),
                 PopulationSize = 500, NumberOfGenerations = 50, FitnessFunction = c('RMSE','RMSE','RMSE','RMSE','RMSE'), 
                 MaxDepthIni = 2, MaxDepthRun = 3, TournamentSize = 4, RoundingFactor = 6, 
                 VariationProbs = c(Pc = 0.7, PrSelected = 0.05, PmTree = 0.5, PmSeparation = 0.2, PmNode = 0.5, PmConstant = 0.7),
                 TerminalNodeProbs = 0.7, PunishOneNodeIndividuals = FALSE, SimpleOutput = TRUE, YacasSimplification = FALSE, 
                 DEoptimization = FALSE, DEiterMax = 200, DEfitness = 'SameAsGPfit',
                 DynamicPopulation = FALSE, CTtypesOccurence = FALSE, n_jobs=4
){
  
  if(DEoptimization) library('DEoptim')
  
  if(YacasSimplification){
    
    yacTest = system('which yacas', intern = TRUE)
    if(length(yacTest) == 0) stop('Yacas is not installed!')
    
  }
  library('parallel')
  library('nsga2R')
  library('fuse')
  
  P_E_funs = c("TANK", "MI", "MII", "MIII", "MIV", "MV", "MVI", "MVII", "MVIII", "MIX", "MX","MXI","MXII", "CT","FUSE","DCT","CT_new","MM_29")
  
  SetGlobalVariables(FunctionSet, IndependentVariables, ConstantRange, StrictFunctionsArgumentsRange, RoundingFactor, DataSet, DependentVariable, 
                     FitnessFunction, Weights, PunishOneNodeIndividuals, MaxDepthRun)
  
  constantMutationFactors = seq(2, 0.2, length = NumberOfGenerations)
  
  #######################
  #create individuals in first population
  population = FirstGeneration(MaxDepthIni, PopulationSize)
  
  #######################
  #Apply DDS algorithm to part of initial population
  DDS_indArray = array(c(NA), dim = c(64, 20))
  DDS_indArray[1,1] <- "CT_new"
  DDS_indArray[2,1] <- "WR"
  DDS_indArray[3,1] <- "IR"
  DDS_indArray[4,1] <- "RR"
  DDS_indArray[5,1] <- "UR"
  DDS_indArray[6,1] <- "FR"
  DDS_indArray[7,1] <- "SR"
  DDS_indArray[8,1] <- "CR"
  DDS_indArray[1,2:20] <-as.character(c(2:20))
  DDS_indArray[2,2:6] <- as.character(c(21:25))
  DDS_indArray[3,2:5] <- as.character(c(26:29))
  DDS_indArray[4,2:4] <- as.character(c(30:32))
  DDS_indArray[5,2:8] <- as.character(c(33:39))
  DDS_indArray[6,2:8] <- as.character(c(40:46))
  DDS_indArray[7,2:6] <- as.character(c(47:51))
  DDS_indArray[8,2:14] <- as.character(c(52:64))
  DDS_indArray[9:64,1] <- as.character(c(1:56))
  
  DDS_IND <- list(IndArray = DDS_indArray, Equation = NA, Fitness = c(rep(NA,length(FitnessFunction))), DriveVec = c(rep("Function",8),rep("Terminal",56)),
                  RestrictedRows = c(2:64), IndLength = 64, Changed = TRUE, Front = NA, Crowding_distance = NA,Sim =c())
  
  no_DDS_IND <- n_jobs*2 
  DDS_IND_list <- list()
  for (i in 1:no_DDS_IND){
    DDS_IND_list[[i]] <-runif(56,0,1)
  }
  
  Sep_DDS_IND_list_1 <- mclapply(as.list(DDS_IND_list),DDS,mc.cores = n_jobs)
  #Sep_DDS_IND_list_1 <- lapply(DDS_IND_list,DDS)
  for (i in 1:no_DDS_IND){
    DDS_IND[['IndArray']][9:64,1] <-as.character(Sep_DDS_IND_list_1[[i]])
    DDS_IND_list[[i]] <- DDS_IND
  }
  
  population[1:no_DDS_IND] <- DDS_IND_list
  ###################
  #Remove unsuitable CT_new individuals
  
  for (r in 1:PopulationSize) {
    if ((population[[r]][['IndArray']][1,1]=="CT_new") && (population[[r]][['IndLength']] <= 20)){
      population[[r]] <- CreateIndividual(2,"full",F)
    } else {
      population[[r]]<- population[[r]]
    }
  }
  
  ####################
  #calculate fitness for 1st generation
  #population = lapply(population, ComputeIndividual)
  population <- mclapply(as.list(population),ComputeIndividual, mc.cores = n_jobs)
  #population <- lapply(population,ComputeIndividual)
  names(population) = paste0('I', 1:PopulationSize)
  
  ###################
  #need to calculate front and crowding distances
  # calculating front 
  front_matrix <-t(sapply(population, function(x) x[['Fitness']]))
  
  front_list <-fastNonDominatedSorting(front_matrix)
  
  for (i in 1: length(front_list)){
    for (k in 1:length(front_list[[i]])){
      population[[front_list[[i]][k]]][['Front']] <- i 
    }
  }
  ##################################
  # calculating crowding distance
  
  rnkIndex <- integer(PopulationSize)
     i <- 1
     while (i <= length(front_list)) {
        rnkIndex[front_list[[i]]] <- i
        i <- i + 1
      }
  front_matrix <- cbind(front_matrix,rnkIndex)
  
  
  objrange <- apply(front_matrix[,1:length(FitnessFunction),drop=F],2,max) - apply(front_matrix[,1:length(FitnessFunction),drop=F],2,min)
  cr_dis <- crowdingDist4frnt(front_matrix,front_list,objrange)
  front_matrix <-cbind(front_matrix,apply(cr_dis,1,sum))
  
  for (i in 1:PopulationSize){
    population[[i]][['Crowding_distance']] <- sum(cr_dis[i,])
  }
  
  #####################
  isFitMulti = FitnessFunction %in% MultiObjectiveFitnesses
  if(isFitMulti){
    
    generationFitnesses = MultiObjective(population)
    
  } else{
    
    generationFitnesses = sapply(population, function(x) x[['Fitness']])
    generationFitnesses = generationFitnesses[1,]
  }
  
  ##########################
  #printing & saving results 
  
  runResults = vector("list",PopulationSize)
  printGenerationResults = c(0, floor(seq(0.1, 0.9, 0.1) * NumberOfGenerations), NumberOfGenerations)
  cat('Generation:', printGenerationResults[1], 'Completed,  ')

  ########################
  
  tenthDEopt = FALSE
  if(tenthDEopt){
    
    tenthSeq = seq(10, NumberOfGenerations, 10) + 1
    tenthDEoptResults = data.frame(Generation = tenthSeq, Fit = numeric(length(tenthSeq)), OptModel = character(length(tenthSeq)), stringsAsFactors = FALSE)
    
  }
  
  if(CTtypesOccurence) {
    
    # List for storage of equations from each generation
    eqStorage = vector(mode = 'list', NumberOfGenerations + 1)
    eqStorage = lapply(eqStorage, function(x){character(PopulationSize)})
    
    eqStorage[[1]] = sapply(population, '[[', 'Equation')
    
  }
  ################# 
  #main cycle
  
  for(q in 1:NumberOfGenerations){
    
    goodFitnesses = generationFitnesses[!is.na(generationFitnesses)]
    
    PopulationCheck(goodFitnesses, TournamentSize, q)
    
    newPopulation = vector('list', PopulationSize)
    
    # creating a mating pool
    front_matrix <- cbind(1:PopulationSize,front_matrix)
    matingPool <- tournamentSelection(front_matrix,PopulationSize,TournamentSize)
    
    noOfNewIndividuals = min(sum(rbinom(PopulationSize, 1, VariationProbs[['Pc']])), PopulationSize)
    noOfNewIndividuals = ceiling(noOfNewIndividuals/n_jobs)*n_jobs
    
    if(noOfNewIndividuals > 0){
      parents = matrix(NA, nrow = noOfNewIndividuals, ncol = 2)
      colnames(parents) = c('Parent1', 'Parent2')
      for (i in 1:noOfNewIndividuals){
        ran_ind_cr <- round(runif(2,1,PopulationSize))
        parents[i,1] <- matingPool[ran_ind_cr[1],][1]
        parents[i,2] <- matingPool[ran_ind_cr[2],][1]
      }
      
      positionRange = 1:noOfNewIndividuals
      newPopulation[positionRange] = Crossover(population[parents[,1]], population[parents[,2]], TerminalNodeProbs)
      
    }
    
    noOfMutants = PopulationSize - noOfNewIndividuals
    Mutating_ind <- c(rep(NA,noOfMutants))
    for (i in 1:noOfMutants){
      ran_ind_mu <- round(runif(1,1,PopulationSize))
      Mutating_ind[i] <- matingPool[ran_ind_mu,1]
    }
    if(noOfNewIndividuals < PopulationSize){
      
      rangeOfMutatingPopulation = (noOfNewIndividuals + 1):PopulationSize
      newPopulation[rangeOfMutatingPopulation] = Mutation(population[Mutating_ind], VariationProbs, constantMutationFactors[q])  
    }  
    
    ####################
    
    accPopulation <- c(population,newPopulation)
    
    names(accPopulation) = paste0('I', 1:length(accPopulation))
    
    # for (i in 1:length(accPopulation)){
    #   cat(i,'-',length(accPopulation[[i]][['Fitness']]),',')
    # }
    
    #creating front matrix
    front_matrix <-t(sapply(accPopulation, function(x) x[['Fitness']]))
    
    front_list <-fastNonDominatedSorting(front_matrix)
    for (i in 1: length(front_list)){
      for (k in 1:length(front_list[[i]])){
        accPopulation[[front_list[[i]][k]]][['Front']] <- i 
      }
    }
    
    #crowding_distance
    rnkIndex <- integer(PopulationSize*2)
    i <- 1
    while (i <= length(front_list)) {
      rnkIndex[front_list[[i]]] <- i
      i <- i + 1
    }
    front_matrix <- cbind(front_matrix,rnkIndex)
    
    
    objrange <- apply(front_matrix[,1:length(FitnessFunction),drop=F],2,max) - apply(front_matrix[,1:length(FitnessFunction),drop=F],2,min)
    cr_dis <- crowdingDist4frnt(front_matrix,front_list,objrange)
    front_matrix <-cbind(front_matrix,apply(cr_dis,1,sum))
    
    for (i in seq_along(accPopulation)){
      accPopulation[[i]][['Crowding_distance']] <- sum(cr_dis[i,])
    }
    
    #############
    #NSGA selection
    
    m <- 1
    Ind_cunt <- 0
    if (length(front_list[[1]])>= PopulationSize){
      margin_front <- 1
    } else {
      while(Ind_cunt +length(front_list[[m]]) < PopulationSize){
        Ind_cunt <- Ind_cunt + length(front_list[[m]])
        m <- m + 1
        margin_front <- m
      }
    }
    
    nextPopulation = list()
    
    if (margin_front ==1){
      margin_list <- accPopulation[front_list[[margin_front]]]
      crw_dis <- c(rep(0,length(front_list[[margin_front]])))
      for (i in 1:length(front_list[[margin_front]])){
        crw_dis[i] <- margin_list[[i]][['Crowding_distance']]
      }
      remain_ind <- data.frame(index = c(1:length(front_list[[margin_front]])),crw_dis = crw_dis)
      sort_remain_ind <- remain_ind[order(remain_ind[['crw_dis']],decreasing = T),]
      nextPopulation <- accPopulation[sort_remain_ind[1:PopulationSize,1]]
      population <- nextPopulation
    } else {
      len_front <- 0
      for (i in 1:(margin_front-1)){
        nextPopulation <- c(nextPopulation,accPopulation[front_list[[i]]])
        len_front <- len_front + length(front_list[[i]])
      }
      if (length(nextPopulation) == PopulationSize){
        population = nextPopulation
      } else {
        margin_list <- accPopulation[front_list[[margin_front]]]
        crw_dis <- c(rep(0,length(front_list[[margin_front]])))
        for (i in 1:length(front_list[[margin_front]])){
          crw_dis[i] <- margin_list[[i]][['Crowding_distance']]
        }
        
        remain_ind <- data.frame(index = c(1:length(front_list[[margin_front]])),crw_dis = crw_dis)
        sort_remain_ind <- remain_ind[order(remain_ind[['crw_dis']],decreasing = T),]
        nextPopulation <- c(nextPopulation,accPopulation[sort_remain_ind[1:(PopulationSize- Ind_cunt),1]])
        population <- nextPopulation
      }
    }
    
    names(population) = paste0('I', 1:PopulationSize)
    
    # need to calculate front and crowding distances
    # calculating front
    front_matrix <-t(sapply(population, function(x) x[['Fitness']]))
    front_list <-fastNonDominatedSorting(front_matrix)
    for (i in 1: length(front_list)){
      for (k in 1:length(front_list[[i]])){
        population[[front_list[[i]][k]]][['Front']] <- i 
      }
    }
    
    #################################
    # calculating crowding distance
    rnkIndex <- integer(PopulationSize)
    i <- 1
    while (i <= length(front_list)) {
      rnkIndex[front_list[[i]]] <- i
      i <- i + 1
    }
    front_matrix <- cbind(front_matrix,rnkIndex)
    
    
    objrange <- apply(front_matrix[,1:length(FitnessFunction),drop=F],2,max) - apply(front_matrix[,1:length(FitnessFunction),drop=F],2,min)
    #objrange <-c((max(front_matrix[,1])-min(front_matrix[,1])),(max(front_matrix[,2])-min(front_matrix[,2])),(max(front_matrix[,3])-min(front_matrix[,3])))
    cr_dis <- crowdingDist4frnt(front_matrix,front_list,objrange)
    front_matrix <-cbind(front_matrix,apply(cr_dis,1,sum))
    
    for (i in 1:PopulationSize){
      population[[i]][['Crowding_distance']] <- sum(cr_dis[i,])
    }
    
    ########################
    if(isFitMulti){
      
      generationFitnesses = MultiObjective(population)
      
    } else{
      
      generationFitnesses = sapply(population, function(x) x[['Fitness']])
      
    }
    
    runResults = SaveRunResults(population, front_list, runResults, q, NumberOfGenerations)
    if (q == NumberOfGenerations){
      runResults <- runResults[1:length(front_list[[1]])]
    } else {
      runResults <- runResults
    }
    if(any(q == printGenerationResults)) cat('Generation:', q, 'Completed,  ')    
  } 
  runResults
}

