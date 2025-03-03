# Clear Global Environment
rm(list = ls())

# Setting working directory
setwd('C:/Naila/Singapore/Research/Vidura codes/Final_Codes/Lumped/MARRMOT')
path <- 'C:/Naila/Singapore/Research/Vidura codes/Final_Codes/Lumped/MARRMOT'
# Reading input data
luxHue = read.table('./Syn_data/Blackwater.txt', sep=',',header = TRUE)
luxHue_train =luxHue[8493:10684,]

# set CHMGP folder 
setwd('./CHMGP')
source('./CHMGP.R', local = TRUE)
setwd('../')

# Import relevant libraries 
library('parallel')
library('nsga2R')
library('fuse')
library('matlabr')


###############
# Setting CHMGP parameters

IndependentRuns <- 10
PopulationSize <- 500
NumberOfGenerations <- 50
n_jobs <- 1  #24
specFunct = c('MARRMot')
criterion = c('VE0','KG10','NS0','logNS0')

constantMutationFactors = seq(2, 0.2, length = NumberOfGenerations)
VariationProbs = c(Pc = 0.7, PrSelected = 0.05, PmTree = 0.5, PmSeparation = 0.3, PmNode = 0.3, PmConstant = 0.7)
TournamentSize <- 4
TerminalNodeProbs <- 0.9
dataNameTraining = 'luxHue_train'
toSearch = 'Q'
depVari = paste0(toSearch)
dataSetTrain = get(dataNameTraining)

setwd('./Matlab_fn')
write.csv(dataSetTrain,"Input_Variables.csv")
setwd('../')

Date=as.Date(dataSetTrain[['Date']])  #Date=strptime(dataSetTrain[['Date']],format = '%d/%m/%Y %H:%M')
DependentVariable = depVari
IndependentVariables = c('P','E','T')
FunctionSet = c(specFunct, '+', '-', '*', '/')
ConstantRange = c(0, 1)
StrictFunctionsArgumentsRange = c(lower = 0, upper = 1)
FitnessFunction = criterion
MaxDepthIni <- 1
MaxDepthRun = 3
RoundingFactor = 3
results = vector(mode = 'list', IndependentRuns)
namesRes = paste0('run_', 1:IndependentRuns)
names(results) = namesRes
DataSet <- dataSetTrain
Weights = NULL
PunishOneNodeIndividuals = FALSE
SimpleOutput = TRUE
YacasSimplification = FALSE
DEoptimization = FALSE
DEiterMax = 200
DEfitness = 'SameAsGPfit'
DynamicPopulation = FALSE
CTtypesOccurence = FALSE

# Link Globle variables
SetGlobalVariables(FunctionSet, IndependentVariables, ConstantRange, StrictFunctionsArgumentsRange, RoundingFactor, DataSet, DependentVariable, 
                   FitnessFunction, Weights, PunishOneNodeIndividuals, MaxDepthRun)

################
# Saving console content

sink("Blackwater_Marrmot_Lumped_1", append=FALSE, split=FALSE)

################
# Running several Independent runs

ptm = proc.time()

for(p in namesRes){
  
  print(p)
  
  ################
  # Initialization
  
  population = FirstGeneration(MaxDepthIni, PopulationSize)
  
  # Ind_index <- 1:length(population)
  # Marrmot_Ind_index <-Ind_index[sapply(population, function(x) x[['IndArray']][1,1]=="MARRMot")]
  # Other_Ind_index <- Ind_index[sapply(population, function(x) x[['IndArray']][1,1]!="MARRMot")]
  # pre_index_order <-c(Marrmot_Ind_index,Other_Ind_index)
  
  Marrmot_Ind <-population[sapply(population, function(x) x[['IndArray']][1,1]=="MARRMot")]
  
  Other_Ind <-population[sapply(population, function(x) x[['IndArray']][1,1]!="MARRMot")]
  
  # if(length(Other_Ind)>0){
  #   for(vv in 1:length(Other_Ind)){
  #     Other_Ind[[vv]]<- ComputeIndividual(Other_Ind[[vv]])
  #   }
  # }
  
  if(length(Other_Ind)>0){
    Other_Ind <- mclapply(as.list(Other_Ind),ComputeIndividual, mc.cores = n_jobs)
  }
  
  Output_script_vec <- as.character(1:length(Marrmot_Ind))
  
  MakeEquation_new = function(individual){
    ReservoirSet <-c("WR","IR","RR","UR","FR","SR","CR")
    IA = individual[['IndArray']]
    if(individual[['IndLength']] == 1){
      return(paste0('y = ', IA[, 1]))
    }
    eq_vec = IA[, 1]
    for (i in individual[['IndLength']]:1){
      if( any(FunctionSet == IA[i, 1]) | any(ReservoirSet == IA[i, 1])){
        for(j in 2:ncol(IA)){
          if( !is.na(IA[i,j]) ){
            arg_row = as.numeric(IA[i,j])
            if( j == 2 && is.na(IA[i, (j + 1)]) ){
              if(!any(IA[i, 1] == SpecialFunctions)){
                eq_vec[i] = paste0(eq_vec[i], '(', eq_vec[arg_row])
              }else{
                eq_vec[i] = paste0(eq_vec[i], eq_vec[arg_row])
              }
            }
            if( j == 2 && !is.na(IA[i, (j + 1)]) ){
              if( any(PrefixFunctions == IA[i, 1]) ){
                first = paste0(eq_vec[i],'(', eq_vec[arg_row], ',')
              } else{
                first = paste0('(', eq_vec[arg_row])
              }
            }
            if(j > 2){
              if( any(PrefixFunctions == IA[i, 1]) ){
                eq_vec[i] = paste0(first, eq_vec[arg_row])
                if(j >= 3) first = paste0(eq_vec[i], ',')
              } else{
                eq_vec[i] = paste0(first, eq_vec[i], eq_vec[arg_row])
              }
            }
          }
        }
        eq_vec[i] = paste0(eq_vec[i], ')')
      }
    } 
    return(eq_vec)
  }
  
  Matlab_code <- c()
  
  if(length(Marrmot_Ind)>0){
    
    for(gg in 1:length(Marrmot_Ind)){
      for(ml in 1:nrow(Marrmot_Ind[[gg]][['IndArray']])){
        if(Marrmot_Ind[[gg]][['IndArray']][ml,1] %in% c("P","E","T")){
          Marrmot_Ind[[gg]][['IndArray']][ml,1] <- runif(1)
        }
      }
      M_Eq <- MakeEquation_new(Marrmot_Ind[[gg]])
      Parameters <- c()
      for(hh in 1:33){
        Parameters[hh]<-max(min(eval(parse(text=noquote(M_Eq[[as.numeric(Marrmot_Ind[[gg]][['IndArray']][1,(hh+1)])]]))),1),0)
      }
      Output_script_name <- Output_script_vec[gg]
      Marrmot_Ind[[gg]][['MARRMot_Ind']]<- Output_script_name
      Model_to_use <- Code_name(Parameters[1],Parameters[2],Parameters[3],Parameters[4],Parameters[5],Parameters[6])
      Matlab_code <- c(Matlab_code,Code_generator(Model_to_use,Parameters[7:length(Parameters)],Output_script_name))
    }
    
    TO_MAT <- array(dim = c(length(Marrmot_Ind),1))
    TO_MAT[,1]<-Matlab_code
    write.csv(TO_MAT,paste(path,"/Matlab_fn/TO_MAT.csv",sep = ''))
    
    MAT_CODE <-c(paste("cd '",path,"/Matlab_fn","';",sep = ""),"A=readcell('TO_MAT.csv')","B=string();",paste("for i =2:",(length(Marrmot_Ind)+1),";B(i-1)=string(A(i,2));end",sep = ""),paste("parfor i=1:",length(Marrmot_Ind),";CommitEval(B(1,i));end",sep = ""))
  
    run_matlab_code(MAT_CODE)
    file_location <- paste(path,"/Matlab_fn",sep = '')
    
    Read_Marrmot <- function(Individual){
      Matlab_output <-read.csv(paste(file_location,"/Output_Q_",Individual[['MARRMot_Ind']],".csv",sep = ''))
      Matlab_output <- t(Matlab_output)
      sim <- Matlab_output[,1][paste("Q_",1:nrow(DataSet),sep = "")]
      obs <- DataSet$Q
      for(mn in 1:length(FitnessFunction)){
        Individual$Fitness[mn] <- FitnessComputation(sim,obs,FitnessFunction[mn])
      }
      Individual$Sim <- sim
      Individual$Equation <- MakeEquation(Individual)
      Individual <-Individual[1:(length(Individual)-1)]
      
      Individual
      
    }
    
    Marrmot_Ind <- mclapply(as.list(Marrmot_Ind),Read_Marrmot, mc.cores = n_jobs)
    
    # for(x in 1:length(Marrmot_Ind)){
    #   Matlab_output <-read.csv(paste(file_location,"/Output_Q_",Output_script_vec[x],".csv",sep = ''))
    #   #unlink(paste(file_location,"/Output_Q_",Output_script_vec[x],".csv",sep = ''))
    #   Matlab_output <- t(Matlab_output)
    #   sim <- Matlab_output[,1][paste("Q_",1:nrow(DataSet),sep = "")]
    #   obs <- DataSet$Q
    #   
    #   for(mn in 1:length(FitnessFunction)){
    #     Marrmot_Ind[[x]]$Fitness[mn] <- FitnessComputation(sim,obs,FitnessFunction[mn])
    #   }
    #   Marrmot_Ind[[x]]$Sim <- sim
    #   Marrmot_Ind[[x]]$Equation <- MakeEquation(Marrmot_Ind[[x]])
    # }
  
  }
  
  population <-c(Marrmot_Ind,Other_Ind)
  
  # Index_dataframe <-data.frame(pre_index_order,1:length(population))
  # names(Index_dataframe)<-c("pre","post")
  # Index_dataframe <-Index_dataframe[order(Index_dataframe$pre),]
  # population <-population[Index_dataframe$post]
  
  #population <- mclapply(as.list(population),ComputeIndividual, mc.cores = n_jobs)
  
  names(population) = paste0('I', 1:PopulationSize)
  
  # calculating fronts
  front_matrix <-t(sapply(population, function(x) x[['Fitness']]))
  
  front_list <-fastNonDominatedSorting(front_matrix)
  for (i in 1: length(front_list)){
    for (k in 1:length(front_list[[i]])){
      population[[front_list[[i]][k]]][['Front']] <- i
    }
  }
  
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
  
  printGenerationResults = c(0, floor(seq(0.1, 0.9, 0.1) * NumberOfGenerations), NumberOfGenerations)
  cat('Generation:', printGenerationResults[1], 'Completed,  ')
  
  ##################
  # Main cycle
  
  for (q in 1:NumberOfGenerations){
    
    newPopulation = vector('list', PopulationSize)
    
    # creating a mating pool
    front_matrix <- cbind(1:PopulationSize,front_matrix)
    matingPool <- tournamentSelection(front_matrix,PopulationSize,TournamentSize)
    
    noOfNewIndividuals = min(sum(rbinom(PopulationSize, 1, VariationProbs[['Pc']])), PopulationSize)
    noOfNewIndividuals = ceiling(noOfNewIndividuals/n_jobs)*n_jobs
    
    # Crossover
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
    
    # Mutation
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
    
    # to avoid subscript out of bounds problem in some generations, let's try following code
    for (i in 1:length(newPopulation)){
      if (length(newPopulation[[i]]) <10){
        newPopulation[[i]] <- population[[round(runif(1,1,PopulationSize))]]
      } else if ((length(newPopulation[[i]]) == 10) & (length(newPopulation[[i]][['Sim']]) ==1)){
        newPopulation[[i]] <- population[[round(runif(1,1,PopulationSize))]]
      } else {
        newPopulation[[i]] <- newPopulation[[i]]
      }
    }
    
    # Creating Parent + Child Population
    accPopulation <- c(population,newPopulation)
    
    names(accPopulation) = paste0('I', 1:length(accPopulation))
    
    # Creating front matrix for Parent + Child Population
    front_matrix <-t(sapply(accPopulation, function(x) x[['Fitness']]))
    
    front_list <-fastNonDominatedSorting(front_matrix)
    
    for (i in 1: length(front_list)){
      for (k in 1:length(front_list[[i]])){
        accPopulation[[front_list[[i]][k]]][['Front']] <- i
      }
    }
    
    # Crowding_distance for Parent + Child Population
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
    
    # NSGA Selection
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
    
    # Calculating front for selected population
    front_matrix <-t(sapply(population, function(x) x[['Fitness']]))
    front_list <-fastNonDominatedSorting(front_matrix)
    
    for (i in 1: length(front_list)){
      for (k in 1:length(front_list[[i]])){
        population[[front_list[[i]][k]]][['Front']] <- i
      }
    }
    
    # Calculating crowding distance for selected population
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
    if(any(q == printGenerationResults)) cat('Generation:', q, 'Completed,  ')
  }
  
  # Running Time
  etmh = (proc.time() - ptm)['elapsed']/3600
  etmd = (proc.time() - ptm)['elapsed']/86400
  message(paste(p, 'Elapsed time:', etmh, 'hours, (ie.',etmd, 'days)'))
  
  cat('Elapsed time:',p, etmh,'hours')
  
  ################
  # Saving Results
  
  runResults = vector("list",PopulationSize)
  runResults = SaveRunResults(population, front_list, runResults, q, NumberOfGenerations)
  runResults <- runResults[front_list[[1]]]
  results[[p]] <- runResults
  saveRDS(results, paste0('Backup_',dataNameTraining, '_Blackwater_Marrmot_Lumped_1.rds'))
}
saveRDS(results, paste0(dataNameTraining, '_Blackwater_Marrmot_Lumped_1.rds'))
sink()