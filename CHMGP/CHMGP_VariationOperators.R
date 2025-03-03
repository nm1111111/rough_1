PopulationCheck = function(goodFitnesses, tournamentSize, generation){
  
  if(length(goodFitnesses) < (2 * tournamentSize)){
    
    warning(paste0('PopulationCheck: Too many NA fitnesses in population! Variability operators ineffective. Generation: ', generation))
    
  }
  
}

######## need to a change here to add crowding distance and front into cosideration of selecting the winner ############
TournamentSelection = function(goodFitnesses, tournamentSize, oneWinner = FALSE){
  
  if(oneWinner){
    
    competitors = Sample(goodFitnesses, tournamentSize)
    return(names(competitors[which.min(competitors)]))
    
  }
  
  noCompetitors = 2 * tournamentSize
  
  competitors = Sample(goodFitnesses, noCompetitors)
  
  groupA = head(competitors, tournamentSize)
  groupB = tail(competitors, tournamentSize)
  
  winners = c(groupA[which.min(groupA)], groupB[which.min(groupB)])
  
  names(winners)
  
}

NodeSelection = function(driveVector, restrictedRowsVector, terminalNodeProbability){
  
  if(length(driveVector) == 1) return(list(targetNode = 1, restriction = FALSE))
  
  # Vector of probabilities for terminal type determination (Koza's recommendation is P(function) = 0.1, P(terminal) = 0.9)
  # terminalNodeProbability = 0.9   
  selectedTerminal = as.logical(rbinom(1, 1, terminalNodeProbability))
  
  if(selectedTerminal){  
    
    targetNode = Sample(which(driveVector == 'Terminal' ), 1)
    
  } else{
    
    targetNode = Sample(which(driveVector == 'Function' ), 1)
    
  }
  
  
  if(targetNode %in% restrictedRowsVector) return(list(targetNode = targetNode, restriction = TRUE))
  
  list(targetNode = targetNode, restriction = FALSE)
  
}

RecursiveNodePointers = function(treeArray, checkRow, nodesPointers = c()){
  
  actualPointers = numeric(0)
  
  if(length(checkRow) > 0){
    
    for(i in 1:length(checkRow)){
      
      actualPointers = c(actualPointers, as.numeric(treeArray[checkRow[i], !is.na(treeArray[checkRow[i], ])]))
      
    }
    
    nodesPointers = c(nodesPointers, actualPointers)
    
    checkRow = actualPointers
    
    RecursiveNodePointers(treeArray, checkRow, nodesPointers)
    
  } else {
    
    return(nodesPointers)
    
  }
  
}

RecursiveRestrictedNodes = function(treeArray, checkRow, nodes, nodesPointers = c()){
  
  actualPointers = numeric(0)
  
  if(length(checkRow) > 0){
    
    for(i in 1:length(checkRow)){
      
      node = nodes[checkRow[i]]
      noOfNotRestrictedNodes = NonStrictArguments[node]
      
      identifiedPointers = !is.na(treeArray[ checkRow[i], ])
      
      identifiedPointers[1:noOfNotRestrictedNodes] = FALSE
      
      #       print(identifiedPointers)
      
      #       print(nodes[as.numeric(treeArray[checkRow[i], 1:noOfNotRestrictedNodes])] %in%  c('+', '-'))
      
      againRestricted = nodes[as.numeric(treeArray[checkRow[i], 1:noOfNotRestrictedNodes])] %in% StandardFunctions 
      
      identifiedPointers[againRestricted] = TRUE
      
      #       print(identifiedPointers)
      
      actualPointers = c(actualPointers, as.numeric(treeArray[checkRow[i], identifiedPointers]))
      
      #       notRestrictedNodes = NonStrictArguments[nodes[checkRow[i]]]
      #       identifiedPointers = !is.na(treeArray[ checkRow[i], ])
      #       identifiedPointers[1:notRestrictedNodes] = FALSE
      
      #       actualPointers = c(actualPointers, as.numeric(treeArray[checkRow[i], identifiedPointers]))
      
    }
    
    nodesPointers = c(nodesPointers, actualPointers)
    
    checkRow = actualPointers
    
    RecursiveNodePointers(treeArray, checkRow, nodesPointers)
    
  } else {
    
    return(nodesPointers)
    
  }
  
}

ChangeIndividual = function(newIndArray, newDriveVector, newLength, change){
  
 
  
  strictFpositions = which(newIndArray[, 1] %in% StrictFunctions)
  
  if(length(strictFpositions) > 0){
    
    restrictedRows = unique(RecursiveRestrictedNodes(newIndArray[, 2:MaxNoColumns], strictFpositions, newIndArray[, 1]))
    
    #     if(newLength > 10) browser()
    
    newIndArray[restrictedRows[newIndArray[restrictedRows, 1] %in% IndependentVariables], 1] = runif(1, ConstantRange[1], ConstantRange[2])
    
  } else {
    
    restrictedRows = NA
    
  }
  
  # Osetrit mista s promennymi tam kde nemaji byt - nahrada za konstanty
  
  #   argsOfSf = NonStrictArguments[newIndArray[strictFpositions, 1]]
  #   
  #   newIndArray[strictFpositions , (argsOfSf+1)]
  
  
  #   restrictedRows = as.numeric(as.vector(newIndArray[newIndArray[ , 1] %in% StrictFunctions, 3:ncol(newIndArray)]))
  #   if(length(restrictedRows) == 0) restrictedRows = NA
  #   if(length(restrictedRows) > 1 & min(restrictedRows, na.rm = TRUE) > 2) browser()
  
  list(IndArray = newIndArray, Equation = NA, Fitness = c(rep(NA,length(FitnessFunction))), DriveVec = newDriveVector, RestrictedRows = restrictedRows, 
       IndLength = newLength, Changed = change, Front = NA, Crowding_distance = NA,Sim =c())
  
}

AssignPointers = function(pointersArray, arrayLength){
  
  transposedPA = t(pointersArray)
  transposedPA[which(!is.na(transposedPA), arr.ind = TRUE)] = as.character(2:arrayLength)
  t(transposedPA)
  
}

SepareTrees = function(parent, cNode){
  
  if(!cNode[['variationAllowed']]){
    
    return(NA)
    
  }
  
  cn = cNode[['targetNode']]
  
  if(cn == 1){
    
    return(parent)
    
  }
  
  if(parent[['DriveVec']][cn] == 'Terminal'){
    
    newDriveVector = parent[['DriveVec']][cn] 
    
    newIndArray = array(parent[['IndArray']][cn, ], dim = c(1, MaxNoColumns))
    
    return(ChangeIndividual(newIndArray, newDriveVector, 1, TRUE))
    
  }
  
  if(parent[['DriveVec']][cn] == 'Function'){
    
    # Determination of separated nodes by pointers of node
    separatedRows = RecursiveNodePointers(parent[['IndArray']][ , 2:MaxNoColumns], cn, cn)
    
  }
  
  # Array and driving vector of separated tree
  newLength = length(separatedRows)
  
  newDriveVector = parent[['DriveVec']][separatedRows]  
  newIndArray = array(parent[['IndArray']][separatedRows, ], dim = c(newLength, MaxNoColumns))
  
  newIndArray[, 2:MaxNoColumns] = AssignPointers(newIndArray[, 2:MaxNoColumns], newLength)
  
  ChangeIndividual(newIndArray, newDriveVector, newLength, TRUE)  
  
}

TerminalToFunction = function(parent, cNode, subtree){
  
  # Determination of deleted subtree identified by pointers of deleting node
  removeRows = RecursiveNodePointers(parent[['IndArray']][ , 2:MaxNoColumns], cNode)
  
  newLength = parent[['IndLength']] - length(removeRows)
  
  newDriveVector = parent[['DriveVec']][-removeRows] 
  newDriveVector[cNode] = 'Terminal'
  
  newIndArray = array(parent[['IndArray']][-removeRows, ], dim = c(newLength, MaxNoColumns))
  newIndArray[cNode, ] = subtree[['IndArray']][1, ]
  
  newIndArray[, 2:MaxNoColumns] = AssignPointers(newIndArray[, 2:MaxNoColumns], newLength)
  
  ChangeIndividual(newIndArray, newDriveVector, newLength, TRUE)
  
}

ChangePointers = function(subtreeArray, lastPointer, startRow, stopRow){
  
  changedArray = subtreeArray[startRow:stopRow, 2:MaxNoColumns]
  
  transposedChA = t(changedArray)
  pointersPositions = which(!is.na(transposedChA), arr.ind = TRUE)
  newPointers = as.character((lastPointer + 1):(nrow(pointersPositions) + lastPointer))
  transposedChA[pointersPositions] = newPointers
  subtreeArray[startRow:stopRow, 2:MaxNoColumns] = t(transposedChA)  
  
  subtreeArray
  
}

FunctionToTerminal = function(parent, cNode, subtree){
  
  newLength = parent[['IndLength']] + subtree[['IndLength']] - 1
  
  newIndArray = array(NA, dim = c(newLength, MaxNoColumns))
  
  newDriveVector = rep(NA, newLength)
  
  # Filling new array and driving vector
  newIndArray[1:parent[['IndLength']], ] = parent[['IndArray']]
  newDriveVector[1:parent[['IndLength']]] = parent[['DriveVec']]
  
  subtreeCopy = subtree
  
  # Change of pointers
  lastMaxPointer = max(as.numeric(parent[['IndArray']][1:(cNode-1), 2:MaxNoColumns]), na.rm = TRUE)
  subtree[['IndArray']] = ChangePointers(subtree[['IndArray']], lastMaxPointer, 1, subtree[['IndLength']])
  
  # Addition of subtree root node
  newIndArray[cNode, ] = subtree[['IndArray']][1, ]
  newDriveVector[cNode] = 'Function'
  
  # Last maximal value of pointers
  lastMaxPointer = max(as.numeric(newIndArray[1:cNode, 2:MaxNoColumns]), na.rm = TRUE)
  
  # Change of original pointers by new function addition
  if(any(!is.na(newIndArray[(cNode + 1):newLength, 2:MaxNoColumns]))){
    
    newIndArray = ChangePointers(newIndArray, lastMaxPointer, (cNode + 1), newLength)
    
  }
  
  flowLength = parent[['IndLength']]
  
  for(i in 1:subtree[['IndLength']]){
    
    j = 2
    while(j <= MaxNoColumns && !is.na(subtreeCopy[['IndArray']][i, j])){
      
      subtreePointer = as.numeric(subtree[['IndArray']][i, j])
      subtreeCopyPointer = as.numeric(subtreeCopy[['IndArray']][i, j]) 
      
      if(!is.na(newDriveVector[subtreePointer])){
        
        flowLength = flowLength + 1
        
        newIndArray[(subtreePointer + 1):flowLength, ] = newIndArray[subtreePointer:(flowLength - 1), ]
        newDriveVector[(subtreePointer + 1):flowLength] = newDriveVector[subtreePointer:(flowLength - 1)]     
        
      }
      
      newIndArray[subtreePointer, ] = subtree[['IndArray']][subtreeCopyPointer, ]
      newDriveVector[subtreePointer] = subtree[['DriveVec']][subtreeCopyPointer] 
      
      if(any(subtree[['DriveVec']][subtreeCopyPointer:subtree[['IndLength']]] == 'Function', na.rm = TRUE)){
        
        lastMaxPointer = max(as.numeric(newIndArray[1:(subtreePointer - 1), 2:MaxNoColumns]), na.rm = TRUE)
        newIndArray = ChangePointers(newIndArray, lastMaxPointer, subtreePointer, newLength)
        subtree[['IndArray']] = ChangePointers(subtree[['IndArray']], lastMaxPointer, subtreeCopyPointer, subtree[['IndLength']])
        
      }
      
      j = j + 1
      
    }
    
  }
  
  ChangeIndividual(newIndArray, newDriveVector, newLength, TRUE)
  
}

AddSubtree = function(parent, cNode, subtree){
  
  # In the case, when the crossover node is not allowed, the parent is returned
  # Funguje, ale jestli radsi nezavolat mutaci na dany uzel, aby prohledavani nestagnovalo???
  if(!cNode[['variationAllowed']]){
    
    return(parent)
    
  }
  
  cn = cNode[['targetNode']]
  
  if(cn == 1){
    
    return(subtree)
    
  }
  
  parentConnectionNode = parent[['DriveVec']][cn]
  subtreeConnectionNode = subtree[['DriveVec']][1]
  
  # Case 1: crossover node in parent is terminal and subtree root is terminal
  if(parentConnectionNode == 'Terminal' & subtreeConnectionNode == 'Terminal'){
    
    parent[['IndArray']][cn, ] = subtree[['IndArray']][1, ]
    newIndividual = ChangeIndividual(parent[['IndArray']], parent[['DriveVec']], parent[['IndLength']], TRUE)
    
  }
  
  # Case 2: crossover node in parent is function and subtree root is terminal  
  if(parentConnectionNode == 'Function' & subtreeConnectionNode == 'Terminal'){
    
    newIndividual = TerminalToFunction(parent, cn, subtree)
    
  }
  
  # Case 3: crossover node in parent is terminal and subtree root is function
  if(parentConnectionNode == 'Terminal' & subtreeConnectionNode == 'Function'){
    
    newIndividual = FunctionToTerminal(parent, cn, subtree)
    
  }  
  
  # Case 4: crossover node in parent is function and subtree root is function
  if(parentConnectionNode == 'Function' & subtreeConnectionNode == 'Function'){
    
    terminalIndividual = CreateIndividual(0, 'grow', restricted = FALSE)
    newIndividual = TerminalToFunction(parent, cn, terminalIndividual)
    newIndividual = FunctionToTerminal(newIndividual, cn, subtree)
    
  } 
  
  newIndividual
  
}

NodeDepth = function(individual, nodePosition){
  
  if(nodePosition == 1) return(0)
  
  pointers = individual[['IndArray']][, 2:MaxNoColumns]
  
  oldEndPos = 1
  endPos = 1
  depth = 0
  
  while(endPos < nodePosition){
    
    endPos = max(as.numeric(pointers[1:endPos, ]), na.rm = T)
    depth = depth + 1
  }
  
  depth
  
}

ShrinkTree = function(individual, depthOfIndividual){
  
  while(depthOfIndividual > MaxDepthRun){
    
    #     shrinkNode = Sample(which(individual[['DriveVec']] == 'Function'), 1)
    shrinkNode = list(targetNode = Sample(which(individual[['DriveVec']] == 'Function'), 1), restriction = FALSE, variationAllowed = TRUE)
    
    if(shrinkNode[['targetNode']] %in% individual[['RestrictedRows']]){
      
      shrinkNode[['restriction']] = TRUE
      shrinkNode[['variationAllowed']] = FALSE
      
    }
    
    # Creation of terminal subtree with respect to restricted nodes
    if(shrinkNode[['variationAllowed']]){
      
      terminalIndividual = CreateIndividual(0, 'grow', restricted = FALSE)
      
    } else{
      
      terminalIndividual = CreateIndividual(0, 'grow', restricted = TRUE)
      #       shrinkNode[['variationAllowed']] = TRUE
      
    }
    #     terminalIndividual = CreateIndividual(0, 'grow', restricted = FALSE)
    
    individual = TerminalToFunction(individual, shrinkNode[['targetNode']], terminalIndividual)
    
    depthOfIndividual = NodeDepth(individual, individual[['IndLength']])
    
  }
  
  individual
  
}

DepthControl = function(individuals){
  
  changedIndividuals = sapply(individuals, function(x) x[['Changed']])
  
  changedIndividualsPositions = which(changedIndividuals)
  
  checkedIndividuals = individuals[changedIndividualsPositions]
  
  lengths = sapply(checkedIndividuals, function(x) x[['IndLength']])
  
  depths = mapply(function(ind, indLength) NodeDepth(ind, indLength), checkedIndividuals, lengths) 
  
  overMaxDepthPositions = which(depths > MaxDepthRun)
  
  checkedIndividuals = mapply(function(ind, depth) 
    ShrinkTree(ind, depth), 
    checkedIndividuals, depths, SIMPLIFY = FALSE)
  
  individuals[changedIndividualsPositions] = checkedIndividuals
  
  individuals
  
}

Crossover = function(parents1, parents2, terminalNodeProbability){
  # browser()
  noIndividuals = length(parents1)
  
  # Driving vectors of parents
  driveVectors1 = lapply(parents1, function(x) x[['DriveVec']])
  driveVectors2 = lapply(parents2, function(x) x[['DriveVec']])
  
  # Restricted rows vectors
  resRowVectors1 = lapply(parents1, function(x) x[['RestrictedRows']])
  resRowVectors2 = lapply(parents2, function(x) x[['RestrictedRows']])
  
  crossNodeIdentification1 = mapply(function(driveVec, resRowVec) 
    NodeSelection(driveVec, resRowVec, terminalNodeProbability), 
    driveVectors1, resRowVectors1, SIMPLIFY = FALSE)
  
  crossNodeIdentification2 = mapply(function(driveVec, resRowVec) 
    NodeSelection(driveVec, resRowVec, terminalNodeProbability), 
    driveVectors2, resRowVectors2, SIMPLIFY = FALSE) 
  
  crossNodes1 = data.frame(targetNode = sapply(crossNodeIdentification1, function(x) x[['targetNode']]),
                           restriction = sapply(crossNodeIdentification1, function(x) x[['restriction']]))  
  
  crossNodes2 = data.frame(targetNode = sapply(crossNodeIdentification2, function(x) x[['targetNode']]),
                           restriction = sapply(crossNodeIdentification2, function(x) x[['restriction']])) 
  
  #################################################
  if ("CT_new" %in% FunctionSet){
    CT_vector1 <- lapply(parents1, function(x) x[['IndArray']][1,1])
    CT_vector2 <- lapply(parents2, function(x) x[['IndArray']][1,1])
    
    CT_elements1 <- grep("CT_new", CT_vector1)
    CT_elements2 <- grep("CT_new", CT_vector2)
    
    for (i in 1:noIndividuals){
      if ((CT_vector1[[i]]=="CT_new") && (crossNodes1[i,1] %in% c(2:8)) && (CT_vector2[[i]]=="CT_new")){
        crossNodes2[i,1] <- crossNodes1[i,1]
        crossNodes2[i,2] <- TRUE
      } else if ((CT_vector1[[i]]=="CT_new") && (crossNodes1[i,1] %in% c(2:8)) && (CT_vector2[[i]] !="CT_new")){
        crossNodes2[i,1] <- crossNodes1[i,1]
        crossNodes2[i,2] <- TRUE
        parents2[[i]] <- parents2[[CT_elements2[round(runif(1,1,length(CT_elements2)))]]]
      } 
    }
    for (i in 1:noIndividuals){
      if ((CT_vector2[[i]]=="CT_new") && (crossNodes2[i,1] %in% c(2:8)) && (CT_vector1[[i]]=="CT_new")){
        crossNodes1[i,1] <- crossNodes2[i,1]
        crossNodes1[i,2] <- TRUE
      } else if ((CT_vector2[[i]]=="CT_new") && (crossNodes2[i,1] %in% c(2:8)) && (CT_vector1[[i]] !="CT_new")){
        crossNodes1[i,1] <- crossNodes2[i,1]
        crossNodes1[i,2] <- TRUE
        parents1[[i]] <- parents1[[CT_elements1[round(runif(1,1,length(CT_elements1)))]]]
      } 
    }
  } else {
    crossNodes1 <- crossNodes1
    crossNodes2 <- crossNodes2
  }
  
  #changing target node if it is 1 in marmmot
  for (u in 1:nrow(crossNodes1)){
    if(crossNodes1[u,1]==1 && parents1[[u]][['IndArray']][1,1]=="MARRMot"){
      crossNodes1[u,1]=round(runif(1,2,parents1[[u]][['IndLength']]))
      crossNodes1[u,2]=TRUE
    }
    if(crossNodes2[u,1]==1 && parents2[[u]][['IndArray']][1,1]=="MARRMot"){
      crossNodes2[u,1]=round(runif(1,2,parents2[[u]][['IndLength']]))
      crossNodes2[u,2]=TRUE
    }
  }
  
  
  variationAllowed = !xor(crossNodes1[['restriction']], crossNodes2[['restriction']])
  crossNodes1 = cbind(crossNodes1, variationAllowed)  
  crossNodes2 = cbind(crossNodes2, variationAllowed)
  
  crossNodes1 = split(crossNodes1, seq(nrow(crossNodes1))) 
  crossNodes2 = split(crossNodes2, seq(nrow(crossNodes2))) 
  
  # Separation of crossovering subtrees
  separedTrees1 = mapply(function(parent, cNode) SepareTrees(parent, cNode), parents1, crossNodes1, SIMPLIFY = FALSE)  
  separedTrees2 = mapply(function(parent, cNode) SepareTrees(parent, cNode), parents2, crossNodes2, SIMPLIFY = FALSE)  
  
  # Offsprings creation
  offsprings1 = mapply(function(parent, cNode, separedTree) AddSubtree(parent, cNode, separedTree), 
                       parents1, crossNodes1, separedTrees2, SIMPLIFY = FALSE)
  
  offsprings2 = mapply(function(parent, cNode, separedTree) AddSubtree(parent, cNode, separedTree), 
                       parents2, crossNodes2, separedTrees1, SIMPLIFY = FALSE)
  
  
  # Kontroly, smazat                     
  #   lapply(offsprings1, function(x) {
  #     a = x$IndArray
  #     if(a[1,1] %in% StrictFunctions & any(tail(a[-2,1]) %in% c('P','ET', 'CRI'))) browser() 
  #   }) 
  #   
  #   lapply(offsprings2, function(x) {
  #     a = x$IndArray
  #     if(a[1,1] %in% StrictFunctions & any(tail(a[-2,1]) %in% c('P','ET', 'CRI'))) browser() 
  #   })    
  #
  
  
  names(offsprings1) = names(offsprings2) = paste0('Off', 1:noIndividuals)
  
  # Depth control
  offsprings1 = DepthControl(offsprings1)
  offsprings2 = DepthControl(offsprings2)
  #offsprings1 <- mclapply(as.list(offsprings1), ComputeIndividual, mc.cores = n_jobs)
  #offsprings2 <- mclapply(as.list(offsprings2), ComputeIndividual, mc.cores = n_jobs)
  
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
  
  #calculating offspring1
  # Ind_index <- 1:length(offsprings1)
  # Marrmot_Ind_index <-Ind_index[sapply(offsprings1, function(x) x[['IndArray']][1,1]=="MARRMot")]
  # Other_Ind_index <- Ind_index[sapply(offsprings1, function(x) x[['IndArray']][1,1]!="MARRMot")]
  # pre_index_order <-c(Marrmot_Ind_index,Other_Ind_index)
  
  Marrmot_Ind <-offsprings1[sapply(offsprings1, function(x) x[['IndArray']][1,1]=="MARRMot")]
  
  Other_Ind <-offsprings1[sapply(offsprings1, function(x) x[['IndArray']][1,1]!="MARRMot")]
  
  # if(length(Other_Ind)>0){
  #   for(vv in 1:length(Other_Ind)){
  #     Other_Ind[[vv]]<- ComputeIndividual(Other_Ind[[vv]])
  #   }
  # }
  
  if(length(Other_Ind)>0){
    Other_Ind <- mclapply(as.list(Other_Ind),ComputeIndividual, mc.cores = n_jobs)
  }
  
  Output_script_vec <- as.character(1:length(Marrmot_Ind))
  
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
  
  offsprings1 <-c(Marrmot_Ind,Other_Ind)
  
  # Index_dataframe <-data.frame(pre_index_order,1:length(offsprings1))
  # names(Index_dataframe)<-c("pre","post")
  # Index_dataframe <-Index_dataframe[order(Index_dataframe$pre),]
  # offsprings1 <-offsprings1[Index_dataframe$post]
  
  #calculating offspring2
  # Ind_index <- 1:length(offsprings2)
  # Marrmot_Ind_index <-Ind_index[sapply(offsprings2, function(x) x[['IndArray']][1,1]=="MARRMot")]
  # Other_Ind_index <- Ind_index[sapply(offsprings2, function(x) x[['IndArray']][1,1]!="MARRMot")]
  # pre_index_order <-c(Marrmot_Ind_index,Other_Ind_index)
  
  Marrmot_Ind <-offsprings2[sapply(offsprings2, function(x) x[['IndArray']][1,1]=="MARRMot")]
  
  Other_Ind <-offsprings2[sapply(offsprings2, function(x) x[['IndArray']][1,1]!="MARRMot")]
  
  # if(length(Other_Ind)>0){
  #   for(vv in 1:length(Other_Ind)){
  #     Other_Ind[[vv]]<- ComputeIndividual(Other_Ind[[vv]])
  #   }
  # }
  
  if(length(Other_Ind)>0){
    Other_Ind <- mclapply(as.list(Other_Ind),ComputeIndividual, mc.cores = n_jobs)
  }
  
  Output_script_vec <- as.character(1:length(Marrmot_Ind))
  
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
  
  offsprings2 <-c(Marrmot_Ind,Other_Ind)
  
  # Index_dataframe <-data.frame(pre_index_order,1:length(offsprings2))
  # names(Index_dataframe)<-c("pre","post")
  # Index_dataframe <-Index_dataframe[order(Index_dataframe$pre),]
  # offsprings2 <-offsprings2[Index_dataframe$post]
  
  
  
  # to avoid subscript out of bounds problem in some generations, let's try following code
  for (i in 1:length(offsprings1)){
    if (length(offsprings1[[i]]) <10){
      offsprings1[[i]] <- population[[round(runif(1,1,PopulationSize))]]
    } else if ((length(offsprings1[[i]]) == 10) & (length(offsprings1[[i]][['Sim']]) ==1)){
      offsprings1[[i]] <- population[[round(runif(1,1,PopulationSize))]]
    } else {
      offsprings1[[i]] <- offsprings1[[i]]
    }
  }
  
  for (i in 1:length(offsprings2)){
    if (length(offsprings2[[i]]) <10){
      offsprings2[[i]] <- population[[round(runif(1,1,PopulationSize))]]
    } else if ((length(offsprings2[[i]]) == 10) & (length(offsprings2[[i]][['Sim']]) ==1)){
      offsprings2[[i]] <- population[[round(runif(1,1,PopulationSize))]]
    } else {
      offsprings2[[i]] <- offsprings2[[i]]
    }
  }
  
  
  #offsprings1<-lapply(offsprings1,ComputeIndividual)
  #offsprings2<-lapply(offsprings2,ComputeIndividual)
  # cl <- makeCluster(n_jobs, 'FORK')
  # 
  # sub_population1 = split(offsprings1,c(1:n_jobs))
  # sub_population2 = split(offsprings2,c(1:n_jobs))
  # # Computation of fitness for both offspring groups (when the tree is changed)
  # ParComputeInd = function(sub_population) lapply(sub_population, ComputeIndividual)
  # sub_population1 = parLapply(cl, sub_population1, ParComputeInd)
  # sub_population2 = parLapply(cl, sub_population2, ParComputeInd)
  # stopCluster(cl)
  # offsprings1 = unlist(sub_population1, recursive = FALSE, use.names = FALSE)
  # offsprings2 = unlist(sub_population1, recursive = FALSE, use.names = FALSE)
  names(offsprings1) = names(offsprings2) = paste0('Off', 1:noIndividuals)
  #offsprings1 = mclapply(offsprings1, function(x) if(x[['Changed']]) ComputeIndividual(x) else x,mc.cores = getOption("mc.cores", 3L))
  #offsprings2 = mclapply(offsprings2, function(x) if(x[['Changed']]) ComputeIndividual(x) else x,mc.cores = getOption("mc.cores", 3L))
  #   browser()
  # isFitMulti = FitnessFunction %in% MultiObjectiveFitnesses
  # if(isFitMulti){
  #   
  #   parFit1 = MultiObjective(parents1)
  #   parFit2 = MultiObjective(parents2)
  #   
  #   parents1 = mapply(function(x, fit) {x[['Fitness']] = fit; x}, parents1, parFit1, SIMPLIFY = FALSE)
  #   parents2 = mapply(function(x, fit) {x[['Fitness']] = fit; x}, parents2, parFit2, SIMPLIFY = FALSE)    
  #   
  #   offFit1 = MultiObjective(offsprings1)
  #   offFit2 = MultiObjective(offsprings2)
  #   
  #   offsprings1 = mapply(function(x, fit) {x[['Fitness']] = fit; x}, offsprings1, offFit1, SIMPLIFY = FALSE)
  #   offsprings2 = mapply(function(x, fit) {x[['Fitness']] = fit; x}, offsprings2, offFit2, SIMPLIFY = FALSE)
  #   
  # }
  #
  nn <-round(runif(1,1,length(FitnessFunction)))
   betterParents = mapply(function(parent1, parent2) list(parent1, parent2)[[which.min(c(parent1[['Fitness']][nn], parent2[['Fitness']][nn]))]], 
                          parents1, parents2, SIMPLIFY = FALSE)
   
   # Determination of best individual i.e. best ind. from one row ofoffspring1, offspring2, betterParent (all together = family, without poorer parent)
   bestFamilyMembers = mapply(function(offspring1, offspring2, betterParent) {
     if(is.na(offspring1[['Fitness']][nn]) & is.na(offspring2[['Fitness']][nn])){
       return(betterParent);
     } 
     list(offspring1, offspring2)[[which.min(c(offspring1[['Fitness']][nn], offspring2[['Fitness']][nn]))]]
   },
   offsprings1, offsprings2, betterParents, SIMPLIFY = FALSE)
   
  # # The fitness must be in form of vector of all multi fitnesses, fit/rank is computed for whole population in main cycle                      
  # if(isFitMulti){
  #   
  #   bestFamilyMembers = lapply(bestFamilyMembers, function(x) ComputeIndividual(x))
  #   
  # }
  # 
   bestFamilyMembers
   
}

TreeMutation = function(individual){
  
  # Determination of mutation node 
  if (individual[['IndArray']][1,1]=="CT_new"){
    mutationNode = list(targetNode = Sample(9:individual[['IndLength']], 1), restriction = FALSE, variationAllowed = TRUE)
  } else {
  mutationNode = list(targetNode = Sample(1:individual[['IndLength']], 1), restriction = FALSE, variationAllowed = TRUE)
  }
  if(mutationNode[['targetNode']] %in% individual[['RestrictedRows']]){
    
    mutationNode[['restriction']] = TRUE
    mutationNode[['variationAllowed']] = FALSE
    
  }
  
  # Determination of mutation node depth
  mutationNodeDepth = NodeDepth(individual, mutationNode[['targetNode']])
  
  # Max. depth of new mutated subtree
  mutationTreeDepth = MaxDepthRun - mutationNodeDepth
  
  # Creation of mutation subtree
  if(mutationNode[['variationAllowed']]){
    
    mutationTree = CreateIndividual(mutationTreeDepth, 'grow', restricted = TRUE)
     #mutationTree = CreateIndividual(mutationTreeDepth, 'grow', restricted = FALSE)
    
  } else{
    
    mutationTree = CreateIndividual(mutationTreeDepth, 'grow', restricted = TRUE)
    mutationNode[['variationAllowed']] = TRUE
    
  }
  
  # Inserting of mutation subtree to mutating individual
  individual = AddSubtree(individual, mutationNode, mutationTree)
  
  # Mutant must be computed
  individual = ChangeIndividual(individual[['IndArray']], individual[['DriveVec']], individual[['IndLength']], TRUE)
  
  individual
  
}

TreeDependentMutatingIndividuals = function(mutatingIndividuals, noMutatingIndividuals, mutationType, variationProbs){
  
  if(mutationType == 'Separation'){
    
    # Individuals lengther than 1 are interesting for separation (else is just replication)
    possibleIndividuals = sapply(mutatingIndividuals, function(x) x[['IndLength']] > 1)
    P = variationProbs[['PmSeparation']]
    
  }
  
  if(mutationType == 'Constant'){
    
    possibleIndividuals = sapply(mutatingIndividuals, function(x) any(!is.na(suppressWarnings(as.numeric(x[['IndArray']][, 1])))))
    P = variationProbs[['PmConstant']]
    
  }
  
  possibleIndividuals[1] = FALSE # First position is elite individual (for sure, minimal probability, that position 1 of population will be in mutatingIndividuals)
  possibleMutants = sum(possibleIndividuals)
  noMutants = min(sum(rbinom(noMutatingIndividuals - 1, 1, P)), possibleMutants)
  MutantsPositions = Sample(which(possibleIndividuals), noMutants)
  
  MutantsPositions
  
}

SeparationMutation = function(individual){
  
  # Determination of mutation node
  if (individual[['IndArray']][1,1]=="CT_new"){
    mutationNode = c(targetNode = Sample(9:individual[['IndLength']], 1),
                     restriction = FALSE,
                     variationAllowed = TRUE)
  } else {
  mutationNode = c(targetNode = Sample(2:individual[['IndLength']], 1),
                   restriction = FALSE,
                   variationAllowed = TRUE)
  }
  
  individual = SepareTrees(individual, mutationNode)
  
  individual = ChangeIndividual(individual[['IndArray']], individual[['DriveVec']], individual[['IndLength']], TRUE)
  
  individual
  
}

NodeMutation = function(individual){
  
  usedArities = FunctionsDefinition[FunctionSet, ]
  
  if (individual[['IndArray']][1,1]=="CT_new") {
    noNodes = Sample(1:(individual[['IndLength']]-8), 1)
    mutationNodes = Sample(9:individual[['IndLength']], noNodes)
  } else {
  noNodes = Sample(1:individual[['IndLength']], 1)
  mutationNodes = Sample(1:individual[['IndLength']], noNodes)
  }
  
  functionsPositions = which(individual[['DriveVec']] == 'Function')
  mutatingFunctionsPositions = functionsPositions[functionsPositions %in% mutationNodes] 
  
  mutatingFunctions = individual[['IndArray']][mutatingFunctionsPositions, 1]  
  
  noNewFunctionNodes = length(mutatingFunctionsPositions)
  
  if(noNewFunctionNodes > 0){
    
    for(i in 1:noNewFunctionNodes){
      
      # Functions with same arities can be in node mutation process (possible problem with functions of same arity, but different strictness)
      possibleFunctions = usedArities[usedArities[ ,'Arity'] == usedArities[mutatingFunctions[i], 'Arity'], 'Functions']    
      individual[['IndArray']][mutatingFunctionsPositions[i], 1] = Sample(possibleFunctions, 1)
      
    }
    
  }
  
  terminalsPositions = which(individual[['DriveVec']] == 'Terminal')
  
  # mutating terminals positions
  mtp = terminalsPositions[terminalsPositions %in% mutationNodes] 
  
  noNewTerminalNodes = length(mtp)
  
  if(noNewTerminalNodes > 0){
    
    possibleConstants = runif(noNewTerminalNodes, ConstantRange[1], ConstantRange[2])
    
    possibleVariables = Sample(IndependentVariables, noNewTerminalNodes, replace = TRUE)
    
    variablePositions = mtp[which(individual[['IndArray']][mtp, 1]  %in% IndependentVariables)]
    constantPositions = mtp[which(!individual[['IndArray']][mtp, 1]  %in% IndependentVariables)]
    
    # number of new Terminal nodes
    nTn = length(variablePositions)
    
    # number of new Constant nodes
    nCn = length(constantPositions)
    
    if(nTn > 0){
      
      individual[['IndArray']][variablePositions, 1] = Sample(c(possibleConstants, possibleVariables), nTn)  
      
    }
    
    if(nCn > 0){
      
      individual[['IndArray']][constantPositions, 1] = Sample(possibleConstants, nCn)  
      
    }
    
  }
  
  individual = ChangeIndividual(individual[['IndArray']], individual[['DriveVec']], individual[['IndLength']], TRUE)
  
  individual
  
}

ConstantMutation = function(individual, constantMutationFactor){
  
  constantsPositions = which(individual[['DriveVec']] == 'Terminal' & !is.na(suppressWarnings(as.numeric(individual[['IndArray']][, 1]))))
  noNodes = Sample(1:length(constantsPositions), 1)
  
  # Selecton of nodes with constants. With repetitions = more mutations
  mutationNodesPositions = Sample(constantsPositions, noNodes)
  
  for(i in mutationNodesPositions){
    
    mutatingConstant = as.numeric(individual[['IndArray']][i, 1])
    actualRange = sort(c(mutatingConstant * (1 + constantMutationFactor), mutatingConstant * (1 - constantMutationFactor)))
    individual[['IndArray']][i, 1] = as.character(runif(1, actualRange[1], actualRange[2]))
    
  } 
  
  individual = ChangeIndividual(individual[['IndArray']], individual[['DriveVec']], individual[['IndLength']], TRUE)
  
  individual
  
}

Mutation = function(mutatingIndividuals, variationProbs, constantMutationFactor){
  
  noMutatingIndividuals = length(mutatingIndividuals)
  #   functionSet = computationVariables[['FunctionSet']]
  
  # Separation - individuals determination
  separationMutantsPositions = TreeDependentMutatingIndividuals(mutatingIndividuals, noMutatingIndividuals, 'Separation', variationProbs)  
  # Separation        
  mutatingIndividuals[separationMutantsPositions] = lapply(mutatingIndividuals[separationMutantsPositions], SeparationMutation)
  
  # Subtree mutation - individuals determination 
  noTreeMutants = sum(rbinom(noMutatingIndividuals, 1, variationProbs[['PmTree']]))
  treeMutantsPositions = Sample(1:noMutatingIndividuals, noTreeMutants)
  
  # Subtree mutation
  mutatingIndividuals[treeMutantsPositions] = lapply(mutatingIndividuals[treeMutantsPositions], TreeMutation)
  # Subtree mutation - control and reparation of tree depth 
  mutatingIndividuals = DepthControl(mutatingIndividuals)
  
  # Node mutation - individuals determination
  noNodeMutants = sum(rbinom(noMutatingIndividuals - 1, 1, variationProbs[['PmNode']]))
  nodeMutantsPositions = Sample(1:noMutatingIndividuals, noNodeMutants)
  
  # Node mutation
  mutatingIndividuals[nodeMutantsPositions] = lapply(mutatingIndividuals[nodeMutantsPositions], NodeMutation)
  
  # Constant mutation - individuals determination
  constantMutantsPositions = TreeDependentMutatingIndividuals(mutatingIndividuals, noMutatingIndividuals, 'Constant', variationProbs)   
  # Constant mutation
  mutatingIndividuals[constantMutantsPositions] = lapply(mutatingIndividuals[constantMutantsPositions], ConstantMutation, 
                                                         constantMutationFactor)
  
  mutantsPositions = unique(c(separationMutantsPositions, treeMutantsPositions, nodeMutantsPositions, constantMutantsPositions))
  mutantsPositions = paste0('I', mutantsPositions)
  # Computation of mutants
  delete.NULLs  <-  function(x.list){   # delele null/empty entries in a list
    x.list[unlist(lapply(x.list, length) != 0)]
  }
  names(mutatingIndividuals) = paste0('I', 1:noMutatingIndividuals)
  
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
  
  #calculating mutatingIndividuals
  # Ind_index <- 1:length(mutatingIndividuals)
  # Marrmot_Ind_index <-Ind_index[sapply(mutatingIndividuals, function(x) x[['IndArray']][1,1]=="MARRMot")]
  # Other_Ind_index <- Ind_index[sapply(mutatingIndividuals, function(x) x[['IndArray']][1,1]!="MARRMot")]
  # pre_index_order <-c(Marrmot_Ind_index,Other_Ind_index)
  
  Marrmot_Ind <-mutatingIndividuals[sapply(mutatingIndividuals, function(x) x[['IndArray']][1,1]=="MARRMot")]
  
  Other_Ind <-mutatingIndividuals[sapply(mutatingIndividuals, function(x) x[['IndArray']][1,1]!="MARRMot")]
  
  # if(length(Other_Ind)>0){
  #   for(vv in 1:length(Other_Ind)){
  #     Other_Ind[[vv]]<- ComputeIndividual(Other_Ind[[vv]])
  #   }
  # }
  
  if(length(Other_Ind)>0){
    Other_Ind <- mclapply(as.list(Other_Ind),ComputeIndividual, mc.cores = n_jobs)
  }
  
  Output_script_vec <- as.character(1:length(Marrmot_Ind))
  
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
  
  mutatingIndividuals <-c(Marrmot_Ind,Other_Ind)
  
  # Index_dataframe <-data.frame(pre_index_order,1:length(mutatingIndividuals))
  # names(Index_dataframe)<-c("pre","post")
  # Index_dataframe <-Index_dataframe[order(Index_dataframe$pre),]
  # mutatingIndividuals <-mutatingIndividuals[Index_dataframe$post]
  
  #mutatingIndividuals <- mclapply(as.list(mutatingIndividuals),ComputeIndividual,mc.cores = n_jobs)
  #mutatingIndividuals<- lapply(mutatingIndividuals,ComputeIndividual)
  # sub_mutationInd = split(mutatingIndividuals, c(1:n_jobs))
  # ParComputeMutInd = function(n) {
  #   sub_mutationInd[[n]] = delete.NULLs(sub_mutationInd[[n]][mutantsPositions])
  #   sub_mutationInd[[n]] = lapply(sub_mutationInd[[n]], ComputeIndividual)
  #   sub_mutationInd[[n]]
  # }
  # cl <- makeCluster(n_jobs,"FORK")
  # sub_mutationInd = parLapply(cl, 1:n_jobs, ParComputeMutInd)
  # stopCluster(cl)
  # 
  # #mutatingIndividuals[mutantsPositions] = mclapply(mutatingIndividuals[mutantsPositions], ComputeIndividual,mc.cores = getOption("mc.cores", 3L))
  # mutatingIndividuals[mutantsPositions] = unlist(sub_mutationInd, recursive = FALSE, use.names = FALSE)
  mutatingIndividuals   
  
}  

controlIndividual = function(x, iv, co) {
  
  a = x[['RestrictedRows']]
  if(!is.null(a) && !is.na(a[1]) && min(a > 3))  {print(x); stop(co)}
  #   if(a[1,1] %in% StrictFunctions & any(tail(a[-1,1]) %in% iv)) {print(x); stop(co)}
  
}

controlIndividualArray = function(x, iv, co) {
  
  a = x[['IndArray']]
  #   if(!is.null(a) && !is.na(a[1]) && min(a > 2))  {print(x); stop(co)}
  if(a[1,1] %in% StrictFunctions & any(tail(a[-2,1]) %in% iv)) {print(x); stop(co)}
  
}
