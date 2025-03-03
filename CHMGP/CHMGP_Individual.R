NodeType = function(method, functionSet, parentFun = NA, terminal = FALSE, restrictedNode = FALSE){
  
  if(terminal && restrictedNode) nodeType = "Constant"
  
  if(terminal && !restrictedNode){
    
    if(is.na(parentFun)) nodeType = Sample(c("Constant", "Variable"), 1)
    nodeType = ifelse(parentFun %in% StrictFunctions, "Variable", Sample(c("Constant", "Variable"), 1))
    
  } 
  
  if(!terminal && restrictedNode){
    # No special functions, because these functions returns a vector!
    # Random choice of function or constant in grow method
    possibleNodes = functionSet[!functionSet %in% SpecialFunctions]
    if(length(possibleNodes) == 0) possibleNodes = c("Constant", "Variable")
    
    if(method == "grow"){
      
      nodeType = Sample(c("Function", "Constant"), 1)
      if(nodeType == "Function") nodeType = Sample(possibleNodes, 1)
      
    }
    
    # Random choice of function or constant in full method
    if(method == "full"){
      
      nodeType = Sample(possibleNodes, 1)
      
    }
    
  }  
  
  if(!terminal && !restrictedNode){
    
    # Boosting probability of Strict function occurence
    stricFuncBoost = rbinom(1,1,0.5)
    if(stricFuncBoost & any(functionSet %in% StrictFunctions)) fn = Sample(functionSet[functionSet %in% StrictFunctions], 1) else fn = Sample(functionSet, 1)
    
    # Random choice of function or terminal in grow method
    if(method == "grow"){
      
      nodeType = Sample(c("Function", "Terminal"), 1)
      
      if(parentFun %in% StrictFunctions) {
        
        nodeType = ifelse(nodeType == "Function", fn, "Variable")      
        
      } else{
        
        nodeType = ifelse(nodeType == "Function", fn, Sample(c("Constant", "Variable"), 1))
        
      }
      
    }
    
    # Random choice of function or constant in full method
    if(method == "full"){
      
      nodeType = fn
      
    }
    
  }
  
  nodeType
  
}

NodeCompletion = function(new_row, nodeType, usedArities, functionSet){
  
  tuple = rep(NA, MaxNoColumns)
  tuple[1] = nodeType
  
  if(any(functionSet == nodeType)) {
    
    if(any(names(VariableArityRanges) == tuple[1])){
      
      funcsArityRanges = VariableArityRanges[[tuple[1]]]
      ## Random number of arities, sequence is created by specific range for TANK, see variableArityRanges in GlobalVar
      n_arit = sample(seq(funcsArityRanges[1], funcsArityRanges[2], funcsArityRanges[1]), 1)
      
    } else {
      
      n_arit = usedArities[which(functionSet == tuple[1])]
      
    }
    
    tuple[2:(n_arit + 1)] = (new_row + 1):(new_row + n_arit)
    
    new_row = new_row + n_arit
    
  }
  
  return(list(tuple = tuple, new_row = new_row))
  
}

InsertTerminals = function(indArray){
  
  constantsPositions = which(indArray[ ,1] == 'Constant')
  lcp = length(constantsPositions)
  
  if(lcp > 0) {
    
    indArray[constantsPositions, 1] = runif(lcp, ConstantRange[1], ConstantRange[2])
    
  }
  
  variablesPositions = which(indArray[ ,1] == 'Variable')
  lvp = length(variablesPositions)
  
  if(lvp > 0) {
    
    indArray[variablesPositions, 1] = Sample(IndependentVariables, lvp, replace = TRUE)
    
  }
  
  indArray
  
}

CreateIndividual = function(maxDepth, method, restricted){
  
  # Arities of used functions
  usedArities = FunctionsDefinition[FunctionSet, 'Arity']
  
  ##########
  ReservoirSet <-c("WR","IR","RR","UR","FR","SR","CR")
  AddusedArities = FunctionsDefinition[ReservoirSet, 'Arity']
  
  #########
  # Test of one node individual
  singleNodeIndividual = ifelse(maxDepth == 0 | is.null(usedArities), TRUE, FALSE)
  
  if(restricted){  
    
    # vector of nonrestricted functions
    funSet = FunctionSet[!FunctionSet %in% StrictFunctions]
    ##### this might be the reason for subscript out of problem, this will added to all codes ###
    usedArities = usedArities[!FunctionSet %in% StrictFunctions]
    
  } else{
    
    funSet = FunctionSet
    usedArities = usedArities
    
  }
  
  if(length(funSet) == 0) singleNodeIndividual = TRUE
  
  # Max. number of rows in individual array
  noRow = ifelse(singleNodeIndividual, 1, sum(max(usedArities)^(0:maxDepth)))
  
  # Individual array preparation
  
  indArray = array(c(NA), dim = c(noRow, MaxNoColumns))
  
  # First node
  nodeType = NodeType(method, funSet, parentFun = NA, terminal = singleNodeIndividual)
  nodeComplete = NodeCompletion(1, nodeType, usedArities, funSet)
  
  indArray[1, ] = nodeComplete[['tuple']]
  
  # Finishing one node individual
  if(singleNodeIndividual){
    
    indArray = InsertTerminals(indArray) 
    
    return(list(IndArray = indArray, Equation = NA, Fitness = c(rep(NA,length(FitnessFunction))), DriveVec = 'Terminal',  IndLength = 1, Changed = TRUE,Front = NA, Crowding_distance = NA,Sim =c())) 
    
  }
  
  # Initial settings of actual depth, row, maximal row number in given depth and determination of pointers of actual node
  actualDepth = 1
  actualRow = 1
  maxRowLevel = 1
  newRow = nodeComplete[['new_row']]
  
  # Vector of array rows which will have restriction due to occurence of strict function (positions of restricted arguments of this strict function)
  transferredRestrictions = numeric()
  
  ########################################
  
  if (indArray[1,1]=="CT_new"){
    while(actualDepth <= maxDepth) {
      if (actualDepth==1){
        actualColumn = 2
        while(actualColumn <= MaxNoColumns && !is.na(indArray[actualRow, actualColumn])){
          actFun = indArray[actualRow, 1]
          actNSA = NonStrictArguments[actFun]
          restrictedNode = (actFun %in% StrictFunctions && actualColumn > (1 + actNSA)) | actualRow %in% transferredRestrictions
          if (actualColumn==2){
            nodeType <- "WR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else if (actualColumn==3) {
            nodeType <- "IR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else if (actualColumn==4) {
            nodeType <- "RR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else if (actualColumn==5) {
            nodeType <- "UR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else if (actualColumn==6) {
            nodeType <- "FR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else if (actualColumn==7) {
            nodeType <- "SR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else if (actualColumn==8) {
            nodeType <- "CR"
            nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
            
          } else {
            #actFun = indArray[actualRow, 1]
            #actNSA = NonStrictArguments[actFun]
            #restrictedNode = (actFun %in% StrictFunctions && actualColumn > (1 + actNSA)) | actualRow %in% transferredRestrictions
            
            if(actualDepth < maxDepth && !restrictedNode){
              
              nodeType = NodeType(method, funSet, parentFun = actFun, terminal = FALSE, restrictedNode = FALSE)
              
            } else if(actualDepth < maxDepth && restrictedNode) {
              
              nodeType = NodeType(method, funSet, parentFun = actFun, terminal = FALSE, restrictedNode = TRUE)
              
            } else if(actualDepth == maxDepth && !restrictedNode){
              
              nodeType = NodeType(method, funSet, parentFun = actFun, terminal = TRUE, restrictedNode = FALSE)
              
            } else if(actualDepth == maxDepth && restrictedNode){
              
              nodeType = NodeType(method, funSet, parentFun = actFun, terminal = TRUE, restrictedNode = TRUE)
              
            }
            nodeComplete = NodeCompletion(newRow, nodeType, usedArities, funSet)
          }
          #nodeComplete = NodeCompletion(newRow, nodeType, AddusedArities, ReservoirSet)
          
          targetRow = as.numeric(indArray[actualRow, actualColumn])
          
          indArray[targetRow, ] = nodeComplete[['tuple']]
          
          if(restrictedNode) transferredRestrictions = c(transferredRestrictions, targetRow)
          
          newRow = nodeComplete[['new_row']]
          
          actualColumn = actualColumn + 1
        }
      } else {
        actualColumn = 2
        while(actualColumn <= MaxNoColumns && !is.na(indArray[actualRow, actualColumn])){
          actFun = indArray[actualRow, 1]
          actNSA = NonStrictArguments[actFun]
          restrictedNode = (actFun %in% StrictFunctions && actualColumn > (1 + actNSA)) | actualRow %in% transferredRestrictions
          
          if(actualDepth < maxDepth && !restrictedNode){
            
            nodeType = NodeType(method, funSet, parentFun = actFun, terminal = FALSE, restrictedNode = FALSE)
            
          } else if(actualDepth < maxDepth && restrictedNode) {
            
            nodeType = NodeType(method, funSet, parentFun = actFun, terminal = FALSE, restrictedNode = TRUE)
            
          } else if(actualDepth == maxDepth && !restrictedNode){
            
            nodeType = NodeType(method, funSet, parentFun = actFun, terminal = TRUE, restrictedNode = FALSE)
            
          } else if(actualDepth == maxDepth && restrictedNode){
            
            nodeType = NodeType(method, funSet, parentFun = actFun, terminal = TRUE, restrictedNode = TRUE)
            
          }
          
          
          nodeComplete = NodeCompletion(newRow, nodeType, usedArities, funSet)
          
          targetRow = as.numeric(indArray[actualRow, actualColumn])
          
          indArray[targetRow, ] = nodeComplete[['tuple']]
          
          if(restrictedNode) transferredRestrictions = c(transferredRestrictions, targetRow)
          
          newRow = nodeComplete[['new_row']]
          
          actualColumn = actualColumn + 1
          
        }
      }
      if(actualRow == maxRowLevel){
        
        if(any(!is.na(as.numeric(indArray[1:actualRow, 2:MaxNoColumns])))){ ## as.numeric pryc
          
          maxRowLevel = max(as.numeric(indArray[1:actualRow, 2:MaxNoColumns]), na.rm = TRUE) ## as.numeric pryc
          
        }
        
        actualDepth = actualDepth + 1
        
      }   
      
      actualRow = actualRow + 1
      
      # Termination of main cycle
      if(is.na(indArray[actualRow, 1])) {
        
        actualDepth = maxDepth + 1
        
      }
      
    }
  } else {
  
  # Construction of the individual array
  while(actualDepth <= maxDepth){
    
    actualColumn = 2
    
    # Construction of target node, eg. on which the actual node points out
    while(actualColumn <= MaxNoColumns && !is.na(indArray[actualRow, actualColumn])){
      
      # some simplification of variables from global variables
      actFun = indArray[actualRow, 1]
      # actNSA = NonStrictArguments[which(actFun == StrictFunctions)]
      actNSA = NonStrictArguments[actFun]
      
      restrictedNode = (actFun %in% StrictFunctions && actualColumn > (1 + actNSA)) | actualRow %in% transferredRestrictions
      
      if(actualDepth < maxDepth && !restrictedNode){
        
        nodeType = NodeType(method, funSet[-1], parentFun = actFun, terminal = FALSE, restrictedNode = FALSE)
        
      } else if(actualDepth < maxDepth && restrictedNode) {
        
        nodeType = NodeType(method, funSet, parentFun = actFun, terminal = FALSE, restrictedNode = TRUE)
        
      } else if(actualDepth == maxDepth && !restrictedNode){
        
        nodeType = NodeType(method, funSet[-1], parentFun = actFun, terminal = TRUE, restrictedNode = FALSE)
        
      } else if(actualDepth == maxDepth && restrictedNode){
        
        nodeType = NodeType(method, funSet, parentFun = actFun, terminal = TRUE, restrictedNode = TRUE)
        
      }
      
      
      nodeComplete = NodeCompletion(newRow, nodeType, usedArities, funSet)
      
      targetRow = as.numeric(indArray[actualRow, actualColumn])
      
      indArray[targetRow, ] = nodeComplete[['tuple']]
      
      if(restrictedNode) transferredRestrictions = c(transferredRestrictions, targetRow)
      
      newRow = nodeComplete[['new_row']]
      
      actualColumn = actualColumn + 1
      
    }
    
    # Change of actual row
    if(actualRow == maxRowLevel){
      
      if(any(!is.na(as.numeric(indArray[1:actualRow, 2:MaxNoColumns])))){ ## as.numeric pryc
        
        maxRowLevel = max(as.numeric(indArray[1:actualRow, 2:MaxNoColumns]), na.rm = TRUE) ## as.numeric pryc
        
      }
      
      actualDepth = actualDepth + 1
      
    }   
    
    actualRow = actualRow + 1
    
    # Termination of main cycle
    if(is.na(indArray[actualRow, 1])) {
      
      actualDepth = maxDepth + 1
      
    }
    
  }
  }
  
  # Shortening of individual
  IndRows = max(which(!is.na(indArray[, 1])))
  indArray = head(indArray, IndRows)
  
  # Driving vector for next use
  driveVec = indArray[, 1]
  driveVec[driveVec == 'Variable' | driveVec == 'Constant'] = 'Terminal'
  driveVec[driveVec != 'Terminal'] = 'Function'
  
  if(restricted){
    
    # Changing all terminals to constants
    indArray[driveVec == 'Terminal', 1] = 'Constant'
    
  }
  
  # Terminals adding
  indArray = InsertTerminals(indArray) 
  
  if(length(transferredRestrictions) == 0) transferredRestrictions = NA
  
  list(IndArray = indArray, Equation = NA, Fitness = c(rep(NA,length(FitnessFunction))), DriveVec = driveVec, RestrictedRows = transferredRestrictions, IndLength = IndRows, Changed = TRUE, Front = NA, Crowding_distance = NA,Sim =c())
  
}
