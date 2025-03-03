DEoptimCHMGP = function(DEoptimization, DEiterMax, DEfitness, yacasSimp, bestEqSimp, bestInd, fitMulti){

  positions = function(target, regexString, firstCoeff){
  
    pos = gregexpr(regexString, target)[[1]]
    posDF = data.frame(position = as.vector(pos), length = attributes(pos)[["match.length"]] - firstCoeff) # -1 is here because of last character which is comma and must not be in selection
    if(posDF[["position"]][1] == -1) posDF = NULL
    posDF
  
  }
  
  if(DEoptimization){
  
    # DE optim on simplified equations (Yacas simplification must be on)
    if(!yacasSimp) stop('Use DE with Yacas Simplification, i.e. YacasSimplification must be TRUE')
    
    ## Add exit command, if no constant is present in best individual ???
    
    cat('DEoptim optimization of numeric values...\n\n')
    
    # Position of first arguments 
    regex4DEargsFirsts  = '\\([0-9]+\\.[0-9]+,|\\([0-9]+,|\\(-[0-9]+\\.[0-9]+,|\\(-[0-9]+,'
    posArgsDfFirsts = positions(bestEqSimp, regex4DEargsFirsts, 1)
    
    # Position of other arguments 
    regex4DEargs = ',[0-9]+\\.[0-9]+|,[0-9]+'
    posArgsDf = positions(bestEqSimp, regex4DEargs, 0)
    
    # Position of non argument constants
    regex4DEconst = '[-*/+ (^][0-9]+\\.[0-9]+|[-*/+ (^][0-9]+'
    posConstDf = positions(bestEqSimp, regex4DEconst, 0)
    
    allConst = unique(rbind(posArgsDfFirsts, posArgsDf, posConstDf)) # unique, because position of non argument constants finds the first argument of SF functions as a constant, dont know how to set it in regex
    
    if(!is.null(allConst)){
    
      allConst = allConst[order(allConst[["position"]]),]
      
      allConst[['type']] = allConst[["position"]] %in% c(as.vector(posArgsDf[['position']]), as.vector(posArgsDfFirsts[['position']]))
      allConst[['type']] = ifelse(allConst[['type']], 'argument', 'constant')  
      
      # Separation of numVals
      nrAllConst = nrow(allConst)
      constRaw = character(nrAllConst) 
      
      for(i in seq_len(nrAllConst)){
      
        from = allConst[['position']][i]
        to = from + allConst[['length']][i] - 1
        if(i == 1 & allConst[['type']][i] == 'argument') constRaw[i] = substr(bestEqSimp, from, to + 1) else constRaw[i] = substr(bestEqSimp, from, to)
        
      } 
      
      numVals = sub('[*/ (,^]', '', constRaw) # 
      numVals = sub('\\(|,', '', numVals) # must be done twice because of strict SF functions first argument 
      numVals = as.numeric(numVals)
      names(numVals) = allConst[['type']]
            
      # Parameter vector
      parvec = character(length(constRaw))
      for (i in seq_along(constRaw)) {
      
        firstSbstr = substring(constRaw[i], 1, 1)
        parvec[i] = paste0(firstSbstr, 'par[', i, ']')
        if(firstSbstr == '(' & allConst[['type']][i] == 'argument') parvec[i] = paste0(parvec[i], ',')
      
      }
      
      # Replacement of constant values by par[x] string, where x = 1,2,..., no. numVals
      parEq = bestEqSimp
      for(i in seq_along(parvec)) parEq = sub(paste0(regex4DEargsFirsts,'|',regex4DEargs,'|',regex4DEconst), parvec[i], parEq)

      fitMeasure = ifelse(DEfitness == 'SameAsGPfit', FitnessFunction, DEfitness)
      
      # Objective function
      optFun = function(par){ 
      
        computedValues = suppressWarnings(try(eval(parse(text = parEq), envir = DataSet)))
        if(class(computedValues) == 'try-error' | any(!is.finite(computedValues))) computedValues = NA           
        fitOpt = FitnessComputation(computedValues, DataSet[[DependentVariable]], fitMeasure)
      
#         if(fitMulti) fitOpt = sum(scale(fitOpt, center = FALSE, scale = FALSE)) ## THIS IS ONLY TEMPORARY, MUST BE PROBABLY DONE IN BETTER WAY - MULTIOBJECTIVE PROBLEM ??? !!!
        
        if(!is.finite(fitOpt)) fitOpt = .Machine$double.xmax
        
        fitOpt
      
      }
      
      lowBound = upBound = numeric(length(parvec))
      
      sar = StrictFunctionsArgumentsRange
      
      lowBound[names(numVals) == 'argument'] = sar[['lower']]
      upBound[names(numVals) == 'argument'] = sar[['upper']]
      
      lowBound[names(numVals) == 'constant'] = numVals[names(numVals) == 'constant'] - 0.5 * abs(numVals['constant'])
      upBound[names(numVals) == 'constant'] = numVals[names(numVals) == 'constant'] + 0.5 * abs(numVals['constant'])
     
      # Initial population for DE optimization        
      nVal2opt = length(numVals)
      NPdeOpt = 10*length(lowBound)
      
      initPop = matrix(numVals + runif(nVal2opt * NPdeOpt, -0.3, 0.3), nrow = NPdeOpt, ncol = nVal2opt, byrow = TRUE)
      initPop[,names(numVals) == 'argument'][initPop[,names(numVals) == 'argument'] < sar[['lower']]] = sar[['lower']] # only the arguments of functions must be in given range
      initPop[,names(numVals) == 'argument'][initPop[,names(numVals) == 'argument'] > sar[['upper']]] = sar[['upper']]
      
      # DE optimization
      DE = DEoptim(optFun, lower = lowBound, upper = upBound, DEoptim.control(itermax = DEiterMax, trace = FALSE, initialpop = initPop))

      # Initial fitness computation
      # initialFitness = ifelse(fitMulti, sum(scale(bestInd[['Fitness']], center = FALSE, scale = FALSE)), bestInd[['Fitness']])      
      compVals = suppressWarnings(try(eval(parse(text = bestEqSimp), envir = DataSet)))
      if(class(compVals) == 'try-error' | any(!is.finite(compVals))) compVals = NA   
      initialFitness = FitnessComputation(compVals, DataSet[[DependentVariable]], fitMeasure)
      
      cat(paste0('Best equation before optimization (',fitMeasure,': ',round(initialFitness, 6),'): ', bestEqSimp, '\n'))
      
      if(DE[[1]][['bestval']] < initialFitness){
      
        for(i in seq_along(parvec)){
        
          regPar = paste0('par\\[', i, '\\]' )
          parEq = sub(regPar, round(DE[[1]][['bestmem']][i], RoundingFactor) , parEq)
        
        }
        
        cat(paste0('Best equation after optimization  (',fitMeasure,': ',round(DE[[1]][['bestval']] , 6),'): ', parEq, '\n'))      
        
        optimizedEquation = parEq        
              
      } else{
      
        optimizedEquation = bestEqSimp
        message('No improvement of efficiency by DEoptim')
      
      }
      
    } else{
    
      optimizedEquation = bestEqSimp
        
    }
  
  } else{
    
    optimizedEquation = bestEqSimp
        
  }
      
  optimizedEquation    
    
}
