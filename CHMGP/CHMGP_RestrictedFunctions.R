# Outflow - function for runoff response and infiltration (threshold = 0)
Outflow = function(storage, threshold, K){

  if(storage > threshold){
  
    out = K * (storage - threshold)
  
  } else {
  
    out = 0
    
  }
  
  out
  
}

# Mutable state of runoff total amount of whole serie of tanks
Qglobal = function(){

  Qg = 0
  
  function(changeQ){
  
    Qg <<- Qg + changeQ
    
    Qg
   
  }

}


# Input - outflow (I-Q) block
IQblock = function(input, parVec){

#   print(input)

  t1 =  parVec[['t1']]
  a1 =  parVec[['a1']] 
  t2 =  parVec[['t2']]
  a2 =  parVec[['a2']] 
  infilPar = parVec[['a23']]

  n = length(input)
  S = runoff = infil = numeric(n)
  
#   S[1] = 1  
  S[1] = 20
    
  runoff[1] = Outflow(S[1], t1, a1) + Outflow(S[1], t2, a2)
  infil[1] = Outflow(S[1], 0, infilPar)
  totaloutflow = runoff[1] + infil[1]
  
  is1stReservoir = QGLOBAL(0) == 0 # a little weird, possible source of errors

  fOutput = .Fortran("TankBlockCycle", S = as.double(S), x = as.double(input), n = as.integer(n), is1stReservoir = as.integer(is1stReservoir), 
                     Q = as.double(runoff), totQ = as.double(totaloutflow), E = as.double(DataSet[['E']]), 
                     infil = as.double(infil), infilPar = as.double(infilPar), 
                     a1 = as.double(a1), a2 = as.double(a2), 
                     t1 = as.double(t1), t2 = as.double(t2))
                      
  QGLOBAL(fOutput[['Q']])

  fOutput[['infil']]
  
}

ControlAndChanges = function(parlst, len){
# browser()
  # control of proper length (maybe not necessary, solved in creation of individual array???), for block with 5 parameters
  if(!len %% 5 == 0) print('Problem with parameters number in TANK function')
  stopifnot(len %% 5 == 0)
  
  # control of proper length (maybe not necessary after strict functions)
  nonScalars = sapply(parlst, function(x) length(x) != 1)  
  # change of vectors to scalars with value 0 (maybe not necessary after strict functions)
  parlst[nonScalars] = 0
  
  # absolute value of parameters
#   parlst = lapply(parlst, abs)

#   # positions without threshold values - need to be set in interval (0, 1), thresholds have range (0, Inf) # IMPORTANT, let it commented only not remove
#   noTpositons = cumsum(c(rep(c(2,2,1), len / 5)))

#   parlst[noTpositons] = lapply(parlst[noTpositons], function(x) {x[x > 1] = 1; x})  


  # changing parameters to given range
  if(exists('StrictFunctionsArgumentsRange')) SFAR = StrictFunctionsArgumentsRange else SFAR = c(lower = 0, upper = 1)

  parlst = lapply(parlst, function(x) {x[x > SFAR[['upper']]] = SFAR[['upper']]; x})  
  parlst = lapply(parlst, function(x) {x[x < SFAR[['lower']]] = SFAR[['lower']]; x})  
  
  unlist(parlst)

}


TANK = function(...){
  
  # It is OK, to have it here, TankBlockCycle function from f90 file is prepared to work with it.
  RItank = DataSet[['P']] - DataSet[['E']]
#   RItank = c(2.00000, 1.56000, 1.43200, 1.52240, 1.76568)
  
#   if(length(data) == 1) return(data)
  
  parList = list(...)
  
  lenParList = length(parList)
    
  allPars = ControlAndChanges(parList, lenParList)
   
  if(any(sapply(parList, function(x) any(is.na(x))))) return(NA)
  
  if(any(sapply(parList, function(x) any(length(x) > 1 )))) return(NA)
  
  tankPars = lapply(as.data.frame(matrix(allPars, ncol = lenParList / 5)), function(x) x)
  
  tankPars = lapply (tankPars, function(x) {names(x) = c('t1','a1','t2','a2','a23'); x[c('t1','t2')] = ScaleParams(x[c('t1','t2')]); x})
#   bltxt = 'data'
  bltxt = 'RItank'
  
  infilPars = sapply(tankPars, '[[', 'a23')
  seqForTANKeq = seq(min(which(infilPars == 0), length(tankPars)))

  for(i in seqForTANKeq){
  
    bltxt = paste0('IQblock(', bltxt, ', tankPars[[', i, ']])')    
    
  }
 
  QGLOBAL <<- Qglobal()
  
  lastInfil = eval(parse(text = bltxt))
  
#   print(lastInfil)
    
  QGLOBAL(0)  

}



## RESERVOIR FUNCTIONS

CreateVector = function(par, n){
  
  if(length(par) == 1) par = rep(par, n)
  par
  
}

## R2T: Reservoir function two tanks with threshold

R2T = function(x, t1, a1, t2, a2, a23, a3){

  if(any(is.na(c(x, t1, a1, t2, a2, a23, a3)))) return(NA)
  
  scaledPars = ScaleParams(c(t1 = t1, t2 = t2))
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x)
  if(n == 1) return(x)
  
  reservoirParsList = list(t1 = t1, a1 = a1, t2 = t2, a2 = a2, a23 = a23, a3 = a3)
  
  reservoirParsList = lapply(reservoirParsList, abs)
  
  if(exists('StrictFunctionsArgumentsRange')) SFAR = StrictFunctionsArgumentsRange else SFAR = c(lower = 0, upper = 1)

  reservoirParsList[c(1:6)] = lapply(reservoirParsList[c(1:6)], function(x) {x[x > SFAR[['upper']]] = SFAR[['upper']]; x})
  reservoirParsList[c(1:6)] = lapply(reservoirParsList[c(1:6)], function(x) {x[x < SFAR[['lower']]] = SFAR[['lower']]; x})
  
  # Create vectors - to be deleted after strict functions will work
  reservoirParsList = lapply(reservoirParsList, CreateVector, n)
  
  toAssign = names(reservoirParsList)
  
  for(i in toAssign){
  
    assign(i, reservoirParsList[[i]])
  
  }
  
  output = numeric(n)

  fOutput = .Fortran("R2T", x = as.double(x), n = as.integer(n), 
                     t1 = as.double(t1), a1 = as.double(a1), 
                     t2 = as.double(t2), a2 = as.double(a2), 
                     a23 = as.double(a23), a3 = as.double(a3),
                     output = as.double(output))
                       
  fOutput[['output']]

}


## R4T

R4T = function(x, t1, a1, t2, a2, a23, t3, a3, a34, t4, a4, a45, a5){

  if(any(is.na(c(x, t1, a1, t2, a2, a23, t3, a3, a34, t4, a4, a45, a5)))) return(NA)
  
  scaledPars = ScaleParams(c(t1 = t1, t2 = t2, t3 = t3, t4 = t4))
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x)
  if(n == 1) return(x)
  
  reservoirParsList = list(t1 = t1, a1 = a1, t2 = t2, a2 = a2, a23 = a23, 
                           t3 = t3, a3 = a3, a34 = a34, 
                           t4 = t4, a4 = a4, a45 = a45, 
                           a5 = a5 )
  
  reservoirParsList = lapply(reservoirParsList, abs)
#   browser()
  if(exists('StrictFunctionsArgumentsRange')) SFAR = StrictFunctionsArgumentsRange else SFAR = c(lower = 0, upper = 1)

  reservoirParsList[c(1:12)] = lapply(reservoirParsList[c(1:12)], function(x) {x[x > SFAR[['upper']]] = SFAR[['upper']]; x})
  reservoirParsList[c(1:12)] = lapply(reservoirParsList[c(1:12)], function(x) {x[x < SFAR[['lower']]] = SFAR[['lower']]; x})
  
  reservoirParsList = lapply(reservoirParsList, CreateVector, n)
  
  toAssign = names(reservoirParsList)
  
  for(i in toAssign){
  
    assign(i, reservoirParsList[[i]])
  
  }
   
  output = numeric(n)

  fOutput = .Fortran("R4T", x = as.double(x), n = as.integer(n), 
                     t1  = as.double(t1), a1  = as.double(a1), t2  = as.double(t2), a2  = as.double(a2), a23 = as.double(a23), 
                     t3  = as.double(t3), a3  = as.double(a3), a34 = as.double(a34), 
                     t4  = as.double(t4), a4  = as.double(a4), a45 = as.double(a45), 
                     a5  = as.double(a5), 
                     output = as.double(output))
                       
  fOutput[['output']]

}
 
