# InitOutflow - function for runoff response and infiltration (threshold = 0)
InitOutflow = function(storage, threshold, K){

  if(storage > threshold){
  
    out = K * (storage - threshold)
  
  } else {
  
    out = 0
    
  }
  
  out
  
}


isRB = function(x) grepl('RB', x)


isRJOIN = function(x) grepl('RJOIN', x)


howManyRB = function(x) {

  rb = gregexpr('RB', x)[[1]]
  if(rb[1] != -1) y = length(rb) else y = 0
  y
  
}


######### RJOIN FUNCTION #########

RJOIN <- function(x1, x2, x3, x4){
  
  rjoinArgs = length(formals(RJOIN))
  argg = tail(as.character(match.call()), rjoinArgs)
#   if(any(grepl('RJOIN', argg))) stop('Bad call of function RJOIN - there are not four arguments. Vojta')

  # No RB function is presented, RJOIN returns NA
  isRBxs = sapply(argg, isRB, USE.NAMES = FALSE)  
  if(all(!isRBxs)) return(NA)
  
  # No RJOIN function is allowed to be in another RJOIN function.
  isRJOINs = sapply(argg, isRJOIN, USE.NAMES = FALSE) 
  if(any(isRJOINs)) return(NA)    
  
  # Check for serial or parallel structure (This and condition below can be merged into one condition, separate commands are simpler to understand)
  if(sum(isRBxs) == 1) { 
  
    rSerial = TRUE
    rParall = FALSE
    
  } else { 
  
    rSerial = FALSE
    rParall = TRUE
  
  }

  # Parallel tank assign and other settings
  assign('rParall', rParall, envir = parent.env(environment()))
  
  # GOUTF Global OUTflow
  assign('GOUTF', 0, envir = parent.env(environment()))
  
  # Assigning isInRJOIN variable which is needed to not evaluate RB selfstanding
  assign('isInRJOIN', TRUE, envir = parent.env(environment())) 
  
  # This counter is used for parallel and serial aswell - it serves as a counter of depth in series of reservoirs.
  assign('LFRBcounter', 1, envir = parent.env(environment()))   
   
  # Calculation of RB in branches
  nrRB = sapply(argg, howManyRB, USE.NAMES = FALSE)
  
  brSeq = seq_len(rjoinArgs)
  brNams = paste0('branch', brSeq)
  names(nrRB) = names(argg) = brNams  
  
  rjoinOutParall = vector('list', rjoinArgs)
  names(rjoinOutParall) = brNams

  brON = which(isRBxs)
  
#   browser()
  
  ## PARALLEL STRUCTURE
  if(rParall){

    ## Computation of first branch
    
    activeBrNam = brNams[[ brON[1] ]]
    activeBrNrRB = nrRB[[activeBrNam]]     
    assign('LastRBinSeries', activeBrNrRB, envir = parent.env(environment())) 
    
    LFRBout <<- as.data.frame(matrix(NA, nrow = NROW(DataSet), ncol = activeBrNrRB))
     
    assign('firstBranch', TRUE, envir = parent.env(environment()))

    branchCalculation = eval(parse(text = argg[[activeBrNam]]), envir = DataSet)
    if(is.na(branchCalculation[1])) {rm(isInRJOIN, envir = parent.env(environment())); return(NA)}   
    
    rjoinOutParall[[activeBrNam]] = GOUTF
   
    # Restarting GOUTF and is1stRB for next branch, if there would be more branches, 
    # this line must be everywhere to penultimate rjoinOutParall assigning!
    assign('GOUTF', 0, envir = parent.env(environment()))
    rm(is1stRB, envir = parent.env(environment())) 
    
    ## Computation of next branches
    
    for(i in tail(seq_along(brON), -1)){ # positions in brON withou first position
    
      previousBrNam = brNams[[ brON[i-1] ]] 
      activeBrNam = brNams[[ brON[i] ]] 
      
      previousBrNrRB = nrRB[[previousBrNam]]
      activeBrNrRB = nrRB[[activeBrNam]]      
      assign('LastRBinSeries', activeBrNrRB, envir = parent.env(environment())) 
      
      LFRBin  <<- LFRBout
      LFRBout <<- as.data.frame(matrix(NA, nrow = NROW(DataSet), ncol = activeBrNrRB))
      
      # Case of more RB in first branch than in second, pbnf (Previou Branch No-target Flow)
      if(previousBrNrRB > activeBrNrRB) pbnf = rowSums(as.data.frame(LFRBin[ , (activeBrNrRB+1):previousBrNrRB ])) else pbnf = 0
      
      assign('LFRBcounter', 1, envir = parent.env(environment()))    
      assign('firstBranch', FALSE, envir = parent.env(environment()))
      
      branchCalculation = eval(parse(text = argg[[activeBrNam]]), envir = DataSet)
      if(is.na(branchCalculation[1])) {rm(isInRJOIN, envir = parent.env(environment())); return(NA)}     
      
      rjoinOutParall[[activeBrNam]] = GOUTF + pbnf 
    
    }
    
    # If there is parallel structure, overall output is GOUTF from last branch.
    rjoinOut = rjoinOutParall[[tail(brON, 1)]] 
    
  }
  
  ## SERIAL STRUCTURE
  
  if(rSerial){
  
    activeBrNam = brNams[[ brON[1] ]]
    activeBrNrRB = nrRB[[activeBrNam]]    
    assign('LastRBinSeries', activeBrNrRB, envir = parent.env(environment())) 
    
    # If the RB is serial structure, only one argument is returned to evaluation, the other one is omitted
    branchCalculation = eval(parse(text = argg[activeBrNam]), envir = DataSet)
    if(is.na(branchCalculation[1])) return(NA)
    
    rjoinOut = GOUTF
    
  }
  
  # Cleaning, probably useless
  rm(isInRJOIN, envir = parent.env(environment()))
  rm(GOUTF, envir = parent.env(environment()))
  rm(is1stRB, envir = parent.env(environment()))
  
  rjoinOut

}


######### RB FUNCTION #########

RB <- function(x, t1, a1, t2, a2, i1){

  # Checking whether is RB within the RJOIN function
  if(!exists('isInRJOIN')) return(NA)
  if(!isInRJOIN) return(NA)
  
  # This is the MOST CURIOUS THINK. Without this command (or any command with x) here RJOIN with RBs does not work same as TANK.
  RBinflow = x
  
  # Following condition has been a pretty taught task. 
  # Now it works on the basis that in RJOIN no is1stRB is defined, thus first RB assign is1stRB value TRUE. 
  # In any following RB is1stRB alreadyexists, thus is1stRB is always set to FALSE  
  if(!exists('is1stRB')) assign('is1stRB', TRUE, envir = parent.env(environment())) else is1stRB = FALSE  

  if(rParall && !firstBranch){
  
    if(is1stRB) {
 
      RBinflow = LFRBin[[LFRBcounter]] 
      
    } else { 
    
      # This case is for less RBs in previous branch, oposite case is solved in RJOIN
      if(ncol(LFRBin) < LFRBcounter) RBinflow = x else RBinflow = x + LFRBin[[LFRBcounter]]
      
    }

  }
  
  # Checking parameters of RB
  if(any(is.na(c(RBinflow, t1, a1, t2, a2, i1)))) return(NA) 

  # Checking of non variable inputs in RBinflow
  if(length(RBinflow) <= 1) return(NA)
  
  # Scaling parameters of RB
  modPars = c(t1, a1, t2, a2, i1)
  names(modPars) = tail(names(formals(RB)), 5)

  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)

  # Preparation of vectors for reservoir storage (rS) lateral flow (rOutQ) and infiltration from reservoir (rOutInfil)
  n = length(RBinflow)
  rS = rOutQ = rOutInfil = numeric(n)
  
  rS[1] = 20

  rOutQ[1] = InitOutflow(rS[1], t1, a1) + InitOutflow(rS[1], t2, a2)
  rOutInfil[1] = InitOutflow(rS[1], 0, i1)
  rAllOutflow = rOutQ[1] + rOutInfil[1]

  fOutput = .Fortran("TankBlockCycle", S = as.double(rS), x = as.double(RBinflow), n = as.integer(n), is1stReservoir = as.integer(is1stRB), 
                     Q = as.double(rOutQ), totQ = as.double(rAllOutflow), E = as.double(DataSet[['E']]), 
                     infil = as.double(rOutInfil), infilPar = as.double(i1), 
                     a1 = as.double(a1), a2 = as.double(a2), 
                     t1 = as.double(t1), t2 = as.double(t2))   
                      
  if(rParall) {
  
    # Outflow from RB unit. When it is last RB in series, and it is not deeper than 2, the outflow is increased by infiltration from this last unit.
    RBoutflow = fOutput$Q
    
    if(LFRBcounter == LastRBinSeries & LastRBinSeries <= 2) RBoutflow = fOutput$Q + fOutput$infil

    LFRBout[[LFRBcounter]] = RBoutflow

    assign('LFRBout', LFRBout, envir = parent.env(environment()))
    
  }
  
  # Outflow from RB unit. When it is last RB in series, and it is not deeper than 2, the outflow is increased by infiltration from this last unit.
  RBoutflow = fOutput$Q

  if(LFRBcounter == LastRBinSeries & LastRBinSeries <= 2) RBoutflow = fOutput$Q + fOutput$infil
  
  assign('LFRBcounter', LFRBcounter + 1, envir = parent.env(environment())) 
    
  assign('GOUTF', GOUTF + RBoutflow, envir = parent.env(environment()))

  fOutput$infil
  
}
