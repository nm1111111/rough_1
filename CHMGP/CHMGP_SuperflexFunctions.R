ScaleParams = function(parPrcg){
  
  parPrcg[parPrcg < 0] = 0
  parPrcg[parPrcg > 1] = 1    
  mapply(function(prcg, rng){prcg * (rng[['max']] - rng[['min']]) + rng[['min']]}, parPrcg, SFPR[names(parPrcg)])
  
}

ScaleParamsManualy = function(parPrcg, functionNam){
  
  names(parPrcg) = names(formals(functionNam))
  parPrcg[parPrcg < 0] = 0
  parPrcg[parPrcg > 1] = 1    
  mapply(function(prcg, rng){prcg * (rng[['max']] - rng[['min']]) + rng[['min']]}, parPrcg, SFPR[names(parPrcg)])
  
}


## MI: Fast Reservoir 
MI = function(alpha_Qq_FR, K_Qq_FR, Ce, m_E_FR){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  modPars = c(alpha_Qq_FR, K_Qq_FR, Ce, m_E_FR)
  names(modPars) = names(formals(MI))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MI", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     alpha_Qq_FR = as.double(alpha_Qq_FR), K_Qq_FR  = as.double(K_Qq_FR), Ce  = as.double(Ce), m_E_FR= as.double(m_E_FR), dT = as.double(1),
                     output = as.double(output))
  
  fOutput[['output']]
  
}


## MII: Fast Reservoir 
MII = function(Ce, Beta_Qq_UR, Smax_UR, K_Qb_UR, Beta_E_UR, SiniFr_UR){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  modPars = c(Ce, Beta_Qq_UR, Smax_UR, K_Qb_UR, Beta_E_UR, SiniFr_UR)
  names(modPars) = names(formals(MII))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MII", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     Ce  = as.double(Ce), Beta_Qq_UR= as.double(Beta_Qq_UR), Smax_UR = as.double(Smax_UR), K_Qb_UR = as.double(K_Qb_UR), 
                     Beta_E_UR=as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), dT = as.double(1),
                     output = as.double(output))
  
  fOutput[['output']]
  
}

## M3: Unsaturated and Fast Reservoir with four parameters
MIII = function(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  modPars = c(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR)
  names(modPars) = names(formals(MIII))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MIII", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     alpha_Qq_FR = as.double(alpha_Qq_FR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), Smax_UR = as.double(Smax_UR), Beta_Qq_UR = as.double(Beta_Qq_UR), 
                     Beta_E_UR = as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), dT=as.double(1),
                     output = as.double(output))                  
  
  fOutput[['output']]
  
}


##M4: Unsaturated and Fast Reservoir with five parameters
MIV = function(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  modPars = c(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR)
  names(modPars) = names(formals(MIV))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MIV", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     alpha_Qq_FR = as.double(alpha_Qq_FR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), Smax_UR = as.double(Smax_UR), Beta_Qq_UR = as.double(Beta_Qq_UR), 
                     Beta_E_UR = as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), dT=as.double(1),
                     output = as.double(output))
  
  fOutput[['output']]
  
}

##M5: Unsaturated and Fast Reservoir with six parameters
MV = function(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, Tlag){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  modPars = c(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, Tlag)
  names(modPars) = names(formals(MV))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MV", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     alpha_Qq_FR = as.double(alpha_Qq_FR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), Smax_UR = as.double(Smax_UR), Beta_Qq_UR = as.double(Beta_Qq_UR), 
                     Beta_E_UR = as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), dT=as.double(1), Tlag=as.double(Tlag),
                     output = as.double(output))
  
  fOutput[['output']]
  
}


##M6: Interception, Unsaturated and Fast Reservoir with eight parameters
MVI = function( alpha_Qq_FR, Beta_Qq_UR, K_Qq_FR, Ce, Smax_UR, Smax_IR, m_QE_IR, Beta_E_UR, SiniFr_UR, Tlag){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(alpha_Qq_FR, Beta_Qq_UR, K_Qq_FR, Ce, Smax_UR, Smax_IR, m_QE_IR, Beta_E_UR, SiniFr_UR, Tlag)
  names(modPars) = names(formals(MVI))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars) # head is due to Tlag
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MVI", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     alpha_Qq_FR=as.double(alpha_Qq_FR), Beta_Qq_UR=as.double(Beta_Qq_UR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), Smax_UR=as.double(Smax_UR), 
                     Smax_IR=as.double(Smax_IR), m_QE_IR = as.double(m_QE_IR), Beta_E_UR= as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), dT=as.double(1), Tlag=as.double(Tlag),
                     output = as.double(output))
  
  fOutput[['output']]
  
}


##M7: Unsaturated, Riparian and Fast Reservoir with seven parameters
MVII = function(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, D_R, K_Qq_RR, Tlag){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(alpha_Qq_FR, K_Qq_FR, Ce, Smax_UR, Beta_Qq_UR, Beta_E_UR, SiniFr_UR, D_R, K_Qq_RR, Tlag)
  names(modPars) = names(formals(MVII))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MVII", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     alpha_Qq_FR = as.double(alpha_Qq_FR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), Smax_UR = as.double(Smax_UR), Beta_Qq_UR = as.double(Beta_Qq_UR), 
                     Beta_E_UR = as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), D_R=as.double(D_R), K_Qq_RR=as.double(K_Qq_RR), dT=as.double(1), Tlag=as.double(Tlag),
                     output = as.double(output)) 
  #   ctrlVar1 = ctrlVar2 = numeric(n)                     
  #   fOutput = .Fortran("MVII", P = as.double(x1), E = as.double(x2), n = as.integer(n),
  #                      alpha_Qq_FR = as.double(alpha_Qq_FR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), Smax_UR = as.double(Smax_UR), Beta_Qq_UR = as.double(Beta_Qq_UR), 
  #                      Beta_E_UR = as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), D_R=as.double(D_R), K_Qq_RR=as.double(K_Qq_RR), dT=as.double(1), Tlag=as.double(Tlag),
  #                      output = as.double(output), ctrlVar1 = as.double(ctrlVar1), ctrlVar2 = as.double(ctrlVar2))
  
  #   ControlVariable1 <<- fOutput[['ctrlVar1']]        
  #   ControlVariable2 <<- fOutput[['ctrlVar2']]
  
  fOutput[['output']]
  
}


##M8: Fast and slow Reservoirs
MVIII = function( K_Qq_FR, Ce, K_Qq_SR, D_S, m_E_FR){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(K_Qq_FR, Ce, K_Qq_SR, D_S, m_E_FR)
  names(modPars) = names(formals(MVIII))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MVIII", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), K_Qq_SR= as.double(K_Qq_SR), D_S=as.double(D_S), m_E_FR=as.double(m_E_FR),dT=as.double(1),
                     output = as.double(output))
  
  fOutput[['output']]
  
}


##M9: Unsaturated, Fast and slow Reservoirs with six parameters
MIX = function(Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR)
  names(modPars) = names(formals(MIX))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MIX", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     Beta_Qq_UR=as.double(Beta_Qq_UR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), K_Qq_SR= as.double(K_Qq_SR), 
                     D_S=as.double(D_S), Smax_UR=as.double(Smax_UR), Beta_E_UR= as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), 
                     dT=as.double(1), output = as.double(output))
  
  fOutput[['output']]
  
}


##M10: Unsaturated, Fast and slow Reservoirs with five parameters
MX = function(K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, Tlag){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, Tlag)
  names(modPars) = names(formals(MX))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MX", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), K_Qq_SR= as.double(K_Qq_SR), 
                     D_S=as.double(D_S), Smax_UR=as.double(Smax_UR), Beta_E_UR= as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), 
                     dT= as.double(1), Tlag=as.double(Tlag), output = as.double(output))
  
  fOutput[['output']]
  
}

##M11: Unsaturated, Fast and slow Reservoirs with six parameters
MXI = function(Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, Tlag){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, Smax_UR, Beta_E_UR, SiniFr_UR, Tlag)
  names(modPars) = names(formals(MXI))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MXI", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     Beta_Qq_UR=as.double(Beta_Qq_UR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce), K_Qq_SR= as.double(K_Qq_SR), 
                     D_S=as.double(D_S), Smax_UR=as.double(Smax_UR), Beta_E_UR= as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), 
                     dT= as.double(1), Tlag=as.double(Tlag), output = as.double(output))
  
  fOutput[['output']]
  
}


##M12: Interception, Unsaturated, Fast and slow Reservoirs with eight parameters
MXII = function(Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, SiniFr_UR, Smax_UR, Smax_IR, m_QE_IR, Beta_E_UR, Tlag){
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  #   modPars = unlist(as.list(match.call())[-(1:3)])
  modPars = c(Beta_Qq_UR, K_Qq_FR, Ce, K_Qq_SR, D_S, SiniFr_UR, Smax_UR, Smax_IR, m_QE_IR, Beta_E_UR, Tlag)
  names(modPars) = names(formals(MXII))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  n = length(x1)
  # if(n == 1) return(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("MXII", P = as.double(x1), E = as.double(x2), n = as.integer(n),
                     Beta_Qq_UR=as.double(Beta_Qq_UR), K_Qq_FR= as.double(K_Qq_FR), Ce  = as.double(Ce),
                     K_Qq_SR= as.double(K_Qq_SR), D_S=as.double(D_S), SiniFr_UR=as.double(SiniFr_UR), Smax_UR=as.double(Smax_UR), 
                     Smax_IR=as.double(Smax_IR), m_QE_IR = as.double(m_QE_IR), Beta_E_UR= as.double(Beta_E_UR), 
                     dT=as.double(1), Tlag=as.double(Tlag), output = as.double(output))
  
  fOutput[['output']]
  
}



### combined_tank - CT
CT = function(Ce, Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR, K_Qq_FR, K_Qb_FR,  m_E_FR, alpha_Qq_FR, Beta_Qq_UR, Smax_UR, Beta_E_UR,
              SiniFr_UR, K_Qb_UR, mu_Qq_UR, Smax_IR, m_QE_IR, K_Qq_RR, Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR, K_Qq_SR, m_E_SR, alpha_Qq_SR, P_ED_max, m_P_ED, 
              D_S, D_I, D_F, D_R, D_C, Tlag, Lag_RR, Lag_FR, Lag_SR, option_i, option_u, option_f, option_s, option_c, w_res, i_res, r_res, u_res, f_res,s_res,c_res)
{
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  x3 = DataSet[['T']]
  
  modPars = c(Ce, Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR, K_Qq_FR, K_Qb_FR,  m_E_FR, alpha_Qq_FR, Beta_Qq_UR, Smax_UR, Beta_E_UR,
              SiniFr_UR, K_Qb_UR, mu_Qq_UR, Smax_IR, m_QE_IR, K_Qq_RR, Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR, K_Qq_SR, m_E_SR, alpha_Qq_SR, P_ED_max, m_P_ED, 
              D_S, D_I, D_F, D_R, D_C, Tlag, Lag_RR, Lag_FR, Lag_SR, option_i, option_u, option_f, option_s, option_c, w_res, i_res, r_res, u_res, f_res,s_res,c_res)
  
  names(modPars) = names(formals(CT))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  # Scaling parameters
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  # Scaling option and res parameters
  optresNam = tail(names(modPars), 15)
  sapply(optresNam, function(nm) assign(nm, round(get(nm)), envir = parent.env(environment())))
  
  ress = c(w = w_res, i = i_res, r = r_res, u = u_res, f = f_res, s = s_res,c = c_res)
  
  # Testing if function is from meaningfull combination
  #   print(0)
  testMeaningfulness = sapply(MeaningfullCTcombinations, function(x) identical(unlist(x),ress))
  if(!any(testMeaningfulness)) return(NA)
  
  # Turning on and off the arguments in dependency of res values
  nam2zero = unlist(DepArgsOfCT[!ress])
  sapply(nam2zero, function(nm) assign(nm, 0, envir = parent.env(environment())))
  #   print(1)
  
  n = length(x1)
  
  output = numeric(n)
  
  fOutput = .Fortran("combined_tank", P = as.double(x1), E = as.double(x2), T = as.double(x3), n = as.integer(n),
                     Ce = as.double(Ce),Cp_WR = as.double(Cp_WR), m_Q_WR = as.double(m_Q_WR), Kq_WR = as.double(Kq_WR), Tp_WR = as.double(Tp_WR), Tm_WR = as.double(Tm_WR), K_Qq_FR= as.double(K_Qq_FR), K_Qb_FR= as.double(K_Qb_FR), m_E_FR=as.double(m_E_FR), alpha_Qq_FR=as.double(alpha_Qq_FR), Beta_Qq_UR=as.double(Beta_Qq_UR),
                     Smax_UR=as.double(Smax_UR), Beta_E_UR= as.double(Beta_E_UR), SiniFr_UR=as.double(SiniFr_UR), K_Qb_UR= as.double(K_Qb_UR), mu_Qq_UR=as.double(mu_Qq_UR), Smax_IR=as.double(Smax_IR),
                     m_QE_IR = as.double(m_QE_IR), K_Qq_RR=as.double(K_Qq_RR),Smax_CR=as.double(Smax_CR), Umax_uCR=as.double(Umax_uCR), Smin_uCR=as.double(Smin_uCR), Beta_Qq_uCR=as.double(Beta_Qq_uCR), mu_Qq_uCR=as.double(mu_Qq_uCR), K_Qb_uCR=as.double(K_Qb_uCR), Sevmax_CR=as.double(Sevmax_CR), Beta_E_CR=as.double(Beta_E_CR), Beta_Qq_sCR=as.double(Beta_Qq_sCR), K_Qb_sCR=as.double(K_Qb_sCR), K_Qd_sCR=as.double(K_Qd_sCR), K_Qq_SR= as.double(K_Qq_SR), m_E_SR= as.double(m_E_SR), alpha_Qq_SR= as.double(alpha_Qq_SR),
                     P_ED_max= as.double(P_ED_max), m_P_ED= as.double(m_P_ED), 
                     D_S=as.double(D_S), D_I=as.double(D_I), D_F=as.double(D_F), D_R=as.double(D_R),D_C=as.double(D_C),
                     Tlag=as.double(Tlag), dT=as.double(1),Lag_RR=as.integer(Lag_RR),Lag_FR=as.integer(Lag_FR),Lag_SR=as.integer(Lag_SR),
                     option_i=as.integer(option_i),option_u=as.integer(option_u),option_f=as.integer(option_f), option_s=as.integer(option_s),option_c=as.integer(option_c), 
                     w_res=as.integer(w_res), i_res=as.integer(i_res), r_res=as.integer(r_res), u_res=as.integer(u_res), f_res=as.integer(f_res),  s_res=as.integer(s_res),c_res=as.integer(c_res),
                     output = as.double(output))
  
  fOutput[['output']]
  
}

### FUSE

FUSE = function(rferr_add,rferr_mlt, maxwatr_1 ,maxwatr_2, fracten, frchzne, fprimqb, rtfrac1, percrte, percexp,sacpmlt, sacpexp, percfrac, iflwrte, 
                baserte, qb_powr, qb_prms, qbrate_2a, qbrate_2b, sareamax, axv_bexp, loglamb, tishape, timedelay,rferr, arch1, arch2, qsurf, qperc, esoil, qintf, q_tdh)
{
  
  x1 = DataSet[['P']]
  x2 = DataSet[['E']]
  
  modPars = c(rferr_add,rferr_mlt, maxwatr_1 ,maxwatr_2, fracten, frchzne, fprimqb, rtfrac1, percrte, percexp,sacpmlt, sacpexp, percfrac, iflwrte, 
              baserte, qb_powr, qb_prms, qbrate_2a, qbrate_2b, sareamax, axv_bexp, loglamb, tishape, timedelay,rferr, arch1, arch2, qsurf, qperc, esoil, qintf, q_tdh)
  
  names(modPars) = names(formals(FUSE))
  
  if(any(is.na(c(modPars)))) return(NA)
  
  scaledPars = ScaleParams(modPars)
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
  
  optresNam = tail(names(modPars), 8)
  sapply(optresNam, function(nm) assign(nm, round(get(nm)), envir = parent.env(environment())))
  
  n = length(x1)
  
  output = numeric(n)
  
  library(fuse)
  data("modlist")
  DATA <- data.frame(P=x1, E=x2)
  mid <-which(modlist[['rferr']]==rferr & modlist[['arch1']]==arch1 & modlist[['arch2']]==arch2 & modlist[['qsurf']]==qsurf & modlist[['qperc']]==qperc & modlist[['esoil']]==esoil & modlist[['qintf']]==qintf & modlist[['q_tdh']]==q_tdh)
  parameters <- data.frame(rferr_add= rferr_add,rferr_mlt=rferr_mlt, maxwatr_1=maxwatr_1 ,maxwatr_2=maxwatr_2, fracten=fracten, frchzne=frchzne, fprimqb=fprimqb, rtfrac1=rtfrac1, percrte=percrte, percexp=percexp,sacpmlt=sacpmlt, sacpexp=sacpexp, percfrac=percfrac, iflwrte=iflwrte,baserte=baserte, qb_powr=qb_powr, qb_prms=qb_prms, qbrate_2a=qbrate_2a, qbrate_2b=qbrate_2b, sareamax=sareamax, axv_bexp=axv_bexp, loglamb=loglamb, tishape=tishape, timedelay=timedelay)
  output <- fuse(DATA,mid,1,parameters)
  return(output)
  
}

####### DCT #####
DCT <- function(Ce1,Cp_WR1, m_Q_WR1, Kq_WR1, Tp_WR1, Tm_WR1, K_Qq_FR1, K_Qb_FR1,  m_E_FR1, alpha_Qq_FR1, Beta_Qq_UR1, Smax_UR1, Beta_E_UR1,
                SiniFr_UR1, K_Qb_UR1, mu_Qq_UR1, Smax_IR1, m_QE_IR1, K_Qq_RR1, Smax_CR1, Umax_uCR1, Smin_uCR1, Beta_Qq_uCR1, mu_Qq_uCR1, K_Qb_uCR1, Sevmax_CR1, Beta_E_CR1, Beta_Qq_sCR1, K_Qb_sCR1, K_Qd_sCR1, K_Qq_SR1, m_E_SR1, alpha_Qq_SR1, P_ED_max1, m_P_ED1,
                D_S1, D_I1, D_F1, D_R1, D_C1, Tlag1, Lag_RR1, Lag_FR1, Lag_SR1, option_i1, option_u1, option_f1, option_s1, option_c1,w_res1, i_res1, r_res1, u_res1, f_res1,s_res1,c_res1,Ce2,Cp_WR2, m_Q_WR2, Kq_WR2, Tp_WR2, Tm_WR2, K_Qq_FR2, K_Qb_FR2,  m_E_FR2, alpha_Qq_FR2, Beta_Qq_UR2, Smax_UR2, Beta_E_UR2,
                SiniFr_UR2, K_Qb_UR2, mu_Qq_UR2, Smax_IR2, m_QE_IR2, K_Qq_RR2, Smax_CR2, Umax_uCR2, Smin_uCR2, Beta_Qq_uCR2, mu_Qq_uCR2, K_Qb_uCR2, Sevmax_CR2, Beta_E_CR2, Beta_Qq_sCR2, K_Qb_sCR2, K_Qd_sCR2, K_Qq_SR2, m_E_SR2, alpha_Qq_SR2, P_ED_max2, m_P_ED2,
                D_S2, D_I2, D_F2, D_R2, D_C2, Tlag2, Lag_RR2, Lag_FR2, Lag_SR2, option_i2, option_u2, option_f2, option_s2, option_c2,w_res2, i_res2, r_res2, u_res2, f_res2,s_res2,c_res2,Ce3,Cp_WR3, m_Q_WR3, Kq_WR3, Tp_WR3, Tm_WR3, K_Qq_FR3, K_Qb_FR3,  m_E_FR3, alpha_Qq_FR3, Beta_Qq_UR3, Smax_UR3, Beta_E_UR3,
                SiniFr_UR3, K_Qb_UR3, mu_Qq_UR3, Smax_IR3, m_QE_IR3, K_Qq_RR3, Smax_CR3, Umax_uCR3, Smin_uCR3, Beta_Qq_uCR3, mu_Qq_uCR3, K_Qb_uCR3, Sevmax_CR3, Beta_E_CR3, Beta_Qq_sCR3, K_Qb_sCR3, K_Qd_sCR3, K_Qq_SR3, m_E_SR3, alpha_Qq_SR3, P_ED_max, m_P_ED,
                D_S3, D_I3, D_F3, D_R3, D_C3, Tlag3, Lag_RR3, Lag_FR3, Lag_SR3, option_i3, option_u3, option_f3, option_s3, option_c3,w_res3, i_res3, r_res3, u_res3, f_res3,s_res3,c_res3,R_t1,R_t2,R_t3,R_1,R_2,R_3){

  modPars = c(Ce1,Cp_WR1, m_Q_WR1, Kq_WR1, Tp_WR1, Tm_WR1, K_Qq_FR1, K_Qb_FR1,  m_E_FR1, alpha_Qq_FR1, Beta_Qq_UR1, Smax_UR1, Beta_E_UR1,
              SiniFr_UR1, K_Qb_UR1, mu_Qq_UR1, Smax_IR1, m_QE_IR1, K_Qq_RR1, Smax_CR1, Umax_uCR1, Smin_uCR1, Beta_Qq_uCR1, mu_Qq_uCR1, K_Qb_uCR1, Sevmax_CR1, Beta_E_CR1, Beta_Qq_sCR1, K_Qb_sCR1, K_Qd_sCR1, K_Qq_SR1, m_E_SR1, alpha_Qq_SR1, P_ED_max1, m_P_ED1,
              D_S1, D_I1, D_F1, D_R1, D_C1, Tlag1, Lag_RR1, Lag_FR1, Lag_SR1, option_i1, option_u1, option_f1, option_s1, option_c1,w_res1, i_res1, r_res1, u_res1, f_res1,s_res1,c_res1,Ce2,Cp_WR2, m_Q_WR2, Kq_WR2, Tp_WR2, Tm_WR2, K_Qq_FR2, K_Qb_FR2,  m_E_FR2, alpha_Qq_FR2, Beta_Qq_UR2, Smax_UR2, Beta_E_UR2,
              SiniFr_UR2, K_Qb_UR2, mu_Qq_UR2, Smax_IR2, m_QE_IR2, K_Qq_RR2, Smax_CR2, Umax_uCR2, Smin_uCR2, Beta_Qq_uCR2, mu_Qq_uCR2, K_Qb_uCR2, Sevmax_CR2, Beta_E_CR2, Beta_Qq_sCR2, K_Qb_sCR2, K_Qd_sCR2, K_Qq_SR2, m_E_SR2, alpha_Qq_SR2, P_ED_max2, m_P_ED2,
              D_S2, D_I2, D_F2, D_R2, D_C2, Tlag2, Lag_RR2, Lag_FR2, Lag_SR2, option_i2, option_u2, option_f2, option_s2, option_c2,w_res2, i_res2, r_res2, u_res2, f_res2,s_res2,c_res2,Ce3,Cp_WR3, m_Q_WR3, Kq_WR3, Tp_WR3, Tm_WR3, K_Qq_FR3, K_Qb_FR3,  m_E_FR3, alpha_Qq_FR3, Beta_Qq_UR3, Smax_UR3, Beta_E_UR3,
              SiniFr_UR3, K_Qb_UR3, mu_Qq_UR3, Smax_IR3, m_QE_IR3, K_Qq_RR3, Smax_CR3, Umax_uCR3, Smin_uCR3, Beta_Qq_uCR3, mu_Qq_uCR3, K_Qb_uCR3, Sevmax_CR3, Beta_E_CR3, Beta_Qq_sCR3, K_Qb_sCR3, K_Qd_sCR3, K_Qq_SR3, m_E_SR3, alpha_Qq_SR3, P_ED_max, m_P_ED,
              D_S3, D_I3, D_F3, D_R3, D_C3, Tlag3, Lag_RR3, Lag_FR3, Lag_SR3, option_i3, option_u3, option_f3, option_s3, option_c3,w_res3, i_res3, r_res3, u_res3, f_res3,s_res3,c_res3,R_t1,R_t2,R_t3,R_1,R_2,R_3)

  names(modPars) = names(formals(DCT))

  if(any(is.na(c(modPars)))) return(NA)

  # Scaling parameters
  scaledPars = ScaleParams(modPars[169:174])
  mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)

  # Scaling option and res parameters
  optresNam = tail(names(modPars), 3)
  sapply(optresNam, function(nm) assign(nm, round(get(nm)), envir = parent.env(environment())))

  flow_1 <- CT1(Ce1,Cp_WR1, m_Q_WR1, Kq_WR1, Tp_WR1, Tm_WR1, K_Qq_FR1, K_Qb_FR1,  m_E_FR1, alpha_Qq_FR1, Beta_Qq_UR1, Smax_UR1, Beta_E_UR1,
                SiniFr_UR1, K_Qb_UR1, mu_Qq_UR1, Smax_IR1, m_QE_IR1, K_Qq_RR1, Smax_CR1, Umax_uCR1, Smin_uCR1, Beta_Qq_uCR1, mu_Qq_uCR1, K_Qb_uCR1, Sevmax_CR1, Beta_E_CR1, Beta_Qq_sCR1, K_Qb_sCR1, K_Qd_sCR1, K_Qq_SR1, m_E_SR1, alpha_Qq_SR1, P_ED_max1, m_P_ED1,
                D_S1, D_I1, D_F1, D_R1, D_C1, Tlag1, Lag_RR1, Lag_FR1, Lag_SR1, option_i1, option_u1, option_f1, option_s1, option_c1,w_res1, i_res1, r_res1, u_res1, f_res1,s_res1,c_res1)

  flow_2 <- CT2(Ce2,Cp_WR2, m_Q_WR2, Kq_WR2, Tp_WR2, Tm_WR2, K_Qq_FR2, K_Qb_FR2,  m_E_FR2, alpha_Qq_FR2, Beta_Qq_UR2, Smax_UR2, Beta_E_UR2,
                SiniFr_UR2, K_Qb_UR2, mu_Qq_UR2, Smax_IR2, m_QE_IR2, K_Qq_RR2, Smax_CR2, Umax_uCR2, Smin_uCR2, Beta_Qq_uCR2, mu_Qq_uCR2, K_Qb_uCR2, Sevmax_CR2, Beta_E_CR2, Beta_Qq_sCR2, K_Qb_sCR2, K_Qd_sCR2, K_Qq_SR2, m_E_SR2, alpha_Qq_SR2, P_ED_max2, m_P_ED2,
                D_S2, D_I2, D_F2, D_R2, D_C2, Tlag2, Lag_RR2, Lag_FR2, Lag_SR2, option_i2, option_u2, option_f2, option_s2, option_c2,w_res2, i_res2, r_res2, u_res2, f_res2,s_res2,c_res2)

  flow_3 <- CT3(Ce3,Cp_WR3, m_Q_WR3, Kq_WR3, Tp_WR3, Tm_WR3, K_Qq_FR3, K_Qb_FR3,  m_E_FR3, alpha_Qq_FR3, Beta_Qq_UR3, Smax_UR3, Beta_E_UR3,
                SiniFr_UR3, K_Qb_UR3, mu_Qq_UR3, Smax_IR3, m_QE_IR3, K_Qq_RR3, Smax_CR3, Umax_uCR3, Smin_uCR3, Beta_Qq_uCR3, mu_Qq_uCR3, K_Qb_uCR3, Sevmax_CR3, Beta_E_CR3, Beta_Qq_sCR3, K_Qb_sCR3, K_Qd_sCR3, K_Qq_SR3, m_E_SR3, alpha_Qq_SR3, P_ED_max, m_P_ED,
                D_S3, D_I3, D_F3, D_R3, D_C3, Tlag3, Lag_RR3, Lag_FR3, Lag_SR3, option_i3, option_u3, option_f3, option_s3, option_c3,w_res3, i_res3, r_res3, u_res3, f_res3,s_res3,c_res3)

  n = length(flow_1)

  output = numeric(n)

  routed_flow_1 <- fuserouting.sim(flow_1,R_1,1,R_t1)
  routed_flow_2 <- fuserouting.sim(flow_2,R_2,1,R_t2)
  routed_flow_3 <- fuserouting.sim(flow_3,R_3,1,R_t3)

  output <- routed_flow_1 + routed_flow_2 + routed_flow_3

  return(output)

}

#########################
WR <- function(Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR){
  modPars = c(Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR)
  names(modPars) = names(formals(WR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR))
}

IR <- function(Smax_IR, m_QE_IR, D_I, option_i){
  modPars = c(Smax_IR, m_QE_IR, D_I, option_i)
  names(modPars) = names(formals(IR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(Smax_IR, m_QE_IR, D_I, option_i))
}

RR <- function(K_Qq_RR, D_R, Lag_RR){
  modPars = c(K_Qq_RR, D_R, Lag_RR)
  names(modPars) = names(formals(RR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(K_Qq_RR, D_R, Lag_RR))
}

UR <- function(Beta_Qq_UR, Smax_UR, Beta_E_UR, SiniFr_UR, K_Qb_UR, mu_Qq_UR, option_u){
  modPars = c(Beta_Qq_UR, Smax_UR, Beta_E_UR, SiniFr_UR, K_Qb_UR, mu_Qq_UR, option_u)
  names(modPars) = names(formals(UR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(Beta_Qq_UR, Smax_UR, Beta_E_UR, SiniFr_UR, K_Qb_UR, mu_Qq_UR, option_u))
}

FR <- function(K_Qq_FR, K_Qb_FR,  m_E_FR, alpha_Qq_FR, D_F, Lag_FR, option_f){
  modPars = c(K_Qq_FR, K_Qb_FR,  m_E_FR, alpha_Qq_FR, D_F, Lag_FR, option_f)
  names(modPars) = names(formals(FR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(K_Qq_FR, K_Qb_FR,  m_E_FR, alpha_Qq_FR, D_F, Lag_FR, option_f))
}

SR <- function(K_Qq_SR, m_E_SR, alpha_Qq_SR, Lag_SR, option_s ){
  modPars = c(K_Qq_SR, m_E_SR, alpha_Qq_SR, Lag_SR, option_s )
  names(modPars) = names(formals(SR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(K_Qq_SR, m_E_SR, alpha_Qq_SR, Lag_SR, option_s))
}

CR <- function(Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR, D_C, option_c){
  modPars = c(Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR, D_C, option_c)
  names(modPars) = names(formals(CR))
  if(any(is.na(c(modPars)))) return(NA)
  
  return(c(Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR, D_C, option_c)) 
}

CT_new <- function(WR, IR, RR, UR, FR, SR, CR, Ce, P_ED_max, m_P_ED, Tlag, D_S, w_res, i_res, r_res, u_res, f_res, s_res, c_res){
  
  modPars = c(WR, IR, RR, UR, FR, SR, CR, Ce, P_ED_max, m_P_ED, Tlag, D_S, w_res, i_res, r_res, u_res, f_res, s_res, c_res)
  names(modPars) = names(formals(CT_new))
  if(any(is.na(c(modPars)))) return(NA)
  
  WR_aug <- WR#(Cp_WR, m_Q_WR, Kq_WR, Tp_WR, Tm_WR, w_res)
  IR_aug <- IR#(Smax_IR, m_QE_IR, D_I, option_i, i_res)
  RR_aug <- RR#(K_Qq_RR, D_R, Lag_RR, r_res)
  UR_aug <- UR#(Beta_Qq_UR, Smax_UR, Beta_E_UR, SiniFr_UR, K_Qb_UR, mu_Qq_UR, option_u, u_res)
  FR_aug <- FR#(K_Qq_FR, K_Qb_FR,  m_E_FR, alpha_Qq_FR, D_F, Lag_FR, option_f, f_res)
  SR_aug <- SR#(K_Qq_SR, m_E_SR, alpha_Qq_SR, Lag_SR, option_s, s_res)
  CR_aug <- CR#(Smax_CR, Umax_uCR, Smin_uCR, Beta_Qq_uCR, mu_Qq_uCR, K_Qb_uCR, Sevmax_CR, Beta_E_CR, Beta_Qq_sCR, K_Qb_sCR, K_Qd_sCR, D_C, option_c, c_res)
  
  output <- CT(Ce, WR_aug[1],WR_aug[2],WR_aug[3],WR_aug[4],WR_aug[5], FR_aug[1],FR_aug[2],FR_aug[3],FR_aug[4], UR_aug[1],UR_aug[2],UR_aug[3],UR_aug[4],UR_aug[5],UR_aug[6], IR_aug[1],IR_aug[2], RR_aug[1], CR_aug[1],
               CR_aug[2],CR_aug[3],CR_aug[4],CR_aug[5],CR_aug[6],CR_aug[7],CR_aug[8],CR_aug[9],CR_aug[10],CR_aug[11], SR_aug[1],SR_aug[2],SR_aug[3], P_ED_max, m_P_ED, D_S, IR_aug[3], FR_aug[5], RR_aug[2], 
               CR_aug[12], Tlag, RR_aug[3], FR_aug[6], SR_aug[4], IR_aug[4], UR_aug[7], FR_aug[7], SR_aug[5], CR_aug[13], w_res, i_res, r_res, u_res, f_res, s_res, c_res)
  
  return(output)
}


#################  MARRMoT Models ###########################

# Model 29 - HYMOD

# MM_29 <- function(M29_smax,M29_b,M29_a,M29_kf,M29_ks,M29_s1,M29_s2,M29_s3,M29_s4,M29_s5){
#   modPars = c(M29_smax,M29_b,M29_a,M29_kf,M29_ks,M29_s1,M29_s2,M29_s3,M29_s4,M29_s5)
#   names(modPars) = names(formals(MM_29))
#   if(any(is.na(c(modPars)))) return(NA)
#   scaledPars = ScaleParams(modPars)
#   mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(scaledPars), scaledPars)
#   
#   M29_theta <- c(M29_smax,M29_b,M29_a,M29_kf,M29_ks)
#   M29_s0 <- c(M29_s1*M29_smax,M29_s2,M29_s3,M29_s4,M29_s5)
#   M29_theta_mat <- rvec_to_matlab(M29_theta)
#   M29_s0_mat <- rvec_to_matlab(M29_s0)
#   
#   file_location <- paste(path,"/Matlab_fn",sep = '')
#   output_script_name <- as.character(round(runif(1,1,1e6)))
#   
#   code <- c(paste("cd '",noquote(file_location),"';",sep = ''),"Input_Var = readmatrix('Input_Variables.csv');","input_climatology.precip = Input_Var(:,3);","input_climatology.temp   = Input_Var(:,5);",
#             "input_climatology.pet    = Input_Var(:,4);","input_climatology.delta_t  = 1;","model = 'm_29_hymod_5p_5s';",paste("input_theta =",noquote(M29_theta_mat)),paste("input_s0 =",noquote(M29_s0_mat)),
#             "input_solver.name              = 'createOdeApprox_IE';","input_solver.resnorm_tolerance = 0.1;","input_solver.resnorm_maxiter   = 6;",
#             "[output_ex,...                                                              
#  output_in,...                                                              
#  output_ss,....                                                             
#  output_waterbalance] = ...                                                 
#                     feval(model,...                                         
#                           input_climatology,...                             
#                           input_s0,...                                      
#                           input_theta,...                                   
#                           input_solver);",paste("writetable(struct2table(output_ex), '",noquote(file_location),"/output_Q_MM29_",noquote(output_script_name),".csv","')",sep = ''))
#   
#   run_matlab_code(code)
#   matlab_out <-read.csv(paste(file_location,"/output_Q_MM29_",output_script_name,".csv",sep = ''))
#   unlink(paste(file_location,"/output_Q_MM29_",output_script_name,".csv",sep = ''))
#   matlab_out <- t(matlab_out)
#   sim <- as.vector(matlab_out[((nrow(matlab_out)/2)+1):nrow(matlab_out),])
#   return(sim)
# }

#write.csv(M29_parameters,file.path(file_location,"Input_Parameters_MM29.csv"))
# if(have_matlab()){
#   get_matlab(try_defaults = TRUE, desktop = FALSE, splash = FALSE, display = FALSE, wait = TRUE)
# }
# system(paste("matlab -wait -nodesktop -nosplash -nodisplay -r \"run('",file_location,"/MM_29.m'); exit\"",sep = ''))

#M29_parameters <- cbind(M29_theta,M29_s0)
#paste("matlab -wait -nodesktop -nosplash -nodisplay -r \"run('",file_location,"/MM_29.m'); exit\"",sep = '')

# code <- c("cd 'C:/Users/E0149661/Desktop/R/Code_Server/MARMoT/Matlab_fn';","Input_Var = readmatrix('Input_Variables.csv');","input_climatology.precip = Input_Var(:,3);","input_climatology.temp   = Input_Var(:,5);",
#           "input_climatology.pet    = Input_Var(:,4);","input_climatology.delta_t  = 1;","model = 'm_29_hymod_5p_5s';","input_theta     = [ 35;3.7;0.4;0.25;0.01];","input_s0       = [15;7;3;8;22];",
#           "input_solver.name              = 'createOdeApprox_IE';","input_solver.resnorm_tolerance = 0.1;","input_solver.resnorm_maxiter   = 6;",
#           "[output_ex,...                                                              
#  output_in,...                                                              
#  output_ss,....                                                             
#  output_waterbalance] = ...                                                 
#                     feval(model,...                                         
#                           input_climatology,...                             
#                           input_s0,...                                      
#                           input_theta,...                                   
#                           input_solver);","writetable(struct2table(output_ex),'output_Q_MM29.csv')")


