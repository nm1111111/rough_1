CHMGPpredict = function(model, data, testing = FALSE){
  
  if(testing){
    SFPR = SetSFPR()
    assign('SFPR', SFPR, envir = parent.env(environment())) 
    assign('DataSet', data, envir = parent.env(environment()))  # this is due to variables were set as strict
    DepArgsOfCT = DepArgsOfCTfun()
    assign('DepArgsOfCT', DepArgsOfCT, envir = parent.env(environment()))
    MeaningfullCTcombinations = MeaningfulCTcombsFun()
    assign('MeaningfullCTcombinations', MeaningfullCTcombinations, envir = parent.env(environment()))   
  }
  
  computedValues = try(eval(parse(text = model), envir = data))
  
  if(class(computedValues) == 'try-error' | any(!is.finite(computedValues))){ 
    
    computedValues = rep(NA,nrow(data))
    
  }
  
  ## Threshold of Q, for simulation of Q only!!!
  computedValues[computedValues < 0] = 0
  
  computedValues
  
}
FcriteriaPrecomputation = function(obs, sim, fitF){
  
  segmentFlows = function(obs, sim, p1, p2 = NA){
    
    seg_1 = quantile(obs, 1 - p1, na.rm = TRUE, names = FALSE)
    if(!is.na(p2)) seg_2 = quantile(obs, 1 - p2, na.rm = TRUE, names = FALSE)
    
    if(!is.na(p2)) index = which(obs >= seg_1 & obs <= seg_2) else index = which(obs >= seg_1)
    
    if(!is.na(p2)) {
      
      return(list(Qo_seg = obs[index], Qs_seg = sim[index], seg1 = seg_1, seg2 = seg_2))
      
    } else {
      
      return(list(Qo_seg = obs[index], Qs_seg = sim[index], seg1 = seg_1))
      
    }
    
  }
  
  Qoh = Qsh = Qom = Qsm = Qoi = Qsi = lQo_i1 = lQo_i2 = lQs_i1 = lQs_i2 = Qol = Qsl = lQol = lQsl = lQoL = lQsL = NA
  
  # High flows
  if(fitF == 'FHV' | fitF == 'FDC'){ 
    
    p_high = 0.02 # percentage of time that a flow is equalled or exceeded
    QsegsH = segmentFlows(obs, sim, p_high)
    Qoh = QsegsH[['Qo_seg']]
    Qsh = QsegsH[['Qs_seg']] 
    
  }
  
  # Medium flows
  if(fitF == 'FMV' | fitF == 'FDC'){ 
    
    p_medium1 = 0.2
    p_medium2 = 0.02
    QsegsM = segmentFlows(obs, sim, p_medium1, p_medium2)
    Qom = QsegsM[['Qo_seg']]
    Qsm = QsegsM[['Qs_seg']] 
    
  }
  
  # Intermediate flows  
  if(fitF == 'FMS' |fitF == 'FDC'){ 
    
    p_intermediate1 = 0.7 #m1
    p_intermediate2 = 0.2 #m2   
    QsegsI = segmentFlows(obs, sim, p_intermediate1, p_intermediate2)
    Qoi = QsegsI[['Qo_seg']]
    Qsi = QsegsI[['Qs_seg']]  
    
    if(fitF == 'FMS' |fitF == 'FDC'){
      
      Qo_i1 = quantile(obs, 1-p_intermediate1, na.rm = TRUE, names = FALSE)
      Qo_i2 = quantile(obs, 1-p_intermediate2, na.rm = TRUE, names = FALSE)
      Qs_i1 = quantile(sim, 1-p_intermediate1, na.rm = TRUE, names = FALSE)
      Qs_i2 = quantile(sim, 1-p_intermediate2, na.rm = TRUE, names = FALSE)
      
      ## Safe logarithm of Qo, Qs
      lQo_i1 = Qlog(Qo_i1)
      lQo_i2 = Qlog(Qo_i2)
      lQs_i1 = Qlog(Qs_i1)
      lQs_i2 = Qlog(Qs_i2)
      
    }
    
  }
  
  # Low flows  
  if(fitF == 'FLV' | fitF == 'FDC'){ 
    
    p_low1 = 0.9 ## WHY 1 ??? #l 
    p_low2 = 0.7 #L  
    if(p_low2 > p_low1 | (1 - p_low2) < 0.1) stop('Bad settings of FLV criteria (p_low1, p_low2)')
    QsegsL = segmentFlows(obs, sim, p_low1, p_low2)    
    Qol = QsegsL[['Qo_seg']]
    Qsl = QsegsL[['Qs_seg']]   
    
    ## Safe logarithm of Qol,L, Qsl,
    if(fitF == 'FLV' | fitF == 'FDC' ){
      
      QoL = QsegsL[['seg1']]
      QsL = quantile(sim, 1 - p_low1, na.rm = TRUE, names = FALSE)
      if(QsL == 0) QsL = quantile(sim[sim > 0], 0.1, na.rm = TRUE, names = FALSE)
      
      lQol = Qlog(Qol)
      lQsl = Qlog(Qsl)
      lQoL = Qlog(QoL)
      lQsL = Qlog(QsL)
      
    }
    
  }
  
  list(Qoh = Qoh, Qsh = Qsh, Qom = Qom, Qsm = Qsm, Qoi = Qoi, Qsi = Qsi, 
       lQo_i1 = lQo_i1, lQo_i2 = lQo_i2, lQs_i1 = lQs_i1, lQs_i2 = lQs_i2,
       Qol = Qol, Qsl = Qsl, 
       lQol = lQol, lQsl = lQsl, lQoL = lQoL, lQsL = lQsL)
  
}

## CED computation function
CEDcomputation = function(Qobs, Qsim){
  
  # empty vector
  vector.is.empty <- function(x) return(length(x) ==0 )
  
  Qobs = sort(Qobs, decreasing = TRUE)
  Qsim = sort(Qsim, decreasing = TRUE)
  
  Qsim_min=min(Qsim)
  Qsim_max=max(Qsim)
  
  Qsim2 =quantile(Qsim,0.98) # 2 percentile of simulated flows
  Qsim20=quantile(Qsim,0.80) # 20 percentile of simulated flows
  Qsim70=quantile(Qsim,0.30) # 70 percentile of simulated flows
  
  unscaled_entropy=function(start_range, end_range, Bins, tol) 
  {
    
    
    if(Bins==1)
    {
      Hu_obs = 0
      Hu_sim = 0
      
    } else {
      
      ftshistQobs=binposQobs=ftshistQsim=binposQsim=numeric()
      
      Bins_Qobs  = floor(min (Bins, (Qobs[ceiling(start_range*lenQo)]-Qobs[ceiling(end_range*lenQo)])/(2*tol)))
      Bins_Qsim  = floor(min (Bins, (Qsim[ceiling(start_range*lenQo)]-Qsim[ceiling(end_range*lenQo)])/(2*tol)))
      
      if(Bins_Qobs==0)
      {
        ftshistQobs=ftshistQobs
        binposQobs=binposQobs
      }else{
        a=min(Qobs[ceiling(start_range*lenQo):ceiling(end_range*lenQo)])
        b=max(Qobs[ceiling(start_range*lenQo):ceiling(end_range*lenQo)])
        w1=hist(Qobs[ceiling(start_range*lenQo):ceiling(end_range*lenQo)],seq(a,b,length.out=Bins_Qobs+1),plot=FALSE)
        ftshistQobs=w1[['counts']]
        binposQobs =w1[['mids']]
      }
      
      if(Bins_Qsim==0)
      {
        ftshistQsim=ftshistQsim
        binposQsim=binposQsim 
      }else{
        c=min(Qsim[ceiling(start_range*lenQo):ceiling(end_range*lenQo)])
        d=max(Qsim[ceiling(start_range*lenQo):ceiling(end_range*lenQo)])
        w2=hist(Qsim[ceiling(start_range*lenQo):ceiling(end_range*lenQo)], seq(c,d,length.out=Bins_Qsim+1),plot=FALSE)
        ftshistQsim=w2[['counts']]
        binposQsim=w2[['mids']]
      }
      
      
      Hu_obs = (-sum((ftshistQobs/sum(ftshistQobs))*log2(ftshistQobs/sum(ftshistQobs)),na.rm=TRUE)) / log2(Bins_Qobs)
      if (is.nan(Hu_obs)){Hu_obs=0}
      if (vector.is.empty(ftshistQsim)){ # In case all ftshistQsim values are 0
        Hu_sim = 0
      } else {
        Hu_sim = (-sum((ftshistQsim/sum(ftshistQsim))*log2(ftshistQsim/sum(ftshistQsim)),na.rm=TRUE)) / log2(Bins_Qsim)
      }
      if (is.nan(Hu_sim)){Hu_sim=0}
    }
    list(Hu_obs, Hu_sim)
  }
  
  
  scaled_entropy=function(quant1,quant2,quant1_sim,quant2_sim,start_range,end_range,Bins,tol) 
  {
    MAX = quant1+tol
    MIN = quant2-tol
    Bins = floor(min (Bins, (MAX-MIN)/(2*tol)))
    
    a=min(c(MIN,MAX))
    b=max(c(MIN,MAX))
    w1=hist(c(MIN,MAX),seq(a,b,length.out=Bins+1),plot=FALSE)
    ftshistQobs=w1[['counts']]
    binposQobs=w1[['mids']]
    
    binposQobsN = binposQobs - (binposQobs[2]-binposQobs[1])/2
    binposQobsN[length(binposQobsN)+1] = binposQobsN[length(binposQobsN)] + (binposQobs[2]-binposQobs[1])
    
    if (quant2_sim<MIN)
    {
      binposQobsN = c(quant2_sim,binposQobsN)
      Bins = Bins + 1
    }
    
    if (MAX < max(quant1_sim,MAX))
    {
      binposQobsN = c(binposQobsN,max(quant1_sim,MAX)+tol)
      Bins = Bins + 1
    }
    
    ftshistQobs=length(Bins)
    ftshistQsim=length(Bins)
    for (m in 1 : Bins)
    {
      ftshistQobs[m] = length(which((Qobs[ceiling(start_range*lenQo):ceiling(end_range*lenQo)] >=binposQobsN[m])&(Qobs[ceiling(start_range*lenQo):ceiling(end_range*lenQo)]<binposQobsN[m+1])))
      ftshistQsim[m] = length(which((Qsim>=binposQobsN[m])&(Qsim<binposQobsN[m+1])))
    }
    
    J=which(ftshistQobs<=0)
    if (vector.is.empty(I))
    {
      ftshistQobs = ftshistQobs
    } else {
      ftshistQobs = ftshistQobs[-J]
    }
    
    I=which(ftshistQsim<=0)
    if (vector.is.empty(I))
    {
      ftshistQsim = ftshistQsim
    } else {
      ftshistQsim  = ftshistQsim[-I]
    }
    
    Hs_obs = (-sum((ftshistQobs/sum(ftshistQobs))*log2(ftshistQobs/sum(ftshistQobs)),na.rm=TRUE)) / log2(Bins)
    if (is.nan(Hs_obs)){Hs_obs=0}
    if (vector.is.empty(ftshistQsim)) #In case all ftshistQsim values are 0
    {
      Hs_sim = 0
    } else {
      Hs_sim = (-sum((ftshistQsim/sum(ftshistQsim))*log2(ftshistQsim/sum(ftshistQsim)),na.rm=TRUE)) / log2(Bins)
    }
    if (is.nan(Hs_sim)){Hs_sim=0}
    list(Hs_obs, Hs_sim)  
  }
  
  
  # equal probable binning
  #scaled entropy
  scaled_entropy_low=function(quant1,quant2,quant1_sim,Bins,tol)
  {
    ftshistQobs=ftshistQsim=binposQobsN=as.numeric(vector())
    MAX = quant1 + tol
    MIN = max(0,quant2- tol)
    Bins = floor(min (Bins, (MAX-MIN)/(2*tol)))
    
    L_F = Qobs[ceiling(0.7*lenQo):length(Qobs)]
    SQobs = sort(L_F [L_F >=0])
    SQobs[1] = MIN
    SQobs[length(SQobs)] = MAX
    temp=floor(length(L_F)/Bins)
    temp1=length(SQobs)
    temp2=seq(1,temp1,by=temp)   
    if(Bins<=1){
      binposQobsN =binposQobsN  
    } else{
      binposQobsN = SQobs[temp2]
      binposQobsN = t(binposQobsN)
    }
    
    if (MIN > 0 && (min(Qsim)<MIN)){
      binposQobsN = c(0,binposQobsN)
      Bins = Bins + 1
    }
    
    if (MAX < max(quant1_sim,MAX)){
      binposQobsN = c(binposQobsN,max(Qsim70,MAX)+tol)
      Bins = Bins + 1
    }
    
    if (Bins+1!=length(binposQobsN)) {
      Bins
      length(binposQobsN)
    }
    
    
    for (m in 1 : Bins)
    {
      ftshistQobs[m] = length(which((L_F>=binposQobsN[m])&(L_F<binposQobsN[m+1])))
      ftshistQsim[m] = length(which((Qsim>=binposQobsN[m])&(Qsim<binposQobsN[m+1])))
    }
    
    J=which(ftshistQobs<=0)
    if (vector.is.empty(I))
    {
      ftshistQobs = ftshistQobs
    } else {
      ftshistQobs = ftshistQobs[-J]
    }
    
    I=which(ftshistQsim<=0)
    if (vector.is.empty(I))
    {
      ftshistQsim = ftshistQsim
    } else {
      ftshistQsim  = ftshistQsim[-I]
    }
    
    Hs_obs = (-sum((ftshistQobs/sum(ftshistQobs))*log2(ftshistQobs/sum(ftshistQobs)),na.rm=TRUE)) / log2(Bins)
    if (is.nan(Hs_obs)){Hs_obs=0}
    if (vector.is.empty(ftshistQsim)){ # In case all ftshistQsim values are 0
      Hs_sim = 0
    } else {
      Hs_sim = (-sum((ftshistQsim/sum(ftshistQsim))*log2(ftshistQsim/sum(ftshistQsim)),na.rm=TRUE)) / log2(Bins)
    }
    if (is.nan(Hs_sim)){Hs_sim=0}
    list(Hs_obs, Hs_sim) 
    
  }
  
  #SUSE
  SUSE = function(Hu_sim, Hu_obs, Hs_sim, Hs_obs){
    
    unscaled  = abs(Hu_sim - Hu_obs)
    scaled    = abs(Hs_sim - Hs_obs)
    criterion = max(unscaled,scaled)
    criterion
    
  }
  
  
  #computation
  #unscaled entropy
  
  Hu_obs_high=unscaled_entropy(1/lenQo,0.02,OptBins_high,0.001)[[1]]
  Hu_sim_high=unscaled_entropy(1/lenQo,0.02,OptBins_high,0.001)[[2]]
  
  Hu_obs_medium=unscaled_entropy(0.02,0.2,OptBins_medium,0.001)[[1]]
  Hu_sim_medium=unscaled_entropy(0.02,0.2,OptBins_medium,0.001)[[2]]
  
  Hu_obs_intermediate=unscaled_entropy(0.2,0.7,OptBins_intermediate,0.001)[[1]]
  Hu_sim_intermediate=unscaled_entropy(0.2,0.7,OptBins_intermediate,0.001)[[2]]
  
  Hu_obs_low=unscaled_entropy(0.7,1,OptBins_low,0.001)[[1]]
  Hu_sim_low=unscaled_entropy(0.7,1,OptBins_low,0.001)[[2]]
  
  #scaled entropy
  
  Hs_obs_high=scaled_entropy(Qobs_max,Qobs2,Qsim_max,Qsim2,1/lenQo,0.02,OptBins_high,0.001)[[1]]
  Hs_sim_high=scaled_entropy(Qobs_max,Qobs2,Qsim_max,Qsim2,1/lenQo,0.02,OptBins_high,0.001)[[2]]
  
  
  Hs_obs_medium=scaled_entropy(Qobs2,Qobs20,Qsim2,Qsim20,0.02,0.2,OptBins_medium,0.001)[[1]]
  Hs_sim_medium=scaled_entropy(Qobs2,Qobs20,Qsim2,Qsim20,0.02,0.2,OptBins_medium,0.001)[[2]]
  
  Hs_obs_intermediate=scaled_entropy(Qobs20,Qobs70,Qsim20,Qsim70,0.2,0.7,OptBins_intermediate,0.001)[[1]]
  Hs_sim_intermediate=scaled_entropy(Qobs20,Qobs70,Qsim20,Qsim70,0.2,0.7,OptBins_intermediate,0.001)[[2]]
  
  Hs_obs_low=scaled_entropy_low(Qobs70,Qobs_min,Qsim70,OptBins_low,0.001)[[1]]
  Hs_sim_low=scaled_entropy_low(Qobs70,Qobs_min,Qsim70,OptBins_low,0.001)[[2]]
  
  # SUSE
  SUSE_high         = SUSE(Hu_sim_high, Hu_obs_high, Hs_sim_high, Hs_obs_high)
  SUSE_medium       = SUSE(Hu_sim_medium, Hu_obs_medium, Hs_sim_medium, Hs_obs_medium)
  SUSE_intermediate = SUSE(Hu_sim_intermediate, Hu_obs_intermediate, Hs_sim_intermediate, Hs_obs_intermediate)
  SUSE_low          = SUSE(Hu_sim_low, Hu_obs_low, Hs_sim_low, Hs_obs_low)
  
  
  #CED
  CED = max(SUSE_high, SUSE_medium, SUSE_intermediate, SUSE_low)
  
  CED 
  
}


## CED computation function whole set
CEDcomputation_all = function(Qobs, Qsim){
  
  # empty vector
  vector.is.empty <- function(x) return(length(x) ==0 )
  
  unscaled_entropy = function(observed,simulated,NoBins){
    
    a = observed[observed >= 0]
    w1 = hist(a, seq(min(a), max(a), length.out = NoBins + 2), plot=FALSE) 
    ftshistobs = w1[['counts']]
    binposobs  = w1[['mids']]
    
    b = simulated[simulated>=0]
    w2 = hist(b, seq(min(b), max(b), length.out = NoBins + 2), plot=FALSE) 
    ftshistsim = w2[['counts']]
    binpossim  = w2[['mids']]
    
    J=which(ftshistobs<=0)
    if (vector.is.empty(I))
    {
      ftshistobs = ftshistobs
    } else {
      ftshistobs = ftshistobs[-J]
    }
    
    I=which(ftshistsim<=0)
    if (vector.is.empty(I))
    {
      ftshistsim = ftshistsim
    } else {
      ftshistsim  = ftshistsim[-I]
    }
    
    Hu_obs = (-sum((ftshistobs/sum(ftshistobs)) * log2(ftshistobs / sum(ftshistobs)),na.rm=TRUE)) / log2(NoBins)
    Hu_sim = (-sum((ftshistsim/sum(ftshistsim)) * log2(ftshistsim / sum(ftshistsim)),na.rm=TRUE)) / log2(NoBins)
    
    
    list(Hu_obs,Hu_sim)
    
  }
  
  scaled_entropy = function(observed,simulated,NoBins){
    
    ftshistobs = ftshistsim = as.numeric(vector())
    
    MAX = max(max(observed), max(simulated))
    a = observed[observed >= 0]
    b = simulated[simulated >= 0]
    MIN = min(min(a),min(b))
    
    vec = c(MIN,MAX)
    w3 = hist(vec, seq(min(vec), max(vec), length.out=NoBins+2), plot=FALSE)
    ftshistobs = w3[['counts']]
    binposobs = w3[['mids']]
    binposobsN = binposobs - (binposobs[2] - binposobs[1]) / 2
    end = length(binposobsN)
    binposobsN[end] = binposobsN[end] * 1.00001
    
    for (m in 1 : length(binposobsN)-1){
      
      ftshistobs[m] = length(which((observed >= binposobsN[m]) & (observed < binposobsN[m+1])))
      ftshistsim[m] = length(which((simulated >= binposobsN[m]) & (simulated < binposobsN[m+1])))
      
    }
    
    J=which(ftshistobs<=0)
    if (vector.is.empty(I))
    {
      ftshistobs = ftshistobs
    } else {
      ftshistobs = ftshistobs[-J]
    }
    
    I=which(ftshistsim<=0)
    if (vector.is.empty(I))
    {
      ftshistsim = ftshistsim
    } else {
      ftshistsim  = ftshistsim[-I]
    }
    
    #scaled entropy
    Hs_obs = (-sum((ftshistobs / sum(ftshistobs)) * log2(ftshistobs / sum(ftshistobs)), na.rm=TRUE)) / log2(NoBins)
    Hs_sim = (-sum((ftshistsim / sum(ftshistsim)) * log2(ftshistsim / sum(ftshistsim)), na.rm=TRUE)) / log2(NoBins)
    
    list(Hs_obs,Hs_sim)
    
  }
  
  #SUSE
  SUSE = function(Hu_sim, Hu_obs, Hs_sim, Hs_obs){
    
    unscaled  = abs(Hu_sim - Hu_obs)
    scaled    = abs(Hs_sim - Hs_obs)
    criterion = max(unscaled,scaled)
    criterion
    
  }
  
  #computation
  #unscaled entropy
  Hu_obs        = unscaled_entropy(Qobs, Qsim, OptBins)[[1]]
  Hu_sim        = unscaled_entropy(Qobs, Qsim, OptBins)[[2]]
  
  
  #scaled entropy
  Hs_obs        = scaled_entropy(Qobs, Qsim, OptBins)[[1]]
  Hs_sim        = scaled_entropy(Qobs, Qsim, OptBins)[[2]]
  
  
  CED        = SUSE(Hu_sim, Hu_obs, Hs_sim, Hs_obs)
  
  CED
  
}

#%BiasTlag

Find_Max_CCF<- function(a,b)
{
  d <- ccf(a, b, plot = FALSE)
  cor = d[['acf']][,,1]
  lag = d[['lag']][,,1]
  res = data.frame(cor,lag)
  res_max = res[which.max(res[['cor']]),]
  return(res_max)
} 

## FitnessComputation: Computation of fitness

FitnessComputation = function(computedValues, originalValues, fitnessFunction = FitnessFunction, testing = FALSE, PforTlag = NA){
  
  if(all(is.na(computedValues))) return(NA)
  
  ## FOR COMPUTATION OF WEI, ORIGINAL Q DATA HAS ZEROES IN IT!!!
  if((fitnessFunction == 'MARE' | fitnessFunction == 'Vis_2'|fitnessFunction == 'rel_d0'|fitnessFunction == 'md0') & any(originalValues == 0)){
    
    originalValues = originalValues + quantile(originalValues[originalValues > 0], 0.1)
    
  }
  
  ## Setting of Global fitness values, dependent on testing
  if(testing) FitGlobs = FitnessGlobal(originalValues, assignIt = TRUE)
  
  ## Threshold of Q, for simulation of Q only!!!
  computedValues[computedValues < 0] = 0
  
  ## Constant result modification
  if(length(computedValues) == 1) computedValues = rep(computedValues, length(originalValues))
  
  computedValuesLow = computedValues[LowValPositions]
  computedValuesHigh = computedValues[HighValPositions]
  computedValuesMedium = computedValues[MediumValPositions]
  
  difference = originalValues - computedValues
  differenceLow = OriginalValuesLow - computedValuesLow
  differenceHigh = OriginalValuesHigh - computedValuesHigh
  differenceMedium = OriginalValuesMedium - computedValuesMedium
  
  if(fitnessFunction == 'logNS0' | fitnessFunction == 'Mai0' | fitnessFunction == 'Vis_1'|fitnessFunction == 'Shafii'){
    
    # logComputedValues = Qlog(computedValues)
    # logDifference = logOriginalValues - logComputedValues
    logComputedValues = Qlog(computedValues + quantile(originalValues[originalValues > 0], 0.1))
    #     logOriginalValues = Qlog(originalValues + quantile(originalValues[originalValues > 0], 0.1))
    logDifference = logOriginalValues - logComputedValues
    
    
  }
  
  if(fitnessFunction == 'FHV' | fitnessFunction == 'FMV' | fitnessFunction == 'FMS' | 
     fitnessFunction == 'FLV' | fitnessFunction == 'FDC'){
    
    fcp = FcriteriaPrecomputation(originalValues, computedValues, fitnessFunction)
    mapply(function(nm, val) assign(nm, val, envir = parent.env(environment())), names(fcp), fcp)
    
  }
  
  # P definition for Tlag computation
  if(!testing & fitnessFunction == 'Tlag') {
    
    if(!'P' %in% IndependentVariables) stop('Tlag fitness function needs a P variable (precipitation) in set of independent variables')
    P = DataSet[['P']]
    
  }
  if(testing & fitnessFunction == 'Tlag') P = PforTlag
  
  if(fitnessFunction == 'Shafii'){
    SumP=sum(DataSet[['P']]) #?????????
    #orrOriginalValues=sum(DataSet[['Q']])/SumP  #?????????????
    orrOriginalValues=sum(originalValues)/SumP
    orrlogOriginalValues=sum(logOriginalValues)/SumP  #??????????
    orrComputedValues= sum(computedValues)/SumP  #1
    orrlogComputedValues= sum(logComputedValues)/SumP  #2
    indexm1=which.min(abs(computedValues-(as.numeric(quantile(computedValues, probs = 0.2, na.rm = TRUE)))))
    computedm1=computedValues[indexm1]
    indexm2=which.min(abs(computedValues-(as.numeric(quantile(computedValues, probs = 0.7, na.rm = TRUE)))))
    computedm2=computedValues[indexm2]
    fdcmComputedValues=Qlog10(computedm2)-Qlog10(computedm1) #3
    fdchComputedValues=sum(computedValuesHigh) #4 
    computedl2=as.numeric(quantile(computedValues, probs = 0, na.rm = TRUE))
    fdclComputedValues=-1*(sum(Qlog10(computedValuesLow)-Qlog10(computedl2))) #5
    MeanComputedValues=mean(computedValues) #6
    VarComputedValues=sqrt(sum((computedValues-MeanComputedValues)^2)/(orlength-1))#7
    MedianComputedValues=median(computedValues) #8
    PeakComputedValues=max(computedValues) #9
    lag_num_Computed=numeric(orlength-1)
    for (i in 1:(orlength-1))
    {
      lag_num_Computed[i]=(computedValues[i]-MeanComputedValues)*(computedValues[i+1]-MeanComputedValues)
    }
    Lag_ac_Computed=sum(lag_num_Computed)/sum((computedValues-MeanComputedValues)^2) #10
    MeanlogComputedValues=mean(logComputedValues) #11
    VarlogComputedValues=sqrt(sum((logComputedValues-MeanlogComputedValues)^2)/(orlength-1)) #12
    Month_1=months(as.Date(DataSet[['Date']])) # ??????????????
    MaxMonthlyOriginal=max(aggregate( DataSet[['Q']] ~ Month_1 , DataSet , mean )[2]) #??????????????
    MaxMonthlyComputed=max(aggregate( computedValues ~ Month_1 , DataSet , mean )[2]) #13
    accept_tresh = .05 #acceptability treshold, choose [0.05,0.1 or 0.2] -> Shafii 2015
    sign_1   = 1-as.integer(as.logical((orrOriginalValues*(1-accept_tresh) <= orrComputedValues & orrOriginalValues*(1+accept_tresh) >= orrComputedValues)))
    sign_2   = 1-as.integer(as.logical((orrlogOriginalValues*(1-accept_tresh) <= orrlogComputedValues & orrlogOriginalValues*(1+accept_tresh) >= orrlogComputedValues)))
    sign_3   = 1-as.integer(as.logical((fdcmOriginalValues*(1-accept_tresh) <= fdcmComputedValues & fdcmOriginalValues*(1+accept_tresh) >= fdcmComputedValues)))
    sign_4   = 1-as.integer(as.logical((fdchOriginalValues*(1-accept_tresh) <= fdchComputedValues & fdchOriginalValues*(1+accept_tresh) >= fdchComputedValues)))
    sign_5   = 1-as.integer(as.logical((fdclOriginalValues*(1-accept_tresh) <= fdclComputedValues & fdclOriginalValues*(1+accept_tresh) >= fdclComputedValues)))
    sign_6   = 1-as.integer(as.logical((MeanOriginalValues*(1-accept_tresh) <= MeanComputedValues & MeanOriginalValues*(1+accept_tresh) >= MeanComputedValues)))
    sign_7   = 1-as.integer(as.logical((VarOriginalValues*(1-accept_tresh) <= VarComputedValues & VarOriginalValues*(1+accept_tresh) >= VarComputedValues)))
    sign_8   = 1-as.integer(as.logical((MedianOriginalValues*(1-accept_tresh) <= MedianComputedValues & MedianOriginalValues*(1+accept_tresh) >= MedianComputedValues)))
    sign_9   = 1-as.integer(as.logical((PeakOriginalValues*(1-accept_tresh) <= PeakComputedValues & PeakOriginalValues*(1+accept_tresh) >= PeakComputedValues)))
    sign_10  = 1-as.integer(as.logical((Lag_ac_Original*(1-accept_tresh) <= Lag_ac_Computed & Lag_ac_Original*(1+accept_tresh) >= Lag_ac_Computed)))
    sign_11  = 1-as.integer(as.logical((MeanlogOriginalValues*(1-accept_tresh) <= MeanlogComputedValues & MeanlogOriginalValues*(1+accept_tresh) >= MeanlogComputedValues)))
    sign_12  = 1-as.integer(as.logical((VarlogOriginalValues*(1-accept_tresh) <= VarlogComputedValues & VarlogOriginalValues*(1+accept_tresh) >= VarlogComputedValues)))
    sign_13  = 1-as.integer(as.logical((MaxMonthlyOriginal*(1-accept_tresh) <= MaxMonthlyComputed & MaxMonthlyOriginal*(1+accept_tresh) >= MaxMonthlyComputed)))
    
  }
  beta = mean(computedValues)/MeanOriginalValues
  #index_cor =  which(originalValues! = computedValues)
  sdo = sd(originalValues)
  sdc = sd(computedValues)
  if(is.na(sdo) | sdo == 0) sdo = .Machine[['double.xmin']] 
  if(is.na(sdc) | sdc == 0) sdc = .Machine[['double.xmin']] 
  alpha = sdc/sdo
  gamma = (sdc/mean(computedValues))/(sdo/MeanOriginalValues)
  
  r      = function()  cov(originalValues, computedValues) / (sdo * sdc)
  r0     = function()  1 - r() ## Wrong criteria - everything must go to zero. H(r) = -1 ,1 ??? !!!
  RMSE   = function()  sqrt(mean((difference)^2))
  NS0    = function()  sum(difference^2) / sum(NS0denominator^2) #14
  logNS0 = function()  sum(logDifference^2) / sum(logNS0denominator^2) #15
  RSq0   = function()  1 - r()^2
  PI0    = function()  sum((difference[2:length(originalValues)])^2) / sum((originalValues[2:length(originalValues)] - originalValues[1:(length(originalValues)-1)])^2)
  MAE    = function()  abs(mean(difference))
  MARE   = function()  ((1/length(originalValues))*(sum(abs(difference)/originalValues)))
  KG     = function(x) sqrt((r() - 1)^2 + (x - 1)^2 + (beta - 1)^2)
  VE0    = function()  sum(abs(difference)) / sum(originalValues)
  MNS0   = function() sum(abs(difference))/sum(abs(originalValues-MeanOriginalValues))
  RSD0   = function() 1-alpha
  CED0   = function() CEDcomputation(originalValues, computedValues)
  
  #yilmaz(2008): %Bias RR, FHV, FLV, FMM , FMS and Tlag
  fitnessDefinition = c(
    MARE       = 'fit = MARE()', #Mean Absolute Relative Error
    MAE        = 'fit = MAE()', # Mean Absolute Error
    RMSE       = 'fit = RMSE()', # Root Mean Square Error
    NS0        = 'fit = NS0()', #14
    logNS0     = 'fit = logNS0()', #15
    KG10       = 'fit = KG(alpha)', # KGE_2009
    KG20       = 'fit = KG(gamma)', # KGE_2012
    KGE1       = 'fit = 1-KG(alpha)', # KGE_2009
    KGE2       = 'fit = 1-KG(gamma)', # KGE_2012
    MNS0       = 'fit = MNS0()',# Modified NSE
    RSD0       = 'fit = RSD0()', # Ratio of Standard Deviation
    RSq0       = 'fit = RSq0()', # Coefficient of Determination
    PI0        = 'fit = PI0()', # Persistence Index
    VE0        = 'fit = VE0()', # Volumetric Efficiency
    r0         = 'fit = r0()', # Pearson correlation coefficient
    FHV        = 'fit = ((sum(Qsh-Qoh))/(sum(Qoh)))*100',  # remove abs : abs(((sum(Qsh-Qoh))/(sum(Qoh)))*100) %BiasFHV : Yilmaz
    FMV        = 'fit = ((sum(Qsm-Qom))/(sum(Qom)))*100', #remove abs %BiasFMV : Ley
    FMS        = 'fit = (((lQs_i1-lQs_i2)-(lQo_i1-lQo_i2))/(lQo_i1-lQo_i2))*100', #remove abs %BiasFMS : Yilmaz
    FLV        = 'fit = -1*((sum(lQsl-lQsL)-sum(lQol-lQoL))/(sum(lQol-lQoL)))*100', #remove abs %BiasFLV : Yilmaz
    FMM        = 'fit = ((log(median(computedValues))-log(median(originalValues)))/log(median(originalValues)))*100', #%BiasFMM : Yilmaz
    Tlag       = 'fit = ((Find_Max_CCF(computedValues,P)[["lag"]]-Find_Max_CCF(originalValues,P[["lag"]])/(Find_Max_CCF(originalValues,P)[["lag"]]))*100', #%BiasTlag
    RR         = 'fit = (sum(difference)/sum(originalValues))*100', #%BiasRR
    CED        = 'fit = CEDcomputation(originalValues, computedValues)',
    SUSE       = 'fit = CEDcomputation_all(originalValues, computedValues)',
    lowRMSE    = 'fit = sqrt(mean(differenceLow^2))',
    highRMSE   = 'fit = sqrt(mean(differenceHigh^2))',
    mediumRMSE = 'fit = sqrt(mean(differenceMedium^2))',
    pbias      = 'fit = (sum(difference)/sum(originalValues))*100', # Percentage Bias 
    md0        = 'fit = sum(abs(difference))/sum(abs(computedValues-MeanOriginalValues),(abs(originalValues-MeanOriginalValues)))', # modified index of agreement
    
    Borsanyi   = 'fit = sqrt((r0())^2+(NS0())^2+(VE0())^2+(KG(alpha))^2)',
    Mai0       = 'fit = sqrt(NS0()^2 + logNS0()^2)',
    Vis_1      = 'fit = sqrt(NS0()^2 + logNS0()^2 + VE0()^2)',  #Vis et al.,2015
    NSE        = 'fit = 1-NS0()',
    r          = 'fit = r()',
    Vis_3      = 'fit = sqrt((r0())^2+ (VE0())^2)', #Vis et al.,2015
    Price      = 'fit = sqrt(NS0()^2 + (MNS0())^2 + RSD0()^2)',
    CED_new    = 'fit = sqrt((CED0())^2 +(KG(alpha))^2)',
    rel_d0     = 'fit = sum((difference/originalValues)^2)/sum((sum(abs(computedValues-MeanOriginalValues),(abs(originalValues-MeanOriginalValues)))/MeanOriginalValues)^2)', #relative index of agreement
    Shafii     = 'fit = NS0()+logNS0()+sign_1+sign_2+sign_3+sign_4+sign_5+sign_6+sign_7+sign_8+sign_9+sign_10+sign_11+sign_12+sign_13' #Shafii 2015
    #AIC       = 'fit = (length(originalValues)*log(RMSE()))+(2*p)', #p: number of free parameters???
    #BIC       = 'fit = (length(originalValues)*log(RMSE()))+(p*log(length(originalValues)))',
    #C2M       = 'fit = NSE/(2-NSE)' #Mathevet et al.,2006,
  )
  
  if(fitnessFunction == 'Multi_Madsen'){
    
    fitOut = data.frame(
      lowRMSE    = eval(parse(text = fitnessDefinition['lowRMSE'])), 
      RMSE       = eval(parse(text = fitnessDefinition['RMSE'])),
      highRMSE   = eval(parse(text = fitnessDefinition['highRMSE'])), 
      MAE        = eval(parse(text = fitnessDefinition['MAE']))
    )
    
  } 
  
  else if (fitnessFunction == 'Dawson'){
    #     Dawson     = 'fit = sqrt((RMSE())^2+(RSq0())^2+(PI0())^2+(MAE())^2)', 
    fitOut = data.frame(
      RMSE       = eval(parse(text = fitnessDefinition['RMSE'])),
      RSq0       = eval(parse(text = fitnessDefinition['RSq0'])),
      PI0        = eval(parse(text = fitnessDefinition['PI0'])),
      MAE        = eval(parse(text = fitnessDefinition['MAE']))
    )
    
  } else if (fitnessFunction == 'Vis_2'){
    #     Vis_2      = 'fit = sqrt((r0())^2+(NS0())^2+ (VE0())^2 +(MARE())^2)', #Vis et al.,2015
    fitOut = data.frame(
      r0      = eval(parse(text = fitnessDefinition['r0'])),
      NS0     = eval(parse(text = fitnessDefinition['NS0'])),
      VE0     = eval(parse(text = fitnessDefinition['VE0'])),
      MARE    = eval(parse(text = fitnessDefinition['MARE']))
    )  
    
  } 
  
  else {
    
    fitOut = eval(parse(text = fitnessDefinition[fitnessFunction]))
    
  }
  
  if(fitnessFunction %in% MultiObjectiveFitnesses){
    
    fitOut = sum(scale(fitOut, center = FALSE, scale = FALSE)) ## THIS IS ONLY TEMPORARY, MUST BE PROBABLY DONE IN BETTER WAY - MULTIOBJECTIVE PROBLEM ??? !!!
    
  }
  
  fitOut
  
}

## Multi-objective criteria computation (Eucledian Distance also can be used!!!???)

MultiObjective = function(population, changeCheck = FALSE){
  
  ## BALANCED FITNESS  
  fitList = lapply(population, function(x) x[['Fitness']])
  fitDF = do.call('rbind', fitList)
  
  colMins = apply(fitDF, 2, min, na.rm = TRUE)
  
  A = max(colMins) - colMins
  
  balancedFitness = apply(fitDF, 1, function(x) sqrt(sum((x + A)^2)))
  
  balancedFitness 
  
}


## MakeEquation: Equation compilation from individual array

MakeEquation = function(individual){
  ReservoirSet <-c("WR","IR","RR","UR","FR","SR","CR")
  IA = individual[['IndArray']]
  
  # Vyreseni jedincu slozenych pouze z terminalu (nebo funkce s aritou 0 - zatim nejsou)
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
        
      } # j cycle end
      
      eq_vec[i] = paste0(eq_vec[i], ')')
      
    }
    
  } # i cycle end
  
  paste0('y = ', eq_vec[1])
  
}

## ComputeIndividual: Computation of individual
# individual - individual for computation
# compVar - variables for computation

ComputeIndividual = function(individual){
  
  if(PunishOneNodeIndividuals){
    
    if(individual[['IndLength']] == 1){
      
      individual[['Fitness']] = Inf
      individual[['Changed']] = FALSE
      return(individual)
      
    }
    
  }
  
  # Transfer of individual array to equation
  individual[['Equation']] = MakeEquation(individual)
 
  # Equation evaluation
  nrDataSet = nrow(DataSet)
  computedValues = suppressWarnings(try(eval(parse(text = individual[['Equation']]), envir = DataSet), silent = TRUE))
  
  # Identification of nonsenses
  if(class(computedValues) == 'try-error' | any(!is.finite(computedValues))){ 
    
    computedValues = NA 
    
  }
  
  # Extension of constant results to length of training set
  if(length(computedValues) == 1 & !is.na(computedValues[1])){
    
    computedValues = rep(computedValues, nrDataSet)
    
  }
  individual[['Sim']] = computedValues
  # Fitness computation
  originalValues = DataSet[[DependentVariable]]
  for (i in 1:length(FitnessFunction)){
  individual[['Fitness']][i] = FitnessComputation(computedValues, originalValues,fitnessFunction = FitnessFunction[i]) 
  }
  for (i in 1:length(FitnessFunction)){
  if (is.na(individual[['Fitness']][i])){
    individual[['Fitness']]<- rep(10000,length(FitnessFunction))
  }
  if (is.infinite(individual[['Fitness']][i])){
    individual[['Fitness']] <- rep(10000,length(FitnessFunction))
  }
  }
  
  individual[['Changed']] = FALSE
  
  individual
  
}
