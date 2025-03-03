DDS <- function(x,n) {
  r <- 0.2
  n <- 10000
  sBest <- x
  sCur <- x
  NS_In <- NSO(x)
  if (is.na(NS_In)){
    NS_In <- 100000
  } else {
    NS_In <- NS_In
  }
  CostBest <- NS_In
  dimen <- length(x)
  x_min <- 0
  x_max <- 1
  x_range <- x_max -x_min
  k <- 0
  for(i in 1:n) {
    for(j in 1:dimen) {
      if (runif(1) <(1-(log10(i)/log10(n)))) {
      k <- k +1
      sCur[j] <- sCur[j] + rnorm(1)*r*x_range
        if(sCur[j]<x_min){
        sCur[j] <- x_min + (x_min - sCur[j])
          if(sCur[j]>x_max){
          sCur[j] <-x_min
          }
        }
       else if(sCur[j]> x_max){
       sCur[j] <- x_max -(sCur[j]-x_max)
          if(sCur[j]< x_min){
          sCur[j] <- x_max
          }
        }
      }
    }
    if(k==0){
      index = round(runif(1,1,dimen))
      sCur[index] <- sCur[index] + rnorm(1)*r*x_range
      if(sCur[index]<x_min){
        sCur[index] <- x_min + (x_min - sCur[index])
        if(sCur[index]>x_max){
          sCur[index] <-x_min
        }
      }
      else if(sCur[index]> x_max){
        sCur[index] <- x_max -(sCur[index]-x_max)
        if(sCur[index]< x_min){
          sCur[index] <- x_max
        }
      }
    }
    k <- 0
    NS <- NSO(sCur)
    if (is.na(NS)){
      NS <- 100000
    } else {
      NS <- NS
    }
    if(NS < CostBest){
      sBest <- sCur
      CostBest <- NS
    } else {
      sBest <- sBest
      CostBest <- CostBest
    }
  }
  return(sBest)
}

NSO <- function(x){
  sim <- CT_new(WR(x[1],x[2],x[3],x[4],x[5]),IR(x[6],x[7],x[8],x[9]),RR(x[10],x[11],x[12]),UR(x[13],x[14],x[15],x[16],x[17],x[18],x[19]),FR(x[20],x[21],x[22],x[23],x[24],x[25],x[26]),SR(x[27],x[28],x[29],x[30],x[31]),CR(x[32],x[33],x[34],x[35],x[36],x[37],x[38],x[39],x[40],x[41],x[42],x[43],x[44]),x[45],x[46],x[47],x[48],x[49],x[50],x[51],x[52],x[53],x[54],x[55],x[56])
  obs <- DataSet[,5]
  NSO_value <- FitnessComputation(sim,obs,"NS0")
  if (is.na(NSO_value)){
    NSO_value <- 100000
  } else {
    NSO_value <- NSO_value
  }
  return(NSO_value)
}

