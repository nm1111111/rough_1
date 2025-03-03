## Heaviside function - unit step

Hstep = function(x){
  
  as.numeric(x >= 0)
  
}


## Gr: greater than

Gr = function(x1, x2){
  
  as.numeric(x1 > x2)
  
}


## GrEq: greater or equal

GrEq = function(x1, x2){
  
  as.numeric(x1 >= x2)
  
}


## Eq: equal

Eq = function(x1, x2){
  
  as.numeric(x1 == x2)
  
}

Pdiv = function(x1, x2){
  
  if(any(x2 == 0)){
    
    x2[which(x2 == 0)] = 1
    
  } 
  
  x1 / x2 
  
}

Psqrt = function(x){
  
  sqrt(abs(x))
  
}


## DLY: delay function
# x - vector or scalar for shift
# shift - vector or scalar of window length (zero is actual-present value)

DLY = function(x, shift){
  
  if(any(is.na(c(x, shift)))) return(NA)
  
  n = length(x)
  m = length(shift) 
  output = numeric(n)
  
  if(n == 1 | (m == 1 && shift == 0)) return(x)
  
  fOutput = .Fortran("SlidingWindow", x = as.double(x), n = as.integer(n), shift = as.double(shift), m = as.integer(m),
                     output = as.double(output), fun = as.character('DLY'))
  
  fOutput[['output']]
  
}


## SSUM: Shifted sum, sum of sliding window
# x - vector or scalar for shift
# shift - vector or scalar of window length (zero is actual-present value) 

SSUM = function(x, shift){
  
  if(any(is.na(c(x, shift)))) return(NA)
  
  n = length(x)
  m = length(shift) 
  output = numeric(n)
  
  if(n == 1 | (m == 1 && shift == 0)) return(x)
  
  fOutput = .Fortran("SlidingWindow", x = as.double(x), n = as.integer(n), shift = as.double(shift), m = as.integer(m),
                     output = as.double(output), fun = as.character('SUM'))
  
  fOutput[['output']]
  
}


## SMA: Simple Moving Average
# x - vector or scalar for shift
# shift - vector or scalar of window length (zero is actual-present value) 

SMA = function(x, shift){
  
  if(any(is.na(c(x, shift)))) return(NA)
  
  n = length(x)
  m = length(shift) 
  output = numeric(n)
  
  if(n == 1 | (m == 1 && shift == 0)) return(x)
  
  
  fOutput = .Fortran("SlidingWindow", x = as.double(x), n = as.integer(n), shift = as.double(shift), m = as.integer(m),
                     output = as.double(output), fun = as.character('SMA'))
  
  fOutput[['output']]
  
}

## RES:Reservoir function, function of simple reservoir (linear in storage response)
# x - input data for reservoir
# K - parametr nadrze
RES = function(x, K){
  
  if(any(is.na(c(x, K)))) return(NA)
  
  #   K[K == 0] = 1e-9
  K[K < 1] = 1
  
  n = length(x)
  m = length(K)
  output = numeric(n)
  
  if(n == 1) return(x)
  
  K = abs(K)
  
  if(m == 1) K = rep(K, n)
  
  fOutput = .Fortran("Reservoir", x = as.double(x), n = as.integer(n), a = as.double(K), output = as.double(output))
  
  fOutput[['output']]
  
}
