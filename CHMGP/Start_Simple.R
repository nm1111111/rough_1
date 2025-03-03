rm(list = ls())

setwd('/home/vojta/Dropbox/PRACE/Jayashree/CHMGPdevelopement/DATA')
nDat = 20 #50 is too less to test CED
data = read.table('Wollefsbach240h.dat',skip =7, col.names = c('year','month','day','hour','min','QC','P','E','Q','T'))
dates = apply(data[,1:3], 1, function(x) paste0(x, collapse = '-'))
testing = data.frame(Date = as.Date(dates), P = data$P, E = data$E, Q = data$Q)[1:nDat,]

setwd('/home/vojta/Dropbox/PRACE/Jayashree/CHMGPdevelopement/CHMGP_strict_SF')
source('CHMGP.R', local = TRUE)

SFPR = SetSFPR()

# testing[["Qtest"]] = MI(0.5, 0.5, 0.5, 0.5)


## TESTING TANKS
# testing = data.frame(P = 1:5, E = (1:5)/10)

testing[["RI"]] = testing[["P"]] - testing[["E"]]
DataSet = testing

## OK
testing[['Qtest']] = RJOIN(99, 99, RB(testing[["RI"]], 0.1, 0.1, 0.9, 0.9, 0.1), 99)
## OK
# testing[['Qtest']] = RJOIN(RB(RB(testing[["RI"]], 0, 0.7, 0.9, 0.1, 0.6),0.2, 0.3, 0.1, 0.6, 0), 99)
## OK
# testing[['Qtest']] = RJOIN(RB(testing[["RI"]],0.2, 0.3, 0.1, 0.6, 0), RB(testing[["RI"]], 0, 0.7, 0.9, 0.1, 0.6))
## OK
# testing[['Qtest']] = RJOIN(RB(RB(testing[["RI"]], 0, 0.7, 0.9, 0.1, 0.6),0.2, 0.3, 0.1, 0.6, 0), RB(testing[["RI"]], 0, 0.7, 0.9, 0.1, 0.6))
## OK
# testing[['Qtest']] = RJOIN(RB(testing[["RI"]], 0, 0.7, 0.9, 0.1, 0.6), RB(RB(testing[["RI"]], 0, 0.7, 0.9, 0.1, 0.6),0.2, 0.3, 0.1, 0.6, 0))

DataSet = testing

# str(DataSet)
# plot(testing[['Qtest']], type = 'l')


## OK
# TANK(0.1, 0.1, 0.1, 0.1, 0.1)
# RJOIN(99, 99, 99, RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1))

## OK
# TANK(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1) # This TANK works with infiltration for second reservoir, not a case of RJOIN
# RJOIN(testing[["P"]], 99, RB( RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1),  0.1, 0.1, 0.1, 0.1, 0.1), 99)

## OK, ONLY THIS TEST CAN BE USED FOR RJOIN WITH DEPTH IDENTIFICATION
# TANK(0.1, 0.1, 0.1, 0.1, 0.1,  0.1, 0.1, 0.1, 0.1, 0.1,  0.1, 0.1, 0.1, 0.1, 0.1) # This TANK works with infiltration for second reservoir, not a case of RJOIN
# RJOIN(99, 99, RB( RB( RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1),  0.1, 0.1, 0.1, 0.1, 0.1),  0.1, 0.1, 0.1, 0.1, 0.1), 99)

##Parallel tanks implementation
## OK
# RJOIN(RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1), 99 , 99, RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1))
# RJOIN(RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1), RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1), RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1), RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1))


## OK 
# RJOIN(99, 99, RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1), RB(RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1),0.1, 0.1, 0.1, 0.1, 0.1))
# RJOIN(RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1), RB(RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1),0.1, 0.1, 0.1, 0.1, 0.1), RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1), RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1))

## OK 
# RJOIN(99, RB(RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1),0.1, 0.1, 0.1, 0.1, 0.1), 99, RB(testing[["RI"]],0.1, 0.1, 0.1, 0.1, 0.1))

# options(warn = 2, error=recover) # warn = 2 means all warnings are turned into errors
# options(error=recover)
# options(error=NULL)

# setwd('/home/vojta/Dropbox/PRACE/Jayashree/CHMGPdevelopement/CHMGP_strict_SF')
# source('CHMGP.R', local = TRUE)

useModels = c('RJOIN', 'RB')

GPresult = CHMGP(testing,

                  DependentVariable = 'Qtest',
                  IndependentVariables = c('RI'), 
                  FunctionSet = c(useModels), # , '+', '-', '*', '/'
                  ConstantRange = c(0, 1), 
                  StrictFunctionsArgumentsRange = c(lower = 0, upper = 1), 
                  PopulationSize = 500, 
                  NumberOfGenerations = 100, 
                  FitnessFunction = 'NS0',
                  MaxDepthIni = 2, 
                  MaxDepthRun = 3, 
                  TournamentSize = 4, 
                  RoundingFactor = 3,
                  VariationProbs = c(Pc = 0.7, PrSelected = 0.05, PmTree = 0.5, PmSeparation = 0.3, PmNode = 0.3, PmConstant = 0.7),
                  SimpleOutput = FALSE,
                  YacasSimplification = TRUE,                  
                  DEoptimization = FALSE,
                  DEiterMax = 100,
                  DEfitness = 'SameAsGPfit', #'SameAsGPfit'
                  CTtypesOccurence = FALSE
                  
)

GPresult



## Testing Simplification
# options(error=recover)
# 
# setwd('/home/vojta/Dropbox/PRACE/Jayashree/CHMGPdevelopement/CHMGP_strict_SF')
# source('CHMGP.R', local = TRUE)
# 
# SimplifyCT("y = CT(0.982,0.061,0,0.084,0,1,0.037,0.513,0.136,0.061,0.419,0.038,0.476,0.686,0.038,0.038,0.061,0.583,0.333,1,0,0.455,0.072,0,0.038,0.829,0.806,0.011,0.24,1.000,1.000,1.000,1.000,1.000)")

 
# setwd('/home/vojta/Plocha/pokus')
# CTtypesPlot(GPresult[['CTtypes']])

# SimplifyCT("y = CT(0.858,0.965,0.81,0.679,0.397,0.278,0.443,0.412,0.44,0.741,0.43,0.293,0.329,0.246,0.517,0.105,0.867,0.155,0.504,0.951,0.934,0.883,0.448,0.068,0.289,0.92,0.687,0.363,0.184,0.000,1.000,1.000,0.000,0.000)")
