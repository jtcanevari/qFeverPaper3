#-----------------------------------------------------------
# run Dutch version Rcpp model. 
#-----------------------------------------------------------

library("Rcpp") 

#-----------------------------------------------------------
# source the model
Rcpp::sourceCpp("C:\\Users\\jcanevari\\Dropbox\\qFeverPaper3\\code\\DutchVersion.cpp")
Rcpp::sourceCpp("C:\\Temp\\DutchVersion.cpp")

#-----------------------------------------------------------
# select starting conditions and parameters

xstart <- c(SNP = 425, SP = 424, INP1 = 0, INP2 = 0, INP3 = 0, INP4 = 0, IP = 1, IP2 = 0,
            JNP = 0, JP = 0, RNP = 0, RP = 0, Y = 150, A = 0, K = 0, KI = 0, E = 0)

# params <- c(muE = 1/20, beta = 1, gamma = 4/14.4, epsilon_p = 1, epsilon_f = 10^(-6)*1000/365, 
            # phi = 0.95, mu = (1+(22/30)*0.95)/(3.1*365), alpha =0.7, fI = 0.75, fJ = 0.25, p = 0.5)

params <- c(phi = 0.95, mu = (1+(22/30)*0.7)/(3.1*365), fI = 0.75, fJ = 0.25, gamma = 4/14.4,
            alpha = 0.2, beta = 1.8, p = 0.5, epsilon.p = 1, epsilon.f = 10^(-6)*1000/365, muE = 1/20)
#-----------------------------------------------------------
# run and store
nsims = 1
system.time(
modelOut <- lapply(1:nsims, function(i) CaneGillespie(103,365*15, xstart, params))
)
#-----------------------------------------------------------