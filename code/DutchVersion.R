#este es el modelo como lo publico Bontje. las ff son ligeramente distintas

library(glogis); library(parallel)

#------------------------------------------------------------#
# modulus function (turns day 366 to calendar day 1)
#------------------------------------------------------------#
mod <- function(a,b,c)
{ R <- a%%b
limit <- function (x){ifelse(x >= c, x, while(x < c) {x <- x + b; return(x)})}
result <- unlist(lapply(R, function(x) limit(x)))
return(result)
}

#------------------------------------------------------------#
#forcing functions
#------------------------------------------------------------#

tp <- 253 #day of parturition
Tg <- 150
tc <- tp - Tg #day conception
Tl <- tp- 4*7 #from here on if infected move to IP2
# ta <- tp - 25 # average day of abortion
ta <- 50 # day abortions start

#create hazard functions for a logistic dist fitted to Autumn 2015 kiddings.
HF <- function(x) {
  exp(dglogis(x, location = tp, scale = 4.987, shape = 1.125, log = TRUE) - pglogis(x, location = tp, scale = 4.987, shape = 1.125, log = TRUE, lower=FALSE))
}

HF2 <- function(x) {
  exp(dglogis(x, location = tc, scale = 4.987, shape = 1.125, log = TRUE) - pglogis(x, location = tc, scale = 4.987, shape = 1.125, log = TRUE, lower=FALSE))
}

# HF3 <- function(x) {
#   exp(dglogis(x, location = ta, scale = 6.645, shape = 1.068, log = TRUE) - pglogis(x, location = ta, scale = 6.645, shape = 1.068, log = TRUE, lower=FALSE))
# }

Qp <- function(t) {ifelse(mod(t,365,1) > (tp+40) | mod(t,365,1) < (tp-30) ,0,HF(mod(t,365,1)))}
Qc <- function(t) {ifelse(mod(t,365,1) > (tc+40) | mod(t,365,1) < (tc-30),0,HF2(mod(t,365,1)))}
Qin <- function(t) {ifelse(mod(t,365,1)<tp+60|mod(t,365,1)>tp+90,0,0.214)}
Q4 <- function(t) {ifelse(mod(t,365,1)>Tl & mod(t,365,1)<tp+40,0,1)}
Qa <- function(t) {ifelse(mod(t,365,1) < (ta-(11.44*3)) | mod(t,365,1) > (ta+(11.44*3)),0,HF3(mod(t,365,1)))} #sd=11.44
#------------------------------------------------------------#
# Gillespie
#------------------------------------------------------------#

# First we define our one-step function.
SIR.onestep <- function (x, params) { #function to calculate one step of stochastic SIR
  t    <- x[1] #local variable for time
  SNP  <- x[2] #local variable for susceptibles
  SP   <- x[3]
  INP1 <- x[4] #local variable for infecteds
  INP2 <- x[5] #local variable for infecteds
  INP3 <- x[6] #local variable for infecteds
  INP4 <- x[7] #local variable for infecteds
  IP   <- x[8]
  IP2  <- x[9]
  JNP  <- x[10]
  JP   <- x[11]
  RNP  <- x[12] #local variable for recovereds
  RP   <- x[13]
  Y    <- x[14]
  A    <- x[15]
  K    <- x[16] #local variable for environmental contamination
  KI   <- x[17] #kidding when infected. not a compartment, just to keep track
  E    <- x[18] #environmnet
  N <- SNP+SP+INP1+INP2+INP3+INP4+IP+IP2+JNP+JP+RNP+RP #total population size (subject to demographic change)
  
  with( #use with as in deterministic model to simplify code
    as.list(params),
    {
      rates <- c(Y*Qin(t),                      #births
                 (phi*(SNP+SP)-SP)*Qc(t),       #conceptions
                 (phi*(INP1+IP/4)-IP/4)*Qc(t),
                 (phi*(INP2+IP/4)-IP/4)*Qc(t),
                 (phi*(INP3+IP/4)-IP/4)*Qc(t),
                 (phi*(INP4+IP/4)-IP/4)*Qc(t),
                 (phi*(RNP+RP)-RP)*Qc(t),
                 (JNP - (1-phi)*(JP + JNP))*Qc(t), #DEJA JNP Y VA A JP
                 SP*Qp(t),                         #kiddings
                 (1-alpha)*IP*Qp(t),
                 alpha*IP*Qp(t),
                 (1-alpha)*IP2*Qp(t),
                 alpha*IP2*Qp(t),
                 JP*Qp(t),                       #DEJA JP Y VA A RNP
                 RP*Qp(t),
                 mu*SNP,                           #deaths
                 mu*INP1,
                 mu*INP2,
                 mu*INP3,
                 mu*INP4,
                 mu*JNP,                           #DEJA JNP Y VA A Y
                 mu*RNP,
                 (beta*SNP*E)/N,                   #transmission terms
                 Q4(t)*beta*SP*E/N,
                 beta*(1-Q4(t))*SP*E/N,
                 gamma*INP1,
                 gamma*INP2,
                 gamma*INP3,
                 (1-p)*gamma*INP4,                #recoveries
                 p*gamma*INP4,
                 alpha*fI*IP*Qa(t),   #abortions
                 (1-alpha)*fI*IP*Qa(t),
                 fJ*JP*Qa(t)         #DEJA JP Y VA A RNP
      )
      rates <- ifelse(rates > 0, rates, exp(-10))     #otherwise negative rates sometimes
      
      changes <- matrix(
        c(1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0, # 1 Y*Qin(t)
          -1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # 2 (phi*(SNP+SP)-SP)*Qc(t)
          0,0,-1,0,0,0,1,0,0,0,0,0,0,0,0,0, # 3 (phi*(INP1+IP/4)-IP/4)*Qc(t)
          0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0,0, # 4 (phi*(INP2+IP/4)-IP/4)*Qc(t)
          0,0,0,0,-1,0,1,0,0,0,0,0,0,0,0,0, # 5 (phi*(INP3+IP/4)-IP/4)*Qc(t)
          0,0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0, # 6 (phi*(INP4+IP/4)-IP/4)*Qc(t)
          0,0,0,0,0,0,0,0,0,0,-1,1,0,0,0,0, # 7 (phi*(RNP+RP)-RP)*Qc(t)
          0,0,0,0,0,0,0,0,-1,1,0,0,0,0,0,0, # 8 (JNP-(1-phi)*(JP+JNP))*Qc(t)
          1,-1,0,0,0,0,0,0,0,0,0,0,0,0,1,0, # 9 SP*Qp.A(t)
          0,0,0,0,0,0,-1,0,0,0,1,0,0,0,1,1, # 10 (1-alpha)*IP*Qp(t)
          0,0,0,0,0,0,-1,0,1,0,0,0,0,0,1,1, # 11 alpha*IP*Qp(t)
          0,0,0,0,0,0,0,-1,0,0,1,0,0,0,1,1, # 12 (1-alpha)*IP2*Qp(t)
          0,0,0,0,0,0,0,-1,1,0,0,0,0,0,1,1, # 13 alpha*IP2*Qp(t)
          0,0,0,0,0,0,0,0,0,-1,1,0,0,0,1,1, # 14 JP*Qp(t)  
          0,0,0,0,0,0,0,0,0,0,1,-1,0,0,1,0, # 15 RP*Qp(t)
          -1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0, # 16 mu*SNP
          0,0,-1,0,0,0,0,0,0,0,0,0,1,0,0,0, # 17 mu*INP1
          0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0,0, # 18 mu*INP2
          0,0,0,0,-1,0,0,0,0,0,0,0,1,0,0,0, # 19 mu*INP3
          0,0,0,0,0,-1,0,0,0,0,0,0,1,0,0,0, # 20 mu*INP4
          0,0,0,0,0,0,0,0,-1,0,0,0,1,0,0,0, # 21 mu*JNP
          0,0,0,0,0,0,0,0,0,0,-1,0,1,0,0,0, # 22 mu*RNP
          -1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0, # 23 (beta*SNP*E)/N
          0,-1,0,0,0,0,1,0,0,0,0,0,0,0,0,0, # 24 Q4(t)*beta*SP*E/N
          0,-1,0,0,0,0,0,1,0,0,0,0,0,0,0,0, # 25 beta*(1-Q4(t))*SP*E/N
          0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0,0, # 26 gamma*INP1
          0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0,0, # 27 gamma*INP2
          0,0,0,0,-1,1,0,0,0,0,0,0,0,0,0,0, # 28 gamma*INP3
          1,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0, # 29 (1-p)*gamma*INP4
          0,0,0,0,0,-1,0,0,0,0,1,0,0,0,0,0, # 30 p*gamma*INP4
          0,0,0,0,0,0,-1,0,1,0,0,0,0,1,0,0, # 31 alpha*fI*IP*Qa(t)
          0,0,0,0,0,0,-1,0,0,0,1,0,0,1,0,0, # 32 (1-alpha)*fI*IP*Qa(t)
          0,0,0,0,0,0,0,0,0,-1,1,0,0,1,0,0  # 33 fJ*JP*Qa(t)
        ),
        ncol=16, byrow=TRUE)
      
      tau <- rexp(n=1,rate=sum(rates)) # exponential waiting time
      
      if(tau <= 14){
        U <- runif(1) #uniform random deviate 
        m <- min(which(cumsum(rates)>=U*sum(rates)))
        tx <- x[2:17] + changes[m,]
        
        if(all(tx >= 0)){
          x <- x[2:17] <-tx
          if (E>0) {E = E - E*muE*tau} #when E >0 it'll dacay at specified rate
          if (E<0) {E <- 0} #but can never be <0
          if (m %in% c(10:14,31:33)) {E <- E + epsilon.p}
        }
        
        if(any(tx < 0)){
        tau = 0.2
        x <- x[2:17]
        if (E>0) {E = E - E*muE*tau} #when E >0 it'll dacay at specified rate
        if (E<0) {E <- 0} #but can never be <0
        }
      }
      
      #if rexp(sum(rates)) is a long step force not to be longer than 14d  
      if(tau > 14){
        tau = 0.2
        x <- x[2:17]
        if (E>0) {E = E - E*muE*tau} #when E >0 it'll dacay at specified rate
        if (E<0) {E <- 0} #but can never be <0  
      }
      
      E <- E + (INP1+INP2+INP3+INP4+IP+IP2+JNP+JP)*epsilon.f*tau
      t <- t + tau
    
      return(out <- c(t, x, E))
    }
  )
}

# Next we write our function for simulating a whole epidemic.
ModelA <- function (x, params) {
  output <- list() 
  output[[1]] <- x #first record of output is initial condition
  k=1
  stop<- params$stop
  while (as.logical(output[[k]][1]<stop)){
    x <- SIR.onestep(x,params)
    output[[k+1]] <- x
    k=k+1
    if(sum(x[c(4:11,18)]) < exp(-10)){break}
  }
  out <- do.call(rbind,output) #return output
  colnames(out) <- c("time","SNP","SP","INP1","INP2","INP3","INP4","IP","IP2","JNP","JP","RNP","RP","Y","A","K","KI","E") 
  out
}

#------------------------------------------------------------#
# params
#------------------------------------------------------------#
params <- data.frame(
  #demographic parameters
  phi = 0.95, # fraction of animals to conceive
  mu = (1+(22/30)*0.7)/(3.1*365), # daily death rate. Median age = 3.1, phi=0.7,Lav=3.1*365
  fI = 0.75, # fraction of infected animals in IP(infected preg) to abort
  fJ = 0.25, # fraction of infected animals in JP(persistently inf preg) to abort (Berri 2007)
  gamma = 4/14.4, #Infectious period of infected non preg (infperiod = 14.4, ninf = 4)
  alpha = 0.2, #fraction persistently infected (we have a lower number ~ 0.2)
  beta = 1.8, #Transmission rate per equivalent Coxiella excretion (day ^-1 scaled unit^-1)
  p = 0.5, #prob INP individual of clearing infection and gain long-term immunity, no data
  #shedding parameters
  epsilon.p = 1,     # C. burnetii excretion in partus material (per parturition)
  epsilon.f = 10^(-6)*1000/365, # C. burnetii excretion in feces and urine per day
  #C.b. survival parameters
  muE = 1/20, #mortality rate for Cb in manure. Bontje et al 2016.
  # Additional parameters for the stochastic simulations
  start = 103, #avergae day of conception
  stop = 365*15
)
# 
# #------------------------------------------------------------#
# # run model
# #------------------------------------------------------------#

#initial conditions
xstart <- c(time=params$start,SNP=425, SP=424,INP1=0,INP2=0,INP3=0,INP4=0,IP=1,IP2=0,JNP=0,JP=0,RNP=0,RP=0,Y=150,A=0,K=0,KI=0,E=0) 

# out <- as.data.frame(ModelA(xstart,params))

#run multiple iterations
# how many sims?
nsims=1

#run in parallel in pc cores
# Initiate cluster
cl <- makeCluster(3)
#need to export the global environment to the clusters
clusterExport(cl, names(.GlobalEnv))
clusterEvalQ(cl, c(library(glogis)))

#run model (takes about 2 mins to run 3 iterations. much faster in C++!
modelOut <- parLapply(cl,1:nsims, function(i) ModelA(xstart,params))

stopCluster(cl)
