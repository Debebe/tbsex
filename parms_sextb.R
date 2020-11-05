
#===global parameters====
globalPars <- list(beta = 7.5,react_mid = 0.0031, fast_mid= 0.135,stab=1.34, RRprog= 2, RRcdr= 0.52, 
                   psi = 0.79,relapse_mid =0.019,rho=0.50,tau=1/0.5,theta=0.73, pm=0.43,pf=0.38)  
#== intermediate parameters====
#localPars <-list(CDR_mid =0.61,life_expect= 75.2,infect_dur= 3)#V
localPars <-list(CDR_mid =0.50,life_expect= 61.373,infect_dur= 3)#Uganda


#==wrap up params used by odin
genPars<-function(globalPars,localPars){
  list2env(c(globalPars, localPars),envir = environment())
  mu    <- 1/life_expect
  recov <- 0.5/infect_dur
  mutb  <- 0.5/infect_dur
  list(mu=mu,mutb=mutb, recov=recov,cdr_mid=CDR_mid) 
}

pars <-genPars(globalPars, localPars)  
pms<- c(globalPars, list(mu=pars$mu, mutb=pars$mutb, 
                         recov=pars$recov, cdr_mid=pars$cdr_mid))
