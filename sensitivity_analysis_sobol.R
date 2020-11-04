

library(odin)
library(data.table)
library(ggplot2)
library(reticulate)


SexTB <- odin::odin({
  nsex <- user(2)
  n    <- nsex
  #=====create beta matrix=====
  
  comb_matrix[1,1] <- (1-(1-rho)/2)
  comb_matrix[2,2] <- (1-(1-rho)/2)
  comb_matrix[1,2] <- (1-rho)/2
  comb_matrix[2,1] <- (1-rho)/2
  
  Beta[,]              <- beta*comb_matrix[i,j]
  
  dim(comb_matrix)     <- c(n,n)
  
  
  lambda_prod[ , ] <- Beta[i, j] * prop_infec[j]*I[j]/N[j] 
  lambda[]         <- sum(lambda_prod[i, ])             
  output(lambda) <- TRUE                                
  output(comb_matrix) <- TRUE                                 
  
  #=====derivatives=====
  N[]  <- S[i] + E[i] + L[i] + I[i]+ T[i] + R[i]
  totN <- sum(N[])
  output(N[])  <- TRUE
  output(totN) <- TRUE
  
  contacts <- user()
  rho = 2*contacts-1
  # differential progression
  react_mid     <-user() 
  fast_mid      <-user()
  cdr_mid       <-user()
  relapse_mid   <-user()
  
  react[1]   <-2*RRprog*react_mid/(1+RRprog)
  react[2]   <-2*react_mid/(1+RRprog)
  fast[1]    <-2*RRprog*fast_mid/(1+RRprog)
  fast[2]    <-2*fast_mid/(1+RRprog)
  
  relapse[1]    <-2*RRprog*relapse_mid/(1+RRprog)
  relapse[2]    <-2*relapse_mid/(1+RRprog)
  
  #cdr
  RRcdr    <- user()
  cdr_mid_rate <-cdr_mid*(recov + mutb +mu)/(1-cdr_mid);
  
  cdr[1]   <- 2*cdr_mid_rate*RRcdr/(1+RRcdr)
  cdr[2]   <- 2*cdr_mid_rate/(1+RRcdr);
  
  dim(react) <-n
  dim(fast)  <-n
  dim(relapse)  <-n
  RRprog     <- user()
  
  dim(cdr) <- n
  prop_infec[1]<- pm 
  prop_infec[2]<- pf
  pm <-user()
  pf <-user()
  #====deriv=====
  deriv(S[]) <- I[i] * mutb  + N[i] * mu + T[i]*(1-theta)*tau- S[i] * (lambda[i] + mu) 
  deriv(E[]) <- S[i] * lambda[i] + L[i]*(lambda[i]*(1-psi))- (mu+stab +fast[i])*E[i]
  deriv(L[]) <- E[i] * stab + I[i] * recov   - (mu + react[i] +lambda[i]*(1-psi))*L[i]
  deriv(I[]) <- E[i] * fast[i] + L[i]*react[i] + R[i]*relapse[i]-I[i]*(mu+mutb+recov + cdr[i])
  deriv(T[]) <- I[i] * cdr[i] - T[i] * (tau + mu) 
  deriv(R[]) <- theta * (tau)*T[i] - R[i] * ( relapse[i] + mu) 
  
  initial(S[]) <- 1-1e-6
  initial(I[]) <- 1e-6
  initial(L[]) <- 0
  initial(E[]) <- 0
  initial(T[]) <- 0
  initial(R[]) <- 0
  
  
  prev[]           <- 1e5*I[i]/N[i]
  output(prev[])   <- TRUE
  output(prev_mid) <-  sum(prev[])/2
  output(MFratio)  <-  prev[1]/prev[2]
  
  beta        <- user()
  recov       <- user()             
  mu          <- user()              
  mutb        <- user()             
  stab        <- user()
  psi         <- user()             
  tau         <- user()          # inverse of treatment duration
  theta       <- user()          # proportion of treated cured
  
  # dimensions
  dim(Beta)         <- c(n, n)
  dim(lambda_prod)  <- c(n, n)
  dim(lambda) <- n
  dim(prop_infec)<-n
  dim(S)      <- n
  dim(E)      <- n
  dim(L)      <- n
  dim(I)      <- n
  dim(T)      <- n
  dim(R)      <- n
  dim(N)      <- n
  dim(prev)   <- n
},target="c")


use_python("/Library/Frameworks/Python.framework/Versions/3.8/bin/python3")
reticulate::py_config()
reticulate::py_available()

#---------------------------------------
# import some functions from SALib
SA       = import('SALib')
saltelli = import("SALib.sample.saltelli")
sobol    = import('SALib.analyze.sobol')
#-----------------------------------------

# Load priors
source('~/Dropbox/SexTB/Clean/new/priorsRevised_new.R')  

priors$rrcdr_a=-0.298
priors$rrcdr_b=0.2026 
priors$rrprog_a <- -0.298
priors$rrprog_b <- 4*0.2026 
priors$rho_a=-0.57
priors$rho_b=0.085
priors$stab_s <- 0.34/5



# Generate samples
problem = list(
  num_vars = 12L,
  names = list('beta', 'contacts','reactivation', 'fast-progression', 
               'stabilisation', 'CDR','progression-RR', 'CDR-RR', 
               'protection', 'relapse', 'TB-CFR', 'recovery'),
  bounds=list(list(0, 1),list(0, 1),list(0, 1),
              list(0, 1),list(0, 1),list(0, 1),
              list(0, 1),list(0, 1),list(0, 1),
              list(0, 1),list(0, 1),list(0, 1)))

params = as.data.frame(saltelli$sample(problem, 3000L))
names(params) = problem$names

# ---timer------
start_time <- Sys.time()
#------------
t <- seq(0, 5000, length.out = 5000) 
k <- J <- 0
outputs <- outprev <- outmfr <- preve <- mfre <- list()
for(i in 1:nrow(params)) {
      k <- k+1
  outputs[[k]] <- data.frame(SexTB(
    beta      = qlnorm(params[i,1], priors$beta_m, priors$beta_s),
    contacts  = qlnorm(params[i,2], priors$rho_a, priors$rho_b),
    react_mid = qlnorm(params[i,3], priors$react_m, priors$react_s),
    fast_mid  = qlnorm(params[i,4], priors$fast_m, priors$fast_s),
    stab      = qlnorm(params[i,5], priors$stab_m, priors$stab_s),
    cdr_mid   = qbeta(params[i,6],  priors$cdr_m, priors$cdr_s),
    RRprog    = qlnorm(params[i,7], priors$rrprog_a, priors$rrprog_b),
    RRcdr     = qlnorm(params[i,8], priors$rrcdr_a, priors$rrcdr_b),
    psi       = qbeta(params[i,9], priors$psi_a, priors$psi_b),
    relapse_mid   = qlnorm(params[i,10], priors$relapse_m, priors$relapse_s),
    mutb      = qlnorm(params[i,11], priors$mutb_m, priors$mutb_s),
    recov     = qlnorm(params[i,12], priors$recov_m, priors$recov_s),
    mu      = 1/62,           
    pm      = 0.43,           
    pf      = 0.38,           
    tau     = 2,
    theta   = 0.73            
  )$run(t));
  
  outprev[[k]] = tail(outputs[[k]]$prev_mid,1)
  outmfr[[k]]  = tail(outputs[[k]]$MFratio,1)
  preve[[k]]   = outputs[[k]]$prev_mid[50]
  mfre[[k]]    = outputs[[k]]$MFratio[50]
  
  if(k==1000){
    print(J)
    J <- J+1
    save(outprev,file=paste0("tmp/outprev",J,".Rdata"))
    save(outmfr,file=paste0("tmp/outmfr",J,".Rdata"))
    save(preve,file=paste0("tmp/preve",J,".Rdata"))
    save(mfre,file=paste0("tmp/mfre",J,".Rdata"))
    k <- 0
    outputs <- outprev <- outmfr <- preve <- mfre <- list()
  }
}

## any left-overs from loop
if(length(outprev)>0){
  J <- J+1
  save(outprev,file=paste0("tmp/outprev",J,".Rdata"))
  save(outmfr,file=paste0("tmp/outmfr",J,".Rdata"))
  save(preve,file=paste0("tmp/preve",J,".Rdata"))
  save(mfre,file=paste0("tmp/mfre",J,".Rdata"))
}

end_time<-Sys.time()



time_diff <-end_time-start_time



library(data.table)

equil_prev <- list()
equil_mfr  <- list()
early_prev   <- list()
early_mfr    <- list()
for(i in 1:78){
  load(paste0("outprev",i,".Rdata"))
  load(paste0("outmfr",i,".Rdata"))
  load(paste0("preve",i,".Rdata"))
  load(paste0("mfre",i,".Rdata"))
  
  equil_prev[[i]] <- outprev
  equil_mfr[[i]]  <- outmfr
  early_prev[[i]] <- preve
  early_mfr[[i]] <- mfre
}

str(outprevL)
equilPrev <- unlist(equil_prev)
equilMfr <- unlist(equil_mfr)
earlyPrev <- unlist(early_prev)
earlyMfr <- unlist(equil_prev)

## analyse model
equil_prev <- SA$analyze$sobol$analyze(problem, as.array(equilPrev),
                                       calc_second_order=TRUE, print_to_console=TRUE)
equil_mfr <- SA$analyze$sobol$analyze(problem, as.array(equilMfr),
                                      calc_second_order=TRUE, print_to_console=TRUE)



SB$parameter <- factor(SB$parameter, levels = c("beta","contacts","reactivation","fast-progression",
                                                "stabilisation", 'CDR', 'progression-RR','CDR-RR',
                                                'protection', 'relapse','TB-CFR', 'recovery'),
                   ordered = TRUE, labels=c(expression(beta),expression(rho),"reactivation","fast-progression",
                                            "stabilisation", 'CDR', 'progression-RR','CDR-RR',
                                            'protection', 'relapse','TB-CFR', 'recovery'))



library(dplyr)
library(ggplot2)
SB<-SB %>% mutate(
  conf.upper = value + confint,
  conf.lower = if_else(confint > value, 0, value - confint))

ggplot(SB) + aes(x = parameter, y = value, fill = index,
                 ymin = conf.lower, ymax = conf.upper) + 
  geom_col(position = "dodge", width = 0.5, colour= 'navy') + 
  geom_errorbar(position = position_dodge(0.5), width = 0.25) +
  facet_wrap(~output, scales = "free_y", ncol = 2) +
  theme(legend.position = "top") + 
  # scale_fill_brewer(NULL) + theme_bw() + theme(text = element_text(size = 14))+
  #scale_fill_manual(values=c("#E69F00", "#66CC99"))+ theme_bw() + theme(text = element_text(size = 13))+
  scale_fill_brewer(palette = 'Dark2') + theme_bw() + theme(text = element_text(size = 14))+
  coord_flip()+scale_x_discrete(labels = c(expression(beta),
                                           expression(rho),
                                           expression(nu), 
                                           expression(epsilon),
                                           expression(kappa),
                                           expression(delta),
                                           expression(alpha),
                                           expression(pi),
                                           expression(psi),
                                           expression(omega),
                                           expression(mu[t]), 
                                           expression(gamma)))+
  ylab("indices") + xlab(NULL)+ggtitle('Sensitivity of model outputs to input parameters')
