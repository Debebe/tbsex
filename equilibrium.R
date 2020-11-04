

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
  lambda[]         <- sum(lambda_prod[i, ])             # rowSums
  output(lambda) <- TRUE                                

  #=====derivatives=====
  N[]  <- S[i] + E[i] + L[i] + I[i]+ T[i] + R[i]
  totN <- sum(N[])
  output(N[])  <- TRUE
  output(totN) <- TRUE
  
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
  rho         <- user()
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


library(data.table)
library(rstan)
load('~/Dropbox/SexTB/Clean/Data/sampsU.Rdata')   # load R0 & MF function
source('~/Dropbox/SexTB/Clean/Codes/R0.R')        # load R0 & MF function

posterior <- rstan::extract(samps, permuted = FALSE, inc_warmup=FALSE)  
dim(posterior)
DT= as.data.table(posterior[,1,])

globalPars <- list(beta = mean(DT$beta),
                   react_mid = mean(DT$react),
                   fast_mid= mean(DT$fast),
                   stab=mean(DT$stab),
                   RRprog= mean(DT$rrprog),
                   RRcdr= mean(DT$rrcdr),
                   psi = mean(DT$psi),      #partial protection of previous infection
                   relapse_mid =mean(DT$relapse),
                   rho=mean(DT$rho),
                   tau=1/0.5,      # duration of treatment 
                   theta=0.73,
                   pm=0.43,
                   pf=0.38)  


#== intermediate parameters=====

localPars <-list(CDRaverage =mean(DT$cdrprop),
                 life_expect= 61.373,
                 infect_dur= 3)

#==wrap up params used by odin
genPars<-function(globalPars,localPars){
  list2env(c(globalPars, localPars),envir = environment())
  
  mu    <- 1/life_expect
  recov <- mean(DT$recov)
  mutb  <- mean(DT$mutb)
  list(mu=mu,mutb=mutb, recov=recov,cdr_mid=CDRaverage) 
}

pars <-genPars(globalPars, localPars)  

pms<- c(globalPars, list(mu=pars$mu, mutb=pars$mutb, 
                         recov=pars$recov, cdr_mid=pars$cdr_mid))

par <- pms

mm <-pms
mm$RRcdr<-1
mm$RRprog<-1
mm$rho <-0   
mm$pf <-mm$pm

#================================
model      <- SexTB(user = mm)
y          <- model$run(t= seq(0, 4000,len = 4000)) 
out        <- data.frame(y) 
new_prev   <-round(tail(out$prev_mid, 1))
#================================

#-------prepare params for R0-------
par$rho <-0          # 0  is random mixing

cdr_rate  <- pms$cdr_mid *(pms$mutb +pms$recov + pms$mu)/(1-pms$cdr_mid)

female_params <- list(
  beta =par$beta,
  react= 2*par$react_mid/(1+par$RRprog),
  fast = 2*par$fast_mid/(1+par$RRprog), 
  relapse = 2*par$relapse_mid/(1+par$RRprog), 
  
  stab = par$stab, 
  recov= par$recov,
  mutb = par$mutb, 
  psi  = par$psi,     # partial protection 
  mu=par$mu,
  fracinf =par$pf,
  tau=par$tau,
  theta=par$theta,
  cdr =2*cdr_rate/(par$RRcdr +1)) 

male_params         <- female_params
male_params$fast    <- par$RRprog*female_params$fast
male_params$react   <- par$RRprog*female_params$react
male_params$relapse <- par$RRprog*female_params$relapse

male_params$fracinf <- par$pm
male_params$cdr     <- par$RRcdr*female_params$cdr


#-----------handle MF ratio in Reproduction number-------


RR <- function (mf_ratio, r0) {
  #mf_ratio=1
  # produces "RRprog" related to a particular M/F ratio in R 
  solveRRprog <-function(x){
    par$RRcdr = 1
    f   <- female_params
    m   <- male_params
    
    f$react   <- 2*par$react_mid/(1+x)
    f$fast    <- 2*par$fast_mid/(1+x)
    f$relapse <- 2*par$relapse_mid/(1+x)
    m$fast    <- x*f$fast
    m$react   <- x*f$react
    m$relapse   <- x*f$relapse
    
    f$cdr     <-cdr_rate/(1 +par$RRcdr)
    m$cdr     <- par$RRcdr*f$cdr
    
    m$fracinf <-f$fracinf
    Rm  <- singleR0(m)$R0
    Rf  <- singleR0(f)$R0
    R0  <- R0combined(par$rho,Rm,Rf)
    
    #-----
    Rm/Rf-mf_ratio
  }
  (mfrRR <<-uniroot(solveRRprog,lower=0,upper=10)$root)
  
  
  thresh <- function(x,verbose=FALSE){
    par$RRcdr = 1
    
    f   <- female_params
    m   <- male_params
    
    f$react <- 2*par$react_mid/(1+mfrRR)
    f$fast  <- 2*par$fast_mid/(1+mfrRR)
    f$relapse  <- 2*par$relapse_mid/(1+mfrRR)
    m$fast  <- mfrRR*f$fast
    m$react <- mfrRR*f$react
    m$relapse <- mfrRR*f$relapse
    
    f$cdr   <- 2*par$cdr_mid/(1+par$RRcdr)
    m$cdr   <- par$RRcdr*f$cdr
    m$fracinf <-f$fracinf
    
    f$beta <- x
    m$beta <- x
    
    Rm  <- singleR0(m)$R0
    Rf  <- singleR0(f)$R0
    R0  <- R0combined(par$rho,Rm,Rf)
    
    R0-r0
   }
  (equilbeta <<- uniroot(f=thresh,interval = c(1e-5,40))$root) #solve beta with R0=1
  
  #--------------------------------
  # determine prevalence at complete random mixing and RRprog(mfrRR==1)
  pp<-pms
  pp$RRcdr<-1
  pp$rho<-0
  pp$RRprog  <-1.00000   # RR when MF ratio is 
  if (r0==1) {
  pp$beta <-6.405194   # when R0=1 and mfr=1
   
  } else if (r0==2) {
  pp$beta <-11.3208  # when R0=2 and mfr=1
    
  } else {
  pp$beta <-16.9812  # when R0=3 and mfr=1
  }

  
  pp$pm<-pp$pf
  model      <- SexTB(user = pp)
  y          <- model$run(t= seq(0, 4000,len = 4000))
  out        <- data.frame(y)
  (prev_at_rho_0 <<-tail(out$prev_mid, 1))
  #-----------------------------------
  
  # tune beta for prevalence
  func <-function(x){
    pmm <-pp
    pmm$pm<-pmm$pf
    pmm$RRprog<-mfrRR
    pmm$beta <-x
    pmm$rho=0
    model      <- SexTB(user = pmm)
    y          <- model$run(t= seq(0, 4000,len = 4000)) 
    out        <- data.frame(y) 
    new_prev <<-round(tail(out$prev_mid, 1) )
    new_prev-prev_at_rho_0
    
  }
  
  (newbeta <<- uniroot(f=func,interval = c(1e-5,40))$root) #solve beta with R0=1
  
  #---------------------
  rho <- runif(100, 0,1)  # deliberately including 0 and 1
  equil_out <- list(length= length(rho))   
  
  for (i in 1: length(rho)) {
    RRprog=mfrRR
    RRcdr=1
    equil_params<-pms
    equil_params$beta   <- newbeta
    equil_params$RRprog <- RRprog
    equil_params$RRcdr  <- RRcdr
    equil_params$pm     <- equil_params$pf 
    par <- equil_params
    
    f   <- female_params
    m   <- male_params
    f$beta <-equilbeta
    m$beta <-equilbeta
    f$react= 2*par$react_mid/(1+par$RRprog)
    f$fast = 2*par$fast_mid/(1+par$RRprog)    
    f$relapse = 2*par$relapse_mid/(1+par$RRprog)

    
    f$cdr  = 2*par$cdr/(par$RRcdr +1)
    f$fracinf <- par$pf
    m <-f
    m$fast    <- par$RRprog*f$fast
    m$react   <- par$RRprog*f$react
    m$relapse    <- par$RRprog*f$relapse
    m$fracinf <- par$pm
    m$cdr     <- par$RRcdr*f$cdr
    
    par$rho   <- rho[i]
    
    # call model
    model   <- SexTB(user = par)
    y       <- model$run(t= seq(0, 4000,len = 4000)) 
    out     <- data.frame(y) 
    #-----R0----------------------------------------
    Rm  <- singleR0(m)$R0
    Rf  <- singleR0(f)$R0
    R0  <- R0combined(par$rho,Rm,Rf)
    
    equil_out[[i]] <-data.frame(MFR= tail(out$MFratio, 1), 
                                prev=tail(out$prev_mid, 1),
                                pm =tail(out$prev.1., 1),
                                pf =tail(out$prev.2., 1),
                                rho= par$rho,
                                R0= R0,
                                Rm=Rm,
                                Rf=Rf,
                                Rmfr=mf_ratio,
                                id   = i)
    print(i)
  }
  return(equil_out)
}

r0<-2
equil_out1<- RR(1, r0)
equil_out2<- RR(2, r0)
equil_out3<- RR(3, r0)
equil_out4<- RR(4, r0)



library(gridExtra)
library(data.table)
library(ggplot2)
library(patchwork)

df1<-do.call("rbind", equil_out1)
df2<-do.call("rbind", equil_out2)
df3<-do.call("rbind", equil_out3)
df4<-do.call("rbind", equil_out4)

df <-rbind(df1, df2, df3,df4)
df$Rmfr <-round(df$Rmfr)
df$Rmfr <-factor(df$Rmfr)
df <- data.table(df)

# initial MF ratio during exponential growth


rho <- seq(from=0,to=1,len=1e2)          
r0data <- data.table(rho=rep(rho,4),ratio=rep(1:4,each=length(rho))) #data
r0data[,R0val:=R0(rho,ratio*1,1)/R0(0,ratio*1,1)]              #evaluate Rf = 1
r0data[,ratio:=factor(ratio)]                     #make a factor

r0data[,MFrat:=MFinitial(rho,as.numeric(as.character(ratio)),1)] #NOTE this is assuming Prat=1



setnames(df, 'Rmfr', 'Rm:Rf')
setnames(r0data, 'ratio', 'Rm:Rf')

n=12
a <-ggplot(df, aes(rho, R0, group=`Rm:Rf`, lty=`Rm:Rf`)) +geom_line()+ 
  theme_classic() +  xlab("Assortativity") + ylab('R0')  + 
  theme(legend.position = "none") +ggtitle('a.') + 
  theme(text = element_text(size = n,family = 'Times'))

b <- ggplot(r0data[rho<1],aes(rho,log(MFrat),group=`Rm:Rf`, lty=`Rm:Rf`)) +
  geom_line() +theme_classic() +
  xlab('Assortativity') + ylab('logarithm M:F ratio ') +
  theme(legend.position = "none") +ggtitle('b. ') + 
  theme(text = element_text(size = n,family = 'Times'))



c<-ggplot(df, aes(rho, log(MFR), group=`Rm:Rf`, lty=`Rm:Rf`)) +geom_line() + 
  xlab("Assortativity") + ylab('logarithm M:F ratio ')+ 
   theme_classic()+ggtitle("c.")+ theme(text = element_text(size = n,family = 'Times'))


d<- ggplot(df, aes(rho, prev, group=`Rm:Rf`, lty=`Rm:Rf`)) +geom_line() + theme_classic() + theme(legend.position = "none") +
  xlab("Assortativity") + ylab("prevalence /100,000") + ggtitle('d. ') + 
  theme(text = element_text(size = n,family = 'Times')) 


all_figs <- a+b+c +d+
  plot_layout(ncol = 4)

# use common legend at bottom
combined <- all_figs & theme_classic() + 
  theme(plot.title = element_text(hjust = 0))
combined <- combined & theme(legend.position = "bottom") 
combined + plot_layout(guides = "collect")

ggsave('Clean/plots/equil.png', width=7, height = 3, dpi=300)



