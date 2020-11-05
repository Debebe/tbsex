
setwd("~/Dropbox/SexTB")

library(data.table)
library(patchwork)
library(ggplot2)
library(GGally)
library(ggpointdensity)
library(sp)              # to extract polygon features
library(PlaneGeometry)   # to draw ellipse
library("plot3D")
library("scatterplot3d")
library(rstan)
library(Rmisc)
library(dplyr)


load("~/Dropbox/SexTB/Clean/Data/sampsE.Rdata")      # load posterior samples from Stan
source('~/Dropbox/SexTB/Clean/Codes/parms_sextb.R')  # Model parameters
source('~/Dropbox/SexTB/Clean/Codes/Sex_TB_IntervR.R') 

#===== Extract parameters from stan===========
#  this gives an array of iterations, chain and parameters
posterior <- rstan::extract(samps, permuted = FALSE, inc_warmup=FALSE)  
dim(posterior)
DT= as.data.table(posterior[,1,]) # extracting all iterations, chain one and all params

#retaian only relevant params

DT <-DT[, c("psi", "recov", "react","fast","stab","relapse", "mutb",
            "cdrprop", "beta","rho","rrcdr",  "rrprog", "lp__" )]
class(DT)
colnames(DT)
setnames(DT, "react", "react_mid")
setnames(DT, "fast", "fast_mid")
setnames(DT, "cdrprop", "cdr_mid")
setnames(DT, "lp__" , 'lp')

n=200
data <- data.table(psi=sample(DT$psi, n), recov= sample(DT$recov, n), react_mid=sample(DT$react_mid, n),
                  fast_mid=sample(DT$fast_mid, n), stab=sample(DT$stab, n),relapse_mid=sample(DT$relapse, n),
                  mutb= sample(DT$mutb, n), cdr_mid =sample(DT$cdr_mid, n),beta =sample(DT$beta, n), 
                  rho = sample(DT$rho, n),  rrcdr = sample(DT$rrcdr, n), rrprog = sample(DT$rrprog, n))


##Wrapper function to run intervention
out<-list()
PL= pms
PL=list()
PL$RRcdr <- NULL   # remove RRcdr as it is an intervention

run_intervention <-function(data, RR) {
 
  
  for (i in 1:nrow(data)){
    
    PL$beta <- data$beta[i]
    PL$cdr_mid <- data$cdr_mid[i]
    PL$mutb <- data$mutb[i]
    PL$relapse_mid  <- data$relapse[i]
    PL$fast_mid <- data$fast_mid[i]        
    PL$stab <- data$stab[i]     
    PL$react_mid <- data$react_mid[i]
    PL$recov <- data$recov[i]
    PL$psi <- data$psi[i]
    PL$theta <- 0.89
    PL$pf  <- 0.35
    PL$pm  <- 0.45
    PL$rho <-data$rho[i]
    PL$RRprog<-data$rrprog[i]
    #====intervention
    tmax=2030
    interv_year <- 0:tmax
    PL$interv_year <- interv_year 
    
    baseline <- rep(data$rrcdr[i], tmax+1)
    interv= baseline
    interv[2021:tmax+1] <-1

    if (RR==0){
      PL$RRcdr_interv =baseline
    }
    else {
      PL$RRcdr_interv=interv
    }
    
    mod0 <- Sex_TB_Interv(user=PL)
    outl <- mod0$run(t= 0:tmax, len = tmax) #TODO
    
    outl <- data.frame(outl)

    
    out[[i]] <- data.table(MFR= tail(outl$MFratio, 15),
                           prev= tail(outl$prev_mid,15),
                           RRcdr= tail(PL$RRcdr_interv,15),
                           year= tail(outl$t, 15),
                           index=i)
  }
  return(out)
}

# call function
baseline <- run_intervention(data , RR=0)
Interven <- run_intervention(data , RR=1)

# extract data

out_baseline <- do.call('rbind', baseline)
out_intervention <- do.call('rbind', Interven)

# create new variables for plotting

out_baseline[, c('Intervention', 'Region'):= list('No', 'whole')]
out_intervention[, c('Intervention', 'Region'):= list('Yes', 'whole')]

DTT <-rbind(out_baseline, out_intervention)

yr=2030
#percentage reduction in outcomes
pop_attr_frac_prev<- DTT[year==yr, .(PAF=100*(first(prev)-last(prev))/first(prev)), by=.(index, Region)]
pop_attr_frac_mfr<- DTT[year==yr, .(PAF=100*(first(MFR)-last(MFR))/first(MFR)), by=.(index, Region)]

#======confidence interval========
PA_prev <-pop_attr_frac_prev[, .(mean= mean(PAF), low= quantile(PAF, 0.025), hi=quantile(PAF, 0.975)), by = Region]
PA_mfr <-pop_attr_frac_mfr[, .(mean= mean(PAF), low= quantile(PAF, 0.025), hi=quantile(PAF, 0.975)), by = Region]



PA_prev$outcome <- 'prevalence'
PA_mfr$outcome  <- 'M:F ratio'
PA_combined =rbind(PA_prev, PA_mfr)


