


#setwd("~/Dropbox/SexTB/Clean/linux")
# setwd("C:/Users/dshaweno/Dropbox/SexTB/Clean/linux")
# setwd("C:/Users/dshaweno/Dropbox/SexTBD/October")
# 
# setwd("C:/Users/dshaweno/Dropbox/SexTBD/stan_october")
setwd("C:/Users/dshaweno/Dropbox/SexTB/Clean/new")
setwd("~/Dropbox/SexTB/Clean/new")

library(rstan)
library(bayesplot)

set.seed(1234)
source('~/Dropbox/SexTB/Clean/linux/priorsRevised_new.R')  # source priors
source('C:/Users/dshaweno/Dropbox/SexTB/Clean/linux/priorsRevised_new.R')  # source priors
## New prior for rrcdr

# priors$rrprog_a <- 0.49
# priors$rrprog_b <- 0.096

priors$rrcdr_a=-0.298
priors$rrcdr_b=0.2026 #  # determined based on the M:F ratio of (notifm/prevm)/(notiff/prevf)


priors$rrprog_a <- -0.298
priors$rrprog_b <- 4*0.2026 

# priors$rho_a=23.2559
# priors$rho_b=17.61564

priors$rho_a=-0.57
priors$rho_b=0.085

## reduce variance in kappa/stab
priors$stab_s <- 0.34/5

priors$cdr_m <-5.6
priors$cdr_s <-3.5


lm <- stan_model('sextb_mid.stan')




#--------------Uganda(2014/2015- WHO-2016, )
priorsU <- c(priors,list(targetd=401e-5,targetdsig=401e-5/10,
                         #targetd_f=178e-5,targetd_fsig=178e-5/10,
                         #targetd_m=734e-5,targetd_msig=734e-5/10,
                         targetMF=4.1,targetMFsig=4.1/10,
                         targetNf= 79e-5, targetNfsig=79e-5/10,
                         targetNm= 161e-5, targetNmsig=161e-5/10,
                         #targetN= 118e-5, targetNsig=118e-5/10,
                         Qtol=1e-6, pm=0.43,pf=0.38, tau=2, theta=0.73, mu=0.01629))



#--------------Ethiopia
priorsE <- c(priors,list(targetd=277e-5,targetdsig=277e-5/10,
                         # targetd_f=246e-5,targetd_fsig=246e-5/10,
                         # targetd_m=304e-5,targetd_msig=304e-5/10,
                         targetMF=1.2,targetMFsig=1.2/10,
                         #targetN= 129e-5, targetNsig=129e-5/10,
                         targetNf= 169e-5, targetNfsig=169e-5/10,
                         targetNm= 209e-5, targetNmsig=209e-5/10,
                         Qtol=1e-6, pm=0.44,pf=0.35, tau=2, theta=0.89, mu= 0.016129))



niter <- 30000/3; nchains <- 4

sampsU <- sampling(lm, data = priorsU, iter = niter, chains = nchains, verbose = TRUE,
                   control = list(adapt_delta=0.95))
sampsU <- sampling(lm, data = priorsU, iter = niter, chains = nchains, verbose = TRUE)


sampsE <- sampling(lm, data = priorsE, iter = niter, chains = nchains, verbose = TRUE,
                   control = list(adapt_delta=0.95))
sampsE <- sampling(lm, data = priorsE, iter = niter, chains = nchains, verbose = TRUE)

save(sampsE, file= "~/Dropbox/SexTB/Clean/Data/sampsE.Rdata") #Ethiopia
save(sampsU, file= "~/Dropbox/SexTB/Clean/Data/sampsU.Rdata") #Uganda





save(samps, file= ("C:/Users/dshaweno/Dropbox/SexTB/Clean/new/sampsU.Rdata")) #Uganda

samps <-sampsE

save(samps, file= ("C:/Users/dshaweno/Dropbox/SexTB/Clean/new/sampsE.Rdata")) #Uganda
     
dev.off()
to_plot <- c('d', 'n', 'nf','nm', 'Qerr','MF', 'beta','rho', 'contacts','fast','react','cdrprop',  'rrcdr', 'rrprog', 'relapse', 'recov', 'stab', 'mutb','psi')
TP <- traceplot(samps, pars = to_plot)
print(samps, pars = to_plot, digits=5)
posterior <- as.matrix(samps)
dens<-mcmc_dens(posterior,pars=to_plot) 

