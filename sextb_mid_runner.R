

setwd("~/Dropbox/SexTB/Clean/Codes")

library(rstan)
library(bayesplot)


source('priorsRevised_new.R')  # source priors


lm <- stan_model('sextb_mid.stan') # compile model

#--------------Uganda
priorsU <- c(priors,list(targetd=401e-5,targetdsig=401e-5/10,
                         targetMF=4.1,targetMFsig=4.1/10,
                         targetNf= 79e-5, targetNfsig=79e-5/10,
                         targetNm= 161e-5, targetNmsig=161e-5/10,
                         Qtol=1e-6, pm=0.43,pf=0.38, tau=2, theta=0.73, mu=0.01629))



#--------------Ethiopia
priorsE <- c(priors,list(targetd=277e-5,targetdsig=277e-5/10,
                         targetMF=1.2,targetMFsig=1.2/10,
                         targetNf= 169e-5, targetNfsig=169e-5/10,
                         targetNm= 209e-5, targetNmsig=209e-5/10,
                         Qtol=1e-6, pm=0.44,pf=0.35, tau=2, theta=0.89, mu= 0.016129))



niter <- 30000/3; nchains <- 4

sampsU <- sampling(lm, data = priorsU, iter = niter, chains = nchains, verbose = TRUE,
                   control = list(adapt_delta=0.95))    # Uganda

sampsE <- sampling(lm, data = priorsE, iter = niter, chains = nchains, verbose = TRUE,
                   control = list(adapt_delta=0.95))    # Ethiopia

save(sampsE, file= "~/Dropbox/SexTB/Clean/Data/sampsE.Rdata") #Ethiopia
save(sampsU, file= "~/Dropbox/SexTB/Clean/Data/sampsU.Rdata") #Uganda


samps <-sampsE


dev.off()
to_plot <- c('d', 'n', 'nf','nm', 'Qerr','MF', 'beta','rho', 'contacts','fast','react','cdrprop',  'rrcdr', 'rrprog', 'relapse', 'recov', 'stab', 'mutb','psi')
TP <- traceplot(samps, pars = to_plot)
print(samps, pars = to_plot, digits=5)
posterior <- as.matrix(samps)
dens<-mcmc_dens(posterior,pars=to_plot) 

