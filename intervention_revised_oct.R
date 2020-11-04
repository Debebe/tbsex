
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
library(Rmisc) # for CI
library(dplyr)


load("~/Dropbox/SexTB/Clean/Data/sampsU.Rdata")  # load posterior samples from Stan
source('~/Dropbox/SexTB/Clean/Codes/parms_sextb.R')  # Model parameters
source('~/Dropbox/SexTB/Clean/Codes/Sex_TB_IntervR.R') 

#===== Extract parameters from stan===========
#  this gives an array of iterations, chain and parameters
posterior <- rstan::extract(samps, permuted = FALSE, inc_warmup=FALSE)  
dim(posterior)
DT= as.data.table(posterior[,1,]) # extracting all iterations, chain one and all params

#retian only relevant params

DT <-DT[, c("psi", "recov", "react","fast","stab","relapse", "mutb",
            "cdrprop", "beta","rho","rrcdr",  "rrprog", "lp__" )]
class(DT)
colnames(DT)
setnames(DT, "react", "react_mid")
setnames(DT, "fast", "fast_mid")
setnames(DT, "cdrprop", "cdr_mid")
setnames(DT, "lp__" , 'lp')


# sample from different regions
### Reactivation and Differential progression

ggplot(DT, aes(rrcdr, rrprog)) +geom_point() # positive corr

main.r=0.05
minor.r=0.03
# lower region
ellL <- Ellipse$new(center = c(quantile(DT$rrcdr,0.25), quantile(DT$rrprog, 0.25)), rmajor = main.r, rminor = minor.r, alpha=90)
ellpathL <- ellL$path() # ellipse as a path
ellpathL <- rbind(ellpathL, ellpathL[1,]) # the path is not closed; close it

#mid region
ellM <- Ellipse$new(center = c(quantile(DT$rrcdr,0.5), quantile(DT$rrprog, 0.5)), rmajor = main.r, rminor = minor.r, alpha=90)
ellpathM <- ellM$path() # ellipse as a path
ellpathM <- rbind(ellpathM, ellpathM[1,]) # the path is not closed; close it

# upper region
ellU <- Ellipse$new(center = c(quantile(DT$rrcdr,0.75), quantile(DT$rrprog, 0.75)), rmajor = main.r, rminor = minor.r, alpha=90)
ellpathU <- ellU$path() # ellipse as a path
ellpathU <- rbind(ellpathU, ellpathU[1,]) # the path is not closed; close it



#Now plot the data and ellipses based on ellipse parameters above

pd <-ggplot() +
  geom_point(aes(x = rrcdr, y = rrprog), DT) +
  geom_path(aes(x = x, y = y), as.data.frame(ellpathL) , color = "red") +
  geom_path(aes(x = x, y = y), as.data.frame(ellpathM), color = "blue")+
  geom_path(aes(x = x, y = y), as.data.frame(ellpathU), color = "green")

pd

# extract parameters from ellipse

build <- ggplot_build(pd)$data
points <- build[[1]]
ellipseL <- build[[2]]
ellipseM <- build[[3]]
ellipseU <- build[[4]]

head(points, 3)
head(ellipseL,3)

dat <- data.frame(
  points[1:2], 
  in.ellipseL = as.logical(point.in.polygon(points$x, points$y, ellipseL$x, ellipseL$y)),
  in.ellipseM = as.logical(point.in.polygon(points$x, points$y, ellipseM$x, ellipseM$y)),
  in.ellipseU = as.logical(point.in.polygon(points$x, points$y, ellipseU$x, ellipseU$y))
)

## Determine points within ellipses

dat <- data.table(dat)
dat[in.ellipseL==TRUE, c('in_ellipse', 'region'):= list(TRUE, 'L')]
dat[in.ellipseM==TRUE, c('in_ellipse', 'region'):= list(TRUE, 'M')]
dat[in.ellipseU==TRUE, c('in_ellipse', 'region'):= list(TRUE, 'U')]

dsa <- dat[, .(x, y,in_ellipse, region )]
rrcdr_rrprog_sample <-dsa      # save samples with in ellipse

#Visualise points with in selected regions
p1 <-ggplot(dsa) +
  geom_point(aes(x, y,color=region),size=1) + 
  stat_ellipse(aes(x, y,color=region)) +xlab('Reactivation') + theme_classic()+
  ylab('RR-progression')

p1

# Process data for odin
dsa <- data.table(dsa)
names(DT)
data_in_ellipse <-dsa[in_ellipse==TRUE]
length <- dim(data_in_ellipse)[1]
head(data_in_ellipse)

set.seed(123)

# sampling from main database based on size in the ellipse
index <- sample(1:nrow(DT), nrow(data_in_ellipse))
data_3Regions <-DT[index, ]

# replace values in the main database by values from ellipse
data_3Regions$rrcdr     <- data_in_ellipse$x
data_3Regions$rrprog    <- data_in_ellipse$y
data_3Regions$region    <- data_in_ellipse$region
data_3Regions <-data.table(data_3Regions)


DT %>%
  summarize(quant25= quantile(lp,0.25),
            quant50= quantile(lp,0.50),
            quant95= quantile(lp,0.95))

data_3Regions[,.(quant25= quantile(lp,0.25),
                 quant50= quantile(lp,0.50), 
                 quant95=quantile(lp,0.95)), by= 'region']

# densities at each region
all_regions<- DT[, .(rrcdr, rrprog, lp)]   # select columns
all_regions[, region:= 'all']                  # create region variable

sampled <- data_3Regions[, .(rrcdr, rrprog,lp, region)] # select colms
all<- rbind(all_regions,sampled )

rrcdr_rrprog_lp <- all

# look the density of lp posterior
d1 <-ggplot(all,aes(x=lp, fill=region)) + geom_density(alpha=0.5) + 
  xlab('log posterior') + 
  ggtitle("rrcdr and RR-progression") + theme_classic()
d1


#### Wrapper function to run intervention
out<-list()
PL= pms
PL$RRcdr <- NULL   # remove RRcdr as it is intervention

runn <-function(data, group,RR) { #make sure interference with names with in func
  if (group=='L'){
    data <-  data_3Regions[region=='L']
  } else if (group=='M') {
    data <-  data_3Regions[region=='M']
  } else {
    data <-  data_3Regions[region=='U']
  }
  
  
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
    PL$theta <- pms$theta
    PL$pf  <- 0.38
    PL$pm  <- 0.43
    PL$rho <-data$rho[i]
    PL$RRprog<-data$rrprog[i]
    #====intervention
    tmax=2030
    interv_year <- 0:tmax
    PL$interv_year <- interv_year 
    
    baseline <- rep(data$rrcdr[i], tmax+1)
    interv= baseline
    interv[2021:tmax+1] <-1
    #PL$RRcdr_interv<-interv
    
    if (RR==0){
      PL$RRcdr_interv =baseline
    }
    else {
      PL$RRcdr_interv=interv
    }
    
    mod0 <- Sex_TB_Interv(user=PL)
    outl <- mod0$run(t= 0:tmax, len = tmax) #TODO
    
    outl <- data.frame(outl)
    #outl$rr <- PL$RRcdr_interv

    
    out[[i]] <- data.table(MFR= tail(outl$MFratio, 15),
                           prev= tail(outl$prev_mid,15),
                           RRcdr= tail(PL$RRcdr_interv,15),
                           year= tail(outl$t, 15),
                           index=i)
  }
  return(out)
}

# call function
Lower_baseline <- runn(data_3Regions ,group='L', RR=0)
Lower_Interven <- runn(data_3Regions ,group='L', RR=1)
mid_baseline <- runn(data_3Regions ,group='M', RR=0)
mid_Interven <- runn(data_3Regions ,group='M', RR=1)
upper_baseline <- runn(data_3Regions ,group='U', RR=0)
upper_Interven <- runn(data_3Regions ,group='U', RR=1)

# extract data

out_Lb <- do.call('rbind', Lower_baseline)
out_LI <- do.call('rbind', Lower_Interven)

out_Mb <- do.call('rbind', mid_baseline)
out_MI <- do.call('rbind', mid_Interven)

out_Ub <- do.call('rbind', upper_baseline)
out_UI <- do.call('rbind', upper_Interven)

# create new variables for plotting

out_Lb[, c('Intervention', 'Region'):= list('No', 'L')]
out_LI[, c('Intervention', 'Region'):= list('Yes', 'L')]

out_Mb[, c('Intervention', 'Region'):= list('No', 'M')]
out_MI[, c('Intervention', 'Region'):= list('Yes', 'M')]

out_Ub[, c('Intervention', 'Region'):= list('No', 'U')]
out_UI[, c('Intervention', 'Region'):= list('Yes', 'U')]

DTT <-rbind(out_Lb, out_LI, out_Mb, out_MI, out_Ub,out_UI )
tmpprev<- DTT[, .(prev=mean(prev)), by=.(year, Intervention, Region)]
tmpmfr <- DTT[, .(MFR=mean(MFR)), by=.(year, Intervention, Region)]

save(DTT, file= "Clean/Data/rrcdr_rrprog_posterior_regions_ethiopia.Rdata") #Uganda

yr=2030
# pop attributable fraction
prev_2030 <- tmpprev[year==yr]
mfr_2030 <- tmpmfr[year==yr]


PAF_prevalence <- prev_2030 %>%             
  group_by(Region) %>%            
  summarize(PAF = 100*(first(prev)-last(prev))/first(prev))
#summarize(prev2 = 100*(max(prev)-min(prev))/max(prev))

PAF_mfr <- mfr_2030 %>%             
  group_by(Region) %>%            
  summarize(PAF = 100*(first(MFR)-last(MFR))/first(MFR))
#summarize(prev2 = 100*(max(prev)-min(prev))/max(prev))

a<-ggplot(tmpprev[year==yr], aes(x=Region, y=prev, fill= Intervention)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal() + coord_flip() + ylab('Prevalence')

b<-ggplot(tmpmfr[year==yr], aes(x=Region, y=MFR, fill= Intervention)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal() + coord_flip() + ylab('M:F ratio')

a/b  + plot_annotation(
  title = 'Intervention impact at different regions of joint posterior\ 
  distribution of rrcdr and rrprog')

#======PAF intervention=========
#===population attributable fraction====
pop_attr_frac_prev<- DTT[year==yr, .(PAF=100*(first(prev)-last(prev))/first(prev)), by=.(index, Region)]
pop_attr_frac_mfr<- DTT[year==yr, .(PAF=100*(first(MFR)-last(MFR))/first(MFR)), by=.(index, Region)]

#======confidence interval========
PA_prev <-pop_attr_frac_prev[, .(mean= mean(PAF), low= quantile(PAF, 0.025), hi=quantile(PAF, 0.975)), by = Region]
PA_mfr <-pop_attr_frac_mfr[, .(mean= mean(PAF), low= quantile(PAF, 0.025), hi=quantile(PAF, 0.975)), by = Region]



PA_prev$outcome <- 'prevalence'
PA_mfr$outcome  <- 'M:F ratio'
PA_combined =rbind(PA_prev, PA_mfr)

paf1 <- ggplot(PA_combined, aes(x=Region, y=mean, fill=outcome)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=low, ymax=hi), width=.2,position=position_dodge(.9))+
  theme_classic() + scale_fill_brewer(palette="Paired")+
  ylab('Reduction (%)') + coord_flip() + ggtitle('a')


###=====rho and rrprog=====

ggplot(DT, aes(x=rho, y=rrprog)) +
  geom_point(size=1, shape=15)

main.r=0.05
minor.r=0.03

ellL <- Ellipse$new(center = c(quantile(DT$rho,0.25), quantile(DT$rrprog, 0.75)), rmajor = main.r, rminor = minor.r, alpha = 90)
ellpathL <- ellL$path() # ellipse as a path
ellpathL <- rbind(ellpathL, ellpathL[1,]) # the path is not closed; close it


ellM <- Ellipse$new(center = c(quantile(DT$rho,0.50), quantile(DT$rrprog, 0.50)), rmajor = main.r, rminor = minor.r, alpha = 90)
ellpathM <- ellM$path() # ellipse as a path
ellpathM <- rbind(ellpathM, ellpathM[1,]) # the path is not closed; close it


ellU <- Ellipse$new(center = c(quantile(DT$rho,0.75), quantile(DT$rrprog, 0.25)), rmajor = main.r, rminor = minor.r, alpha = 90)
ellpathU <- ellU$path() # ellipse as a path
ellpathU <- rbind(ellpathU, ellpathU[1,]) # the path is not closed; close it


pd <-ggplot() +
  geom_point(aes(x = rho, y = rrprog), DT) +
  geom_path(aes(x = x, y = y), as.data.frame(ellpathL) , color = "red") +
  geom_path(aes(x = x, y = y), as.data.frame(ellpathM), color = "blue")+
  geom_path(aes(x = x, y = y), as.data.frame(ellpathU), color = "green")

pd

#=====extract all points=====

build <- ggplot_build(pd)$data
points <- build[[1]]
ellipseL <- build[[2]]
ellipseM <- build[[3]]
ellipseU <- build[[4]]

head(points, 3)
head(ellipseL, 3)

#### Capture extracted data into a data frame


dat <- data.frame(
  points[1:2], 
  in.ellipseL = as.logical(point.in.polygon(points$x, points$y, ellipseL$x, ellipseL$y)),
  in.ellipseM = as.logical(point.in.polygon(points$x, points$y, ellipseM$x, ellipseM$y)),
  in.ellipseU = as.logical(point.in.polygon(points$x, points$y, ellipseU$x, ellipseU$y))
)

#====Determine points with in ellipses====


dat <- data.table(dat)
dat[in.ellipseL==TRUE, c('in_ellipse', 'region'):= list(TRUE, 'L')]
dat[in.ellipseM==TRUE, c('in_ellipse', 'region'):= list(TRUE, 'M')]
dat[in.ellipseU==TRUE, c('in_ellipse', 'region'):= list(TRUE, 'U')]

# select useful columns
dsa <- dat[, .(x, y,in_ellipse, region )]
rho_rrprog_sample <-dsa

p2 <-ggplot(dsa) +
  geom_point(aes(x, y,color=region),size=1) + 
  stat_ellipse(aes(x, y,color=region)) + xlab('rho') + 
  ylab('rrprog')
p2

#### Process data for odin


dsa <- data.table(dsa)
names(DT)
data_in_ellipse <-dsa[in_ellipse==TRUE]
length <- dim(data_in_ellipse)[1]
head(data_in_ellipse)

set.seed(123)

# sampling from main database based on size in the ellipse
index <- sample(1:nrow(DT), nrow(data_in_ellipse))
data_3Regions <-DT[index, ]

# replace values in the main database by values from ellipse
data_3Regions$rho       <- data_in_ellipse$x
data_3Regions$rrprog    <- data_in_ellipse$y
data_3Regions$region    <- data_in_ellipse$region
data_3Regions <-data.table(data_3Regions)

DT %>%
  summarize(quant25= quantile(lp,0.25),
            quant50= quantile(lp,0.50),
            quant95=quantile(lp,0.95))
#=====summarise regional data=====

data_3Regions[,.(quant25= quantile(lp,0.25),
                 quant50= quantile(lp,0.50), 
                 quant95=quantile(lp,0.95)), by= 'region']

#### Visualise the densities at sampled regions


all_regions<- DT[, .(rho, rrprog, lp)]
all_regions[, region:= 'all']

sampled <- data_3Regions[, .(rho, rrprog,lp, region)]
all<- rbind(all_regions,sampled )

rho_rrprog_lp <- all


# look the density of lp posterior
d2 <-ggplot(all,aes(x=lp, fill=region)) + geom_density(alpha=0.5) + 
  xlab('log posterior') + 
  ggtitle("Social mixing and RR-progression") + theme_classic()
d2

Lower_baseline <- runn(data_3Regions ,group='L', RR=0)
Lower_Interven <- runn(data_3Regions ,group='L', RR=1)
mid_baseline <- runn(data_3Regions ,group='M', RR=0)
mid_Interven <- runn(data_3Regions ,group='M', RR=1)
upper_baseline <- runn(data_3Regions ,group='U', RR=0)
upper_Interven <- runn(data_3Regions ,group='U', RR=1)


out_Lb <- do.call('rbind', Lower_baseline)
out_LI <- do.call('rbind', Lower_Interven)

out_Mb <- do.call('rbind', mid_baseline)
out_MI <- do.call('rbind', mid_Interven)

out_Ub <- do.call('rbind', upper_baseline)
out_UI <- do.call('rbind', upper_Interven)

out_Lb[, c('Intervention', 'Region'):= list('no', 'L')]
out_LI[, c('Intervention', 'Region'):= list('yes', 'L')]

out_Mb[, c('Intervention', 'Region'):= list('no', 'M')]
out_MI[, c('Intervention', 'Region'):= list('yes', 'M')]

out_Ub[, c('Intervention', 'Region'):= list('no', 'U')]
out_UI[, c('Intervention', 'Region'):= list('yes', 'U')]

DTT <-rbind(out_Lb, out_LI, out_Mb, out_MI, out_Ub,out_UI )
tmpprev<- DTT[, .(prev=mean(prev)), by=.(year, Intervention, Region)]
tmpmfr <- DTT[, .(MFR=mean(MFR)), by=.(year, Intervention, Region)]

save(DTT, file= "Clean/Data/rho_rrprog_posterior_regions_uganda.Rdata") #Uganda
save(DTT, file= "Clean/Data/rho_rrprog_posterior_regions_ethiopia.Rdata") #Uganda

c<-ggplot(tmpprev[year==yr], aes(x=Region, y=prev, fill= Intervention)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal() + coord_flip() + ylab('Prevalence')

d<-ggplot(tmpmfr[year==yr], aes(x=Region, y=MFR, fill= Intervention)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal() + coord_flip() + ylab('M:F ratio')

c/d  + plot_annotation(
  title = 'Intervention impact at different regions of joint posterior\ 
 distribution of rho and rho')
#ggsave('FigInt2_rho_rrprog.png')

#========PAF2======
# pop attributable fraction
prev_2030 <- tmpprev[year==yr]
mfr_2030 <- tmpmfr[year==yr]

PAF_prevalence <- prev_2030 %>%             
  group_by(Region) %>%            
  summarize(PAF = 100*(first(prev)-last(prev))/first(prev))
#summarize(prev2 = 100*(max(prev)-min(prev))/max(prev))
PAF_prevalence

PAF_mfr <- mfr_2030 %>%             
  group_by(Region) %>%            
  summarize(PAF = 100*(first(MFR)-last(MFR))/first(MFR))
#summarize(prev2 = 100*(max(prev)-min(prev))/max(prev))
PAF_mfr


#======PAF intervention=========

#===population attributable fraction====
pop_attr_frac_prev<- DTT[year==yr, .(PAF=100*(first(prev)-last(prev))/first(prev)), by=.(index, Region)]
pop_attr_frac_mfr<- DTT[year==yr, .(PAF=100*(first(MFR)-last(MFR))/first(MFR)), by=.(index, Region)]

#======confidence interval========
PA_prev <-pop_attr_frac_prev[, .(mean= mean(PAF), low= quantile(PAF, 0.025), hi=quantile(PAF, 0.975)), by = Region]
PA_mfr <-pop_attr_frac_mfr[, .(mean= mean(PAF), low= quantile(PAF, 0.025), hi=quantile(PAF, 0.975)), by = Region]


PA_prev$outcome <- 'prevalence'
PA_mfr$outcome  <- 'M:F ratio'
PA_combined =rbind(PA_prev, PA_mfr)

paf2 <- ggplot(PA_combined, aes(x=Region, y=mean, fill=outcome)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=low, ymax=hi), width=.2,position=position_dodge(.9))+
  theme_classic() + scale_fill_brewer(palette="Paired")+
  ylab('Reduction (%)') + coord_flip() + ggtitle('b')


#====rho and rrcdr======


ggplot(DT, aes(x=rho, y=rrcdr)) +
  geom_point(size=1, shape=15)


rmain= 0.05
rmin=  0.02
# construct ellipses
ellL <- Ellipse$new(center = c(quantile(DT$rho,0.25), quantile(DT$rrcdr, 0.25)), rmajor = rmain, rminor = rmin, alpha = 90)
ellpathL <- ellL$path() # ellipse as a path
ellpathL <- rbind(ellpathL, ellpathL[1,]) # the path is not closed; close it


ellM <- Ellipse$new(center = c(quantile(DT$rho,0.50), quantile(DT$rrcdr, 0.50)), rmajor = rmain, rminor = rmin, alpha = 90)
ellpathM <- ellM$path() # ellipse as a path
ellpathM <- rbind(ellpathM, ellpathM[1,]) # the path is not closed; close it


ellU <- Ellipse$new(center = c(quantile(DT$rho,0.75), quantile(DT$rrcdr, 0.75)), rmajor = rmain, rminor = rmin, alpha = 90)
ellpathU <- ellU$path() # ellipse as a path
ellpathU <- rbind(ellpathU, ellpathU[1,]) # the path is not closed; close it


pd <-ggplot() +
  geom_point(aes(x = rho, y = rrcdr), DT) +
  geom_path(aes(x = x, y = y), as.data.frame(ellpathL) , color = "red") +
  geom_path(aes(x = x, y = y), as.data.frame(ellpathM), color = "blue")+
  geom_path(aes(x = x, y = y), as.data.frame(ellpathU), color = "green")

pd

#### extract all points

build <- ggplot_build(pd)$data
points <- build[[1]]
ellipseL <- build[[2]]
ellipseM <- build[[3]]
ellipseU <- build[[4]]


head(points,3)
head(ellipseL,3)


## Capture extracted data into a data frame


dat <- data.frame(
  points[1:2], 
  in.ellipseL = as.logical(point.in.polygon(points$x, points$y, ellipseL$x, ellipseL$y)),
  in.ellipseM = as.logical(point.in.polygon(points$x, points$y, ellipseM$x, ellipseM$y)),
  in.ellipseU = as.logical(point.in.polygon(points$x, points$y, ellipseU$x, ellipseU$y))
)


#### Determine points with in ellipses


dat <- data.table(dat)
dat[in.ellipseL==TRUE, c('in_ellipse', 'region'):= list(TRUE, 'L')]
dat[in.ellipseM==TRUE, c('in_ellipse', 'region'):= list(TRUE, 'M')]
dat[in.ellipseU==TRUE, c('in_ellipse', 'region'):= list(TRUE, 'U')]


# select useful columns
dsa <- dat[, .(x, y,in_ellipse, region )]
rho_rrcdr_sample <-dsa

p3<-ggplot(dsa) +
  geom_point(aes(x, y,color=region),size=1) + 
  stat_ellipse(aes(x, y,color=region)) + xlab('Social mixing') + 
  ylab('rrcdr')
p3
#ggsave(p3,file= 'FigSampRegion3.png')

p1 <- p1+ theme_classic()
p2 <- p2+ theme_classic()
p3 <- p3+ theme_classic()
p1 + p2 + p3 + 
  plot_layout(ncol = 2) + 
  plot_annotation('Sampled regions of joint posterior distributions')

ggsave('Clean/plots/SampleRegion.png')
ggsave('Clean/plots/SampleRegion.pdf')



dsa <- data.table(dsa)
names(DT)

data_in_ellipse <-dsa[in_ellipse==TRUE]
length <- dim(data_in_ellipse)[1]
head(data_in_ellipse)


set.seed(123)

# sampling from main database based on size in the ellipse
index <- sample(1:nrow(DT), nrow(data_in_ellipse))
data_3Regions <-DT[index, ]

# replace values in the main database by values from ellipse
data_3Regions$rho     <- data_in_ellipse$x
data_3Regions$rrcdr    <- data_in_ellipse$y
data_3Regions$region   <- data_in_ellipse$region
data_3Regions <-data.table(data_3Regions)


# summarise all data
DT %>%
  summarize(quant25= quantile(lp,0.25),
            quant50= quantile(lp,0.50),
            quant95=quantile(lp,0.95))
# summarise regional data

data_3Regions[,.(quant25= quantile(lp,0.25),
                 quant50= quantile(lp,0.50), 
                 quant95=quantile(lp,0.95)), by= 'region']


# plot and visualise the densities

all_regions<- DT[, .(rho, rrcdr, lp)]
all_regions[, region:= 'all']

sampled <- data_3Regions[, .(rho, rrcdr,lp, region)]

all<- rbind(all_regions,sampled )
rho_rrcdr_lp <- all

save(rrcdr_rrprog_lp, rho_rrprog_lp,rho_rrcdr_lp, file= 'Clean/Data/logPosterior_ETH.Rdata')
save(rrcdr_rrprog_sample, rho_rrprog_sample,rho_rrcdr_sample, file= 'Clean/Data/sample_posterior_ETH.Rdata')


#=====Density  of lp posterior=====


d3 <-ggplot(all,aes(x=lp, fill=region)) + 
  geom_density(alpha=0.5) + 
  xlab('log posterior') + 
  ggtitle("Social mixing and Effective contact") + theme_classic()
d3

d1+d2+d3 +
  plot_layout(ncol = 2) + 
  plot_annotation('Log-posterior at sampled regions of joint posterior distributions')

ggsave('Clean/plots/FigDensityRegion.png')
ggsave('Clean/plots/FigDensityRegion.pdf')



##call warpper function to run intervention


Lower_baseline <- runn(data_3Regions ,group='L', RR=0)
Lower_Interven <- runn(data_3Regions ,group='L', RR=1)
mid_baseline <- runn(data_3Regions ,group='M', RR=0)
mid_Interven <- runn(data_3Regions ,group='M', RR=1)
upper_baseline <- runn(data_3Regions ,group='U', RR=0)
upper_Interven <- runn(data_3Regions ,group='U', RR=1)


out_Lb <- do.call('rbind', Lower_baseline)
out_LI <- do.call('rbind', Lower_Interven)

out_Mb <- do.call('rbind', mid_baseline)
out_MI <- do.call('rbind', mid_Interven)

out_Ub <- do.call('rbind', upper_baseline)
out_UI <- do.call('rbind', upper_Interven)

out_Lb[, c('Intervention', 'Region'):= list('no', 'L')]
out_LI[, c('Intervention', 'Region'):= list('yes', 'L')]

out_Mb[, c('Intervention', 'Region'):= list('no', 'M')]
out_MI[, c('Intervention', 'Region'):= list('yes', 'M')]

out_Ub[, c('Intervention', 'Region'):= list('no', 'U')]
out_UI[, c('Intervention', 'Region'):= list('yes', 'U')]

DTT <-rbind(out_Lb, out_LI, out_Mb, out_MI, out_Ub,out_UI )
tmpprev<- DTT[, .(prev=mean(prev)), by=.(year, Intervention, Region)]
tmpmfr <- DTT[, .(MFR=mean(MFR)), by=.(year, Intervention, Region)]
save(DTT, file= "Clean/Data/rho_rrcdr_posterior_regions_uganda.Rdata") #Uganda
save(DTT, file= "Clean/Data/rho_rrcdr_posterior_regions_ethiopia.Rdata") #Uganda

yr=2030
e<-ggplot(tmpprev[year==yr], aes(x=Region, y=prev, fill= Intervention)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal() + coord_flip() + ylab('Prevalence')

f<-ggplot(tmpmfr[year==yr], aes(x=Region, y=MFR, fill= Intervention)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_minimal() + coord_flip() + ylab('M:F ratio')

e/f  + plot_annotation(
  title = 'Intervention impact at different regions of joint posterior\ 
  distribution of social mixing and rrcdr')

# Population attributable fraction of intervention
prev_2030 <- tmpprev[year==yr]
mfr_2030 <- tmpmfr[year==yr]

PAF_prevalence <- prev_2030 %>%             
  group_by(Region) %>%            
  summarize(PAF = 100*(first(prev)-last(prev))/first(prev))
#summarize(prev2 = 100*(max(prev)-min(prev))/max(prev))
PAF_prevalence

PAF_mfr <- mfr_2030 %>%             
  group_by(Region) %>%            
  summarize(PAF = 100*(first(MFR)-last(MFR))/first(MFR))
#summarize(prev2 = 100*(max(prev)-min(prev))/max(prev))
PAF_mfr

# label by joint distribution type a- react & rrprog, b= rho and rrprog, c=rho & beta
a <-a +ggtitle('a')
c <-c+ggtitle('b')
e <-e+ggtitle('c')


allplots <-a+c+e+b+ d+f+
  plot_layout(ncol = 3) + 
  plot_annotation('')

# use common legend at bottom
combined <- allplots & theme_classic() + 
  theme(plot.title = element_text(hjust = 0))
combined <- combined & theme(legend.position = "bottom") 
combined + plot_layout(guides = "collect")

ggsave('Clean/Plots/interv_allU.png')
ggsave('Clean/PLots/interv_allU.pdf')

#======PAF intervention=========

#===population attributable fraction====
pop_attr_frac_prev<- DTT[year==yr, .(PAF=100*(first(prev)-last(prev))/first(prev)), by=.(index, Region)]
pop_attr_frac_mfr<- DTT[year==yr, .(PAF=100*(first(MFR)-last(MFR))/first(MFR)), by=.(index, Region)]
pop_attr_frac_prev[PAF<0, PAF:=0]

#======confidence interval========
PA_prev <-pop_attr_frac_prev[, .(mean= mean(PAF), low= quantile(PAF, 0.025), hi=quantile(PAF, 0.975)), by = Region]
PA_mfr <-pop_attr_frac_mfr[, .(mean= mean(PAF),  low= quantile(PAF, 0.025), hi=quantile(PAF, 0.975)), by = Region]

PA_prev$outcome <- 'prevalence'
PA_mfr$outcome  <- 'M:F ratio'
PA_combined =rbind(PA_prev, PA_mfr)

# plot
paf3 <- ggplot(PA_combined, aes(x=Region, y=mean, fill=outcome)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=low, ymax=hi), width=.2,position=position_dodge(.9))+
  theme_classic() + scale_fill_brewer(palette="Paired")+
  ylab('Reduction (%)') + coord_flip() + ggtitle('c')



allfigs <-paf1 + paf2 + paf3 + 
  plot_layout(ncol = 3) + 
  plot_annotation()


# use common legend at bottom
combined <- allfigs & theme_classic() + 
  theme(plot.title = element_text(hjust = 0))
combined <- combined & theme(legend.position = "bottom") 
combined + plot_layout(guides = "collect")

ggsave('Clean/Plots/reduction_intervU.png')

#=====dynamics=====
dts <- DTT[, .(prev=mean(prev)), by =.(year, Intervention, Region)]
ggplot(dts, aes(year, prev, col= Intervention)) +geom_line() +facet_wrap(~Region)


# =======proccess pposterior sample======

load("~/Dropbox/SexTB/Clean/Data/sample_posterior.Rdata")

p1<-ggplot(rho_rrcdr_sample) +
  geom_point(aes(x, y,color=region),size=1) + 
  stat_ellipse(aes(x, y,color=region)) + xlab(expression(rho)) + 
  theme_classic() +
  ylab(expression(pi))

p2<-ggplot(rho_rrprog_sample) +
  geom_point(aes(x, y,color=region),size=1) + 
  stat_ellipse(aes(x, y,color=region)) + xlab(expression(rho)) + 
  theme_classic() +
  ylab(expression(alpha))

p3<-ggplot(rrcdr_rrprog_sample) +
  geom_point(aes(x, y,color=region),size=1) + 
  stat_ellipse(aes(x, y,color=region)) + xlab(expression(pi)) + 
  theme_classic() +
  ylab(expression(alpha))

p1+p2+p3


all_figs <- p1+p2+p3+
  plot_layout(ncol = 3)

# use common legend at bottom
combined <- all_figs & theme_classic() + 
  theme(plot.title = element_text(hjust = 0))
combined <- combined & theme(legend.position = "bottom") 
combined + plot_layout(guides = "collect")


ggsave('Clean/plots/samps_UGA.png', height=3, width=4.5, dpi=300)
