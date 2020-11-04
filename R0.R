



#=======initial MF ratio during exponential growth=======
R0 <- function(rho,Rm,Rf){
  am <- (Rm+Rf)/2
  gm <- sqrt(Rf*Rm)
  fac <- (1-(1-rho)/2) + sqrt((1-(1-rho)/2)^2-rho*(gm/am)^2)
  am*fac
}

rho <- seq(from=0,to=1,len=1e2)          #make a load of assortativities
r0data <- data.table(rho=rep(rho,4),ratio=rep(1:4,each=length(rho))) #data
#r0data[,R0val:=R0(rho,ratio*Rf,Rf)]              #evaluate
r0data[,ratio:=factor(ratio)]                     #make a factor
#relative to start


MFinitial <- function(rho,              #assortativity with 1 most assortative
                      phi,              #Rm/Rf
                      Prat=1){          #Prat = Pm/Pf
  rr <- (sqrt(1-16*phi*rho/((1+phi)*(1+rho))^2))
  2*Prat * ((1-rho)/(1+rho)) /(1-phi+rr*(1+phi))
}

r0data[,MFrat:=MFinitial(rho,as.numeric(as.character(ratio))*2,1)]


# =====R0============
singleR0 <- function(L){
  list2env(L,envir=environment())
  T1 <-  (1/(cdr + recov + mu + mutb))    #mean durn 1 disease episode
  Pf <- fast / (fast + mu + stab)                #prob of TB via fast
  Ps <- (stab/(fast+stab+mu)) * (react/(react+mu)) #prob of TB via slow
  Prel <- (cdr * T1) * (theta*tau/(tau+mu)) *
    (relapse/(relapse+mu)) #prob of relapse after tx
  Psc <- (recov * T1) * (react / (react + mu))        #prob of recurrence after self-cure
  P2 <- Prel + Psc                          #prob of additional TB episode
  Ttot <- T1 / (1 - P2)                     #mean TB durn inclusing additional episodes
  R0 <- fracinf * beta * (Pf + Ps) * Ttot   #infectiousness x progression x duration
  list(R0=R0,T=T1, Ttot=Ttot,Pf=Pf,Ps=Ps,Psc=Psc,Prel=Prel)
}


## PJD combining R0
R0combined <- function(rho,Rm,Rf){
  am <- (Rm+Rf)/2
  gm <- sqrt(Rf*Rm)
  fac <- (1-(1-rho)/2) + sqrt((1-(1-rho)/2)^2-rho*(gm/am)^2)
  am*fac
}



eqmQ <- function(L,verbose=FALSE){
  list2env(L,envir=environment())
  
  beta <- beta *  fracinf
  
  C <- beta*mu*(mu + relapse)*(fast*(mu + react) + react*stab)*(mu + tau) - 
    mu*(fast + mu + stab)*(((mu^2) + mutb*react + 
                              mu*(mutb + react + recov))*(mu + relapse)*(mu + tau) + 
                             cdr*(mu + react)*((mu^2) + mu*(relapse + tau) - 
                                                 relapse*tau*(-1 + theta)))
  
  B<- -(beta*mu*(beta*fast*(-1 + psi)*(mu + relapse)*(mu + tau) - 
                   (mu + relapse)*((mu^2)*(-2 + psi) - mutb*react + 
                                     fast*(mu*(-2 + psi) + mutb*(-1 + psi) - react - recov) + 
                                     mu*(mutb*(-2 + psi) - react - 2*recov + psi*recov - stab) - 
                                     mutb*stab - react*stab - recov*stab)*(mu + tau) + 
                   cdr*(-((mu^3)*(-2 + psi)) + react*relapse*stab + 
                          react*relapse*tau + relapse*stab*tau + 
                          (mu^2)*(react - (-2 + psi)*relapse + stab + 2*tau - 
                                    psi*tau) + mu*(stab*tau + react*(relapse + stab + tau) + 
                                                     relapse*(stab + (-2 + psi)*tau*(-1 + theta))) - 
                          react*relapse*tau*theta + react*stab*tau*theta - 
                          relapse*stab*tau*theta + 
                          fast*(-((mu^2)*(-2 + psi)) + 
                                  (-1 + psi)*relapse*tau*(-1 + theta) + 
                                  react*(relapse + tau*theta) + 
                                  mu*(react + 2*relapse - psi*relapse + tau - psi*tau + tau*theta)
                          ))))
  
  A <- (beta^2)*mu*(-1 + psi)*((fast + mu + mutb + recov)*(mu + relapse)*
                                 (mu + tau) + cdr*((mu^2) + mu*(relapse + tau) - 
                                                     relapse*tau*(-1 + theta) + fast*(mu + relapse + tau*theta)))
  ## output
  if(verbose) print(c(A,B,C,d))
  ANS <- A * d^2 + B * d + C * d^0
  ANS
}




