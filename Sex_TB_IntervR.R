Sex_TB_Interv <- odin::odin({
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
  
  #------------------------------
  RRcdr_interv[]     <- user()               
  interv_year[]      <- user() 
  dim(interv_year)   <- user()
  length_interv_year <- length(interv_year)
  dim(RRcdr_interv)  <- length_interv_year
  RRcdr    <-  interpolate(interv_year, RRcdr_interv,'linear')
  
  #------------------------------

  cdr_mid_rate <-cdr_mid*(recov + mutb +mu)/(1-cdr_mid);
  
  cdr[1]   <- 2*cdr_mid_rate*RRcdr/(1+RRcdr)
  cdr[2]   <- 2*cdr_mid_rate/(1+RRcdr);
  dim(cdr) <- n
  
  dim(react) <-n
  dim(fast)  <-n
  dim(relapse) <-1
  RRprog     <- user()
  
  
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
  
  
  #Beta[,]    <- user()              
  beta        <- user()
  
  recov       <- user()             
  mu          <- user()              
  mutb        <- user()             
  stab        <- user()
  psi         <- user()             
  #relapse     <- user()             
  
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

