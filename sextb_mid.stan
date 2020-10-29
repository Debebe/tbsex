functions{
  /* function to generate statistic that is zero for each sex if model @ eqm */
  real Qfun(
            real foi,            /* FOI in particular sex */
            real d,              /* prevalence in particular sex */
            real psi,             
            real stab,           
            real fast,
            real react,
            real recov,
            real relapse,
            real mutb,
            real cdr,
            real theta,
            real tau,  
            real mu  
            
            ){
    real A;                       /* quadratic a */
    real B;                       /* quadratic b */
    real C;   
    real ANS;                       /* quadratic a */
    real beta;                     /* coefficients for transm */
    beta = foi / d;
    /* coefficients for quadratic */

  C = beta*mu*(mu + relapse)*(fast*(mu + react) + react*stab)*(mu + tau) - 
      mu*(fast + mu + stab)*(((mu^2) + mutb*react + 
                              mu*(mutb + react + recov))*(mu + relapse)*(mu + tau) + 
                             cdr*(mu + react)*((mu^2) + mu*(relapse + tau) - 
                                               relapse*tau*(-1 + theta)));

  B= -(beta*mu*(beta*fast*(-1 + psi)*(mu + relapse)*(mu + tau) - 
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
                  ))));

  A = (beta^2)*mu*(-1 + psi)*((fast + mu + mutb + recov)*(mu + relapse)*
      (mu + tau) + cdr*((mu^2) + mu*(relapse + tau) - 
      relapse*tau*(-1 + theta) + fast*(mu + relapse + tau*theta)));
  
  ANS = A * d^2 + B * d + C * d^0;
  return ANS;
  }
}

data{
  real targetd;
  real targetdsig;
  
  // real targetd_f;
  // real targetd_fsig;
  // real targetd_m;
  // real targetd_msig;
  
  real targetMF;
  real targetMFsig;
  real targetNf;          
  real targetNfsig;  
  real targetNm;          
  real targetNmsig;   
  real Qtol;
  real pf;         
  real pm;  
  real tau;                     /* treatment duration */
  real theta; 
  /* parameters */
  real mu;                      /* background mortality */
  /* priors parms */
  real psi_a; real psi_b;  /* protection */
  real recov_m; real recov_s;/* selfcure rate */
  real react_m;
  real react_s;/* slow progression rate */
  real fast_m; real fast_s;/* prob TB from fast */
  real stab_m; real stab_s;/* prob TB from fast */
  real relapse_m; real relapse_s;/* relapse */
  real mutb_m; real mutb_s;/* TB mortality */
  real cdr_m; real cdr_s;/* detect */
  real beta_m; real beta_s;/* infection parameter */
  real rho_a; real rho_b;        
  real rrcdr_a; real rrcdr_b;   
  real rrprog_a; real rrprog_b;  
}

parameters{
  real<lower=0,upper=1> df;                       /* TB prevalence,F */
  real<lower=0,upper=1> dm;                       /* TB prevalence,M */
  real<lower=0,upper=1> psi;      /* LTBI protection */
  real<lower=0> recov;                     /* self-cure rate */
  real<lower=0> react;                     /* slow progression rate */
  real<lower=0> fast;      /* fast progression */
  real<lower=0> stab;      /* stabilisation rate */
  real<lower=0> relapse;   /* relapse rate */
  real<lower=0>  mutb;    /* TB mortality rate */
  real<lower=0, upper=1>cdrprop;                   /* detection prop */
  real<lower=0>  beta;           /* effective contact rate */
  
 // real<lower=0,upper=1> rho        
  //real<lower=0,upper=1> rrcdr;  
  real<lower=0> rrcdr;      
  real<lower=0> rrprog;  
  real<lower=0 ,upper=1>contacts;
}

transformed parameters{
  real Qm;                      /* statistic to check eqm */
  real Qf;                      /* statistic to check eqm */
  real rho;
  real<lower=0> foim;
  real<lower=0> foif;
  real<lower=0> d;              /* target for overall prevalence */
  real<lower=0> n;
  real<lower=0> nf;
  real<lower=0> nm;
  real<lower=0> MF;             /* target for MF prevalence ratio */
  real<lower=0> reactf;         /* reactivation females*/  
  real<lower=0> reactm;         /* reactivation males*/                        
  real<lower=0> fastf;          /* reactivation females*/  
  real<lower=0> fastm;          /* reactivation males*/                       
  real<lower=0> cdrf;
  real<lower=0> cdrm;           /* cdr rate in males*/                         
  real<lower=0> cdr;
  real<lower=0> relapsef;
  real<lower=0> relapsem;  
  
  real<lower=0> pnrf; /*DSA*/
  real<lower=0> pnrm; /*DSA*/
  
  // tau = 2.0;                    /* treatment duration */
  // theta = 0.91;

   cdr=cdrprop*(recov + mutb +mu)/(1-cdrprop);
   reactf= 2*react/(1+rrprog);
   reactm= rrprog*reactf;
   fastf= 2*fast/(1+rrprog);
   fastm= rrprog*fastf;
   cdrf = 2*cdr/(1+rrcdr);
   cdrm = rrcdr*cdrf; 
   relapsef = 2*relapse/(1+rrprog);
   relapsem = rrprog*relapsef; 
   //rho = 2*contacts-1;
   rho = (2*contacts-1);

  
  /* equilibrium surfaces */
  foim = ((1-(1-rho)/2) * dm * pm + (1-rho)/2 * df * pf) * beta;
  foif = ((1-(1-rho)/2) * df * pf + (1-rho)/2 * dm * pm) * beta;
  
   Qf = Qfun(foif, df,psi,stab,fastf,reactf,recov,relapsef,mutb,cdrf,theta,tau,mu);
   Qm = Qfun(foim, dm,psi,stab,fastm,reactm,recov,relapsem,mutb,cdrm,theta,tau,mu);

  d = (df + dm)/2;
  MF = dm / (df + 1e-10);
  //n = (df*cdrf + dm*cdrm)/2;
  nf = df*cdrf;
  nm = dm*cdrm;
  n = (nf + nm)/2; /* average not*/
  pnrf= df/nf;
  pnrm= dm/nm;
}

model{
  /* priors */
  psi   ~ beta(psi_a,psi_b);                    /* protection from LTBI*/
  recov ~ lognormal(recov_m,recov_s);         /* selfcure rate */
  react ~ lognormal(react_m,react_s);         /* slow progression rate */
  fast  ~ lognormal(fast_m,fast_s);            /* prob TB from fast */
  stab  ~ lognormal(stab_m,stab_s);            /* prob TB from fast */
  relapse ~ lognormal(relapse_m,relapse_s);   /* relapse */
  mutb ~ lognormal(mutb_m,mutb_s);            /* TB mortality */
  cdrprop ~ beta(cdr_m,cdr_s);                /* detect prop*/
  beta ~ lognormal(beta_m,beta_s);            /* infection parameter*/
  
  //rho    ~ beta(rho_a, rho_a);
  //contacts ~ normal(rho_a, rho_b);
  //contacts ~ beta(rho_a, rho_b);
  contacts ~ lognormal(rho_a, rho_b);
  rrprog ~ lognormal(rrprog_a,rrprog_b);      /* relative progression parameter*/ /*DSA*/
  rrcdr  ~ lognormal(rrcdr_a,rrcdr_b);        /* relative CDR parameter*/        /*DSA*/

  /* calibration and model defn targets */
  target += -(d-targetd)^2 / (2*targetdsig^2);
  // target += -(df-targetd_f)^2 / (2*targetd_fsig^2);
  // target += -(dm-targetd_m)^2 / (2*targetd_msig^2);

  target += -(MF-targetMF)^2 / (2*targetMFsig^2);
  //target += -(n-targetN)^2 / (2*targetNsig^2);      /* notif*/
  target += -(nf-targetNf)^2 / (2*targetNfsig^2);      /* notif*/
  target += -(nm-targetNm)^2 / (2*targetNmsig^2);      /* notif*/
  target += -Qm^2 / (2*Qtol^2);
  target += -Qf^2 / (2*Qtol^2);
}

generated quantities{
  real Qerr;
  /* model departure z score in total */
  Qerr = (sqrt(Qm^2 / (2*Qtol^2)) + sqrt(Qf^2 / (2*Qtol^2)))/2;
}

