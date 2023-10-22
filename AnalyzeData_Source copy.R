######
###### Source file for NW data analysis
###### All Bayesian models calculate "waic"

###### Basic calculation functions
## pc = probability (risk) in control group
## pa = probability (risk) in active group
OR.cal <- function(pc, pa) {
  OR <- (pa/(1-pa)) / (pc/(1-pc))
  return(OR)
}

RR.cal <- function(pc, pa) {
  RR <- pa/pc
  return(RR)
}

RD.cal <- function(pc, pa) {
  RD <- pa-pc
  return(RD)
}

expit <- function(p) { exp(p)/(1+exp(p)) }
q.025<-function(x){quantile(x, probs = 0.025)}
q.975<-function(x){quantile(x, probs = 0.975)}

#### Save results for naive estimator
## For Data Analysis (DA)
output.pool.DA <- function(mod, measure, const, est, se){
  low <- est - const*se
  up <- est + const*se
  width <- up - low
  tau <- NA
  res <- c(est, se, low, up, width, tau, rep(NA,6))
  names(res) <- paste(measure,".",c("e","se","low","up","width","tau","elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(res)
}
output.pool.RE.DA <- function(mod, measure, const, est, se, tau){
  low <- est - const*se
  up <- est + const*se
  width <- up - low
  tau <- tau
  res <- c(est, se, low, up, width, tau, rep(NA,6))
  names(res) <- paste(measure,".",c("e","se","low","up","width","tau","elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(res)
}

#### Save results for frequentist models (fixed effect)
output.FE.DA <- function(mod,measure,res) {
  if(length(res)>1) {
    out <- c(res$b, res$se, res$ci.lb, res$ci.ub, 
             res$ci.ub-res$ci.lb, NA, rep(NA,13))
  } 
  else if(length(res)<=1) {
    out <- rep(NA,19)
  }
  names(out) <- paste(measure,".",c("e","se","low","up","width","tau","tauc", "taua","corr", "pc", "pa","loc","loa", "elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(out)
}

#### Save results for frequentist models (random effect)
output.RE.DA <- function(mod,measure,res) {
  if(length(res)>1) {
    out <- c(res$b, res$se, res$ci.lb, res$ci.ub, 
             res$ci.ub-res$ci.lb, sqrt(res$tau2), rep(NA,13))
  } 
  else if(length(res)<=1) {
    out <- rep(NA,19)
  }
  names(out) <- paste(measure,".",c("e","se","low","up","width","tau","tauc","taua","corr","pc","pa","loc","loa","elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(out)
}

#### Save results for frequentist model continuous outcome (fixed effect)
output.FE.cont.DA <- function(mod,measure,res) {
  if(length(res)>1) {
    out <- c(res$b, res$se, res$ci.lb, res$ci.ub, 
             res$ci.ub-res$ci.lb, NA, rep(NA,6))
  } 
  else if(length(res)<=1) {
    out <- rep(NA,12)
  }
  names(out) <- paste(measure,".",c("e","se","low","up","width","tau", "elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(out)
}

#### Save results for frequentist model continuous outcome (random effect)
output.RE.cont.DA <- function(mod,measure,res) {
  if(length(res)>1) {
    out <- c(res$b, res$se, res$ci.lb, res$ci.ub, 
             res$ci.ub-res$ci.lb, sqrt(res$tau2), rep(NA,6))
  } 
  else if(length(res)<=1) {
    out <- rep(NA,12)
  }
  names(out) <- paste(measure,".",c("e","se","low","up","width","tau","elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(out)
}



#### Save results for Bayesian models with waic (fixed effect)
####
output.Bayes.cte.DA <- function(mod, measure, mcmc, loglik){
  est <- median(mcmc)
  se <- sd(mcmc)
  low <- q.025(mcmc)
  up <- q.975(mcmc)
  width <- up - low
  tau <- NA
  tauc <- NA
  taua <- NA
  corr <- NA
  pc <- NA
  pa <- NA
  loc <- NA
  loa <- NA
  waic <- waic(loglik)
  waic <- unlist(waic)[c("elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic")]
  res <- c(est, se, low, up, width, tau, tauc, taua, corr, pc, pa,loc,loa, waic)
  names(res) <- paste(measure,".",c("e","se","low","up","width","tau","tauc","taua","corr","pc","pa","loc","loa","elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(res)
}
#### Save rsults for Bayesian models with waic (random effect)
output.Bayes.hte.DA <- function(mod, measure, mcmc, tau.mcmc, loglik){
  est <- median(mcmc)
  se <- sd(mcmc)
  low <- q.025(mcmc)
  up <- q.975(mcmc)
  width <- up - low
  tau <- median(tau.mcmc)
  tauc <- NA
  taua <- NA
  corr <- NA
  pc <- NA
  pa <- NA
  loc <- NA
  loa <- NA
  waic <- waic(loglik)
  waic <- unlist(waic)[c("elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic")]
  res <- c(est, se, low, up, width, tau, tauc, taua,corr, pc, pa,loc,loa, waic)
  names(res) <- paste(measure,".",c("e","se","low","up","width","tau","tauc","taua","corr", "pc","pa","loc","loa","elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(res)
}

output.Bayes.AB.DA <- function(mod, measure, mcmc, tauc.mcmc, taua.mcmc, corr.mcmc, pc.mcmc, pa.mcmc, loc.mcmc, loa.mcmc, loglik){
  est <- median(mcmc)
  se <- sd(mcmc)
  low <- q.025(mcmc)
  up <- q.975(mcmc)
  width <- up - low
  tau <- NA 
  tauc <- median(tauc.mcmc)
  taua <- median(taua.mcmc)
  corr <- median(corr.mcmc)
  pc <- median(pc.mcmc)
  pa <- median(pa.mcmc)
  loc <- median(loc.mcmc)
  loa <- median(loa.mcmc)
  waic <- waic(loglik)
  waic <- unlist(waic)[c("elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic")]
  res <- c(est, se, low, up, width, tau, tauc, taua, corr, pc, pa, loc, loa, waic)
  names(res) <- paste(measure,".",c("e","se","low","up","width","tau", "tauc", "taua","corr", "pc","pa", "loc", "loa", "elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(res)
}

output.Bayes.cte.beta.DA <- function(mod, measure, mcmc, pc.mcmc, pa.mcmc, loglik){
  est <- median(mcmc)
  se <- sd(mcmc)
  low <- q.025(mcmc)
  up <- q.975(mcmc)
  width <- up - low
  tau <- NA
  tauc <- NA
  taua <- NA
  corr <- NA
  pc <- median(pc.mcmc)
  pa <- median(pa.mcmc)
  loc <- NA
  loa <- NA
  waic <- waic(loglik)
  waic <- unlist(waic)[c("elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic")]
  res <- c(est, se, low, up, width, tau, tauc, taua, corr, pc, pa, loc, loa, waic)
  names(res) <- paste(measure,".",c("e","se","low","up","width","tau", "tauc", "taua", "corr", "pc","pa",
                                    "loc", "loa", "elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(res)
}

output.Bayes.beta.DA <- function(mod, measure, mcmc, tauc.mcmc, taua.mcmc, pc.mcmc, pa.mcmc, loglik){
  est <- median(mcmc)
  se <- sd(mcmc)
  low <- q.025(mcmc)
  up <- q.975(mcmc)
  width <- up - low
  tau <- NA
  tauc <- median(tauc.mcmc)
  taua <- median(taua.mcmc)
  corr <- NA
  pc <- median(pc.mcmc)
  pa <- median(pa.mcmc)
  loc <- NA
  loa <- NA
  waic <- waic(loglik)
  waic <- unlist(waic)[c("elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic")]
  res <- c(est, se, low, up, width, tau, tauc, taua, corr, pc, pa, loc, loa, waic)
  names(res) <- paste(measure,".",c("e","se","low","up","width","tau", "tauc", "taua", "corr", "pc", "pa", "loc", "loa", "elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(res)
}

output.Bayes.hte.cont.DA <- function(mod, measure, mcmc, tau.mcmc, delta.mcmc){
  est <- median(mcmc)
  se <- sd(mcmc)
  low <- q.025(mcmc)
  up <- q.975(mcmc)
  width <- up - low
  tau <- median(tau.mcmc)
  waic <- waic(delta.mcmc)
  waic <- unlist(waic)[c("elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic")]
  res <- c(est, se, low, up, width, tau, waic)
  names(res) <- paste(measure,".",c("e","se","low","up","width","tau", "elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(res)
}

output.Bayes.cte.cont.DA <- function(mod, measure, mcmc, delta.mcmc){
  est <- median(mcmc)
  se <- sd(mcmc)
  low <- q.025(mcmc)
  up <- q.975(mcmc)
  width <- up - low
  tau <- NA
  waic <- NA
  waic <- unlist(NA)[c("elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic")]
  res <- c(est, se, low, up, width, tau, waic)
  names(res) <- paste(measure,".",c("e","se","low","up","width","tau", "elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic"),".",mod,sep="")
  return(res)
}

AnalyzeData <- function(DATA){
  ###### Modify data structure
  data <- DATA %>% 
    mutate(study = 1:n(),
           ec_cc1 = ifelse(r2==0 | r1==0, r2+0.5, r2),
           ea_cc1 = ifelse(r2==0 | r1==0, r1+0.5, r1),
           nc_cc1 = ifelse(r2==0 | r1==0, n2+1, n2),
           na_cc1 = ifelse(r2==0 | r1==0, n1+1, n1))
  data.l <- data.frame(rbindlist(list(cbind(data[,c("study","n2","r2")],1), cbind(data[,c("study","n1","r1")],2))))
  colnames(data.l) <- c("study","n","y","t")
  
  ###### Frequentist models ######
  ###### 1. IV, fixed effect (data modification has been implemented in the code; add 0.5 to zero cells)
  #### 5a. LOR
  lor.res.ivfe <- output.FE.DA(mod="ivfe",measure="lor",
                               res=try(rma(ai=ea_cc1, ci=ec_cc1, n1i=na_cc1, n2i=nc_cc1, measure="OR", data=data, method="FE")))
  
  ###### 2. IV, random effect DL (data modification has been implemented in the code; add 0.5 to zero cells)
  #### 6a. LOR
  lor.res.ivre <- output.RE.DA(mod="ivre",measure="lor",
                               res=try(rma(ai=ea_cc1, ci=ec_cc1, n1i=na_cc1, n2i=nc_cc1, measure="OR",data=data,method="DL")))
  ###### Bayesian models ######
  ###### Bayesian models don't use data modification ######
  niter <- 25000
  nburnin <- 5000
  nchains <- 2 
  nthin <- 1
  NS <- nrow(data)
  
  ###### 1. CTE-vague
  #### 1a. LOR 
  init <- list(list(dd=c(NA,0), mu=rep(0,NS)),
               list(dd=c(NA,0.01), mu=rep(0,NS)))
  para <- c("lor","loglik")
  jdata <- list(N=nrow(data.l),NS=NS,s=data.l$study,t=data.l$t,r=data.l$y,n=data.l$n,m1=0,prec1=0.01)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345,progress.bar="text",
         model.file=CTEvague.logit.mod)
  )
  if (1*(length(fit)>1)==0) { 
    lor.res.ctevaguelogit <- rep(NA,19) 
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    lor.mcmc <- fit.mcmc$lor
    loglik.mcmc <- fit.mcmc$loglik
    lor.res.ctevaguelogit <- output.Bayes.cte.DA(mod="ctevaguelogit",measure="lor",mcmc=lor.mcmc, loglik=loglik.mcmc)
  }
  
  ###### 2. HTE-vague
  #### 2a. LOR 
  init <- list(list(lor=0,tau=0.5,mu=rep(0,NS)),
               list(lor=0.1,tau=0.5,mu=rep(0,NS)))
  para <- c("lor","tau","loglik")
  jdata <- list(N=nrow(data.l),NS=NS,s=data.l$study,t=data.l$t,r=data.l$y,n=data.l$n,m1=0,prec1=0.01)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345,progress.bar="text",
         model.file=HTEvague.logit.mod)
  )
  if (1*(length(fit)>1)==0) { 
    lor.res.htevaguelogit <- rep(NA,19) 
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    lor.mcmc <- fit.mcmc$lor
    tau.mcmc <- fit.mcmc$tau
    loglik.mcmc <- fit.mcmc$loglik
    lor.res.htevaguelogit <- output.Bayes.hte.DA(mod="htevaguelogit",measure="lor", mcmc=lor.mcmc, tau.mcmc=tau.mcmc, loglik=loglik.mcmc)
  }
  
  ###### 3. HTE-AB
  #### 3a. LOR 
  init <- list(list(mu=c(0,0), invR=diag(100,2)),
               list(mu=c(0.01, 0.01), invR=diag(100,2)))
  para <- c("lor","loglik", "tau","rho","mu","phat")
  jdata <- list(N=nrow(data.l),NS=NS,s=data.l$study,t=data.l$t,r=data.l$y,n=data.l$n,
                Omega=diag(0.02,2))
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345,progress.bar="text",
         model.file=AB.logit.mod)
  )
  if (1*(length(fit)>1)==0) { 
    lor.res.ablogit <- rep(NA,19) 
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    lor.mcmc <- fit.mcmc$lor
    corr.mcmc <- fit.mcmc$rho
    tauc.mcmc <- fit.mcmc$tau[1]
    taua.mcmc <- fit.mcmc$tau[2]
    pc.mcmc <- fit.mcmc$phat[1]
    pa.mcmc <- fit.mcmc$phat[2]
    loc.mcmc <- fit.mcmc$mu[1]
    loa.mcmc <- fit.mcmc$mu[2]
    loglik.mcmc <- fit.mcmc$loglik
    # AB models' taus are not for variability of lor, so we don't record them.
    lor.res.ablogit<- output.Bayes.AB.DA(mod="ablogit",measure="lor",
                                         mcmc=lor.mcmc,corr.mcmc = corr.mcmc, tauc.mcmc=tauc.mcmc, taua.mcmc=taua.mcmc,
                                         pc.mcmc=pc.mcmc, pa.mcmc=pa.mcmc, loc.mcmc=loc.mcmc,
                                         loa.mcmc=loa.mcmc, loglik=loglik.mcmc)
  }   
  
  ###### 4. CTE-Beta
  #### 4a. LOR
  init <- list(list(pc=0.1, pa=0.1),
               list(pc=0.2, pa=0.2))
  para <- c("lor","loglik","pa","pc")
  jdata <- list(NS=NS,nc=data$n2,na=data$n1,ec=data$r2,ea=data$r1)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345,progress.bar="text",
         model.file=CTEbeta.mod)
  )
  if (1*(length(fit)>1)==0) { 
    lor.res.ctebeta <- rep(NA,19)
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    lor.mcmc <- fit.mcmc$lor
    loglik.mcmc <- fit.mcmc$loglik
    pc.mcmc <- fit.mcmc$pc
    pa.mcmc <-fit.mcmc$pa
    # Beta prior models' taus are not for variability of lor, so we don't record them.
    lor.res.ctebeta<- output.Bayes.cte.beta.DA(mod="ctebeta",measure="lor", mcmc=lor.mcmc, pc.mcmc=pc.mcmc, pa.mcmc=pa.mcmc, loglik=loglik.mcmc)
  }
  
  ###### 5. HTE-Beta
  #### 5a. LOR
  ## sensitive to the prior for Vc and Va!!
  init <- list(list(Uc=0.5, Ua=0.5, Vc=100, Va=100),
               list(Uc=0.51, Ua=0.51, Vc=100.1, Va=100.1))
  para <- c("lor", "loglik","tau_a","tau_c","pa","pc")
  jdata <- list(NS=NS,nc=data$n2,na=data$n1,ec=data$r2,ea=data$r1)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345,progress.bar="text",
         model.file=HTEbeta.mod)
  )
  if (1*(length(fit)>1)==0) { 
    lor.res.htebeta <- rep(NA,19)
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    lor.mcmc <- fit.mcmc$lor
    loglik.mcmc <- fit.mcmc$loglik
    tauc.mcmc <- fit.mcmc$tau_c
    taua.mcmc <- fit.mcmc$tau_a
    pc.mcmc <- fit.mcmc$pc
    pa.mcmc <- fit.mcmc$pa
    # Beta prior models' taus are not for variability of lor, so we don't record them.
    lor.res.htebeta<- output.Bayes.beta.DA(mod="htebeta",measure="lor", mcmc=lor.mcmc, tauc.mcmc, taua.mcmc,pc.mcmc,pa.mcmc, loglik=loglik.mcmc)
  }
  
  ######
  ###### Combine results
  #### 1. LOR
  lor <- data.frame(rbind(lor.res.ivfe, lor.res.ivre, lor.res.ctevaguelogit,
                          lor.res.ctebeta, lor.res.htevaguelogit, lor.res.ablogit, lor.res.htebeta))
  lor$modelname <- c("ivfe", "ivre", "ctevaguelogit",
                     "ctebeta", "htevaguelogit", "ablogit", "htebeta")
  lor$col.group <- factor(c(rep(1,2), rep(2,5))); levels(lor$col.group) <- c("Frequentist", "Bayesian")
  lor$pch.group <- factor(c(1,2,1,1,2,2,2)); levels(lor$pch.group) <- c("CTE model", "HTE model")
 
  colnames(lor) <- c("est","se","low","up","width","tau","tauc","taua","corr", "pc","pa", "loc", "loa",
                                                      "elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic",
                                                      "modelname","col.group","pch.group")

  output <- lor
  return(output)
}

Cal.I_2 <- function(DATA){
  data <- DATA %>% 
    mutate(study = 1:n(),
           ec_cc1 = ifelse(r2==0 | r1==0, r2+0.5, r2),
           ea_cc1 = ifelse(r2==0 | r1==0, r1+0.5, r1),
           nc_cc1 = ifelse(r2==0 | r1==0, n2+1, n2),
           na_cc1 = ifelse(r2==0 | r1==0, n1+1, n1))
  res=try(rma(ai=ea_cc1, ci=ec_cc1, n1i=na_cc1, n2i=nc_cc1, measure="OR", data=data, method="FE"))$I2
  return(res)
}

AnalyzeData.tau50 <- function(DATA){
  data <- DATA %>% 
    mutate(study = 1:n(),
           ec_cc1 = ifelse(r2==0 | r1==0, r2+0.5, r2),
           ea_cc1 = ifelse(r2==0 | r1==0, r1+0.5, r1),
           nc_cc1 = ifelse(r2==0 | r1==0, n2+1, n2),
           na_cc1 = ifelse(r2==0 | r1==0, n1+1, n1))
  data.l <- data.frame(rbindlist( list(cbind(data[,c("study","n2","r2")],1), cbind(data[,c("study","n1","r1")],2))))
  colnames(data.l) <- c("study","n","y","t")
  ###### Bayesian models ######
  ###### Bayesian models don't use data modification ######
  niter <- 25000
  nburnin <- 5000
  nchains <- 2 
  nthin <- 1
  NS <- nrow(data)
  
  ###### HTE-vague
  #### LOR 
  init <- list(list(lor=0,tau=0.5,mu=rep(0,NS)),
               list(lor=0.1,tau=0.5,mu=rep(0,NS)))
  para <- c("lor","tau","loglik")
  jdata <- list(N=nrow(data.l),NS=NS,s=data.l$study,t=data.l$t,r=data.l$y,n=data.l$n,m1=0,prec1=0.01)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345,progress.bar="text",
         model.file=HTEvague.logit.mod.tau50)
  )
  if (1*(length(fit)>1)==0) { 
    lor.res.htevaguelogit <- rep(NA,19) 
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    lor.mcmc <- fit.mcmc$lor
    tau.mcmc <- fit.mcmc$tau
    loglik.mcmc <- fit.mcmc$loglik
    lor.res.htevaguelogit <- output.Bayes.hte.DA(mod="htevaguelogit",measure="lor", mcmc=lor.mcmc, tau.mcmc=tau.mcmc, loglik=loglik.mcmc)
  }
}

AnalyzeData.ab.tau.var.50 <- function(DATA){
  data <- DATA %>% 
    mutate(study = 1:n(),
           ec_cc1 = ifelse(r2==0 | r1==0, r2+0.5, r2),
           ea_cc1 = ifelse(r2==0 | r1==0, r1+0.5, r1),
           nc_cc1 = ifelse(r2==0 | r1==0, n2+1, n2),
           na_cc1 = ifelse(r2==0 | r1==0, n1+1, n1))
  data.l <- data.frame(rbindlist( list(cbind(data[,c("study","n2","r2")],1), cbind(data[,c("study","n1","r1")],2))))
  colnames(data.l) <- c("study","n","y","t")
  ###### Bayesian models ######
  ###### Bayesian models don't use data modification ######
  niter <- 25000
  nburnin <- 5000
  nchains <- 2 
  nthin <- 1
  NS <- nrow(data)
  ###### 3. HTE-AB
  #### 3a. LOR 
  ####Set prior for tau = 10, Omega=diag(0.1,2)
  init <- list(list(mu=c(0,0), invR=diag(100,2)),
               list(mu=c(0.01, 0.01), invR=diag(100,2)))
  para <- c("lor","loglik", "tau","rho","mu","phat")
  jdata <- list(N=nrow(data.l),NS=NS,s=data.l$study,t=data.l$t,r=data.l$y,n=data.l$n,
                Omega=diag(100,2))
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345,progress.bar="text",
         model.file=AB.logit.mod)
  )
  if (1*(length(fit)>1)==0) { 
    lor.res.ablogit <- rep(NA,19) 
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    lor.mcmc <- fit.mcmc$lor
    corr.mcmc <- fit.mcmc$rho
    tauc.mcmc <- fit.mcmc$tau[1]
    taua.mcmc <- fit.mcmc$tau[2]
    pc.mcmc <- fit.mcmc$phat[1]
    pa.mcmc <- fit.mcmc$phat[2]
    loc.mcmc <- fit.mcmc$mu[1]
    loa.mcmc <- fit.mcmc$mu[2]
    loglik.mcmc <- fit.mcmc$loglik
    # AB models' taus are not for variability of lor, so we don't record them.
    lor.res.ablogit<- output.Bayes.AB.DA(mod="ablogit",measure="lor",
                                         mcmc=lor.mcmc,corr.mcmc = corr.mcmc, tauc.mcmc=tauc.mcmc, taua.mcmc=taua.mcmc,
                                         pc.mcmc=pc.mcmc, pa.mcmc=pa.mcmc, loc.mcmc=loc.mcmc,
                                         loa.mcmc=loa.mcmc, loglik=loglik.mcmc)
  }
}

AnalyzeData.cont <- function(DATA, niter, nburnin) {
  data <- DATA
  ###### Frequentist model######
  ####fixed effect model
  md.res.ivfe <- output.FE.cont.DA(mod = "ivre", measure = "md", res = 
                                try(rma(yi = data$y, sei = sqrt(data$se_MD), measure = "MD",
                                        data = data, method = "FE")))
  cohend.res.ivfe <- output.FE.cont.DA(mod = "ivfe", measure = "cohen_d", res = 
                                    try(rma(yi = data$cohen_d, sei = data$se_cohen_d, measure = "SMD",
                                            data = data, method = "FE")))
  hedgeg.res.ivfe <- output.FE.cont.DA(mod = "ivfe", measure = "hedge_g", res = 
                                    try(rma(yi = data$Hedge_g, sei = data$se_hedge_g, measure = "SMD",
                                            data = data, method = "FE")))
  ####random effect model
  md.res.ivre <- output.RE.cont.DA(mod = "ivre", measure = "md", res = 
                                try(rma(yi = data$y, sei = sqrt(data$se_MD), measure = "MD",
                                        data = data, method = "DL")))
  cohend.res.ivre <- output.RE.cont.DA(mod = "ivre", measure = "cohen_d", res = 
                                    try(rma(yi = data$cohen_d, sei = data$se_cohen_d, measure = "SMD",
                                            data = data, method = "DL")))
  hedgeg.res.ivre <- output.RE.cont.DA(mod = "ivre", measure = "hedge_g", res = 
                                    try(rma(yi = data$Hedge_g, sei = data$se_hedge_g, measure = "SMD",
                                            data = data, method = "DL")))
  ###### Bayesian models using JAGS ######
  ###### Bayesian models don't use data modification ######
  ###### MD ######
  niter <- 25000
  nburnin <- 5000
  nchains <- 2 
  nthin <- 1
  NS <- nrow(data)
  y <- data$y
  sigma <- sqrt(data$se_MD)
  #### hte model
  init <- list(list(md=0, tau=0.5),
               list(md=0.0001, tau=0.5))
  para <- c("delta", "md", "tau")
  jdata <- list(NS=NS, y=y, sigma=sigma)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345, progress.bar="text",
         model.file=htecont.mod)
  )
  if (1*(length(fit)>1)==0) { 
    md.res.ctecont <- rep(NA,19) 
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    md.mcmc <- fit.mcmc$md
    tau.mcmc <- fit.mcmc$tau
    delta.mcmc <- fit.mcmc$delta
    md.res.htecont <- output.Bayes.hte.cont.DA(mod="md.htecont", measure="md", mcmc=md.mcmc, tau.mcmc = tau.mcmc, delta.mcmc = delta.mcmc)
  }
  #### cte model
  init <- list(list(md=0),
               list(md=0.0001))
  para <- c("md")
  jdata <- list(NS=NS, y=y, sigma=sigma)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345,progress.bar="text",
         model.file=ctecont.mod)
  )
  if (1*(length(fit)>1)==0) { 
    md.res.ctecont <- rep(NA,19) 
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    md.mcmc <- fit.mcmc$md
    delta.mcmc <- fit.mcmc$delta
    md.res.ctecont <- output.Bayes.cte.cont.DA(mod="md.ctecont",measure="md", mcmc=md.mcmc)
  }
  
  ###### Cohen_d ######
  niter <- 25000
  nburnin <- 5000
  nchains <- 2 
  nthin <- 1
  NS <- nrow(data)
  y <- data$cohen_d
  sigma <- data$se_cohen_d
  #### hte model Cohen_d
  init <- list(list(md=0),
               list(md=0.0001))
  para <- c("delta", "md", "tau")
  jdata <- list(NS=NS, y=y, sigma=sigma)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345, progress.bar="text",
         model.file=htecont.mod)
  )
  if (1*(length(fit)>1)==0) { 
    md.res.ctecont <- rep(NA,19) 
  } else {
      fit.mcmc <- fit$BUGSoutput$sims.list
      md.mcmc <- fit.mcmc$md
      tau.mcmc <- fit.mcmc$tau
      delta.mcmc <- fit.mcmc$delta
      cohend.res.htecont <- output.Bayes.hte.cont.DA(mod="cohend.htecont", measure="smd", mcmc=md.mcmc, tau.mcmc = tau.mcmc, delta.mcmc = delta.mcmc)
  }
  #### cte model Cohen_d
  init <- list(list(md=0),
               list(md=0.0001))
  para <- c("md")
  jdata <- list(NS=NS, y=y, sigma=sigma)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345,progress.bar="text",
         model.file=ctecont.mod)
  )
  if (1*(length(fit)>1)==0) { 
    md.res.ctecont <- rep(NA,19) 
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    md.mcmc <- fit.mcmc$md
    delta.mcmc <- fit.mcmc$delta
    cohend.res.ctecont <- output.Bayes.cte.cont.DA(mod="cohend.ctecont",measure="smd", mcmc=md.mcmc)
  }
  
  ###### Hedge_g ######
  niter <- 25000
  nburnin <- 5000
  nchains <- 2 
  nthin <- 1
  NS <- nrow(data)
  y <- data$Hedge_g
  sigma <- data$se_hedge_g
  #### hte model Hedge_g
  init <- list(list(md=0),
               list(md=0.0001))
  para <- c("delta", "md", "tau")
  jdata <- list(NS=NS, y=y, sigma=sigma)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345, progress.bar="text",
         model.file=htecont.mod)
  )
  if (1*(length(fit)>1)==0) { 
    md.res.ctecont <- rep(NA,19) 
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    md.mcmc <- fit.mcmc$md
    tau.mcmc <- fit.mcmc$tau
    delta.mcmc <- fit.mcmc$delta
    hedgeg.res.htecont <- output.Bayes.hte.cont.DA(mod="hedgeg.htecont", measure="smd", mcmc=md.mcmc, tau.mcmc = tau.mcmc, delta.mcmc = delta.mcmc)
  }
  #### cte model Hedge_g
  init <- list(list(md=0),
               list(md=0.0001))
  para <- c("md")
  jdata <- list(NS=NS, y=y, sigma=sigma)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345,progress.bar="text",
         model.file=ctecont.mod)
  )
  if (1*(length(fit)>1)==0) { 
    md.res.ctecont <- rep(NA,19) 
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    md.mcmc <- fit.mcmc$md
    delta.mcmc <- fit.mcmc$delta
    hedgeg.res.ctecont <- output.Bayes.cte.cont.DA(mod="hedgeg.ctecont",measure="smd", mcmc=md.mcmc)
  }
  
  ####Output
  md <- data.frame(rbind(md.res.ivfe, md.res.ivre, md.res.htecont, md.res.ctecont, 
                          cohend.res.ivfe, cohend.res.ivre, cohend.res.htecont, cohend.res.ctecont,
                          hedgeg.res.ivfe, hedgeg.res.ivre, hedgeg.res.htecont, hedgeg.res.ctecont))
  
  md$modelname <- rep(c("ivfe", "ivre", "htebayesian", "ctebayesian"), 3)
  
  md$col.group <- factor(c(rep(1,2), rep(2,2))); levels(md$col.group) <- c("Frequentist", 
                                                                             "Bayesian")
  colnames(md)<- c("est","se","low","up","width","tau",
                    "elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic",
                    "modelname","col.group")
  return(md)
}

AnalyzeData.MD.tau50 <- function(DATA, niter, nburnin) {
  data <- DATA
  ###### Frequentist model######
  ####fixed effect model
  md.res.ivfe <- output.FE.cont.DA(mod = "ivre", measure = "md", res = 
                                     try(rma(yi = data$y, sei = sqrt(data$s2), measure = "MD",
                                             data = data, method = "FE")))
  ####random effect model
  md.res.ivre <- output.RE.cont.DA(mod = "ivre", measure = "md", res = 
                                     try(rma(yi = data$y, sei = sqrt(data$s2), measure = "MD",
                                             data = data, method = "DL")))
  ###### Bayesian models using JAGS ######
  ###### Bayesian models don't use data modification ######
  ###### MD ######
  niter <- 25000
  nburnin <- 5000
  nchains <- 2 
  nthin <- 1
  NS <- nrow(data)
  y <- data$y
  sigma <- sqrt(data$s2)
  #### hte model
  init <- list(list(md=0, tau=0.5),
               list(md=0.0001, tau=0.5))
  para <- c("delta", "md", "tau")
  jdata <- list(NS=NS, y=y, sigma=sigma)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345, progress.bar="text",
         model.file=htecont.mod)
  )
  if (1*(length(fit)>1)==0) { 
    md.res.htecont <- rep(NA,19) 
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    md.mcmc <- fit.mcmc$md
    tau.mcmc <- fit.mcmc$tau
    delta.mcmc <- fit.mcmc$delta
    md.res.htecont <- output.Bayes.hte.cont.DA(mod="md.htecont", measure="md", mcmc=md.mcmc, tau.mcmc = tau.mcmc, delta.mcmc = delta.mcmc)
  }
  #### cte model
  init <- list(list(md=0),
               list(md=0.0001))
  para <- c("md")
  jdata <- list(NS=NS, y=y, sigma=sigma)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345,progress.bar="text",
         model.file=ctecont.mod)
  )
  if (1*(length(fit)>1)==0) { 
    md.res.ctecont <- rep(NA,19) 
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    md.mcmc <- fit.mcmc$md
    delta.mcmc <- fit.mcmc$delta
    md.res.ctecont <- output.Bayes.cte.cont.DA(mod="md.ctecont",measure="md", mcmc=md.mcmc)
  }
  ####Output
  md <- data.frame(rbind(md.res.ivfe, md.res.ivre, md.res.htecont, md.res.ctecont))
  
  md$modelname <- c("ivfe", "ivre", "htebayesian", "ctebayesian")
  
  md$col.group <- factor(c(rep(1,2), rep(2,2))); levels(md$col.group) <- c("Frequentist", 
                                                                           "Bayesian")
  colnames(md)<- c("est","se","low","up","width","tau",
                   "elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic",
                   "modelname","col.group")
  md$outcome <- rep('MD', 4)
  return(md)
}

AnalyzeData.MD.tau100 <- function(DATA, niter, nburnin) {
  data <- DATA
  ###### Bayesian models using JAGS ######
  ###### Bayesian models don't use data modification ######
  ###### MD ######
  niter <- 25000
  nburnin <- 5000
  nchains <- 2 
  nthin <- 1
  NS <- nrow(data)
  y <- data$y
  sigma <- sqrt(data$s2)
  #### hte model
  init <- list(list(md=0, tau=0.5),
               list(md=0.0001, tau=0.5))
  para <- c("delta", "md", "tau")
  jdata <- list(NS=NS, y=y, sigma=sigma)
  fit <- NULL
  fit <- try(
    jags(data=jdata, inits=init, para,
         n.iter=niter, n.burnin=nburnin, n.chains=nchains, n.thin=nthin,
         DIC=TRUE, jags.seed=12345, progress.bar="text",
         model.file=htecont_tau100.mod)
  )
  if (1*(length(fit)>1)==0) { 
    md.res.htecont.tau100 <- rep(NA,19) 
  } else {
    fit.mcmc <- fit$BUGSoutput$sims.list
    md.mcmc <- fit.mcmc$md
    tau.mcmc <- fit.mcmc$tau
    delta.mcmc <- fit.mcmc$delta
    md.res.htecont.tau100 <- output.Bayes.hte.cont.DA(mod="md.hteconttau100", measure="md", mcmc=md.mcmc, tau.mcmc = tau.mcmc, delta.mcmc = delta.mcmc)
  }
  ####Output
  md <- as.data.frame(t(data.frame(md.res.htecont.tau100)))
  md$modelname <- 'htebayesian.tau100'
  md$col.group <- 'bayesian'
  colnames(md)<- c("est","se","low","up","width","tau",
                   "elpd_waic", "p_waic", "waic", "se_elpd_waic", "se_p_waic", "se_waic",
                   "modelname","col.group")
  md$outcome <- 'MD'
  return(md)
}
  