##CODE FILE FOR "Do 'complex' financial models really lead to complex dynamics? Agent-based models and multifractality".
##Kukacka, J., Kristoufek, L. (2020). Journal of Economic Dynamics and Control, doi: 10.1016/j.jedc.2020.103855
##CODES FOR MF-DFA, MF-DMA AND GHE ARE GIVEN HERE AS WELL AS TWO ILLUSTRATIVE EXAMPLES OF MULTIFRACTAL SERIES.
##QUESTIONS CAN BE DIRECTED TO LK@FSV.CUNI.CZ.
##VERSION: FEBRUARY 14, 2020.
##ATTACHED FILES: mmc2.csv (for BTC), mmc3.csv (for SP).

#######################################################################################
#######################################################################################
#Functions used in the manuscript, needed for estimation of the multifractal spectrum.#
#######################################################################################
#######################################################################################

#fluctuation function for fast linear trend fits for MF-DFA
fluctuation<-function(x){
  t<-1:length(x)
  return(mean((lm(x~t)$residuals)^2))
}

#MF_DFA function returns a matrix of generalized Hurst exponents H(q) for a given range of qs based on MF-DFA
#Function arguments:
# - x - time series of (log-)returns
# - smin, smax - minimum and maximum scales to consider for the scaling law (Eq. 6 of the manuscript)
# - sstep - step length for scales between smin and smax
# - qmin, qmax - minim and maximum moments q for estimation of the multifractal spectrum
# - qstep - step length for moments between qmin and qmax
MF_DFA<-function(x,smin,smax,sstep,qmin,qmax,qstep){
  xx<-cumsum(x-mean(x))
  qcount<-(qmax-qmin)/qstep+1
  F2_s<-matrix(rep(NA,(qcount+1)*((smax-smin)/sstep+1)),ncol=qcount+1)
  H_q<-matrix(rep(NA,2*qcount),ncol=2)
  counter_s<-1
  for(s in seq(smin,smax,by=sstep)){
    cut<-length(xx)-floor(length(xx)/s)*s
    xmat_L_PROBLEM<-FALSE
    xmat_R_PROBLEM<-FALSE
    if(cut==0){
      xmat_L<-matrix(xx,nrow=s,ncol=floor(length(xx)/s))
    } else {
      xmat_L<-matrix(xx[-c((length(xx)-cut+1):length(xx))],nrow=s,ncol=floor(length(xx)/s))
    }
    if (is.element(Inf, xmat_L) || is.element(-Inf, xmat_L) || is.element(NaN, xmat_L)){
      xmat_L_PROBLEM<-TRUE
    } else {
      F2_L<-apply(xmat_L,2,fluctuation)
    }
    if(cut==0){
      F2_R<-NULL
    } else {
      xmat_R<-matrix(xx[-c(1:cut)],nrow=s,ncol=floor(length(xx)/s))
      if (is.element(Inf, xmat_R) || is.element(-Inf, xmat_R) || is.element(NaN, xmat_R)){
        xmat_R_PROBLEM<-TRUE
      } else {
        F2_R<-apply(xmat_R,2,fluctuation)
      }
    }
    if (xmat_L_PROBLEM || xmat_R_PROBLEM){
      F2_s[counter_s,1]<-s
      counter_q<-2
      for(q in seq(qmin,qmax,by=qstep)){
        F2_s[counter_s,counter_q]<-NaN
        counter_q<-counter_q+1
      }
    } else {
      F2<-c(F2_L,F2_R)
      F2_s[counter_s,1]<-s
      counter_q<-2
      for(q in seq(qmin,qmax,by=qstep)){
        if(q==0){
          F2_s[counter_s,counter_q]<-exp(mean(log(F2))/2)
        } else {
          F2_s[counter_s,counter_q]<-mean(F2^(q/2))^(1/q)  
        }
        counter_q<-counter_q+1
      }
    }
    counter_s<-counter_s+1
  }
  counter_q<-1
  for(q in seq(qmin,qmax,by=qstep)){
    H_q[counter_q,1]<-q
    if (is.element(NaN, F2_s)) {
      H_q[counter_q,2]<-NaN
    } else {
      H_q[counter_q,2]<-lm(log(F2_s[,(counter_q+1)])~log(F2_s[,1]))$coefficients[2]
    }
    counter_q<-counter_q+1
  }
  return(H_q)
}

#MF_DFA_strength function returns $\Delta H$ and $\Delta \alpha$ according to Eq. 4 and below of the manuscript.
#Whole spectrum can be obtained and slight tweaking of the function (by uncommenting)
#The parameters hold from the MF_DFA function above
MF_DFA_strength<-function(x,smin,smax,sstep,qmin,qmax,qstep){
  H_q<-MF_DFA(x,smin,smax,sstep,qmin,qmax,qstep)
  delta_H<-max(H_q[,2])-min(H_q[,2])
  qcount<-(qmax-qmin)/qstep+1
  alpha_f<-H_q[-1,]
  alpha_f[,1]<-H_q[-1,][,2]+H_q[-1,][,1]*(H_q[-1,][,2]-H_q[-qcount,][,2])
  alpha_f[,2]<-1+H_q[-1,][,1]*(alpha_f[,1]-H_q[-1,][,2])
  delta_alpha<-max(alpha_f[,1])-min(alpha_f[,1])
  #plot(alpha_f[,1],alpha_f[,2])
  #return(matrix(c(alpha_f[,1],alpha_f[,2]),ncol=2))
  return(c(delta_H,delta_alpha))
}

#MF_DMA function returns a matrix of generalized Hurst exponents H(q) for a given range of qs based on MF-DMA
#Function arguments:
# - x - time series of (log-)returns
# - smin, smax - minimum and maximum scales to consider for the scaling law (below Eq. 6 of the manuscript);
    ### as the centered moving average is used, smin and smax should be odd integers
# - sstep - step length for scales between smin and smax
# - qmin, qmax - minim and maximum moments q for estimation of the multifractal spectrum
# - qstep - step length for moments between qmin and qmax
MF_DMA<-function(x,smin,smax,sstep,qmin,qmax,qstep){
  xx<-cumsum(x-mean(x))
  qcount<-(qmax-qmin)/qstep+1
  F2_s<-matrix(rep(NA,(qcount+1)*((smax-smin)/sstep+1)),ncol=qcount+1)
  H_q<-matrix(rep(NA,2*qcount),ncol=2)
  counter_s<-1
  for(s in seq(smin,smax,by=sstep)){
    ma<-c(rep(1,s))/s
    X_MA<-filter(xx,ma)
    X_MA<-X_MA[is.na(X_MA)==FALSE]
    cut<-length(X_MA)-floor(length(X_MA)/s)*s
    xmat_L_PROBLEM<-FALSE
    xmat_R_PROBLEM<-FALSE
    if(cut==0){
      xmat_L<-matrix(X_MA,nrow=s,ncol=floor(length(X_MA)/s))
    } else {
      xmat_L<-matrix(X_MA[-c((length(X_MA)-cut+1):length(X_MA))],nrow=s,ncol=floor(length(X_MA)/s))
    }
    if (is.element(Inf, xmat_L) || is.element(-Inf, xmat_L) || is.element(NaN, xmat_L)){
      xmat_L_PROBLEM<-TRUE
    } else {
      F2_L<-apply(xmat_L,2,fluctuation)
    }
    if(cut==0){
      F2_R<-NULL
    } else {
      xmat_R<-matrix(X_MA[-c(1:cut)],nrow=s,ncol=floor(length(X_MA)/s))
      if (is.element(Inf, xmat_R) || is.element(-Inf, xmat_R) || is.element(NaN, xmat_R)){
        xmat_R_PROBLEM<-TRUE
      } else {
        F2_R<-apply(xmat_R,2,fluctuation)
      }
    }
    if (xmat_L_PROBLEM || xmat_R_PROBLEM){
      F2_s[counter_s,1]<-s
      counter_q<-2
      for(q in seq(qmin,qmax,by=qstep)){
        F2_s[counter_s,counter_q]<-NaN
        counter_q<-counter_q+1
      }
    } else {
      F2<-c(F2_L,F2_R)
      F2_s[counter_s,1]<-s
      counter_q<-2
      for(q in seq(qmin,qmax,by=qstep)){
        if(q==0){
          F2_s[counter_s,counter_q]<-exp(mean(log(F2))/2)
        } else {
          F2_s[counter_s,counter_q]<-mean(F2^(q/2))^(1/q)  
        }
        counter_q<-counter_q+1
      }
    }
    counter_s<-counter_s+1
  }
  counter_q<-1
  for(q in seq(qmin,qmax,by=qstep)){
    H_q[counter_q,1]<-q
    if (is.element(NaN, F2_s)) {
      H_q[counter_q,2]<-NaN
    } else {
      H_q[counter_q,2]<-lm(log(F2_s[,(counter_q+1)])~log(F2_s[,1]))$coefficients[2]
    }
    counter_q<-counter_q+1
  }
  return(H_q)
}

#MF_DMA_strength function returns $\Delta H$ and $\Delta \alpha$ according to Eq. 4 and below of the manuscript.
#Whole spectrum can be obtained and slight tweaking of the function (by uncommenting)
#The parameters hold from the MF_DMA function above
MF_DMA_strength<-function(x,smin,smax,sstep,qmin,qmax,qstep){
  H_q<-MF_DMA(x,smin,smax,sstep,qmin,qmax,qstep)
  delta_H<-max(H_q[,2])-min(H_q[,2])
  qcount<-(qmax-qmin)/qstep+1
  alpha_f<-H_q[-1,]
  alpha_f[,1]<-H_q[-1,][,2]+H_q[-1,][,1]*(H_q[-1,][,2]-H_q[-qcount,][,2])
  alpha_f[,2]<-1+H_q[-1,][,1]*(alpha_f[,1]-H_q[-1,][,2])
  delta_alpha<-max(alpha_f[,1])-min(alpha_f[,1])
  #plot(alpha_f[,1],alpha_f[,2])
  #return(matrix(c(alpha_f[,1],alpha_f[,2]),ncol=2))
  return(c(delta_H,delta_alpha))
}

#GHE function returns a matrix of generalized Hurst exponents H(q) for a given range of qs based on GHE (Eq. 7 of the manuscript)
#Function arguments:
# - x - time series of (log-)returns
# - tmin - minimum time scale (Eq. 8 of the manuscript);
# - tmax1, tmax2 - a set of maximum time scales used for jackknifing the estimation of H(q) (below Eq. 8 of the manuscript)
# - qmin, qmax - minim and maximum moments q for estimation of the multifractal spectrum
# - qstep - step length for moments between qmin and qmax
GHE<-function(x,tmin,tmax1,tmax2,qmin,qmax,qstep){
  x<-cumsum(x)
  qs<-round((qmax-qmin)/qstep,0)+1
  qH<-matrix(rep(NA,2*qs),ncol=2)
  counter_q<-1
  for(q in seq(qmin,qmax,by=qstep)){
    qH[counter_q,1]<-q
    t_mat<-matrix(rep(NA,2*(tmax2-tmax1+1)),ncol=2)
    counter_tmax<-1
    for(tmax in tmax1:tmax2){
      tau_mat<-matrix(rep(NA,2*(tmax-tmin+1)),ncol=2)
      counter_tau<-1
      for(tau in tmin:tmax){
        tau_mat[counter_tau,1]<-tau
        tau_mat[counter_tau,2]<-mean(abs(diff(x,lag=tau,differences=1))^q)
        counter_tau<-counter_tau+1
      }
      t_mat[counter_tmax,1]<-tmax
      if (is.element(Inf, tau_mat) || is.element(-Inf, tau_mat) || is.element(NaN, tau_mat)){
        t_mat[counter_tmax,2]<-NaN
      } else {
        t_mat[counter_tmax,2]<-lm(log(tau_mat[,2])~log(tau_mat[,1]))$coefficients[2]/q
      }
      counter_tmax<-counter_tmax+1
    }
    qH[counter_q,2]<-mean(t_mat[,2],na.rm=TRUE)
    counter_q<-counter_q+1
  }
  return(qH)
}

#GHE_strength function returns $\Delta H$ and $\Delta \alpha$ according to Eq. 4 and below of the manuscript.
#Whole spectrum can be obtained and slight tweaking of the function (by uncommenting)
#The parameters hold from the GHE function above
GHE_strength<-function(x,tmin,tmax1,tmax2,qmin,qmax,qstep){
  H_q<-GHE(x,tmin,tmax1,tmax2,qmin,qmax,qstep)
  delta_H<-max(H_q[,2])-min(H_q[,2])
  qcount<-(qmax-qmin)/qstep+1
  alpha_f<-H_q[-1,]
  alpha_f[,1]<-H_q[-1,][,2]+H_q[-1,][,1]*(H_q[-1,][,2]-H_q[-qcount,][,2])
  alpha_f[,2]<-1+H_q[-1,][,1]*(alpha_f[,1]-H_q[-1,][,2])
  delta_alpha<-max(alpha_f[,1])-min(alpha_f[,1])
  #plot(alpha_f[,1],alpha_f[,2])
  #return(matrix(c(alpha_f[,1],alpha_f[,2]),ncol=2))
  return(c(delta_H,delta_alpha))
}

############################
############################
#Two illustrative examples.#
############################
############################

#Illustrative example of a multifractal series. Binomial multifractal model (Feder, J. 1988. Fractals, Plenum Press, NY, USA) is used.
#Function MFsim(N,a) from the MFDFA package is used where N is the number of observations and a is the exponent describing the strength of multifractality.
require(MFDFA)
set.seed(1)
a<-0.7
N<-2^14
MF_BN<-MFsim(N,a)
estimates<-MF_DFA(MF_BN,10,1000,10,-3,3,0.1)
#Plot the estimates as circles and theoretical values as a solid curve. There is no theoretical value for q=0.
plot(estimates[,1],estimates[,2],
     xlab="q", ylab="H(q)")
curve(1/x-(log(a^x+(1-a)^x))/(x*log(2)),-3,3,add=T)

#Illustrative example of a multifractal series. Pareto distributed i.i.d. series with parameter $alpha$. 
#Kantelhardt et al. (2002, Physica A 316:87-114) provides more details about the MF properties of power-law distributed series.
require(EnvStats)
set.seed(1)
alpha<-1
MF_Pareto<-rpareto(10000,0.1,alpha)
estimates<-MF_DFA(MF_Pareto,10,1000,10,-3,3,0.1)
#Plot the estimates as circles and theoretical values as a solid curve.
plot(estimates[,1],estimates[,2],
     xlab="q", ylab="H(q)")
curve(1/alpha+0*x,-3,alpha,add=T) #H(q)=1/alpha for q<alpha
curve(1/x,alpha,3,add=T) #H(q)=1/qfor q>alpha

##################################################################################################################################
##################################################################################################################################
#Empirical example as given in the manuscript, here for BTC. Can be used for S&P500 as well. This will be shown in several steps.#
##################################################################################################################################
##################################################################################################################################

#STEP1: Load the data, get the returns.
set.seed(1)

price <- read.csv("mmc2.csv", sep=";") #for BTC
#price <- read.csv("mmc3.csv", sep=";") #for SP
ret<-log(price$Close)-log(price$Open)

#Step2: Setting up empty matrices for eventual results, specifying the number of shuffling repetitions.
rep<-1000 #Number of shufling repetitions
results_H<-matrix(NA,nrow=3,ncol=4)
results_alpha<-matrix(NA,nrow=3,ncol=4)
colnames(results_H)<-c("delta_H","delta_H_shuffle_mean","q005_shuffle","q095_shuffle")
rownames(results_H)<-c("MF-DFA","MF-DMA","GHE")
colnames(results_alpha)<-c("delta_alpha","delta_alpha_shuffle_mean","q005_shuffle","q095_shuffle")
rownames(results_alpha)<-c("MF-DFA","MF-DMA","GHE")
estimates<-matrix(NA,nrow=3,ncol=rep)

#Step3: Estimation on the original data
MF_DFA_original<-MF_DFA_strength(ret,10,200,10,-3,3,0.1)
results_H[1,1]<-MF_DFA_original[1]
results_alpha[1,1]<-MF_DFA_original[2]
MF_DMA_original<-MF_DMA_strength(ret,11,201,10,-3,3,0.1)
results_H[2,1]<-MF_DMA_original[1]
results_alpha[2,1]<-MF_DMA_original[2]
MF_GHE_original<-GHE_strength(ret,1,5,20,0.1,3,0.1)
results_H[3,1]<-MF_GHE_original[1]
results_alpha[3,1]<-MF_GHE_original[2]

#Step4: Generating shuffled series
shuffled<-sapply(1:rep,function(x) sample(ret,length(ret),replace=T))

#Step5: Estimation on the shuffled series in parallel.
require(parallel)
no_cores<-parallel::detectCores()
cl<-parallel::makeCluster(no_cores-1,type = "FORK")

#MF-DFA
MF_DFA_shuffled<-parApply(cl,shuffled,2,function(x) MF_DFA_strength(x,10,200,10,-3,3,0.1))
results_H[1,2]<-mean(MF_DFA_shuffled[1,])
results_H[1,3]<-quantile(MF_DFA_shuffled[1,],0.05)
results_H[1,4]<-quantile(MF_DFA_shuffled[1,],0.95)
results_alpha[1,2]<-mean(MF_DFA_shuffled[2,])
results_alpha[1,3]<-quantile(MF_DFA_shuffled[2,],0.05)
results_alpha[1,4]<-quantile(MF_DFA_shuffled[2,],0.95)

#MF-DMA
MF_DMA_shuffled<-parApply(cl,shuffled,2,function(x) MF_DMA_strength(x,11,201,10,-3,3,0.1))
results_H[2,2]<-mean(MF_DMA_shuffled[1,])
results_H[2,3]<-quantile(MF_DMA_shuffled[1,],0.05)
results_H[2,4]<-quantile(MF_DMA_shuffled[1,],0.95)
results_alpha[2,2]<-mean(MF_DMA_shuffled[2,])
results_alpha[2,3]<-quantile(MF_DMA_shuffled[2,],0.05)
results_alpha[2,4]<-quantile(MF_DMA_shuffled[2,],0.95)

#GHE
GHE_shuffled<-parApply(cl,shuffled,2,function(x) GHE_strength(x,1,5,20,0.1,3,0.1))
results_H[3,2]<-mean(GHE_shuffled[1,])
results_H[3,3]<-quantile(GHE_shuffled[1,],0.05)
results_H[3,4]<-quantile(GHE_shuffled[1,],0.95)
results_alpha[3,2]<-mean(GHE_shuffled[2,])
results_alpha[3,3]<-quantile(GHE_shuffled[2,],0.05)
results_alpha[3,4]<-quantile(GHE_shuffled[2,],0.95)

stopCluster(cl)

#Step6: Getting results
results_H
results_alpha

#Step7: H(q) spectrum
#ret<-sample(ret,length(ret),replace=T) #For shuffled series.
MF_DFA(ret,10,200,10,-3,3,0.1)
MF_DMA(ret,11,201,10,-3,3,0.1)
GHE(ret,1,5,20,0.1,3,0.1)

##################################################################################################################################
##################################################################################################################################
#############A specimen of the Monte Carlo simulation study utilized in the manuscript to produce Fig. 4-type chart.##############
###Computational setup trivialized for personal computers. The original setup for a cluster computation is given in Section 4.1.##
##################################################################################################################################
##################################################################################################################################

#Step1: Data-generating process. Either pick one of the two provided here or specify your own DGP.
#By default, there are two options corresponding to theoretical examples above: 
#Pareto distributed i.i.d. series with $alpha=1.5$ (shape) and binomial MF model with $a=0.7$.
#Can be rewritten into any DGP as needed, $obs$ needs to be adjustable for the number of observations.

#For BNMF model, use:
require(MFDFA)
#For Pareto distributed i.i.d., use:
#require(EnvStats)
DGP<-function(obs=1000){
  MF_BN<-MFDFA::MFsim(obs,0.7) #Note: the DGP is deterministic, ie no stochastic component is added, ie estimates over all repetitions are the same (SD = 0). Not for the shuffled data, of course.
  #rpareto(n=obs,0.1,shape=1.5) #Note: multifractality here is only due to distributional properties that remain unchanged after shuffling.
}

#Step2: Simulation and MF estimation function to be repeated. MF_DMA utilized here. Can be easily rewritten for MF_DFA or GHE.
MF_simulation<-function(obs=1000){
  results<-c(NA,NA)
  ret<-DGP(obs)
  ret_shuffle<-sample(ret,replace=T)
  results[1]<-MF_DMA_strength(ret,11,obs/5+1,10,-3,3,0.1)[2]
  results[2]<-MF_DMA_strength(ret_shuffle,11,obs/5+1,10,-3,3,0.1)[2] #Standard setting of MF_DMA as in the paper.
  return(results)
}

#Step3: Simulation parallelized.
#Number of Monte Carlo repetitions $rep$ needs to be set.
#Set a vector of time series lengths $obs_vec$ to be simulated. Should be divisible by 10.
require(parallel)
no_cores<-parallel::detectCores()
#For Macs/Linux-based:
cl<-parallel::makeCluster(no_cores-1,type = "FORK") #Sets the number of used cores to max minus one. Can be adjusted.
#For Windows:
#cl<-parallel::makeCluster(no_cores-1,type = "PSOCK")
MF_simulation_par<-function(rep=100,obs_vec=c(250,500,1000,5000)){
  results_mat<-matrix(NA,nrow=length(obs_vec),ncol=6)
  colnames(results_mat)<-c("obs","mean","median","SD","q0.05","q0.95")
  results_mat_shuffle<-matrix(NA,nrow=length(obs_vec),ncol=6)
  colnames(results_mat_shuffle)<-c("obs","mean","median","SD","q0.05","q0.95")
  for(i in 1:length(obs_vec)){
    series<-sapply(1:rep,function(x) DGP(obs=obs_vec[i]))
    alphas<-parApply(cl,series,2,function(x) MF_simulation(obs=obs_vec[i]))
    results_mat[i,1]<-obs_vec[i]
    results_mat[i,2]<-mean(alphas[1,])
    results_mat[i,3]<-median(alphas[1,])
    results_mat[i,4]<-sd(alphas[1,])
    results_mat[i,5]<-quantile(alphas[1,],0.05)
    results_mat[i,6]<-quantile(alphas[1,],0.95)
    results_mat_shuffle[i,1]<-obs_vec[i]
    results_mat_shuffle[i,2]<-mean(alphas[2,])
    results_mat_shuffle[i,3]<-median(alphas[2,])
    results_mat_shuffle[i,4]<-sd(alphas[2,])
    results_mat_shuffle[i,5]<-quantile(alphas[2,],0.05)
    results_mat_shuffle[i,6]<-quantile(alphas[2,],0.95)
    print(paste("Computation No.",i,"out of",length(obs_vec),"(number of time series lengths) finished."),quote=F)
  }
  parallel::stopCluster(cl)
  plot(results_mat[,1],results_mat[,2],type="l",ylim=c(min(c(results_mat[,-1],results_mat_shuffle[,-1])),
                                                      max(c(results_mat[,-1],results_mat_shuffle[,-1]))),
        main="Replication of Fig. 4-type chart",
        xlab="No. of observations",ylab="delta alpha")
  lines(results_mat[,1],results_mat[,5],lty="dashed")
  lines(results_mat[,1],results_mat[,6],lty="dashed")
  points(results_mat[,1],results_mat[,3],pch=4)
  points(results_mat[,1],results_mat[,4],pch=1)
  lines(results_mat[,1],results_mat_shuffle[,2],lty="solid",col="gray")
  lines(results_mat[,1],results_mat_shuffle[,5],lty="dashed",col="gray")
  lines(results_mat[,1],results_mat_shuffle[,6],lty="dashed",col="gray")
  points(results_mat[,1],results_mat_shuffle[,3],col="gray",pch=4)
  points(results_mat[,1],results_mat_shuffle[,4],col="gray",pch=1)
  return(list(simulated=results_mat,shuffled=results_mat_shuffle))
}

#FinalStep: An example in a single function for 100 Monte Carlo repetitions and four time series lengths.
#Always rerun Step3 before running FinalStep to refresh the parallel setup.
#The whole procedure can be easily used for any data-generating process. It needs to be specified in Step1. 
#There can be different parameters in the DGP that need to be specified in Step1 and also in the MF_simulation function in Step2.
#The example function returns a trivialized chart-type as in Fig. 4 of the paper together with the list of underlying statistics.
MF_simulation_par(100,c(250,500,1000,5000))