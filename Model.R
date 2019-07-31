###################
## packages
####################

## library for using the Sterling number function
library(gmp)
library(ggplot2)
library(gridExtra)
library(grid)
library(scales)

#####################
## global parameters
####################
Nmax <- 350
nsim <- 10000

###################
## functions
####################

# falling factorial
fallingfac <- function(N,n){
  if (N<=170) {
  exp( lfactorial(N) - lfactorial(N-n))
  } else { # Stirling approximation for large factorials
    exp ( (N*log(N)-N) - ((N-n)*log(N-n)-(N-n)) )
  }
}

fallingfac <- Vectorize(fallingfac, "N")


# unique items distribution (UID, from Mendelsen et al 2016)
prob <- function(n, N, m) {
  fallingfac(N,n)/N^m*as.numeric(Stirling2(m, n))
}

# expectation of UID
expect <- function(N, m) {
  (N^m - (N-1)^m)/(N^(m-1))
}

## variance of UID
vari <- function(N, m) {
  N*(N-1)*(1-2/N)^m + N*(1-1/N)^m - N^2*(1-1/N)^(2*m)
}

## posterior approximation by sampling importance resampling
f <- function(par)
{
  n <- par[1]
  m <- par[2]
  ## mean of prior NBD
  mean <- par[3]
  ## overdispersion of prior NBD
  overdisp <- par[4]
  prior <- par[5]
  
  ## SIR algorithm 
  # 1. sample a bunch of values from a uniform - Nstar - needs to be an
  # approximate of the target distribution 
  Nstar <- round(runif(nsim, n, Nmax))
  
  # 2. calculate importance weights/importance ratios for each Xi
  # these are used as resampling weights to select the sample
  if (prior==1) {
  w <- prob(n, Nstar, m)*dnbinom(Nstar, mu=mean, size=overdisp)
  }
  else {
    w <- prob(n, Nstar, m)
  }
  
  # 3. resample (with replacement) a sample of size nsim from
  #the target distribution with the weights w
  
  samp <- (sample(Nstar, size=nsim, prob=w, replace=T))
  
  ## expectation, bias, variance,  confidence intervals, 
  ## ratio bias/m, ratio bias/N
  
  df <- data.frame(expectN=mean(samp), bias=n-mean(samp),
                    varN = var(samp), percentbias = (n-mean(samp))/mean(samp),
                    mn = m/n,
                    lwr = quantile(samp, probs=c(0.025)), 
                    upr = quantile(samp, probs=c(0.975)))
  df
}


########################################################
## 1. explore mean, variance and bias of n for different  m
########################################################

df <- expand.grid(N=seq(1,100), m=c(5, 10, 20, 40))

df <- within(df, {
  expectN <- expect(N,m)
  varN <- vari(N,m)
  bias <- expectN - N
  mN <- m/N
  percentbias <- bias/N
})

########################################################
## 2. explore posterior of N for different n, m and priors
########################################################

## run SIR simulations
df1 <- expand.grid(n=seq(1,30), m=c(5, 10, 20, 40), mean=c(20, 45), 
                   overdisp=c(0.5, 1), prior=1)
df2 <- expand.grid(n=seq(1,30), m=c(5, 10, 20, 40), mean=c(0), 
                   overdisp=c(0), prior=0)

tmp <- vector("list", nrow(df1))
for (i in 1:nrow(df1)) {
  if (df1[i,1]>df1[i,2]) {
    tmp[[i]] <- NA
  } else {
  tmp[[i]] <- f(as.numeric(df1[i,]))
  }
}

df1 <- cbind(df1, do.call(rbind, tmp))


tmp <- vector("list", nrow(df2))
for (i in 1:nrow(df2)) {
  if (df2[i,1]>df2[i,2]) {
    tmp[[i]] <- NA
  } else {
    tmp[[i]] <- f(as.numeric(df2[i,]))
  }
}

df2 <- cbind(df2, do.call(rbind, tmp))

df0 <- rbind(df1, df2)