## Maximum Likelihood Estimation for Right-Censored Data 
## Accompanies "The Shape of Educational Inequality" by Quarles, Budak & Resnick (Science Advances)
## Author: Christopher Quarles
## Email: chrisquarles@gmail.com
## Last Date Updated: May 13, 2020
## About: This file does maximum likelihood estimation on right-censored data. It was designed to
## work with community college data, specifically the distribution of credits earned by a
## cohort of students. 

## This file uses the VGAM package for its zeta() function.
library(VGAM)


#########################################################################################
##  mlecens - right-censored maximum likelihood estimator
#########################################################################################
## This is the main function for doing the analysis. It takes information about individual
## students and returns information about the student capital of the cohort. 
## Specifically, you need to specify the number of credits each student earned as a vector  
## of positive integers (x) and a vector (yc) which tells which students dropped out without 
## earning a degree or transferring (TRUE means they dropped out). 
## For a specified distribution (dist), the function outputs a set of parameters for that
## distribution (outlined below). If func.out=TRUE, the function will also return the 
## log-likelihood function for the given data and distribution.
## The function can also create the point-wise estimate of the distribution of "success 
## points", i.e. the probability that a student will graduate/transfer at any given credit
## level. These pointwise estimates will be very sensitive to the input data, and should not 
## be used to calculate anything but the general shape of when students finish school.
## 
## x = a vector of the number of credits earned for each student, must be positive integers. 
##        Exclude all individuals who earned 0 credits
## yc = a binary(T/F) vector of the same length of x which describes whether the student dropped out.
##        Equivalently, this is 1-(completion?) for each student.
## dist = the name of the assumed distribution (see below)
## lims = the limits (if necessary) within which to search for the parameter(s). Only used
##        for power_Y. Even then, the default limits should work. For power_Y, the interval 
## func.out = Explains whether mlecens should include the log-likelihood function in
##             its output. Does not work for custom_G.
##
##  The variable dist can be:
##    "exponential_Y" for the distribution Y_k,
##      where p_k=(1-q)q^(k-1), and where the cdf from k to Inf =is equal to q^k.  
##      This uses an analytical method, and outputs the value q.
##    "power_Y" for the distribution Y_k, with p_k=(1/zeta(a))*(1/k^a)
##      In this case, lims needs to be a length 2 vector which gives the endpoints of
##      the range for finding a. Note that the calculations may get challenging and/or
##      less reliable as the lower limit gets closer to 1. Both values must be larger
##      than 1.  This uses a computational method, and outputs the value a.
##    "normal_Y" for the truncated discrete normal distribution Y_k, 
##      with p_k= (1/A) *exp(-(k-mu)^2/(2*sigma^2)), where A is a normalizing factor.
##      This uses a computational method and outputs the two-vector c(mu, sigma).
##      In this case, the output log-likelihood function either takes a two-vector ms,
##      where the first element is mu and the second element is sigma, OR it can take
##      in an nx2 matrix or dataframe, where the first column represents mu values and the
##      second column represents sigma values. 
##    "custom_G" for the distribution G_k,
##      which finds each G_k individually so that the KL divergence of the model w.r.t. the
##      data set is minimized. The output of this is a vector whose k-th element gives the 
##      probability of a randomly chosen student graduating after exactly k credits.
##      The length of this vector is kGmax, the largest credit value for someone who graduated.


mlecens <- function(x, yc, dist="exponential_Y", lims=c(1.1, 4), func.out=FALSE){
  ## Create a vector gc, which tells whether a student graduated/transferred (i.e. didn't drop out)
  yc <- as.logical(yc)
  gc <- !yc
  
  ## Do some variable checks
  if(length(x)!=length(yc) | length(x)!=length(gc)){
    stop("Your data vectors need to be all the same length.")
  }
  
  ## Create some useful variables
  n <- length(x)  # number of data points
  xmax <- max(x)
  if(xmax <=7){stop("max(x) must be greater than 7")} #Somebody needs to earn at least 8 credits for this to work.
  xbar <- mean(x)  # Average number of credits earned
  ## Create vectors for the number of times non-censored and censored values take on the value
  ## of each possible data value. 
  tempfunc <- function(k,vec){sum(vec==k)}
  nY <- sapply(X=1:xmax, FUN=tempfunc, vec=x[yc]) #k-th element is # of ppl who dropped out after k credits
  nG <- sapply(X=1:xmax, FUN=tempfunc, vec=x[gc]) #k-th element is # of ppl who grad/trans after k credits
  if(n != sum(nY+nG)){stop("yc+gc does not equal the one vector. You probably have 
                              credit values that are not positive integers.")}
  
  
  ## Fitting Y to a power law (zeta) distribution
  if(dist=="power_Y"){
    ## We need to use the VGAM zeta() function here.
    ## As of 6/13/18, zeta allows you to put vectors into the shift argument, but the 
    ## output vector is incorrect. The next function, zetashift, fixes that so that it's vectorized
    ## in the shift argument. zetashift is NOT vectorized in both components.
    ## x still needs to be atomic.
    zetashift <- function(x, shift){
      if(length(shift)==1){zeta(x=x,shift=shift)}
      if(length(shift)>1){
        zeta.short <- function(sh){
          zeta(x=x, shift=sh)
        }
        return(sapply(X=shift, FUN=zeta.short))
      }
    }
    ## Power laws can give big values, so we only use the observed values of x, rather than
    ## the complete range 1:xmax
    kvals <- as.integer(names(table(x))) 
    nY <- sapply(X=kvals, FUN=tempfunc, vec=x[yc])
    nG <- sapply(X=kvals, FUN=tempfunc, vec=x[gc])
    if(n != sum(nY + nG)){stop("yc+gc does not equal the one vector. You probably have 
                              credit values that are not positive integers.")}
    ylnmean <- mean(log(x[yc]))
    LLunv <- function(param){ #the log-likelihood function. Doesn't accept vector arguments
      -n*log(zeta(param))-param*sum(nY)*ylnmean + sum(nG*log(zetashift(param,kvals)))
    }
    LL <- function(param){ #the vectorized log-likelihood function
      ifelse(length(param==1),
             return(LLunv(param=param)),
             return(sapply(X=param,FUN=LLunv))
      )
    }
    estimate <- optim(par=2, fn=LLunv, method="L-BFGS-B", lower=lims[1], upper=lims[2], 
                      control=list(fnscale=-1))$par
  }
  

  ## Fitting Y to an exponential (geometric) distribution
  if(dist=="exponential_Y"){
    qhat <- (xbar-1)/(xbar-(sum(nG)/n)) #estimator calculated analytically
    estimate <- qhat
    LLunv <- function(param){ #the log-likelihood function. Doesn't accept vector arguments
      sum(nY)*log(1-param)+n*(xbar-1)*log(param)
    }
    LL <- function(param){ #the vectorized log-likelihood function
      ifelse(length(param==1),
             return(LLunv(param=param)),
             return(sapply(X=param,FUN=LLunv))
      )
    }
  }
  
  
  ## Fitting Y to a truncated discrete normal distribution
  if(dist=="normal_Y"){
    ## This method is trickier, because it requires maximizing a complicated log-likelihood function
    ## for the discrete, truncated normal distribution on a non-rectangular trapezoid.
    ## First, create the unvectorized log-likelihood function, which only accepts vectors of length 
    ## two: ms = (mu, sigma)
    ## The normal_Y log-likelihood function is complicated to compute due to overflow/underflow problems.
    ## It will often return NA. These are areas where the log-likelihood function is very negative anyways,
    ## so these points are away from the maximum.
    LLunv <- function(ms){ 
      m <- ms[1]
      s <- ms[2]
      kmax <- ceiling(m+6*s)
      kmin <- max(floor(m-6*s),1)
      kvals <- kmin:kmax
      sq <- -.5*(( (1:max(kvals,xmax)) - m)/s)^2
      kexp <- exp(sq[kmin:kmax])
      partialsum <- rev(cumsum(rev(kexp)))
      if(kmin==1) {weirdsum <- partialsum}
      if(kmin>=2) {weirdsum <- c(rep(partialsum[1], times=kmin-1),partialsum)}
      output <- -n*log(partialsum[1])+sum(nY*sq[1:xmax])+sum(nG*log(weirdsum)[1:xmax])
      return(output)
    }
    ## Now create the vectorized log-likelihood function. It takes an nx2 matrix/data frame as an argument.
    LL <- function(ms){ 
      if(length(ms)==2){
        return(LLunv(ms))
      }
      if(length(ms)>2){
        ms <- as.matrix(ms)
        if(ncol(ms)!=2){ #error checking
          stop("For mlecens, input must be an nx2 matrix or dataframe.")
        }
        return(apply(X=ms,MARGIN=1,FUN=LL))
      } 
    }
    ## Now rescale the input values so that we can run optim on a rectangle.
    ## The output of f is a value of s, the standard deviation.
    f <- function(mv){
      m <- mv[1]
      v <- mv[2]
      (v*(m+5*xmax)-m+xmax)/6
    }
    ## Create the negative log-likelihood function, rescaled to be on a rectangle
    LLf_neg <- function(mv){
      -LLunv(c(mv[1],f(mv))) #use LLunv, because we're only plugging in two-vectors.
    }
    ## Minimize the negative log-likelihood function (or maximize the log-likelihood function)
    optim1 <- optim(par=c(mean(x[yc]),.5), fn=LLf_neg, method="L-BFGS-B", lower=c(1,0), upper=c(xmax-6,1))
    ## Note: For some initial values, optim throws an error. I'm not
    ## sure why this is, and it's super-sensitive. For instance, if (a, .5) was a bad initial 
    ## value, then (a+.001, .5) and (a-.001, .5) can be fine. 
    ## To deal with this, I included the following loop to retry with different starting values
    ## if the initial initial value didn't converge.
    i <- 1
    convergence_test_val <- optim1$convergence
    while(convergence_test_val != 0){ ## Start robustness loop
      ## If this loop has happened too many times, break out and throw an error.
      if(i >= 10){
        cat("For some reason, optim isn't converging when trying to find the max log-likelihood. \n")
        cat("In optim convergence = ", optim1$convergence, "\n")
        stop(optim1$message)
      }
      ## Add a small random amount to the initial value, and then run optim again.
      ## The 5 is there so that if somehow mean(x[yc]) is small enough to get
      ## negative initial values, we don't get errors.
      new_x0 <- max(mean(x[yc]) + .1*runif(1,-i,i), 5)
      optim1 <- optim(par=c(new_x0,.5), fn=LLf_neg, method="L-BFGS-B", lower=c(1,0), upper=c(xmax-6,1))
      convergence_test_val <- optim1$convergence
      i <- i+1
    } ## End robustness loop
    m1 <- optim1$par[1]
    v1 <- f(optim1$par)
    estimate <- c(m1,v1)
  }
  
  ## Fitting G to a custom distribution that minimizes KL divergence
  if(dist=="custom_G"){
    if(func.out){stop("The method for finding custom_G does not give a log-likelihood function. 
                      Change func.out to FALSE.")}
    kGmax <- max(x[gc]) ## The max credit value for a graduating/transfering student
    n <- n[1:kGmax]
    nY <- nY[1:kGmax]
    nG <- nG[1:kGmax]
    ## We assume that G_k=m_k*lambda, where lambda is a normalizing factor. Using the max
    ## likelihood calculations, we find m inductively.
    m <- vector(mode="numeric", length=kGmax)
    m[kGmax] <- 1
    a <- vector(mode="numeric", length=kGmax)
    a[kGmax] <- nG[kGmax]
    for(k in (kGmax-1):1){
      a[k] <- a[k+1] + nY[k]/sum(m[(k+1):kGmax])
      m[k] <- nG[k]/(a[k])
    }
    ## Now we calculate G.
    G <- vector(mode="numeric", length=xmax)
    if(xmax > kGmax){
      G[(kGmax+1):xmax] <- 0
    }
    G[1:kGmax] <- m/sum(m)
    estimate <- G
    }
  
 ## Return the results 
  if(!func.out){  #if func.out is false, return an integer
    return(estimate)
  } else { # otherwise, return a list
    output <- list(estimate=estimate)
    if(func.out){output$LLfunc <- LL}
    ## Give a warning if the user is asking for the log-likelihood function of
    ## the discrete normal distribution.
    if((func.out == TRUE) & (dist == "normal_Y")){
      warning("The log-likelihood function for the truncated discrete normal distribution may return NA at certain points due to computer rounding. These points are away from the maximum value. This is probably not an issue, unless you see an NA where you don't expect it.")
    }
    return(output)
  }
}

