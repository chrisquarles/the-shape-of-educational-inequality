## R Code which contains analysis for
## "The Shape of Educational Inequality" by Quarles, Budak & Resnick (in Science Advances)
## Author: Christopher Quarles
## Email: chrisquarles@gmail.com
## Last Date Updated: May 13, 2020
##
## This code is provided to enable easier reproduction of the images in the paper.
## Due to privacy concerns, we were unable to share the data. So some notes have
## been provided about how to process your own data to use with these scripts.

source('mlecens.R')

##########################################################################
############### A Note on the Data  ######################################
##########################################################################
## The dataset used in the paper contained data for 5 years worth of 
## entering students at each of 28 colleges. Graphs were made of the whole
## group of students, individual colleges (all 5 years), and individual
## year-college cohorts. For instance Figure S2 had 3 graphs for each
## college. Because that data is unavailable, in this document I've 
## shown how to make a single graph of each type for one such cohort 
## or college. The code assumes that the dataset for a single 
## college or cohort is called coldat, and that the dataset has the
## following variables:
##
##  credits_earned = # of credits earned by a given student
##  droppedout = FALSE if the student graduated or transferred, TRUE 
##                otherwise (can also be 1/0)
##  transferred = TRUE iff the student transferred to a 4-year college (can also be 1/0)
##  transnograd = TRUE iff the student transferred but didn't graduate (can also be 1/0)

## The sample data is in this format, and you can run this code on 
## sample_data.csv. Or you can put your own data in that format and
## load it by editing the following line.

coldat <- read.csv("sample_data.csv")

## Make sure the variables are TRUE/FALSE, rather than 1/0
coldat$droppedout <- as.logical(coldat$droppedout)
coldat$transferred <- as.logical(coldat$transferred)
coldat$transnograd <- as.logical(coldat$transnograd)


###########################################################################
################### Preliminary functions #################################
###########################################################################
## These are functions used to process the data. They need to be run first.

###############################
## custom_G functions
## These functions are the standard d, p, q, and r functions typical for distributions in R.
## Except,in these, we first (a) estimate the distribution of success points (G_k) using x
## and yc from real data, then (b) use that distribution to generate the distribution functions.
## I've found rG very useful, because it can generate success points for a simulated 
## group of student.
## Inputs:
## a = vector of quantiles
## p = vector of probabilities
## n = number of observations
## x,yc = as in mlecens

## The probability mass function, defined on positive integers
dG <- function(a, x, yc){
  G <- mlecens(x=x, yc=yc, dist="custom_G")
  if(max(abs(a%%1)) != 0){stop("a must be a vector of positive integers")}
  if(min(a) <= 0){stop("a must be a vector of positive integers")}
  output <- vector(mode="numeric", length=length(a)) #create the output vector of 0's
  # We leave elements of output as 0 if those lements of a are larger than the max G
  kmax <- length(G)
  aok <- (a <= kmax)
  output[aok] <- G[a[aok]]
  return(output)
}

## The cumulative probability function, defined on positive integers
pG <- function(a ,x, yc){
  G <- mlecens(x=x, yc=yc, dist="custom_G")
  cumG <- cumsum(G)
  if(max(abs(a%%1)) != 0){stop("a must be a vector of integers")}
  if(min(a) <= 0){stop("a must be a vector of positive integers")}
  output <- rep(1, times=length(a)) #create the output vector of 1's
  ## This will be the default for values which are so large that the cdf is 1 there.
  ## We calculate those values of a that are smaller than the max value
  kmax <- length(cumG)
  aok <- (a <= kmax)
  output[aok] <- cumG[a[aok]]
  return(output)
}

## qG is the quantile function (inverse cdf) for the distribution defined by custom_G
qG <- function(p, x, yc){
  G <- mlecens(x=x, yc=yc, dist="custom_G")
  cumG <- cumsum(G)
  basefun <- function(b){
    min(which(b <= cumG))
  }
  return(sapply(p, FUN=basefun))
}

## This function generates a random dataset for g_i using the estimated distribution of G.
rG <- function(n, x, yc){
  quantilevec <- runif(n=n)
  output <- qG(p=quantilevec,x=x,yc=yc)
  return(output)
}

################################
## smoothG
## This function takes the output of mlecens with custom_G, and smooths it using 
## splines, so that the distribution is easier to visualize.
smoothG <- function(G){
  myxvals <- 1:length(G)
  myspline <- smooth.spline(myxvals, cumsum(G),spar=.6) #creates a smooth cdf
  smoothedG <- vector(mode="numeric", length = length(G))
  smoothedG[1] <- max(0,myspline$y[1])
  smoothedG[2:length(G)] <- myspline$y[2:length(G)] - myspline$y[1:(length(G)-1)]
  smoothedG <- smoothedG/sum(smoothedG)
  return(smoothedG)
}

###########################################################################
## discrete normal functions 
## These functions are the standard d, p, q, and r functions typical for distributions in R.
## They take the mean and standard deviation of the truncated discrete normal distribution
## as used here, and give the values for that. Note that this function starts at x=1, so
## that P(x<=0) == 0. Note that a truncated normal distribution generated with mean=m will
## not have a mean of m, since we're cutting off the left tail at 0.
## Inputs:
## a, q = vector of quantiles
## p = vector of probabilities
## n = number of observations
## m = nominal mean of the distribution (though this will be different from the true mean)
## sd = nominal sd of the distribution

## The probability mass function, defined on positive integers
ddiscnorm <- function(a, mean, sd){
  ## First, make sure the inputs are discrete
  if(suppressWarnings(!all(round(a)==a))){warning("a must consist of only integer values.")}
  if(!all(a > 0)){warning("a must consist of only positive values.")}
  ## Also, make sure the values aren't too big.
  if(mean + 6*sd > 10^12){stop("This function can't handle values of mean and sd that large.")}
  
  ## Now calculate the upper bound for values to calculate
  max_nonzero_value <- ceiling(mean + 6*sd)
  
  ## The discrete norm is the same as the normal distribution, except it's truncated
  ## at 0, and it's normalized only on the integers. So we calculate the normalizing factor.
  normalization_constant <- 1/sum(dnorm(1:max_nonzero_value, mean = mean, sd=sd))
  
  ## Now calculate the values
  output <- normalization_constant*dnorm(x=a, mean=mean, sd=sd)
  output[a <= 0] <- 0
  output[a > max_nonzero_value] <- 0
  return(output)
}

## The cumulative probability function, defined on positive integers
pdiscnorm <- function(q, mean, sd, lower.tail=TRUE){
  ## First, make sure the inputs are discrete and positive
  if(suppressWarnings(!all(round(q)==q))){warning("q must consist of only integer values.")}
  if(!all(q > 0)){warning("q must consist of only positive values.")}
  ## Also, make sure the values aren't too big.
  if(mean + 6*sd > 10^12){stop("This function can't handle values of mean and sd that large.")}
  
  ## Calculate the full range of pmf values
  pmf_vals <- ddiscnorm(a=1:ceiling(mean+6*sd), mean=mean, sd=sd)
  ## Calculate the full range of cumulative mass function values
  cmf_vals <- cumsum(pmf_vals)
  output <- cmf_vals[q]
  if(lower.tail == TRUE){
    return(output)
  }
  if(lower.tail == FALSE){
    return(1-output)
  }
}

## qdiscnorm is the quantile function (inverse cdf) for the truncated
## discrete normal distribution
## Note: If a percentile == 1, then this function will return (mean+6*sd)
qdiscnorm <- function(p, mean, sd){
  cmf_vals <- pdiscnorm(1:(mean+6*sd), mean=mean, sd=sd)
  basefun <- function(b){
    min(which(b <= cmf_vals))
  }
  return(sapply(p, FUN=basefun))
}

## This function generates a random dataset for the truncated normal distribution
rdiscnorm <- function(n, mean, sd){
  quantilevec <- runif(n=n)
  output <- qdiscnorm(p=quantilevec,mean=mean, sd=sd)
  return(output)
}



###########################################################################
###########################################################################
## Figure 1A & 1B - Models of Student Capital 
##
## These graphs didn't use any real data, since they're for just models.

## Define points on the x-axis
xvals <- 1:100

## Define functions that give PDF for each model
cogabil.pdf <- function(vec){dnorm(vec,mean=50, sd=20)/sum(dnorm(1:10000, mean=50, sd=20))}
finite.pdf <- function(vec){dgeom(vec, prob=.02)}
rich.pdf <- function(vec){dzeta(vec, shape=.1348)}

## Hazard vectors = % of people who drop out/ % of people at k
cogabilpdf <- dnorm(xvals,mean=50, sd=20)/sum(dnorm(seq(1,10000, by=1), mean=50, sd=20))
cogabilcumul <- 1-cumsum(c(0,cogabil.pdf(xvals))) #1-cdf
cogabilhaz <- 1- (cogabilcumul[2:length(cogabilcumul)]/cogabilcumul[1:(length(cogabilcumul)-1)])
finitehaz <- rep(0.02, times=100)
richcumul <- 1-cumsum(c(0,rich.pdf(xvals))) #1-cdf
richhaz <- 1- (richcumul[2:length(richcumul)]/richcumul[1:(length(richcumul)-1)])

### Figure 1A - PDF Plot
pdf("Figure 1A Models Prob.pdf", height = 3.6, width = 4.8)
plot(xvals, cogabil.pdf(xvals), type="l", ylim=c(0, 0.05), col="black", lwd=2, 
     xlab="Student Capital", ylab="Probability Distribution", cex.lab=1, las=1)
lines(xvals, finite.pdf(xvals), col="red", lwd=2)
lines(xvals, rich.pdf(xvals), col="green", lwd=2)
legend(x=30, y=.05, 
       legend=c("Cognitive Ability Model (Normal)",  
                "Finite Resource Model (Exponential)", "Rich-Get-Richer Model (Power Law)"),
       col=c("black", #"blue", 
             "red", "green"), lwd=2, cex = .7)
dev.off()

## Figure 1B - Hazard Plot
pdf("Figure 1B Models Hazard.pdf", height=3.6, width=4.8)
plot(xvals, cogabilhaz, type="l", col="black", ylim=c(0, .14), lwd=2, 
     xlab="Student Capital", ylab="Hazard Rate", cex.lab=1, las=1)
lines(xvals, finitehaz, col="red", lwd=2)
lines(xvals, richhaz, col="green", lwd=2)
legend(x=6, y=.14, 
       legend=c("Cognitive Ability Model (Normal)", 
                "Finite Resource Model (Exponential)", "Rich-Get-Richer Model (Power Law)"),
       col=c("black",  "red", "green"), lwd=2, cex=.7)
dev.off()



###########################################################################
## Figures 2 & S1 - Distributions of credits earned w/graduation & transfer
##
## In the paper, there was one graph for each college. Here, you'll find the
## code to make one such graph.

## kmax is the largest number of credits that you want to plot. Adjust it to fit your needs.
## The bin width (5) was chosen because many classes in WA are 5 credits.
pdf("Figure 2.pdf", height=3, width=4)
kmax <- 265
hist(coldat$credits_earned, breaks=(kmax/5), xlim=c(0,kmax), main=NULL,
     las=1, xlab="")
hist(coldat[!coldat$droppedout,]$credits_earned, breaks=(kmax/5), col="blue", add=T, xlim=c(0,kmax))
hist(coldat[coldat$transferred,]$credits_earned, breaks=(kmax/5), col="green", add=T, xlim=c(0,kmax))
hist(coldat[coldat$transnograd,]$credits_earned, breaks=(kmax/5), col="yellow", add=T, xlim=c(0,kmax))
legend(115, sum(coldat$credits_earned <= 5), legend=c("Dropped Out", "Graduated Only", "Graduated & Transferred", "Transferred Only"), 
       fill=c("white", "blue", "green", "yellow"), bty="n", cex=.7)
title(main="College A", line=0)
title(xlab="Credits Earned", line=2.5)
dev.off()



###########################################################################
## Figure 3 QQ plots for the three parametric models
##
## To create this, we first need to estimate the parameters for the cohort
## and then create synthetic data using those parameters

## Generate synthetic data for the cognitive-ability/normal model
set.seed(-70671)
ms <- mlecens(x=coldat$credits_earned, yc=coldat$droppedout, dist="normal_Y") #estimate parameters
ysim_norm <- rdiscnorm(n=10000, mean=ms[1], sd=ms[2])  # generate synthetic distribution of student capital
gsim_norm <- rG(n=10000, x=coldat$credits_earned, yc=coldat$droppedout) # generate synthetic dist. of graduation point
xsim_norm <- pmin(ysim_norm, gsim_norm)  # Find the # of credits each synthetic student would earn
ycsim_norm <- (ysim_norm < gsim_norm) # Did each student drop out?

## Generate synthetic data for the finite resource/exponential model
set.seed(-70671)
q <- mlecens(x=coldat$credits_earned, yc=coldat$droppedout, dist="exponential_Y") #estimate parameters
ysim_exp <- rgeom(n=10000, prob=1-q) + 1  # generate synthetic distribution of student capital
gsim_exp <- rG(n=10000, x=coldat$credits_earned, yc=coldat$droppedout) # generate synthetic dist. of graduation point
xsim_exp <- pmin(ysim_exp, gsim_exp)  # Find the # of credits each synthetic student would earn
ycsim_exp <- (ysim_exp < gsim_exp) # Did each student drop out?

## Generate synthetic data for the rich-get-richer/power law model
## Note that rzeta, which generates numbers according to a power law, sometimes tries to generate numbers
## very far out in the tail. This can take a long time. Sometimes it throws a warning and rounds to the last
## large number that it found. This doesn't affect our result, because these are very large credit values on the
## order of 10^30 that don't occur in colleges.
set.seed(-70671)
a <- mlecens(x=coldat$credits_earned, yc=coldat$droppedout, dist="power_Y") #estimate parameters
ysim_pow <- suppressWarnings(rzeta(n=1000, shape=a-1)) # generate synthetic distribution of student capital
gsim_pow <- rG(n=10000, x=coldat$credits_earned, yc=coldat$droppedout) # generate synthetic dist. of graduation point
xsim_pow <- pmin(ysim_pow, gsim_pow)  # Find the # of credits each synthetic student would earn
ycsim_pow <- (ysim_pow < gsim_pow) # Did each student drop out?


## Make the graph. Each of the graphs could be generated individually as well.
pdf("Figure 3 QQ plot.pdf", height=4, width=12, useDingbats = FALSE)
par(mfrow=c(1,3))
  ## Normal
qqplot(xsim_norm, coldat$credits_earned, xlab="Normal Synthetic Data", ylab="Real Data", cex.lab=1.3,
       xlim=c(0, 250), ylim=c(0,250))
abline(coef=c(0,1), col="red") # A good fit lies on this line.
  ## Exponential
qqplot(xsim_exp, coldat$credits_earned, xlab="Exponential Synthetic Data", ylab="Real Data", cex.lab=1.3,
       xlim=c(0, 250), ylim=c(0,250))
abline(coef=c(0,1), col="red") # A good fit lies on this line.
  ## Power Law
qqplot(xsim_pow, coldat$credits_earned, xlab="Power Law Synthetic Data", ylab="Real Data", cex.lab=1.3)
abline(coef=c(0,1), col="red") # A good fit lies on this line.
dev.off()



###########################################################################
## Figure 4A - Average Student Capital by Cohort 
##
## This graph was created by calculating
q = mlecens(x=coldat$credits_earned, yc=coldat$droppedout)
mu_s = 1/(1-q)
## for each cohort, and then plotting a histogram of these values.



###########################################################################
## Figures 4B & S2 - Actual vs Simulated Dropout Rates (And correlation)
##
## This graph was created by calculating the following values for each cohort
actual_dropout_rate = mean(coldat$droppedout)
est_dropout_rate = mean(ycsim_exp) ## ycsim_exp is the vector calculated for Fig. 3

## I put each of these types of statistics into their own vector, and
## then created a dotplot with abline(coef=c(0,1), col="red")



###########################################################################
## Figures S3 & S4 
##
## These are histograms of summary statistics by college and cohort, created
## using the hist() function.



###########################################################################
## Figure S5 - Inferred distribution of success points
##
## First calculate the distribution of success points
G = mlecens(x=coldat$credits_earned, yc=coldat$droppedout, dist="custom_G")
## Then make the graph
png("Figure S5 Success Points.png", height=600, width=800)
plot(1:150, G[1:150], pch=19, xlab="Credits", ylab="Probability of Graduation/Transfer",
     cex.lab=1.5, xaxp=c(0,150,15), main="Distribution of Success Points")
dev.off()



###########################################################################
## Figure S6 - Inferred, smoothed distribution of success points
##
## Calculate G as in Figure S5
G <- mlecens(x=coldat$credits_earned, yc=coldat$droppedout, dist="custom_G")
## Then use the smoothG function on this distribution, which smooths with
## splines. Note that it starts at 3, since splines need multiple data points.
png("Figure S6 Smoothed Success Points.png", height=600, width=800)
plot(3:150, smoothG(G)[3:150], pch=19, xlab="Credits", ylab="Probability of Graduation/Transfer",
     cex.lab=1.5, xaxp=c(0,150,15), main="Smoothed Distribution of Success Points",
     type="l", ylim=c(0, .015))
dev.off()



###########################################################################
## Table S1 - AIC for each cohort & model 
##
## The Akaike information criterion (AIC) is given by 
## AIC = 2*(# of parameters) - 2*(log-likelihood of the model).
## The AIC will depend on the sample size, so you can't compare AICs of
## different datasets.
## Here's how it was calculated for the best-fit model of each type.

## Calculate the AIC for the normal distribution on the given college. (Ignore the warning.)
normmle <- mlecens(x=coldat$credits_earned, yc=coldat$droppedout, dist="normal_Y", func.out=TRUE)
AIC_norm <- 2*2 - 2*normmle$LLfunc(normmle$estimate)
## Calculate the AIC for the exponential distribution on the given college.
expmle <- mlecens(x=coldat$credits_earned, yc=coldat$droppedout, dist="exponential_Y", func.out=TRUE)
AIC_exp <- 2*1 - 2*expmle$LLfunc(expmle$estimate)
## Calculate the AIC for the power law distribution on the given college.
powmle <- mlecens(x=coldat$credits_earned, yc=coldat$droppedout, dist="power_Y", func.out=TRUE)
AIC_pow <- 2*1 - 2*powmle$LLfunc(powmle$estimate)

