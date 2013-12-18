# ----- Single-valued calculations ----- #

#' Compute sensitivity, specificity, positive likelihood ratio, negative likelihood ratio for a single 2x2 table
#' @param truePos The number of true positive tests
#' @param totalPos The total number of positives ("sick") in the population
#' @param trueNeg The number of true negatives in the population
#' @param totalNeg The total number of negatives ("well") in the population
#' @return A matrix containing sensitivity, specificity, posLR, negLR results
#' @references Deeks JJ, Altman DG. BMJ. 2004 July 17; 329(7458): 168-169.
confusionStatistics <- function( truePos, totalPos, trueNeg, totalNeg ) {
  n <- length(truePos)
  res <- matrix( NA, ncol=4, nrow=n )
  colnames(res) <- c("sens","spec","posLR","negLR")
  res[,"sens"] <- truePos / totalPos
  res[,"spec"] <- trueNeg / totalNeg
  res[,"posLR"] <- res[,"sens"] / ( 1 - res[,"spec"] )
  res[,"negLR"] <- ( 1 - res[,"sens"] ) / res[,"spec"]
  res
}

# ----- Optimization tools ----- #

#' Find the lowest population probability whose median is consistently one
#' This is the lowest estimate for Sens that is consistently (over 5 runs) most likely to yield a sample estimate of 100/100.
#' @param pr Probability input
#' @param size Number of trials
#' @param R number of bootstrap replications
#' @param nConsistentRuns Number of runs that all have to be identical to return TRUE
#' @param warn Warn if searching outside of the range c(0,1)
#' @return Boolean of length one (TRUE or FALSE)
#' @examples
#' \dontrun{
#' prs <- seq(.990,.995,.0001)
#' bools <- sapply( prs, medianConsistentlyOne, size=truePos, R=R )
#' data.frame( prs, bools )
#' }
medianConsistentlyOne <- function(pr, size, R, nConsistentRuns=5, warn=TRUE) {
  if( 0 > pr | pr > 1) {
    if(warn)  warning("Searching probabilities outside of 0,1. Returning FALSE.")
    return( FALSE )
  } else {
    reps <- replicate( nConsistentRuns, median( rbinom(R, size=size, prob=pr) ) )
    return( all( reps==size ) )
  }
}

#' Optimize a function returning a single numeric value subject to a boolean constraint
#' Utilizes a naive recursive grid search
#' @param f Function to be minimized: takes a single numeric value and returns a single numeric value
#' @param constraint Function of a single variable returning a single boolean value (must be TRUE to be at the optimum)
#' @param bounds A numeric vector of length two which are the upper and lower bounds of the input to try
#' @param nEach Number of points n each round of grid searching to use
#' @param shrink Factor indicating how much (1/shrink) to narrow the search width by each round. Highly recommended that shrink is at least half the size of nEach.
#' @param tol The tolerance (epsilon)
#' @param verbose Whether to display verbose output
#' @param \dots Arguments to pass along to constraint
#' @return The optimized input value (numeric)
sequentialGridSearch <- function( f, constraint, bounds, nEach=40, shrink=10, tol=.Machine$double.eps ^ 0.5, verbose=FALSE, ... ) {
  #! The alabama package or similar might be a better way of doing this in the future.
  if(verbose) cat("Grid searching between",bounds[1],"and",bounds[2],"\n")
  x <- seq( from=bounds[1], to=bounds[2], length.out=nEach )
  fx <- f( x )
  if( any(is.na(fx)) ) stop("NAs produced while evaluating f")
  cx <- constraint( x, ... )
  if( any(is.na(cx)) ) stop("NAs produced while evaluating constraint")
  if( !any(cx) ) stop("No value found while searching between",bounds[1],"and",bounds[2],". Try setting a looser tolerance, a lower shrinkage value, or a higher number for nEach.\n")
  newVal <- x[ which( fx==min( fx[cx] ) ) ] #! Not very efficient
  if(verbose) cat("Newval found as:", newVal,", producing value",f(newVal),"\n")
  if( exists("lastVal") && abs(f(newVal)-f(lastVal))<tol ) { # Successive rounds within tolerance
    if(verbose) cat("Successive rounds within tolerance.  Success!\n")
    return(newVal)
  } else if( sum(cx) > 1 && abs( f(newVal) - min(fx[cx & x!=newVal]) )<tol ) { # Two values this round are within tolerance
    if(verbose) cat("Two values this round within tolerance.  Success!\n")
    return(newVal)
  } else { # No values within tolerance, but at least one value still works--keep recursing!
    if(verbose) cat("No values within tolerance.  Recursing.\n")
    newHalfRange <- ( abs(diff(range(bounds)))/shrink ) / 2
    newBounds <- c( newVal - newHalfRange, newVal + newHalfRange ) # New bounds with narrower range, centered on crude optimum so far
    lastVal <- newVal
    return( sequentialGridSearch( f=f, constraint=constraint, bounds=newBounds, nEach=nEach, shrink=shrink, tol=tol, verbose=verbose, ... ) ) 
  }
}

# ----- Main function and its helpers ----- #

#' Compute the (negative) likelihood ratio with appropriate, bootstrapped confidence intervals
#' @param truePos The number of true positive tests
#' @param totalPos The total number of positives ("sick") in the population
#' @param trueNeg The number of true negatives in the population
#' @param totalNeg The total number of negatives ("well") in the population
#' @param R is the number of replications in each round of the bootstrap (has been tested at 50,000 or greater)
#' @param verbose Whether to display internal operations as they happen
#' @param \dots Arguments to pass along to boot.ci for the BCa confidence intervals
#' @return An object of class lrtest
#' @export BayesianLR.test
#' @examples
#' blrt <- BayesianLR.test( truePos=100, totalPos=100, trueNeg=60, totalNeg=100 )
#' blrt
#' summary(blrt)
#' BayesianLR.test( truePos=98, totalPos=100, trueNeg=60, totalNeg=100 )
#' BayesianLR.test( truePos=60, totalPos=100, trueNeg=100, totalNeg=100 )
#' BayesianLR.test( truePos=60, totalPos=100, trueNeg=99, totalNeg=100 )
BayesianLR.test <- function( truePos, totalPos, trueNeg, totalNeg, R=5*10^4, verbose=FALSE, ... ) {
  # -- Check inputs -- #
  if( R < 5*10^4 ) warning("Setting the number of bootstrap replications to a number lower than 50,000 may lead to unstable results")
  if( totalPos == 0 | totalNeg == 0 ) stop("This package may seem like magic, but not even magic will solve your problem (totalPos or totalNeg = 0).")
  if( trueNeg > totalNeg | truePos > totalPos ) stop("You cannot have more test positive/negative than you have total positive/negative.")
  
  # -- Bootstrap sensitivity and specificity -- #
  cs <- confusionStatistics( truePos=truePos, totalPos=totalPos, trueNeg=trueNeg, totalNeg=totalNeg )
  csExact <- cs
  
  bootmean <- function(x,i)  mean(x[i])
  
  if( truePos == totalPos ) {
    sensb <- drawMaxedOut( n=totalPos, R=R, verbose=verbose )
    cs[,"sens"] <- attr(sensb,"lprb") #! Why is this lprb--that's the lower bound of sensitivity not the central tendency?
  } else {
    sensb <- boot(
      rep( 1:0, c( truePos, totalPos-truePos ) ), 
      bootmean, 
      R=R
    )$t
  }
  
  if( trueNeg == totalNeg ) {
    specb <- drawMaxedOut( n=totalNeg, R=R, verbose=verbose )
    cs[,"spec"] <- attr(specb,"lprb")
  } else {
    specb <- boot(
      rep( 1:0, c( trueNeg, totalNeg-trueNeg ) ), 
      bootmean, 
      R=R
    )$t
  }
  
  # -- Compute pos/neg LRs and their BCa confidence intervals -- #
  negLR <- unname( ( 1 - cs[,"sens"] ) / cs[,"spec"]  )
  negLRexact <- unname( ( 1-csExact[,"sens"] ) / csExact[,"spec"] )
  if( all( specb != 0L ) ) {
    negLR.ci <- bca( ( 1 - sensb) / specb, negLR, ... )$bca[4:5]
  } else {
    negLR.ci <- 1/bca( specb / ( 1 - sensb), 1/negLR, ... )$bca[c(4,5)]
  }
  posLR <- unname( cs[,"sens"] / ( 1 - cs[,"spec"] ) )
  posLRexact <- unname( csExact[,"sens"] / ( 1 - csExact[,"spec"] ) )
  if( all( specb != 1L ) ) {
    posLR.ci <- bca( sensb / ( 1 - specb ), posLR, ... )$bca[4:5]
  } else {
    posLR.ci <- 1/bca( ( 1 - specb ) / sensb, 1/posLR, ... )$bca[c(5,4)] # Reversed because the order inverts when you take the reciprocal
  }
  
  # -- Return lrtest object -- #
  structure( list(
    negLR = negLRexact,
    negLR.ci = negLR.ci,
    posLR = posLRexact,
    posLR.ci = posLR.ci,
    inputs = structure( c( truePos, totalPos, trueNeg, totalNeg ), names=c("truePos","totalPos","trueNeg","totalNeg") ),
    statistics = cs[ , c("sens","spec") ]
  ), 
  class = "lrtest",
  ci.type = "BCa",
  ci.width = .95
  )
}

#' Internal function to draw a set of sensitivities or specificities
#' This is intended for the case where testPos == totalPos or testNeg == totalNeg
#' @param n The total number of positives/negatives in the population
#' @param R is the number of replications in each round of the bootstrap (has been tested at 50,000 or greater)
#' @param verbose Whether to display internal operations as they happen
drawMaxedOut <- function( n, R, verbose ) {
  lprb <- sequentialGridSearch( # lowest probability that consistently produces 1's 
    f=identity, # We just want to minimize pr
    constraint=function(probs,...) vapply( probs, FUN=medianConsistentlyOne, FUN.VALUE=NA, ... ),
    bounds=c(0,1), 
    verbose=verbose,
    size=n, R=R, warn=FALSE,
    shrink=5,
    tol=.0005,
    nEach=80
  )
  res <- rbinom(R, size=n, prob=lprb)/n
  attr( res, "lprb" ) <- lprb
  res
}


#' Internal function to analyze LR bootstrap finding median, and standard and
#' BCa percentile 95% CIs
#' To obtain bca CI on a non-boot result, use a dummy boot
#' and replace t and t0 with the results of interest.
#' @param t The vector to obtain a BCa bootstrap for (e.g. nlr)
#' @param t0 The central value of the vector (e.g. the )
#' @param \dots Pass-alongs to boot.ci
bca <- function( t, t0, ... ) {
  R <- length(t)
  dummy <- rep(1:0,c(5,5)) # Doesn't matter what values are given here, since we're replacing them
  dummyb <- boot(dummy, function(x,i) 1, R=R)
  dummyb$t <- matrix(t,ncol=1)
  dummyb$t0 <- t0
  boot.ci(dummyb, t0=dummyb$t0, t=dummyb$t, type=c("perc", "bca"), ...)
}

# ----- Functions to display the resulting lrtest object ----- #

#' Prints results from the BayesianLR.test
#' @param x The lrtest object created by BayesianLR.test
#' @param \dots Pass-alongs (currently ignored)
#' @return Returns x unaltered
#' @method print lrtest
#' @S3method print lrtest
#' @export print.lrtest
print.lrtest <- function( x, ... ) {
  digits <- 3 # Number of digits to round to for display purposes
  cat("\n")
  cat("Likelihood ratio test of a 2x2 table")
  cat("\n\n")
  cat("data:\n")
  print(x$inputs)
  cat( paste0( "Positive LR: ", round(x$posLR,digits), " (", round(x$posLR.ci[1],digits), " - ", round(x$posLR.ci[2],digits), ")\n" ) )
  cat( paste0( "Negative LR: ", round(x$negLR,digits), " (", round(x$negLR.ci[1],digits), " - ", round(x$negLR.ci[2],digits), ")\n" ) )
  cat( paste0( attr(x,"ci.width")*100, "% confidence intervals computed via ", attr(x,"ci.type"), " bootstrapping.\n" ) )
  cat( "Note: This procedure depends on repeated random sampling.  As such it is subject to some variability in results.  Variability is minimized by large numbers of replications (generally 50,000) [and averaging 5 repeated results], but with small sample sizes or sensitivity or specificity near 0 or 1, variability becomes more pronounced.  This is not an error, it is a function of the nature of the procedure." )
  invisible(x)
}