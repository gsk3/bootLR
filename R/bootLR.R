

#' Find the lowest population probability whose median is consistently one
#' This is the lowest estimate for Sens that is consistently (over 5 runs) most likely to yield a sample estimate of 100/100.
#' @param pr Probability input
#' @param size Number of trials
#' @param R number of bootstrap replications
#' @param nConsistentRuns Number of runs that all have to be identical to return TRUE
#' @param warn Warn if searching outside of the range c(0,1)
#' @return Boolean of length one (TRUE or FALSE)
#' @example
#' prs <- seq(.990,.995,.0001)
#' bools <- sapply( prs, medianConsistentlyOne, size=truePos, R=R )
#' data.frame( prs, bools )
medianConsistentlyOne <- function(pr, size, R, nConsistentRuns=5, warn=TRUE) {
  if( 0 > pr | pr > 1) {
    warning("Searching probabilities outside of 0,1. Returning FALSE.")
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

#' Compute the (negative) likelihood ratio with appropriate, bootstrapped confidence intervals
#' @param truePos The number of true positive tests
#' @param totalPos The total number of positives ("sick") in the population
#' @param trueNeg The number of true negatives in the population
#' @param totalNeg The total number of negatives ("well") in the population
#' @param R is the number of replications in each round of the bootstrap
#' @param verbose Whether to display internal operations as they happen
#' @return An object of class nlrtest
#' @export BayesianLR.test
#' @example
#' truePos <- 100
#' totalPos <- 100
#' trueNeg <- 60
#' totalNeg <- 100
#' nlr <- BayesianLR.test( truePos=truePos, totalPos=totalPos, trueNeg=trueNeg, totalNeg=totalNeg )
#' nlr
#' summary(nlr)
BayesianLR.test <- function( truePos, totalPos, trueNeg, totalNeg, R=5*10^4, verbose=FALSE ) {
  # -- Check inputs -- #
  if( R!=5*10^4 ) warning("Setting the number of bootstrap replications to a number lower than 50,000 may lead to unstable results")
  
  # -- Determine version of algorithm to run -- #
  
  
  # -- Computations -- #
  
  # - Find lowest probability that consistently produces 1's - #
  
  lprb <- sequentialGridSearch( 
    f=identity, # We just want to minimize pr
    constraint=function(probs,...) vapply( probs, FUN=medianConsistentlyOne, FUN.VALUE=NA, ... ),
    bounds=c(0,1), 
    verbose=verbose,
    size=totalPos, R=R, warn=FALSE
  )
  
  
  
}



