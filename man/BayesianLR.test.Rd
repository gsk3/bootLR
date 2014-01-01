\name{BayesianLR.test}
\alias{BayesianLR.test}
\title{Compute the (negative) likelihood ratio with appropriate, bootstrapped confidence intervals}
\usage{
  BayesianLR.test(truePos, totalPos, trueNeg, totalNeg,
    R = 5 * 10^4, verbose = FALSE, ...)
}
\arguments{
  \item{truePos}{The number of true positive tests}

  \item{totalPos}{The total number of positives ("sick") in
  the population}

  \item{trueNeg}{The number of true negatives in the
  population}

  \item{totalNeg}{The total number of negatives ("well") in
  the population}

  \item{R}{is the number of replications in each round of
  the bootstrap (has been tested at 50,000 or greater)}

  \item{verbose}{Whether to display internal operations as
  they happen}

  \item{\dots}{Arguments to pass along to boot.ci for the
  BCa confidence intervals}
}
\value{
  An object of class lrtest
}
\description{
  Compute the (negative) likelihood ratio with appropriate,
  bootstrapped confidence intervals
}
\examples{
blrt <- BayesianLR.test( truePos=100, totalPos=100, trueNeg=60, totalNeg=100 )
blrt
summary(blrt)
BayesianLR.test( truePos=98, totalPos=100, trueNeg=60, totalNeg=100 )
BayesianLR.test( truePos=60, totalPos=100, trueNeg=100, totalNeg=100 )
BayesianLR.test( truePos=60, totalPos=100, trueNeg=99, totalNeg=100 )
}

