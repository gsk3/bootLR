\name{confusionStatistics}
\alias{confusionStatistics}
\title{Compute sensitivity, specificity, positive likelihood ratio, negative likelihood ratio for a single 2x2 table}
\usage{
  confusionStatistics(truePos, totalPos, trueNeg, totalNeg)
}
\arguments{
  \item{truePos}{The number of true positive tests}

  \item{totalPos}{The total number of positives ("sick") in
  the population}

  \item{trueNeg}{The number of true negatives in the
  population}

  \item{totalNeg}{The total number of negatives ("well") in
  the population}
}
\value{
  A matrix containing sensitivity, specificity, posLR,
  negLR results
}
\description{
  Compute sensitivity, specificity, positive likelihood
  ratio, negative likelihood ratio for a single 2x2 table
}
\references{
  Deeks JJ, Altman DG. BMJ. 2004 July 17; 329(7458):
  168-169.
}

