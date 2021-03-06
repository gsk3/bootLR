% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/bootLR.R
\name{bca}
\alias{bca}
\title{Internal function to analyze LR bootstrap finding median, and standard and
BCa percentile 95% CIs.
To obtain bca CI on a non-boot result, use a dummy boot.
and replace t and t0 with the results of interest.}
\usage{
bca(t, t0, ...)
}
\arguments{
\item{t}{The vector to obtain a BCa bootstrap for (e.g. nlr).}

\item{t0}{The central value of the vector (e.g. the ).}

\item{\dots}{Pass-alongs to boot.ci.}
}
\description{
Internal function to analyze LR bootstrap finding median, and standard and
BCa percentile 95% CIs.
To obtain bca CI on a non-boot result, use a dummy boot.
and replace t and t0 with the results of interest.
}

