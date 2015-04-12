library(bootLR)

set.seed(3097959)

# How often does the bootLR 95% CI for the negative LR, LRneg, include the true population value?
# Nominal coverage would mean this occurs at least 95% of the time. 

# ss=sample size for those with disease and those without disease
# total sample = 2*ss 
# psn=population sensitivity
# psp=population specificity
# ns=number of samples drawn from population
SSs = c( 10, 50, 100 )
PSNs = PSPs = c( .6, .98, .99, .999, 1 )
ns = 100
params <- expand.grid( SSs, PSNs, PSPs )
colnames( params ) <- c( "ss", "psn", "psp" )


#' Define a function to calculate nominal CI from a single set of parameters
#' @param out Type of output desired, either "percent" or "raw" for the raw matrix of trials
calcNominalCI <- function( ss, psn, psp, ns, out = "percent" ) {
  
  # PopNLR=population neg LR ratio
  # PopPLR=population pos LR ratio
  PopNLR=(1-psn)/psp
  PopPLR=psn/(1-psp)
  
  # ssn=random sample sensitivities (vector provides a sensitivity result for each sample)
  # ssp=random sample specificity (vector provides a specificity result for each sample)
  
  ssn<-rbinom(ns, size=ss, prob=psn)/ss
  ssn
  
  ssp<-rbinom(ns, size=ss, prob=psp)/ss
  ssp
  
  # m is a matrix of ssn and ssp values for each sample with total ns samples
  m<-cbind(ssn,ssp)
  
  # Perform BootLR for each ssn and ssp pair with resulting Boot LR Negative LR 95% CI (and
  # Boot LR Positive LR 95% CI).  Total four numbers output for each sample.
  

  # Run BLRT on each drawn sens/spec
  res <- apply( m, 1, function(ssboth,ss) {
    ssn <- ssboth[1]
    ssp <- ssboth[2]
    cat("Sample Sensitivity=",ssn,"Specificity=",ssp,"\n")
    blrt <- try( BayesianLR.test( 
      truePos = ssn * ss, 
      totalDzPos = ss, 
      trueNeg = ssp * ss,
      totalDzNeg = ss
    ) )
    if( class(blrt) == "try-error" ) {
      rep( NA, 2 )
    } else {
      c( with( blrt, negLR.ci[1] < PopNLR & PopNLR < negLR.ci[2] ),
      with( blrt, posLR.ci[1] < PopPLR & PopPLR < posLR.ci[2] ) )
    }
  }, ss=ss )
  
  # Check if any didn't work
  if( any(is.na(res)) )  warning( "Number of runs which failed to converge:", sum(is.na(res)) )  
  # Output percentage of population LRs within sample LR CI
  if( out == "percent" ) {
  	apply( res, 1, function(x) sum( x, na.rm=TRUE ) / sum(!is.na( x ) ) )
  # Output instead the whole vector and we can add it up later (allows calculation of how many we need)
  } else if( out == "raw" ) {
  	res
  } else {
  	stop("Improper value of out")
  }
}

# Run it once to test
nominalCI <- calcNominalCI( ss=params[1,"ss"], psn=params[1,"psn"], psp=params[1,"psp"], ns=ns, out = "raw" )

#' Define function to calculate a percentage for the raw matrix
nominalCI_rawToPercent <- function( nominalCI ) {
	apply( nominalCI, 1, function(x) sum( x, na.rm=TRUE ) / sum(!is.na( x ) ) )
}
nominalCI_rawToPercent( nominalCI )
# How much variation is there for smaller samples?
nominalCI_rawToPercent( nominalCI[ , 1:500 ] )
nominalCI_rawToPercent( nominalCI[ , 501:1000 ] )


# Run it over every combination of population parameters
nominalCIs <- apply( params, 1, function( param, ns ) {
  cat("---------------\n")
  cat(param,"\n")
  try( calcNominalCI( ss=param["ss"], psn=param["psn"], psp=param["psp"], ns=ns, out = "percent" ) )
}, ns=ns )
