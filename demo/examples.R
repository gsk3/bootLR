# Example: Producing the 95% CI for a random binomial sample of size 100 with probability
# of one set at 0.6.

p<-rbinom(10000, size=100, prob=.6)/100
q<-quantile(p, c(.025,.975))
q



#Calculate LR upper 95% CI for Sens 100/100 and Spec 60/100

#First, find the lowest population probability whose median is consistently one. 
#This is the lowest estimate for Sens that is consistently most likely to yield a
#sample estimate of 100/100.
#Repeat computation below starting from lower number and going up until median is consistently 1 with 5
#consecutive runs.

p<-rbinom(50000, size=100, prob=.9933)/100
q<-quantile(p, c(.025,.5,.975))
q

#For sample size 100 patients with disease, answer is 0.9933.  Use this value to generate bootstrap
#samples below:



#Bootstrap Neg LR (1-Sens/Spec) where sens=100/100 and spec=60/100

sens<-rbinom(50000, size=100, prob=.9933)/100
spec<-rep(1:0,c(60,40))
specf<- function(spec,i) {return (1/mean(spec[i]))}
specb<- boot(spec, specf, R=50000)
lr<-((1-sens)*specb$t)
#Next four lines are FYI:
#qq<-quantile(sens, c(.025,.5,.975))
#qq
#qqq<-quantile(specb$t, c(.025,.5,.975))
#qqq
q<-quantile(lr, c(.025,.5,.975))
q



# Analyze LR bootstrap finding median, and standard and
# BCa percentile 95% CIs
# To obtain bca CI on a non-boot result, use a dummy boot
# and replace t and t0 with the results of interest.
dummy<-rep(1:0,c(6,4))
dummyf<- function(dummy,i) {return (mean(dummy[i]))}
dummyb<- boot(dummy, dummyf, R=50000)
b$t<-matrix(lr,nrow=50000,byrow=T)
#(1-.9933)/.6 = .0112
b$t0<-.0112
boot.ci(dummyb, t0=b$t0, t=b$t, type=c("perc", "bca"))




#Calculate LR upper 95% CI for Sens 99/100 and Spec 60/100

sens<-rep(1:0,c(99,1))
spec<-rep(1:0,c(60,40))
sensf<- function(sens,i) {return (1 - mean(sens[i]))}
sensb<- boot(sens, sensf, R=10000)
specf<- function(spec,i) {return (1/mean(spec[i]))}
specb<- boot(spec, specf, R=10000)
lr<-(sensb$t*specb$t)
qq<-quantile(sensb$t, c(.025,.5,.975))
qq
qqq<-quantile(specb$t, c(.025,.5,.975))
qqq
q<-quantile(lr, c(.025,.5,.975))
q

# Analyze LR bootstrap finding median, and standard and
# BCa percentile 95% CIs
# To obtain bca CI on a non-boot result, use a dummy boot
# and replace t and t0 with the results of interest.
dummy<-rep(1:0,c(6,4))
dummyf<- function(dummy,i) {return (mean(dummy[i]))}
dummyb<- boot(dummy, dummyf, R=10000)
b$t<-matrix(lr,nrow=10000,byrow=T)
#.01/.6 = .0167
b$t0<-.0167
boot.ci(dummyb, t0=b$t0, t=b$t, type=c("perc", "bca"))

#Note, to see histogram
hist(sensb$t, breaks=20)



#Calculate LR upper 95% CI for Sens 98/100 and Spec 60/100

sens<-rep(1:0,c(98,2))
spec<-rep(1:0,c(60,40))
sensf<- function(sens,i) {return (1 - mean(sens[i]))}
sensb<- boot(sens, sensf, R=10000)
specf<- function(spec,i) {return (1/mean(spec[i]))}
specb<- boot(spec, specf, R=10000)
lr<-(sensb$t*specb$t)
qq<-quantile(sensb$t, c(.025,.5,.975))
qq
qqq<-quantile(specb$t, c(.025,.5,.975))
qqq
q<-quantile(lr, c(.025,.5,.975))
q

# Analyze LR bootstrap finding median, and standard and
# BCa percentile 95% CIs
# To obtain bca CI on a non-boot result, use a dummy boot
# and replace t and t0 with the results of interest.
dummy<-rep(1:0,c(6,4))
dummyf<- function(dummy,i) {return (mean(dummy[i]))}
dummyb<- boot(dummy, dummyf, R=10000)
b$t<-matrix(lr,nrow=10000,byrow=T)
#.02/.6 = .033
b$t0<-.033
boot.ci(dummyb, t0=b$t0, t=b$t, type=c("perc", "bca"))



#Calculate LR upper 95% CI for Sens 97/100 and Spec 60/100

sens<-rep(1:0,c(97,3))
spec<-rep(1:0,c(60,40))
sensf<- function(sens,i) {return (1 - mean(sens[i]))}
sensb<- boot(sens, sensf, R=10000)
specf<- function(spec,i) {return (1/mean(spec[i]))}
specb<- boot(spec, specf, R=10000)
lr<-(sensb$t*specb$t)
qq<-quantile(sensb$t, c(.025,.5,.975))
qq
qqq<-quantile(specb$t, c(.025,.5,.975))
qqq
q<-quantile(lr, c(.025,.5,.975))
q

# Analyze LR bootstrap finding median, and standard and
# BCa percentile 95% CIs
# To obtain bca CI on a non-boot result, use a dummy boot
# and replace t and t0 with the results of interest.
dummy<-rep(1:0,c(6,4))
dummyf<- function(dummy,i) {return (mean(dummy[i]))}
dummyb<- boot(dummy, dummyf, R=10000)
b$t<-matrix(lr,nrow=10000,byrow=T)
#.03/.6 = .05
b$t0<-.05
boot.ci(dummyb, t0=b$t0, t=b$t, type=c("perc", "bca"))




#Calculate LR upper 95% CI for Sens 96/100 and Spec 60/100

sens<-rep(1:0,c(96,4))
spec<-rep(1:0,c(60,40))
sensf<- function(sens,i) {return (1 - mean(sens[i]))}
sensb<- boot(sens, sensf, R=10000)
specf<- function(spec,i) {return (1/mean(spec[i]))}
specb<- boot(spec, specf, R=10000)
lr<-(sensb$t*specb$t)
qq<-quantile(sensb$t, c(.025,.5,.975))
qq
qqq<-quantile(specb$t, c(.025,.5,.975))
qqq
q<-quantile(lr, c(.025,.5,.975))
q

# Analyze LR bootstrap finding median, and standard and
# BCa percentile 95% CIs
# To obtain bca CI on a non-boot result, use a dummy boot
# and replace t and t0 with the results of interest.
dummy<-rep(1:0,c(6,4))
dummyf<- function(dummy,i) {return (mean(dummy[i]))}
dummyb<- boot(dummy, dummyf, R=10000)
b$t<-matrix(lr,nrow=10000,byrow=T)
#.04/.6 = .0667
b$t0<-.0667
boot.ci(dummyb, t0=b$t0, t=b$t, type=c("perc", "bca"))






