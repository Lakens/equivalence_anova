

############################################
#  start of noninf_Fstat  function
############################################

# example:
# aovdata<-simdata(20) 
# aovobject <-aov(aovdata$dat~aovdata$groups)
# noninf_Fstat(aovobject, 0.35)

noninf_Fstat <-  function(aovobject, eq_bound_eta=0.4, alpha = 0.05, within=FALSE, makeplot =FALSE, tol = 1e-09){
 
 

if(within==FALSE){ttt <- summary(aovobject)} 
if(within==TRUE){ttt <- summary(aovobject)$"Error: Within"}
# using the notation as in Okada (2013)

Fstat= ttt[[1]]$F[1]
dfb= ttt[[1]]$Df[1]
dfw= ttt[[1]]$Df[2]
SSb= (ttt[[1]])$"Sum Sq"[1]
SSw= (ttt[[1]])$"Sum Sq"[2]
SSt= SSb + SSw
MSb= (ttt[[1]])$"Mean Sq"[1]
MSw= (ttt[[1]])$"Mean Sq"[2]

if(within==FALSE){N= nobs((aovobject))}
if(within==TRUE){N= (dfw/(dfb+1))+1 }

# we use notation as in Okada (2013): 

eta2_hat =  SSb/SSt
epsilon2_hat = (SSb-dfb*MSw)/SSt
omega2_hat = (SSb-dfb*MSw)/(SSt+MSw)  # Note that there is a typo in eq.5 of Okada (2013)


# compare the three different estimates:
c(eta2_hat, epsilon2_hat, omega2_hat)

  #Then calculate one-sided (1-alpha)% CI for eta2.

  ncp_upper <- uniroot(function(q){ 
  
  								(pf(Fstat, 
                                       dfb, 
                                       dfw, 
                                       ncp = q) - (alpha))},

                       interval = c(0, 1000000),

                       tol = tol)$root

# Convert the CI bounds for the ncp to bounds for eta2:

  LL_CI <- 0
  UL_CI <- ncp_upper/(ncp_upper+N)

# calculate a p-value for the non-inferiority test:

      pval <- uniroot(function(alpha) (
    
    				pf(	Fstat, 
					dfb, 
				    dfw, 
					ncp = (eq_bound_eta * N) / (1 - eq_bound_eta)) - (alpha)),
                    interval = c(0,1),
                    tol = tol)$root


    # add plot if wanted:
if(makeplot==TRUE){
  ## TO DO: add plotting code
}
  
    if(eta2_hat > eq_bound_eta){
    print(paste("Warning: point estimate of", round(eta2_hat,3), "is above eq_bound_eta."))}

  
  if(pval==0){print(paste("warning: pval <=", tol,". Lower tol for additional accuracy."))}

	return(list(eta2_hat= eta2_hat, epsilon2_hat= epsilon2_hat, omega2_hat= omega2_hat, pval=pval, onesidedCI= c(LL_CI, UL_CI), Nobs=N))

}

############################################
#  end of noninf_Fstat  function
############################################


############################################
### Now a function "simdata()" that will be used for
### simulating ANOVA "across-subjects" data with J=3 groups, with 
### a true population signal-to-noise ratio
### of 0.667 and a true population
### eta2 of 0.4 (see Kelly(2007), section 3.2) :
### sigma2w = 1 (see Okada(2013))
### sigma2b = 0.666 (see Okada(2013))

simdata<-function(npergroup=100){

mu1 <- -1; mu2 <- 0; mu3 <- 1; mu <- 0

sigma2w = 1     #(see Okada(2013))
sigma2b = 2/3 #(see Okada(2013))

sigma2t = sigma2w + sigma2b

# We can write eta2 = sigma2b/sigma2t (see  Okada(2013))

eta2 <- sigma2b/sigma2t
eta2

# Or equivalently, eta2 = sigma2b/sigma2t (see  Okada(2013))

eta2 <- 1-sigma2w/sigma2t
eta2

n1 <- npergroup; n2 <- npergroup; n3 <- npergroup

N <- sum(c(n1,n2,n3))

tao1 <- mu1-mu; tao2 <- mu2-mu; tao3 <- mu3-mu

Lambda <-  sum(c(n1*tao1^2,n2*tao2^2,n3*tao3^2))/sigma2w

# Alternatively, Kelly(2007) Equation 42 can be written as:
sigma2p<- sum(c(n1*tao1^2,n2*tao2^2,n3*tao3^2))/N

# so that Lambda=
N*sigma2p/sigma2w

# The signal-to-noise ratio is defined as:
phi2 <- sigma2p/sigma2w
# or:
Lambda/N

# and the proportion of variance in Y accounted for by 
# knowing the level of the factor (or group status in 
# a single factor design) is defined as:
eta2 <- Lambda/(Lambda+N)
eta2

# simulate data:
dat1 <- rnorm(n1,mu1,sqrt(sigma2w))
dat2 <- rnorm(n2,mu2,sqrt(sigma2w))
dat3 <- rnorm(n3,mu3,sqrt(sigma2w))

dat <- c(dat1,dat2,dat3)
groups <- as.factor(c(rep(1,n1),rep(2,n2),rep(3,n3)))

aovdata<-data.frame(dat=dat, groups= groups)

return((aovdata))}

# Examples of using simdata() function: 
#simdata(10)  # n=30 observations
#simdata(20)  # n=60 observations




############################################
# Small simulation study (will take ~5min):
		 # We simulate data from a population
		 # where the true population eta_p^2 value is
		 # equal to 0.4 and the true population ncp = 40.
		 # We define the DELTA eq_bound_eta = 0.3999.
		 # This is a "between-subjects" design.
		
		 
nSim <- 50000
simresults <- apply(cbind(1:nSim),1, 
	function(x){
				if(round(x/37)==(x/37)){print(paste("% complete:",round(x/nSim,2)))}
		aovdata<-simdata(20) 
		aovobj<-aov(aovdata$dat~aovdata$groups)
		result<-noninf_Fstat(aovobj, eq_bound_eta=0.3999)
		return(unlist(result))
				} 
			)


## Is our point estimate, eta2_hat, unbiased for eta2? (true eta2=0.40)
mean(simresults["eta2_hat",], na.rm=TRUE)
median(simresults["eta2_hat",], na.rm=TRUE)
hist(simresults["eta2_hat",])

## Is our point estimate, epsilon2_hat, unbiased for eta2? (true eta2=0.40)
mean(simresults["epsilon2_hat",], na.rm=TRUE)
median(simresults["epsilon2_hat",], na.rm=TRUE)
hist(simresults["epsilon2_hat",])

## Is our point estimate, omega2_hat, unbiased for eta2? (true eta2=0.40)
mean(simresults["omega2_hat",], na.rm=TRUE)
median(simresults["omega2_hat",], na.rm=TRUE)
hist(simresults["omega2_hat",])

## In truth, all three estimators are biased, see Darlington, (1968)
## However, assymptotically (i.e. for very large n) they are all unbiased, see Maxell (1981).

## Does our test has correct type 1 error? (pvalues should be uniformly ditributed)
mean(simresults["pval",]<0.05, na.rm=TRUE)
median(simresults["pval",], na.rm=TRUE)
hist(simresults["pval",])




############################################
############################################

## We would also like to verify if our non-inferiority test for eta2
## works for a "within-subjects" design.

# We will simulate data from a population with the following parameters:
mu1=-1; mu2=0; mu3=1;
sigma2w <- 1
n1<-n2<-n3<-100
COR <- 0.85

# Before generating the data, what are the true population Lambda and pop. eta2 values?:
mu<-mean(c(mu1,mu2,mu3))
tao1 <- mu1-mu; tao2 <- mu2-mu; tao3 <- mu3-mu
(Lambda <-  sum(c(n1*tao1^2,n2*tao2^2,n3*tao3^2))/sigma2w)
(eta2 <- Lambda/(Lambda+N))

# Now here is a function to simulate data from the desired population;
# repeated measures within a subject are highly correlated (cor=0.85):
simdata_within<-function(n1){

mu1=-1; mu2=0; mu3=1;
sigma2w <- 1
n3<-n1
n2<-n1
COR <- 0.85

dat<-rmvnorm(n1, mean=c(mu1,mu2,mu3), sigma=matrix(c(sigma2w,COR,COR,COR, sigma2w, COR,COR, COR, sigma2w),3,3))

design_result <- data.frame(subject=rep(1:n1,3), y=c(dat), a=as.factor(sort(rep(1:3,n1)))) 
design_result}


## Test the function to create a dataset with 100 subjects:

dat<-simdata_within(100)

# alternatively, this could be done with:
#design_result <- ANOVA_design(string = "3w",
#                   n = n1, 
#                   mu = c(mu1, mu2, mu3), 
#                   sd = sigma2w, 
#                   r=0.85, 
#                   p_adjust = "none")


# Verify that the data is consistent with above parameters:
coldata<-cbind((dat[dat$a==1,]$y),(dat[dat$a==2,]$y),(dat[dat$a==3,]$y))
var(coldata)
colMeans(coldata)


# Run an ANOVA anlaysis:
aovobject<- aov(y~a +Error(subject), data= dat)

# Run a non-inferiority test with delta=0.45
# Note that N is the number of observations:
result<-noninf_Fstat(aovobject, eq_bound_eta=0.45, within=TRUE)
result

############################################
# Now a small simulation study (will take ~5min):
		 # We simulate data from a population
		 # where the true population eta_p^2 value is
		 # equal to 0.4 and the true population ncp = 40.
		 # We define the DELTA eq_bound_eta = 0.3999.
		 # This is a "within-subjects" design. 
nSim <- 1000

simresults <- apply(cbind(1:nSim),1, 
	function(x){
				if(round(x/37)==(x/37)){print(paste("% complete:",round(x/nSim,2)))}
			aovdata<-simdata_within(100)
 			   aovobj <-aov(y~a +Error(subject), data= aovdata)
				result<-noninf_Fstat(aovobj, eq_bound_eta=0.3999, within=TRUE)
				return(unlist(result))
				} 
			)


## Is our point estimate, eta2_hat, unbiased for eta2? (true eta2=0.40)
mean(simresults["eta2_hat",], na.rm=TRUE)
median(simresults["eta2_hat",], na.rm=TRUE)
hist(simresults["eta2_hat",])

## Is our point estimate, epsilon2_hat, unbiased for eta2? (true eta2=0.40)
mean(simresults["epsilon2_hat",], na.rm=TRUE)
median(simresults["epsilon2_hat",], na.rm=TRUE)
hist(simresults["epsilon2_hat",])

## Is our point estimate, omega2_hat, unbiased for eta2? (true eta2=0.40)
mean(simresults["omega2_hat",], na.rm=TRUE)
median(simresults["omega2_hat",], na.rm=TRUE)
hist(simresults["omega2_hat",])

## In truth, all three estimators are biased, see Darlington, (1968)
## However, assymptotically (i.e. for very large n) they are all unbiased, see Maxell (1981).

## Does our test has correct type 1 error? (pvalues should be uniformly ditributed)
mean(simresults["pval",]<0.05, na.rm=TRUE)
median(simresults["pval",], na.rm=TRUE)
hist(simresults["pval",])

