
############################################

#  start of noninf_Fstat  function

############################################


noninf_Fstat <-  function(aovobject, eq_bound_eta2=0.4, alpha = 0.05, makeplot =FALSE, tol = 1e-09){


ttt <- summary(aovobject)

# using the notation as in Okada (2013)

Fstat= ttt[[1]]$F[1]
pval <- (ttt[[1]]$"Pr(>F)")[1]
dfb= ttt[[1]]$Df[1]
dfw= ttt[[1]]$Df[2]
SSb= (ttt[[1]])$"Sum Sq"[1]
SSw= (ttt[[1]])$"Sum Sq"[2]
SSt= SSb + SSw
MSb= (ttt[[1]])$"Mean Sq"[1]
MSw= (ttt[[1]])$"Mean Sq"[2]
N= nobs((aovobject))
m <- length(aovobject$coefficients)  # number of groups

## We should have that these are equal:
Fstat
MSb/MSw

pval
1-pf(Fstat, dfb,  dfw)

dfb
m-1

dfw
(N/m-1)*(m)		


# point estimates for eta2, using notation as in Okada (2013): 
eta2_hat =  SSb/SSt
epsilon2_hat = (SSb-dfb*MSw)/SSt
omega2_hat = (SSb-dfb*MSw)/(SSt+MSw)  # (Note that there is a typo in eq.5 of Okada (2013))

# compare the three different point estimates for eta2:
c(eta2_hat, epsilon2_hat, omega2_hat)

# calculate a point estimate for the non-centrality parameter Lambda.
# still need a source for this ???
ncp_hat <- (N-1)*SSb/SSw

# calculate a point estimate for the RMSSE parameter .
# still need a source for this ???
RMSSE_hat <- sqrt(ncp_hat/((m-1)*(N/m)))


# We calculate one-sided (1-alpha)% CI for the non-centrality parameter Lambda.
  ncp_upper <- uniroot(function(q){ 
  								(pf(Fstat, 
                                       dfb, 
                                       dfw, 
                                       ncp = q) - (alpha))},
                       interval = c(0, 1000000),
                       tol = tol)$root


# Convert the CI bounds for the ncp to CI bounds for eta2:
eta2_CI <- ncp_upper/(ncp_upper+N)


# Convert the CI bounds for the ncp to CI bounds for RMSSE (see Steiger2004, eq. 12):
RMSSE_CI <- sqrt(ncp_upper/((m-1)*(N/m)))


# And calculate a p-value for the non-inferiority test, H0: eta2> eq_bound_eta2:

      pval <- uniroot(function(alpha) (
 
    				pf(	Fstat, 
					dfb, 
				    dfw, 
					ncp = (eq_bound_eta2 * N) / (1 - eq_bound_eta2)) - (alpha)),
	                  interval = c(0,1),
                    tol = tol)$root

# add plot if wanted:

if(makeplot==TRUE){

  ## TO DO: add plotting code

}

# print warnings:
    if(eta2_hat > eq_bound_eta2){
    print(paste("Warning: point estimate of", round(eta2_hat,3), "is above eq_bound_eta2."))}

   if(pval==0){print(paste("warning: pval <=", tol,". Lower tol for additional accuracy."))}

	return(list(onesidedCI_ncp= c(0,ncp_upper),  onsidedCI_RMSSE=c(0,RMSSE_CI),  eta2_hat= eta2_hat, epsilon2_hat= epsilon2_hat, omega2_hat= omega2_hat, onesidedCI_eta2= c(0, eta2_CI), Nobs=N, pval=pval, eq_bound_eta2= eq_bound_eta2))

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

p <- 3 # number of groups
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


# we also specify:
n1 <- npergroup; n2 <- npergroup; n3 <- npergroup
N <- sum(c(n1,n2,n3))
tao1 <- mu1-mu; tao2 <- mu2-mu; tao3 <- mu3-mu

## We have the non-centrality parameter:
Lambda <-  sum(c(n1*tao1^2,n2*tao2^2,n3*tao3^2))/sigma2w

# or as in Steiger2004 equation (7) 
alpha1<-tao1;alpha2<-tao2;alpha3<-tao3;
Lambda <-  npergroup*sum(c( (alpha1/sqrt(sigma2w))^2, (alpha2/sqrt(sigma2w))^2, (alpha3/sqrt(sigma2w))^2))


# Alternatively, Kelly(2007) Equation 42 can be written as:
sigma2p<- sum(c(n1*tao1^2,n2*tao2^2,n3*tao3^2))/N

# so that Lambda=
N*sigma2p/sigma2w


# The signal-to-noise ratio is defined as:
phi2 <- sigma2p/sigma2w

# or as in Steiger2004 equation (7):
Lambda/N

# and the proportion of variance in Y accounted for by 
# knowing the level of the factor (or group status in 
# a single factor design) is defined as:

eta2 <- Lambda/(Lambda+N)
eta2

# Finally the RMSSE, as in Steiger2004 eq. 12:
RMSSE=sqrt(Lambda/((p-1)*npergroup))

# simulate data:
dat1 <- rnorm(n1,mu1,sqrt(sigma2w))
dat2 <- rnorm(n2,mu2,sqrt(sigma2w))
dat3 <- rnorm(n3,mu3,sqrt(sigma2w))

dat <- c(dat1,dat2,dat3)

groups <- as.factor(paste("g",c(rep(1,n1),rep(2,n2),rep(3,n3)),sep=""))
aovdata<-data.frame(dat=dat, groups= groups)

print(c(paste("Lambda=",Lambda, "eta2=",eta2,"phi2=", phi2, "RMSSE=",RMSSE )))

return((aovdata))}


############################################

# Small simulation study (will take ~5min):

		 # We simulate data from a population
		 # where the true population eta_p^2 value is
		 # equal to 0.4 and the true population ncp = 40.
		 # We define the DELTA eq_bound_eta2 = 0.3999.
		 # This is a "between-subjects" design.


nSim <- 10000

simresults <- apply(cbind(1:nSim),1, 
	function(x){
		if(round(x/37)==(x/37)){print(paste("% complete:",round(x/nSim,2)))}
		aovdata<-simdata(20) 
		aovobj<-aov(aovdata$dat~aovdata$groups)
		result<-noninf_Fstat(aovobj, eq_bound_eta2=0.3999)
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

# How about coverage for the one-sided 95% CIs:

# how about the coverage for the onesidedCI_ncp?
mean(simresults["onesidedCI_ncp2",]>40, na.rm=TRUE)
# how about the coverage for the onesidedCI_RMSSE?
mean(simresults["onsidedCI_RMSSE2",]>1, na.rm=TRUE)
# how about the coverage for the onesidedCI_eta2?
mean(simresults["onesidedCI_eta22",]>0.40, na.rm=TRUE)



############################################
# Everything looks good.
############################################


