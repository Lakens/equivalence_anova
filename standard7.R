#################################################
#################################################

## Here are functions that give equivalence testing p-values that correspond to the CIs that are calculated in Kelly (2007).

library(MBESS)


# There are eight functions:

# 1. cohen_d()  corresponds to ci.sm()
# 2. cohen_d_twosample()  corresponds to ci.smd()
# 3. pop_sig_to_noise() coresponds to ci.snr()
# 4. pop_prop_of_var() corresponds to ci.pvaf()
# 5. sqrt_pop_sig_to_noise() corresponds to ci.srsnr()
# 6. std_targeted_eff() corresponds to ci.sc()  ### work needs to be done ###
# 7. equiv_R2() corresponds to ci.R2()
# 8. equiv_beta.k() correpsonds to ci.src()  ### work needs to be done ###


#################################################
#################################################

# 1. cohen_d()  corresponds to ci.sm()

## Equivalence (and non-inferiority) test for the 
## standardized mean difference (one group)
## (see Kelley (2007), end of section 2.2)

cohen_d<-function(sm, N, delta_lower, delta_upper){
#sm <- 0.3; N<-150; delta_lower<-(-0.5); delta_upper<-0.5
	
	my_sm <- sm
	my_N <- N

if((delta_upper<sm)|(delta_lower> sm)){	warning("sm outside of margin");return(1)} 
if(!((delta_upper<sm)|(delta_lower> sm))){	


CI<-function(x){unlist((ci.sm(sm=my_sm, N=my_N, conf.level=(1-2*x))[c(1,3)]))}

absmin<-function(x){x[which.min(abs(x))]}

delta_dist<-function(zzz){
			invisible(capture.output(CI_x <- CI(zzz)))
			LargeNumber<-10000
			dist<-LargeNumber*
			as.numeric(!(CI_x[1]>delta_lower &  CI_x[2]<delta_upper)) + 
			absmin(c((CI_x[2]-delta_upper),(delta_lower-CI_x[1])))
			
			return(dist)
		}

pval<-uniroot(delta_dist, lower=0,upper=1)$root
return(pval)
}}

########################################
## Example of using cohen_d():


# The standardized mean is 0.3, with N=43 subjects.  
# The equivalence margin is [-0.5, +0.5].
# Then the equivalence test p-value is:
equiv.pval <- cohen_d(sm=0.3, N=43, delta_lower=-0.5, delta_upper=+0.5)
equiv.pval

# A 90% two-sided confidence interval for this standardized mean does
# not fall entirely within the equiv. margin:
ci.sm(sm=0.3, N=43, conf.level=0.90)
# or:
ci.sm(sm=0.3, N=43, alpha.lower=0.05, alpha.upper=0.05)

# However, a (1-2*equiv.pval)% C.I. two-sided confidence interval
# falls exactly within the margin:
ci.sm(sm=0.3, N=43, conf.level=(1-2*equiv.pval))

## Small simulation study to investigate if there is
## correct type 1 error. (!!!takes 5 minutes!!!):
nSim <- 5000
simresults <- apply(cbind(1:nSim), 1, 
			function(x){
				if(round(x/21)==(x/21)){print(paste("% complete:", round(x/nSim,2)))}
				dat <- rnorm(1)*rnorm(100, 0.5+0.001)
				cohen_d(sm= mean(dat)/sd(dat), 
						N=length(dat), 
						delta_lower=-0.5, 
						delta_upper=+0.5)<0.05
					} 
				)

sum(simresults, na.rm=TRUE)/nSim
# 0.0486

#################################################
#################################################

# 2. cohen_d_twosample()  corresponds to ci.smd()

## Equivalence test for the standardized mean 
## difference for two independent groups
## (see Kelley (2007), section 3.1)

cohen_d_twosample<-function(smd, n.1, n.2, delta_lower, delta_upper){
	
	my_smd <- smd
	my_n.1 <- n.1
	my_n.2 <- n.2


if((delta_upper<smd)|(delta_lower> smd)){	warning("sm outside of margin"); return(1)} 
if(!((delta_upper<smd)|(delta_lower> smd))){	
	
CI<-function(x){unlist((ci.smd(smd=my_smd, n.1= my_n.1, n.2= my_n.2, conf.level=(1-2*x))[c(1,3)]))}


absmin<-function(x){x[which.min(abs(x))]}

delta_dist<-function(zzz){
			invisible(capture.output(CI_x <- CI(zzz)))
			LargeNumber<-10000
			dist<-LargeNumber*
			as.numeric(!(CI_x[1]>delta_lower &  CI_x[2]<delta_upper)) + 
			absmin(c((CI_x[2]-delta_upper),(delta_lower-CI_x[1])))
			
			return(dist)
		}

pval<-uniroot(delta_dist, lower=0,upper=1)$root

return(pval)
}}

########################################
## Example:

# The standardized mean difference is 0.3, 
# with two groups of size n1=43 and n2=183 subjects.  
# The equivalence margin is [-0.5, +0.5].
# Then the equivalence test p-value is:

equiv.pval <- cohen_d_twosample(smd=0.3, n.1=43, n.2=183, delta_lower=-0.5, delta_upper=0.5)
equiv.pval

# 90% C.I. two-sided confidence interval:
ci.smd(smd=0.3, n.1=43, n.2=183, conf.level=0.90)
# or:
ci.smd(smd=0.3, n.1=43, n.2=183, alpha.lower=0.05, alpha.upper=0.05)

# (1-2*equiv.pval)% C.I. two-sided confidence interval:
ci.smd(smd=0.3, n.1=43, n.2=183, conf.level=(1-2*equiv.pval))

## correct type 1 error? (!!!takes five minutes!!!):
nSim <- 5000
simresults <- apply(cbind(1:nSim),1, 
	function(x){
		if(round(x/21)==(x/21)){print(paste("% complete:",round(x/nSim,2)))}
		  n1<-50; n2<-100; mu1<-0; mu2<-0.5+0.001; sig<-1;
		  smp1 <- rnorm(n1, mu1, sig)
		  smp2 <- rnorm(n2, mu2, sig)
		  sig.hat <- sqrt( ((n1-1)*var(smp1)+(n2-1)*var(smp2))/(n1+n2-2) ) 
  
		  cohen_d_twosample(smd=(mean(smp2)-mean(smp1))/sig.hat, 
  					n.1=length(smp1), 
  					n.2=length(smp2), 
					delta_lower=-0.5, 
					delta_upper=+0.5)<0.05
				} 
			)

sum(simresults, na.rm=TRUE)/nSim
# 0.0492
#################################################

####################################
### Now we take a break in our work to 
### construct a function "simdata()" that will be used for
### simulating ANOVA style data with J=3 groups, with 
### a true population signal-to-noise ratio
### of 0.667 and a true population
### eta2 of 0.4 (see Kelly(2007), section 3.2) :

simdata<-function(nfactor=1, equalsizes=TRUE){

sigma2e <- 1

mu1 <- -1
mu2 <- 0
mu3 <- 1
mu <- 0

if(equalsizes){
n1 <- 10*nfactor
n2 <- 10*nfactor
n3 <- 10*nfactor}

if(!equalsizes){
n1 <- 10*nfactor
n2 <- 20*nfactor
n3 <- 30*nfactor}


N <- sum(c(n1,n2,n3))

tao1 <- mu1-mu
tao2 <- mu2-mu
tao3 <- mu3-mu

dat1 <- rnorm(n1,mu1,sqrt(sigma2e))
dat2 <- rnorm(n2,mu2,sqrt(sigma2e))
dat3 <- rnorm(n3,mu3,sqrt(sigma2e))

dat <- c(dat1,dat2,dat3)
groups <- as.factor(c(rep(1,n1),rep(2,n2),rep(3,n3)))

aovdata<-data.frame(dat=dat, groups= groups)

F<-summary(aov(dat ~ as.factor(groups), data= aovdata))[[1]]$F[1]

Lambda <-  sum(c(n1*tao1^2,n2*tao2^2,n3*tao3^2))/sigma2e
Lambda
# Alternatively, Equation 42 can be written as:

sigma2p<- sum(c(n1*tao1^2,n2*tao2^2,n3*tao3^2))/N

N*sigma2p/sigma2e

# The signal-to-noise ratio is defined as:

phi2 <- sigma2p/sigma2e
#print(paste("pop. signal-to-noise ratio =" ,round(phi2,3)))
# or:

Lambda/N

# and the proportion of variance in Y accounted for by 
# knowing the level of the factor (or group status in 
# a single factor design) is defined as:

eta2 <- Lambda/(Lambda+N)
eta2
return((aovdata))}

# Examples of using simdata() funciton: 
simdata(1)

simdata(2)

simdata(1, equalsizes=FALSE)


#################################################
#################################################
## EQUIVALENCE TESTS FOR STANDARDIZED OMNIBUS 
## EFFECTS IN AN ANOVA CONTEXT
## (see Kelley (2007), section 3.2)
#################################################
#################################################

# 3. pop_sig_to_noise() coresponds to ci.snr()

## Non-inferiority test for population signal-to-noise 
## ratio (i.e., phi_p^2) for the pth fixed effects 
## factor in an ANOVA setting

pop_sig_to_noise<-function(F.value, df.1, df.2, N, delta_upper){
# F.value=3.25; df.1=4; df.2=50; N=55; delta_upper=0.5;

	my_F <- F.value
	my_df.1 <- df.1
	my_df.2 <- df.2
	my_N <- N

# First find out what is the estimated signal-to-noise-ratio:
signal_to_noise_hat<-mean(unlist(ci.snr(F.value= my_F, df.1= my_df.1, df.2=my_df.2, N= my_N, alpha.lower=0.49, alpha.upper=0.49)))

#print(signal_to_noise_hat)
if(signal_to_noise_hat> delta_upper){print("sm outside of margin");return(1)}
if(signal_to_noise_hat<= delta_upper){


CI<-function(x){unlist((ci.snr(F.value=my_F, df.1=my_df.1, df.2=my_df.2, N=my_N, alpha.upper=x, alpha.lower=0)[c(1,2)]))}


absmin<-function(x){x[which.min(abs(x))]}

delta_dist<-function(zzz){
			invisible(capture.output(CI_x <- CI(zzz)))
			LargeNumber<-10000
			dist<-LargeNumber*as.numeric(!(CI_x[2]<delta_upper)) + (CI_x[2]-delta_upper)
			
			return(dist)
		}

qq<-0.99; while(is.na(CI(qq)[2])){qq<-qq-0.001}  # This is needed because of potential for ci.snr() to output NAs with large values of "alpha.upper"

pval<-uniroot(delta_dist, lower=0.0001,upper=qq)$root

return(pval)
}}



########################################
## Example:

# Non-inferiority test p-value (H0: phi_p^2>0.5 vs. H1: phi_p^2<=0.5)
noninf.pval <- pop_sig_to_noise(F.value=3.25, df.1=4, df.2=50, N=55, delta_upper=0.5)
noninf.pval

# 95% C.I. one-sided confidence interval:
ci.snr(F.value=3.25, df.1=4, df.2=50, N=55, alpha.lower=0.00, alpha.upper=0.05)

# (1-noninf.pval)% C.I. one-sided confidence interval:
ci.snr(F.value=3.25, df.1=4, df.2=50, N=55, alpha.lower=0.00, alpha.upper=noninf.pval)


## correct type 1 error? (!!could take a few minutes!!)
nSim <- 5000
simresults <- apply(cbind(1:nSim),1, 
	function(x){
				if(round(x/21)==(x/21)){print(paste("% complete:",round(x/nSim,2)))}
		 # simulate data from the null
		 # with true signal-to-noise
		 # equal to 0.667.
		 aovdata<-simdata() 
		ttt<-(summary(aov(aovdata$dat~aovdata$groups)))
		pop_sig_to_noise(F.value= ttt[[1]]$F[1], df.1= ttt[[1]]$Df[1], df.2= ttt[[1]]$Df[2], N=dim(aovdata)[1], delta_upper=0.66)<0.05
				} 
			)

sum(simresults, na.rm=TRUE)/nSim
# 0.0492


# Now an example from: https://rstudio-pubs-static.s3.amazonaws.com/251565_e4ac68427b134454a88b43fb02ea66ef.html
library(AMCP)
library(dplyr)
data(chapter_3_table_3)
chapter_3_table_3
chapter_3_table_3$Condition <- as.factor(chapter_3_table_3$Condition)

Summary.Mood <- 
chapter_3_table_3 %>% 
group_by(Condition) %>%
summarize(
Sum.Y = sum(Rating),
Y.bar = mean(Rating),
Sum.Squared.Errors = sum((Rating-mean(Rating))^2),
n = n())

Summary.Mood
(a <- length(unique(chapter_3_table_3$Condition)))
(N <- length(chapter_3_table_3$Rating))

(E.F <- sum(Summary.Mood$Sum.Squared.Errors))
(E.R <- sum((chapter_3_table_3$Rating - mean(chapter_3_table_3$Rating))^2))
(SS.Between <- E.R - E.F)
(MS.Between <- SS.Between/(a-1))
(MS.Within <- E.F/(N-a))
(df.F <- df.Within <- N-a)
(df.R <- N-1)
(df.Between <- df.R-df.F)
(F_val <- ((E.R - E.F)/(df.R-df.F))/(E.F/df.F))
# Or
((E.R - E.F)/df.Between)/(E.F/df.Within)
# Or, 
MS.Between/MS.Within
summary(aov(Rating ~ as.factor(Condition), data=chapter_3_table_3))


# 90% C.I. two-sided confidence interval for the population signal-to-noise ratio :
ci.snr(F.value=F_val, df.1= df.Between, df.2= df.F, N=N, conf.level=0.90)
# or 95% C.I. one-sided confidence interval for the population signal-to-noise ratio :
ci.snr(F.value=F_val, df.1= df.Between, df.2= df.F, N=N, alpha.lower=0, alpha.upper=0.05)

# Non-inferiority test p-value (H0: phi_p^2>3 vs. H1: phi_p^2<=3)
noninf.pval <- pop_sig_to_noise(F.value= F_val, df.1= df.Between, df.2= df.F, N=N, delta_upper=3)
noninf.pval

# (1-noninf.pval)% C.I. one-sided confidence interval:
ci.snr(F.value=F_val, df.1= df.Between, df.2= df.F, N=N, alpha.lower=0, alpha.upper=noninf.pval)

##########
# Now another example where we simulate a large dataset with
# true pop. signal-to-noise ratio is 0.66667
aovdata<-simdata(1000)
dat<-aovdata$dat
groups<-aovdata$groups
N<-dim(aovdata)[1]

ttt<-(summary(aov(dat~groups)))
 
# 90% two-sided confidence interval: 
ci.snr(F.value=ttt[[1]]$F[1], df.1= ttt[[1]]$Df[1], df.2= ttt[[1]]$Df[2], N=N, conf.level=0.90)

# 95% one-sided confidence interval: 
ci.snr(F.value=ttt[[1]]$F[1], df.1= ttt[[1]]$Df[1], df.2= ttt[[1]]$Df[2], N=N, alpha.lower=0, alpha.upper=0.05)

# Non-inferiority test p-value (H0: phi_p^2>0.65 vs. H1: phi_p^2<=0.65)
noninf.pval <- pop_sig_to_noise(F.value= ttt[[1]]$F[1], df.1= ttt[[1]]$Df[1], df.2= ttt[[1]]$Df[2], N=N, delta_upper=0.68)
noninf.pval


#################################################
#################################################

### !! END FOR NOW.  AFTER THIS POINT: WORK IN PROGRESS !!

#################################################
#################################################

# 4. pop_prop_of_var() corresponds to ci.pvaf()

## Equivalence test for the population proportion 
## of variance accounted for in the dependent 
## variable by knowing group status (i.e., eta_p^2) 
## for the pth fixed effects factor in an ANOVA setting.

pop_prop_of_var<-function(F.value, df.1, df.2, N, delta_upper=0.2, tol=0.001){
# F.value=3.25; df.1=4; df.2=50; N=55; delta_upper=0.5;	
	my_F.value <- F.value
	my_df.1 <- df.1
	my_df.2 <- df.2
	my_N <- N

prop_of_var_hat<-mean(unlist(ci.pvaf(F.value=my_F.value, df.1=my_df.1, df.2=my_df.2, N=my_N, alpha.upper=0.49, alpha.lower=0.49)[c(1,3)]))



#print(prop_of_var_hat)
if(prop_of_var_hat> delta_upper){print("sm outside of margin");return(1)}
if(prop_of_var_hat<= delta_upper){


CI<-function(x){unlist((	
ci.pvaf(F.value=my_F.value, df.1=my_df.1, df.2=my_df.2, N=my_N, alpha.upper=x, alpha.lower=0.00000001))[c(3)])
	}

delta_dist<-function(zzz){
			invisible(capture.output(CI_x <- CI(zzz)))
			LargeNumber<-10000
			dist<-LargeNumber*
			as.numeric(!(CI_x[1]<delta_upper)) + (CI_x[1]-delta_upper)
			
			return(dist)
		}

pval<-uniroot(delta_dist, lower=0.001,upper=0.999)$root

return(pval)}
}

## example:
##########
# true pop. eta_p^2 ratio is 0.4, equal sample sizes
aovdata<-simdata(300)
dat<-aovdata$dat
groups<-aovdata$groups
N<-dim(aovdata)[1]

ttt<-(summary(aov(dat~groups)))
 
# 90% two-sided confidence interval: 
ci.pvaf(F.value=ttt[[1]]$F[1], df.1= ttt[[1]]$Df[1], df.2= ttt[[1]]$Df[2], N=N, conf.level=0.90)

# 95% one-sided confidence interval (!!ERROR!!): 
ci.pvaf(F.value=ttt[[1]]$F[1], df.1= ttt[[1]]$Df[1], df.2= ttt[[1]]$Df[2], N=N, alpha.upper=0.05, alpha.lower=0.00)
# alternative?
ci.pvaf(F.value=ttt[[1]]$F[1], df.1= ttt[[1]]$Df[1], df.2= ttt[[1]]$Df[2], N=N, alpha.upper=0.05, alpha.lower=0.00000001)

# Non-inferiority test p-value (H0: phi_p^2>0.65 vs. H1: phi_p^2<=0.65)
noninf.pval <- pop_prop_of_var(F.value= ttt[[1]]$F[1], df.1= ttt[[1]]$Df[1], df.2= ttt[[1]]$Df[2], N=N, delta_upper=0.48)
noninf.pval

## correct type 1 error?
nSim <- 1000
simresults <- apply(cbind(1:nSim),1, 
	function(x){
				if(round(x/21)==(x/21)){print(paste("% complete:",round(x/nSim,2)))}
		 # simulate data from the null
		 # with true eta_p^2 ratio is
		 # equal to 0.4.
		 aovdata<-simdata(10) 
		ttt<-(summary(aov(aovdata$dat~aovdata$groups)))
		pval<-pop_prop_of_var(F.value= ttt[[1]]$F[1], df.1= ttt[[1]]$Df[1], df.2= ttt[[1]]$Df[2], N=dim(aovdata)[1], delta_upper=0.395)
		return(pval<0.05)
				} 
			)

sum(simresults, na.rm=TRUE)/nSim


#################################################
#################################################

# 5. sqrt_pop_sig_to_noise() corresponds to ci.srsnr()

## Non-inferiority test for the square root of 
## the signal-to-noise ratio for the pth 
## fixed effects factor (i.e., phi_p) in an 
## ANOVA setting.

sqrt_pop_sig_to_noise<-function(F.value, df.1, df.2, N, delta_upper=+0.2, tol=0.001){
	
	my_F.value <- F.value
	my_df.1 <- df.1
	my_df.2 <- df.2
	my_N <- N


CI<-function(x){unlist((ci.srsnr(F.value= my_F.value, df.1=my_df.1, df.2=my_df.2, N=my_N, alpha.upper=x, alpha.lower=0)[c(1,2)]))}

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
within <- CI_x[2]<delta_upper}
if(zzz>0.99){within<-TRUE}
}
return(zzz)
}

## example:


sqrt_pop_sig_to_noise(5.5, 5, 4, 100, 0.15,0.95)

## example:
##########
# true pop. phi_p is sqrt(0.6667)=0.816, equal sample sizes
aovdata<-simdata(300)
dat<-aovdata$dat
groups<-aovdata$groups
N<-dim(aovdata)[1]

ttt<-(summary(aov(dat~groups)))
 
# 90% two-sided confidence interval: 
ci.srsnr(F.value=ttt[[1]]$F[1], df.1= ttt[[1]]$Df[1], df.2= ttt[[1]]$Df[2], N=N, conf.level=0.90)

# 95% one-sided confidence interval: 
ci.srsnr(F.value=ttt[[1]]$F[1], df.1= ttt[[1]]$Df[1], df.2= ttt[[1]]$Df[2], N=N, alpha.upper=0.05, alpha.lower=0.00)

# Non-inferiority test p-value (H0: phi>0.85 vs. H1: phi<=0.85)
noninf.pval <- sqrt_pop_sig_to_noise(F.value= ttt[[1]]$F[1], df.1= ttt[[1]]$Df[1], df.2= ttt[[1]]$Df[2], N=N, delta_upper=0.85)
noninf.pval


## correct type 1 error?
nSim <- 1000
simresults <- apply(cbind(1:nSim),1, 
	function(x){
		if(round(x/21)==(x/21)){print(round(x/nSim,2))}
		 # simulate data from the null
		 # with true true pop. phi_p is sqrt(0.6667)=0.816
		 aovdata<-simdata(10) 
		ttt<-(summary(aov(aovdata$dat~aovdata$groups)))
		pval<-sqrt_pop_sig_to_noise(F.value= ttt[[1]]$F[1], df.1= ttt[[1]]$Df[1], df.2= ttt[[1]]$Df[2], N=dim(aovdata)[1], delta_upper= 0.82)
		return(pval<0.05)
				} 
			)

sum(simresults, na.rm=TRUE)/nSim


#################################################
#################################################

# 6. std_targeted_eff() corresponds to ci.sc()  ### work needs to be done ###

## Equivalence test for standardized targeted effects in ANOVA
## (for the population standardized comparison (i.e., phi) 
##  in an ANOVA setting)
## (see Kelley (2007), section 3.3)

std_targeted_eff<-function(means, s.anova, c.weights, n, N, my_n, delta_lower=0.5, delta_upper=+1.5, tol=0.001){
	
	my_means <- means
	my_s.anova <- s.anova
	my_c.weights <- c.weights
	my_n <- n
	my_N <- N


CI<-function(x){unlist(
	ci.sc(means=my_means, s.anova=my_s.anova, c.weights=my_c.weights, n=my_n, N=my_N, conf.level=(1-2*x))[c(1,3)]    )
	}

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
within <- CI_x[1]>delta_lower &  CI_x[2]<delta_upper
}
if(zzz>0.49){within<-TRUE}
}

return(zzz)
}

## example:
std_targeted_eff(means=c(2, 4, 9, 13), s.anova=.80, c.weights=c(.5, -.5, -.5, .5), n=c(30, 30, 30, 30), N=120, 0.15,0.95)


#################################################
#################################################

# 7. equiv_R2() corresponds to ci.R2()

## Equivalence test for the population 
## squared multiple correlation coeffcient (i.e., P2). 
## (see Kelley (2007), end of section 4.1)



equiv_R2<-function(R2, N, K, Random.Regressors=FALSE, delta_lower, delta_upper, one.sided=FALSE, tol=0.001){

if(one.sided==FALSE){
	if(delta_upper<R2){warning("R2 outside of margin"); return(NA)}
	if(delta_lower> R2){warning("R2 outside of margin"); return(NA)}
CI<-function(x){unlist(	
	ci.R2(R2=my_R2, N=my_N, K=my_K, conf.level=(1-2*x), Random.Regressors= my_Random.Regressors)[c(1,3)]    )}
}


if(one.sided=="upper"){
	if(delta_upper<R2){warning("R2 outside of margin"); return(NA)}
	CI<-function(x){
		unlist((ci.R2(R2=my_R2, N=my_N, K=my_K, alpha.upper=x, alpha.lower=0,Random.Regressors= my_Random.Regressors)[c(1,3)]))}}


if(one.sided=="lower"){
	if(delta_lower>R2){warning("R2 outside of margin"); return(NA)}
	CI<-function(x){
		unlist((ci.R2(R2=my_R2, N=my_N, K=my_K, alpha.upper=0, alpha.lower=x,Random.Regressors= my_Random.Regressors)[c(1,3)]))}}

	
	my_R2 <- R2
	my_N <- N
	my_K <- K
	my_Random.Regressors <- Random.Regressors

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
if(one.sided==FALSE){within <- CI_x[1]>delta_lower &  CI_x[2]<delta_upper}
if(one.sided=="upper"){within <- CI_x[2]<delta_upper}
if(one.sided=="lower"){within <-  CI_x[1]>delta_lower}
}	
}
return(zzz)
}





## simulate linear regression data with true population P2=0.28
simdatalm<-function(nfactor=1, Random.Regressors=TRUE){

if(Random.Regressors==FALSE){
x1 <- seq(5,90,length.out=10*nfactor)
x2 <- c(rep(c(0,1), 5*nfactor))
}

if(Random.Regressors==TRUE){
x1 <- runif(10*nfactor,5,90)
x2 <- rbinom(10*nfactor,1,0.5)
}


b0 <- 17
b1 <- 0.5
b2 <- 0.01

sigma <- 20

eps <- rnorm(10*nfactor,0,sigma)
y <- b0 + b1*x1  + b2*x2  + b3*x3 + eps
(var(y) - var(eps))/var(y)
lmdata<-data.frame(y,x1,x2)
return(lmdata)
}

equiv_R2(R2=0.13, N=30, K=4, delta_upper=0.2, delta_lower=0,Random.Regressors=TRUE)

lmdat<-simdatalm(100)
myR2<-(summary(lm(y~., data=lmdat)))$r.squared
myN<-dim(lmdat)[1]
myK<-dim(lmdat)[2]-1
equiv_R2(R2=myR2, N=myN, K=myK, Random.Regressors=FALSE, delta_upper=0.30)


 
## correct type 1 error? for Random.Regressors=TRUE
nSim <- 1000
simresults <- apply(cbind(1:nSim),1, 
	function(x){
		if(round(x/21)==(x/21)){print(round(x/nSim,2))}
		 # simulate data from the null
		 # with true true P2=0.28
		lmdat<-simdatalm(10, Random.Regressors=TRUE)
myR2<-(summary(lm(y~., data=lmdat)))$r.squared
myN<-dim(lmdat)[1]
myK<-dim(lmdat)[2]-1
pval<-equiv_R2(R2=myR2, N=myN, K=myK, Random.Regressors=TRUE, delta_upper=0.275, one.sided="upper")
	
		return(pval<0.05)
				} 
			)

sum(simresults, na.rm=TRUE)/nSim


## correct type 1 error? for Random.Regressors=FALSE
nSim <- 1000
simresults <- apply(cbind(1:nSim),1, 
	function(x){
		if(round(x/21)==(x/21)){print(round(x/nSim,2))}
		 # simulate data from the null
		 # with true true P2=0.28
		lmdat<-simdatalm(10, Random.Regressors=FALSE)
myR2<-(summary(lm(y~., data=lmdat)))$r.squared
myN<-dim(lmdat)[1]
myK<-dim(lmdat)[2]-1
pval<-equiv_R2(R2=myR2, N=myN, K=myK, Random.Regressors=FALSE, delta_upper=0.275, one.sided="upper")
	
		return(pval<0.05)
				} 
			)

sum(simresults, na.rm=TRUE)/nSim

## NOTE THAT type 1 error is incorrect when Random.Regressors does not match:


### TOO BIG :
nSim <- 1000
simresults <- apply(cbind(1:nSim),1, 
	function(x){
		if(round(x/21)==(x/21)){print(round(x/nSim,2))}
		 # simulate data from the null
		 # with true true P2=0.28
		lmdat<-simdatalm(10, Random.Regressors=TRUE)
myR2<-(summary(lm(y~., data=lmdat)))$r.squared
myN<-dim(lmdat)[1]
myK<-dim(lmdat)[2]-1
pval<-equiv_R2(R2=myR2, N=myN, K=myK, Random.Regressors=FALSE, delta_upper=0.275, one.sided="upper")
	
		return(pval<0.05)
				} 
			)

sum(simresults, na.rm=TRUE)/nSim

### OR too small :
nSim <- 1000
simresults <- apply(cbind(1:nSim),1, 
	function(x){
		if(round(x/21)==(x/21)){print(round(x/nSim,2))}
		 # simulate data from the null
		 # with true true P2=0.28
		lmdat<-simdatalm(10, Random.Regressors=FALSE)
myR2<-(summary(lm(y~., data=lmdat)))$r.squared
myN<-dim(lmdat)[1]
myK<-dim(lmdat)[2]-1
pval<-equiv_R2(R2=myR2, N=myN, K=myK, Random.Regressors=TRUE, delta_upper=0.275, one.sided="upper")
	
		return(pval<0.05)
				} 
			)

sum(simresults, na.rm=TRUE)/nSim



#################################################
#################################################

# 8. equiv_beta.k() correpsonds to ci.src()  ### work needs to be done ###

## Equivalence test for population 
## standardized regression coeffcient
## (see Kelley (2007), end of section 4.2)

equiv_beta.k<-function(beta.k, SE.beta.k , N, K, delta_lower=0, delta_upper=0.55, tol=0.001){
	
if(delta_upper<beta.k){stop("beta.k outside of margin")}
if(delta_lower>beta.k){stop("beta.k outside of margin")}
	
	my_beta.k <- beta.k
	my_SE.beta.k <-SE.beta.k
	my_N <- N
	my_K <- K

CI<-function(x){unlist(
ci.src(beta.k= my_beta.k, SE.beta.k= my_SE.beta.k, N= my_N, K= my_K, conf.level=(1-2*x)    ))[c(1,3)] 
	}

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
within <- CI_x[1]>delta_lower &  CI_x[2]<delta_upper
}
if(zzz>0.489){within<-TRUE}
}

return(zzz)
}

equiv_beta.k(beta.k=0.43, SE.beta.k=0.14, N=100, K=4)

