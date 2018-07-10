#################################################
#################################################

## These are functions that give equivalence testing p-values that correspond to the CIs that are calculated in Kelly (2007).

library(MBESS)

#################################################
#################################################

## Equivalence (and non-inferiority) test for the 
## standardized mean difference (one group)
## (see Kelley (2007), end of section 2.2)

cohen_d<-function(sm, N, delta_lower, delta_upper, one.sided=FALSE, tol=0.001){
	
	my_sm <- sm
	my_N <- N

if(one.sided==FALSE){
	if(delta_upper<sm){warning("sm outside of margin"); return(NA)}
	if(delta_lower> sm){warning("sm outside of margin"); return(NA)}
	CI<-function(x){unlist((ci.sm(sm=my_sm, N=my_N, conf.level=(1-2*x))[c(1,3)]))}}

if(one.sided=="upper"){
	if(delta_upper<sm){warning("sm outside of margin"); return(NA)}
	CI<-function(x){unlist((ci.sm(sm=my_sm, N=my_N, alpha.upper=x, alpha.lower=0)[c(1,3)]))}}

if(one.sided=="lower"){
	if(delta_lower> sm){warning("sm outside of margin"); return(NA)}
	CI<-function(x){unlist((ci.sm(sm=my_sm, N=my_N, alpha.lower=x, alpha.upper=0)[c(1,3)]))}}

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
if(one.sided==FALSE){within <- CI_x[1]>delta_lower &  CI_x[2]<delta_upper}
if(one.sided=="upper"){within <- CI_x[2]<delta_upper}
if(one.sided=="lower"){within <-  CI_x[1]>delta_lower}}
}
return(zzz)
}

########################################
## Example:

# Equivalence test p-value:
equiv.pval <- cohen_d(sm=0.3, N=43, delta_lower=-0.5, delta_upper=+0.5)
equiv.pval

# 90% C.I. two-sided confidence interval:
ci.sm(sm=0.3, N=43, conf.level=0.90)
# or:
ci.sm(sm=0.3, N=43, alpha.lower=0.05, alpha.upper=0.05)
# (1-2*equiv.pval)% C.I. two-sided confidence interval:
ci.sm(sm=0.3, N=43, conf.level=(1-2*equiv.pval))

# Non-inferiority test p-value (H0: sm>0.5 vs. H1: sm<=0.5)
noninf.pval <- cohen_d(sm=0.3, N=43, delta_upper=+0.5, one.sided="upper")
noninf.pval

# 95% C.I. one-sided confidence interval:
ci.sm(sm=0.3, N=43, alpha.lower=0.00, alpha.upper=0.05)

# (1-noninf.pval)% C.I. one-sided confidence interval:
ci.sm(sm=0.3, N=43, conf.level=(1-2*noninf.pval))


## correct type 1 error? (!!!takes 10 minutes!!!):
nSim <- 100
simresults <- apply(cbind(1:nSim), 1, 
			function(x){
				if(round(x/21)==(x/21)){print(round(x/nSim,2))}
				dat <- rnorm(1)*rnorm(100, 0.5+0.001)
				cohen_d(sm= mean(dat)/sd(dat), 
						N=length(dat), 
						delta_lower=-0.5, 
						delta_upper=+0.5)<0.05
					} 
				)

sum(simresults, na.rm=TRUE)/nSim


#################################################
#################################################

## Equivalence test for the standardized mean 
## difference for two independent groups
## (see Kelley (2007), section 3.1)

cohen_d_twosample<-function(smd, n.1, n.2, delta_lower, delta_upper, one.sided=FALSE, tol=0.001){
	
	my_smd <- smd
	my_n.1 <- n.1
	my_n.2 <- n.2

if(one.sided==FALSE){
	if(delta_upper<smd){warning("sm outside of margin"); return(NA)}
	if(delta_lower> smd){warning("sm outside of margin"); return(NA)}
	CI<-function(x){unlist((ci.smd(smd=my_smd, n.1= my_n.1, n.2= my_n.2, conf.level=(1-2*x))[c(1,3)]))}}

if(one.sided=="upper"){
	if(delta_upper<smd){warning("sm outside of margin"); return(NA)}
	CI<-function(x){unlist((ci.smd(smd=my_smd, n.1= my_n.1, n.2= my_n.2, alpha.upper=x, alpha.lower=0)[c(1,3)]))}}

if(one.sided=="lower"){
	if(delta_lower> smd){warning("sm outside of margin"); return(NA)}
	CI<-function(x){unlist((ci.smd(smd=my_smd, n.1= my_n.1, n.2= my_n.2, alpha.lower=x, alpha.upper=0)[c(1,3)]))}}

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
if(one.sided==FALSE){within <- CI_x[1]>delta_lower &  CI_x[2]<delta_upper}
if(one.sided=="upper"){within <- CI_x[2]<delta_upper}
if(one.sided=="lower"){within <-  CI_x[1]>delta_lower}}

}
return(zzz)
}

########################################
## Example:

# Equivalence test p-value:
equiv.pval <- cohen_d_twosample(smd=0.3, n.1=43, n.2=183, delta_lower=-0.5, delta_upper=0.5)
equiv.pval

# 90% C.I. two-sided confidence interval:
ci.smd(smd=0.3, n.1=43, n.2=183, conf.level=0.90)
# or:
ci.smd(smd=0.3, n.1=43, n.2=183, alpha.lower=0.05, alpha.upper=0.05)
# (1-2*equiv.pval)% C.I. two-sided confidence interval:
ci.smd(smd=0.3, n.1=43, n.2=183, conf.level=(1-2*equiv.pval))

# Non-inferiority test p-value (H0: sm>0.5 vs. H1: sm<=0.5)
noninf.pval <- cohen_d_twosample(smd=0.3, n.1=43, n.2=183, delta_upper=+0.5, one.sided="upper")
noninf.pval

# 95% C.I. one-sided confidence interval:
ci.smd(smd=0.3, n.1=43, n.2=183, alpha.lower=0.00, alpha.upper=0.05)

# (1-noninf.pval)% C.I. one-sided confidence interval:
ci.smd(smd=0.3, n.1=43, n.2=183, conf.level=(1-2*noninf.pval))

## correct type 1 error? (!!!takes five minutes!!!):
nSim <- 100
simresults <- apply(cbind(1:nSim),1, 
	function(x){
		if(round(x/21)==(x/21)){print(round(x/nSim,2))}
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
# 0.047

#################################################
#################################################
## EQUIVALENCE TESTS FOR STANDARDIZED OMNIBUS 
## EFFECTS IN AN ANOVA CONTEXT
## (see Kelley (2007), section 3.2)
#################################################
#################################################

## Non-equivalence test for population signal-to-noise 
## ratio (i.e., phi_p^2) for the pth fixed effects 
## factor in an ANOVA setting

pop_sig_to_noise<-function(F.value, df.1, df.2, N, delta_upper, tol=0.001){
	
	my_F <- F.value
	my_df.1 <- df.1
	my_df.2 <- df.2
	my_N <- N

CI<-function(x){unlist((ci.snr(F.value=my_F, df.1=my_df.1, df.2=my_df.2, N=my_N, alpha.upper=x, alpha.lower=0)[c(1,2)]))}

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
within <- CI_x[2]<delta_upper}
}
return(zzz)
}

########################################
## Example:

# Non-inferiority test p-value (H0: phi_p^2>0.5 vs. H1: phi_p^2<=0.5)
noninf.pval <- pop_sig_to_noise(F.value=3.25, df.1=4, df.2=50, N=55, delta_upper=0.5)
noninf.pval

# 95% C.I. one-sided confidence interval:
ci.snr(F.value=3.25, df.1=4, df.2=50, N=55, alpha.lower=0.00, alpha.upper=0.05)

# (1-noninf.pval)% C.I. one-sided confidence interval:
ci.snr(F.value=3.25, df.1=4, df.2=50, N=55, alpha.lower=0.00, alpha.upper=noninf.pval)


## correct type 1 error? (!!!takes five minutes!!!):
nSim <- 100
simresults <- apply(cbind(1:nSim),1, 
	function(x){
		if(round(x/21)==(x/21)){print(round(x/nSim,2))}
		 
		 
		 			pop_sig_to_noise(F.value= , df.1= , df.2= , N= , delta_upper=0.5)<0.05
				} 
			)

sum(simresults, na.rm=TRUE)/nSim


# Example from: https://rstudio-pubs-static.s3.amazonaws.com/251565_e4ac68427b134454a88b43fb02ea66ef.html
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
# 95% C.I. one-sided confidence interval for the population signal-to-noise ratio :
ci.snr(F.value=F_val, df.1= df.Between, df.2= df.F, N=N, alpha.lower=0, alpha.upper=0.05)

# Non-inferiority test p-value (H0: phi_p^2>3 vs. H1: phi_p^2<=3)
noninf.pval <- pop_sig_to_noise(F.value= F_val, df.1= df.Between, df.2= df.F, N=N, delta_upper=3)
noninf.pval

# (1-noninf.pval)% C.I. one-sided confidence interval:
ci.snr(F.value=F_val, df.1= df.Between, df.2= df.F, N=N, alpha.lower=0, alpha.upper=noninf.pval)




#################################################
#################################################

## Equivalence test for the population proportion 
## of variance accounted for in the dependent 
## variable by knowing group status (i.e., eta_p^2) 
## for the pth fixed effects factor in an ANOVA setting.

pop_prop_of_var<-function(F.value, df.1, df.2, N, delta_lower=0, delta_upper=+0.2, tol=0.001){
	
	my_F.value <- F.value
	my_df.1 <- df.1
	my_df.2 <- df.2
	my_N <- N


CI<-function(x){unlist((	
ci.pvaf(F.value=my_F.value, df.1=my_df.1, df.2=my_df.2, N=my_N, conf.level=(1-2*x))[c(1,3)]))
	}

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
within <- CI_x[1]>delta_lower &  CI_x[2]<delta_upper}
if(zzz>0.49){within<-TRUE}
}
return(zzz)
}

## example:

pop_prop_of_var(2.5, 5, 4, 100, 0.0,0.5)



#################################################
#################################################

## Equivalence test for the square root of 
## the signal-to-noise ratio for the pth 
## fixed effects factor (i.e., phi_p) in an 
## ANOVA setting.

sqrt_pop_sig_to_noise<-function(F.value, df.1, df.2, N, delta_lower=-0.2, delta_upper=+0.2, tol=0.001){
	
	my_F.value <- F.value
	my_df.1 <- df.1
	my_df.2 <- df.2
	my_N <- N


CI<-function(x){unlist((ci.srsnr(F.value= my_F.value, df.1=my_df.1, df.2=my_df.2, N=my_N, conf.level=(1-2*x))[c(1,2)]))}

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
within <- CI_x[1]>delta_lower &  CI_x[2]<delta_upper}
if(zzz>0.49){within<-TRUE}
}
return(zzz)
}

## example:


sqrt_pop_sig_to_noise(5.5, 5, 4, 100, 0.15,0.95)


#################################################
#################################################

## Equivalence test for standardized targeted effects in ANOVA
## (for the population standardized comparison (i.e., phi) 
##  in an ANOVA setting)
## (see Kelley (2007), section 3.3)

std_tageted_eff<-function(means, s.anova, c.weights, n, N, my_n, delta_lower=0.5, delta_upper=+1.5, tol=0.001){
	
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
std_tageted_eff(means=c(2, 4, 9, 13), s.anova=.80, c.weights=c(.5, -.5, -.5, .5), n=c(30, 30, 30, 30), N=120, 0.15,0.95)


#################################################
#################################################

## Equivalence test for the population 
## squared multiple correlation coeffcient (i.e., P2). 
## (see Kelley (2007), end of section 4.1)



equiv_R2<-function(R2, N, K, Random.Regressors=FALSE, delta_lower=0, delta_upper=0.25, tol=0.001){

if(delta_upper<R2){stop("R2 outside of margin")}	
if(delta_lower>R2){stop("R2 outside of margin")}	
	
	my_R2 <- R2
	my_N <- N
	my_K <- K
	my_Random.Regressors <- Random.Regressors

CI<-function(x){unlist(
	
	ci.R2(R2=my_R2, N=my_N, K=my_K, conf.level=(1-2*x), Random.Regressors= my_Random.Regressors)[c(1,3)]    )
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

equiv_R2(R2=0.13, N=30, K=4, Random.Regressors=FALSE)


#################################################
#################################################

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

