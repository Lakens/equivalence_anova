library(MBESS)

#################################################
#################################################

## Equivalence test for the standardized 
## mean difference (one group)
## (see Kelley (2007), end of section 2.2)

cohen_d<-function(sm, N, delta_lower=-0.2, delta_upper=+0.2, tol=0.001){
	
	my_sm <- sm
	my_N <- N

CI<-function(x){unlist((ci.sm(sm=my_sm, N=my_N, conf.level=(1-2*x))[c(1,3)]))}

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
within <- CI_x[1]>delta_lower &  CI_x[2]<delta_upper}
if(zzz>0.999){within<-TRUE}
}
return(zzz)
}

## example:
cohen_d(sm=0.3, N=43, delta_lower=-0.5, delta_upper=+0.5)



#################################################
#################################################

## Equivalence test for the standardized mean 
## difference for two independent groups
## (see Kelley (2007), section 3.1)

cohen_d_twosample<-function(smd, N1, N2, delta_lower=-0.2, delta_upper=+0.2, tol=0.001){
	
	my_smd <- smd
	my_N1 <- N1
	my_N2 <- N2

CI<-function(x){unlist((ci.smd(smd=my_smd, n.1= my_N1, n.2= my_N2,conf.level=(1-2*x))[c(1,3)]))}

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
within <- CI_x[1]>delta_lower &  CI_x[2]<delta_upper}
if(zzz>0.999){within<-TRUE}
}
return(zzz)
}


## example:

cohen_d_twosample(0.05,35, 45, -0.5,0.5)

#################################################
#################################################
## EQUIVALENCE TESTS FOR STANDARDIZED OMNIBUS 
## EFFECTS IN AN ANOVA CONTEXT
## (see Kelley (2007), section 3.2)
#################################################
#################################################

## Equivalence test for population signal-to-noise ratio
## (i.e., phi_p^2) for the pth fixed effects factor 
## in an ANOVA setting

pop_sig_to_noise<-function(F.value, df.1, df.2, N, delta_lower=-0.2, delta_upper=+0.2, tol=0.001){
	
	my_F <- F.value
	my_df.1 <- df.1
	my_df.2 <- df.2
	my_N <- N


CI<-function(x){unlist((ci.snr(F.value=my_F, df.1=my_df.1, df.2=my_df.2, N=my_N, conf.level=(1-2*x))[c(1,2)]))}

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
within <- CI_x[1]>delta_lower &  CI_x[2]<delta_upper}
if(zzz>0.999){within<-TRUE}
}
return(zzz)
}

## example:

pop_sig_to_noise(F.value=2.5, df.1=5, df.2=4, N=100, 0.05,0.5)

# fixed error on:
pop_sig_to_noise(F.value=2.5, df.1=5, df.2=4, N=100, 0.25,0.5)

#################################################
#################################################

## Equivalence test for the population proportion 
## of variance accounted for in the dependent 
## variable by knowing group status (i.e., eta_p^2) 
## for the pth fixed effects factor in an ANOVA setting.

pop_prop_of_var<-function(F.value, df.1, df.2, N, delta_lower=0, delta_upper=+0.2, tol=0.001){
	
	my_F <- F.value
	my_df.1 <- df.1
	my_df.2 <- df.2
	my_N <- N


CI<-function(x){unlist((	
ci.pvaf(F.value=my_F, df.1=my_df.1, df.2=my_df.2, N=my_N, conf.level=(1-2*x))[c(1,3)]))
	}

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
within <- CI_x[1]>delta_lower &  CI_x[2]<delta_upper}
if(zzz>0.999){within<-TRUE}
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

sqrt_pop_sig_to_noise<-function(F, df.1, df.2, N, delta_lower=-0.2, delta_upper=+0.2, tol=0.001){
	
	my_F <- F
	my_df.1 <- df.1
	my_df.2 <- df.2
	my_N <- N


CI<-function(x){unlist((ci.srsnr(F.value=F, df.1=my_df.1, df.2=my_df.2, N=my_N, conf.level=(1-2*x))[c(1,2)]))}

zzz<-0.000
within<-FALSE
while(within==FALSE){
zzz <- zzz + tol
invisible(capture.output(CI_x <- CI(zzz)))
if(sum(is.na(CI_x))<1){
within <- CI_x[1]>delta_lower &  CI_x[2]<delta_upper}
if(zzz>0.999){within<-TRUE}
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
if(zzz>0.999){within<-TRUE}
}

return(zzz)
}

## example:
std_tageted_eff(means=c(2, 4, 9, 13), s.anova=.80, c.weights=c(.5, -.5, -.5, .5), n=c(30, 30, 30, 30), N=120, 0.15,0.95)



