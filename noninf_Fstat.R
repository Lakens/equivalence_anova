




noninf_Fstat <-  function(Fstat, df1, df2, N, eq_bound_eta, alpha = 0.05, makeplot =FALSE, tol = 1e-09){
  

# first calculate point-estimate to make sure it is less than eq_bound_eta:

  ncp_hat <- uniroot(function(q) (pf(Fstat, 
                                     df1, 
                                     df2, 
                                     ncp = q) - 0.5),
                     interval = c(0,
                                  Fstat * df1 * 1000),
                     tol = tol)$root

  eta2_hat =  ncp_hat/(ncp_hat+N)
#  eta_pop<-(Fstat*df1)/(Fstat*df1+df2)  # not sure where this is coming from?
#  ncp = N/(1/eta_pop-1)

  
  if(eta2_hat > eq_bound_eta){

    print(paste("Warning: point estimate of", round(eta2_hat,3), "is above eq_bound_eta."))
	return(list(ncp_hat=ncp_hat, eta2_hat=eta2_hat, pval=1))}


  #Then calculate one-side (1-alpha)% CI for eta2.

  ncp_lower <- 0
  

  ncp_upper <- uniroot(function(q) (pf(Fstat, 

                                       df1, 

                                       df2, 

                                       ncp = q) - (alpha)),

                       interval = c(0, ncp_hat*100),

                       tol = tol)$root

# Convert the CI bounds for the ncp to bounds eta2:

  LL_CI <- ncp_lower/(ncp_lower+N)
  UL_CI <- ncp_upper/(ncp_upper+N)

  

  if(eta2_hat <= eq_bound_eta){	

    pval <- uniroot(function(alpha) (
    
    				pf(	Fstat, 
					df1, 
				    df2, 
					ncp = (eq_bound_eta * N) / (1 - eq_bound_eta)) - (alpha)),
                    interval = c(0,1),
                    tol = tol)$root

  }

    #add plot
if(makeplot==TRUE){
    #Eta function

    xmin <- 0

    xmax <- 0.999

    ymax <- 1

    x=seq(xmin,xmax,length=1000)

    ncp = ncp_hat

    #    ncp = N/(1/eta2_hat-1)

    crit_f <- qf(1 - alpha, df1, df2, ncp = 0)

    eta_pop_dist <- function(x) df((x*df2)/(df1-x*df1), df1, df2, ncp = ncp)

    par(bg = "white")

    plot(-10,xlab=substitute(paste(eta[p]^2)), ylab="", axes=FALSE,

         xlim=c(0,xmax),  ylim=c(0, ymax))

    axis(side=1, at=seq(0,xmax, 0.1), labels=seq(0,xmax, 0.1))

    #axis(side=2)

title(main=paste("Equivalence bound = ",round(eq_bound_eta,digits=3),"\nObserved eta = ",round(eta2_hat,digits=3)," \n one-sided ", 100*(1-alpha),"% CI [",round(LL_CI,digits=3),";",round(UL_CI,digits=3),"], p = ",round(pval,digits=3),sep=""), cex.main=1)

    ncp<-0

    eta_pop_crit<-(crit_f*df1)/(crit_f*df1+df2)

    # #Draw null distribution


    #Add Type 2 error rate

    ncp <- ncp_hat #So set ncp to observed value

    curve(eta_pop_dist, 0.00000000001, 0.99999999, n=10000, col="black", lwd=2, add=TRUE)


    abline(v=eq_bound_eta, lty=2)

    abline(v=eta2_hat, lty=2, col="grey")

    points(x=eta2_hat, y=0.7, pch=15, cex=2)

    segments(LL_CI,0.7,UL_CI,0.7, lwd=3)


}
  
  if(pval==0){print(paste("warning: pval <=", tol,". Lower tol for additional accuracy."))}

	return(list(ncp_hat=ncp_hat, eta2_hat=eta2_hat, pval=pval))

}



## Does our test has correct type 1 error? (this takes ~5 min)

		 # We simulate data from the "null"
		 # where the true population eta_p^2 value is
		 # equal to 0.4 and define the DELTA margin = 0.3999.
		 
nSim <- 5000
simresults_pval <- apply(cbind(1:nSim),1, 
	function(x){
				if(round(x/37)==(x/37)){print(paste("% complete:",round(x/nSim,2)))}
		aovdata<-simdata(2) 
		ttt<-(summary(aov(aovdata$dat~aovdata$groups)))
		result<-noninf_Fstat(Fstat= ttt[[1]]$F[1], 
					df1= ttt[[1]]$Df[1], 
					df2= ttt[[1]]$Df[2], 
					N=dim(aovdata)[1], 
					eq_bound_eta =0.3999)
		return(result$pval<0.05)
				} 
			)

sum(simresults_pval, na.rm=TRUE)/nSim



# Small simulation study (will take ~5min):
		 # We simulate data from a population
		 # where the true population eta_p^2 value is
		 # equal to 0.4 and the true population ncp=40.
		 # We define the DELTA eq_bound_eta= 0.3999.
		
		 
nSim <- 10000
simresults <- apply(cbind(1:nSim),1, 
	function(x){
				if(round(x/37)==(x/37)){print(paste("% complete:",round(x/nSim,2)))}
		aovdata<-simdata(2) 
		ttt<-(summary(aov(aovdata$dat~aovdata$groups)))
		result<-noninf_Fstat(Fstat= ttt[[1]]$F[1], 
					df1= ttt[[1]]$Df[1], 
					df2= ttt[[1]]$Df[2], 
					N=dim(aovdata)[1], 
					eq_bound_eta=0.3999)
		return(unlist(result))
				} 
			)

## Is our point estimate for Lambda (i.e. ncp_hat) unbiased?  (true ncp=40)
mean(simresults["ncp_hat",], na.rm=TRUE)
median(simresults["ncp_hat",], na.rm=TRUE)
hist(simresults["ncp_hat",])

## Is our point estimate for eta2 unbiased? (true eta2=0.40)
mean(simresults["eta2_hat",], na.rm=TRUE)
median(simresults["eta2_hat",], na.rm=TRUE)
hist(simresults["eta2_hat",])

## Does our test has correct type 1 error? (pvalues should be uniformly ditributed)
mean(simresults["pval",]<0.05, na.rm=TRUE)
median(simresults["pval",], na.rm=TRUE)
hist(simresults["pval",])