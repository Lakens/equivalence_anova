noninf_Fstat <-  function(Fstat, df1, df2, N, eq_bound_eta, alpha = 0.10, tol = 1e-09){
  
  df1 <- 2
  df2 <- 48
  N <- 50
  
  # first calculate point-estimate to make sure it is less than eq_bound_eta:
  ncp_hat <- uniroot(function(q) (pf(Fstat, 
                                     df1, 
                                     df2, 
                                     ncp = q) - 0.5),
                     interval = c(0,
                                  Fstat * df1 * 1000),
                     tol = tol)$root
  
  eta2_hat =  ncp_hat/(ncp_hat+N)
  eta_pop<-(Fstat*df1)/(Fstat*df1+df2)
  ncp = N/(1/eta_pop-1)
  
  print(c("ncp_hat", "eta2_hat"))
  print(round(c(ncp_hat, eta2_hat),3))
  
  if(eta2_hat > eq_bound_eta){
    print(paste("Warning: point estimate of", eta2_hat, "is above eq_bound_eta."))
    return(1)
  }
  #Then calculate CI around eta.
  ncp_lower <- uniroot(function(q) (pf(Fstat, 
                                       df1, 
                                       df2, 
                                       ncp = q) - (1-alpha/2)),
                       interval = c(0, Fstat * df1),
                       tol = tol)$root
  
  ncp_upper <- uniroot(function(q) (pf(Fstat, 
                                       df1, 
                                       df2, 
                                       ncp = q) - (alpha/2)),
                       interval = c(0, Fstat * df1*100),
                       tol = tol)$root
  #Convert non-central F bounds for the confidence interval to bounds expressed in partial eta-squared
  LL_CI <- ncp_lower/(ncp_lower+N)
  UL_CI <- ncp_upper/(ncp_upper+N)
  
  if(eta2_hat <= eq_bound_eta){	
    pval <- uniroot(function(alpha) (pf(Fstat, 
                                        df1, 
                                        df2, 
                                        ncp = (eq_bound_eta * N) / (1 - eq_bound_eta)) - (alpha)),
                    interval = c(0,
                                 1),
                    tol = tol)$root
  }
    #add plot
    #Eta function
    xmin <- 0
    xmax <- 0.999
    ymax <- 1
    x=seq(xmin,xmax,length=1000)
    ncp = ncp_hat
    #    ncp = N/(1/eta2_hat-1)
    crit_f <- qf(1 - alpha, df1, df2, ncp = 0)
    
    eta_pop_dist <- function(x) df((x*df2)/(df1-x*df1), df1, df2, ncp = ncp)
    par(bg = "aliceblue")
    plot(-10,xlab=substitute(paste(eta[p]^2)), ylab="Density", axes=FALSE,
         main=substitute(paste("distribution for ",eta[p]^2 == eta2_hat,".")), xlim=c(0,xmax),  ylim=c(0, ymax))
    axis(side=1, at=seq(0,xmax, 0.1), labels=seq(0,xmax, 0.1))
    axis(side=2)
    ncp<-0
    eta_pop_crit<-(crit_f*df1)/(crit_f*df1+df2)
    #Draw null distribution
    ncp<-0 #So set ncp to null
    curve(eta_pop_dist, 0.00000000001, 0.99999999, n=10000, col="grey", lwd=2, add=TRUE)
    x=seq(eta_pop_crit,xmax,length=10000) 
    z<-df((x*df2)/(df1-x*df1), df1, df2) #determine upperbounds polygon
    polygon(c(eta_pop_crit,x,xmax),c(0,z,0),col=rgb(1, 0, 0,0.5)) #draw polygon
    #Add Type 2 error rate
    ncp <- ncp_hat #So set ncp to observed value
    curve(eta_pop_dist, 0.00000000001, 0.99999999, n=10000, col="black", lwd=2, add=TRUE)
    y=seq(0.00000000001,eta_pop_crit,length=10000) 
    z<-df((y*df2)/(df1-y*df1), df1, df2, ncp) #determine upperbounds polygon
    polygon(c(y,eta_pop_crit,eta_pop_crit),c(0,z,0),col=rgb(0, 0, 1,0.5))
    segments(eta_pop_crit, 0, eta_pop_crit, ymax-0.03, col= 'black', lwd=2)
    
    
    
    invisible(list(eta2_hat = eta2_hat, TOST_p = pval, alpha = alpha, eq_eqbound_eta = eq_bound_eta, LL_CI_TOST = LL_CI, UL_CI_TOST = UL_CI))
  
}

## a random example:

Fstat <- 1
df1 <- 1
df2 <- 49
N <- 50
eq_bound_eta <- 0.3

res <- noninf_Fstat(Fstat = Fstat, df1 = df1, df2 = df2, N = N, eq_bound_eta = eq_bound_eta)
res
