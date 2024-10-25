# LIBRARIES
library(circular)
library(truncnorm)
library(mvtnorm)

# PROPOSAL DISTRIBUTION METROPOLIS HASTINGS ---------------------------------------------------
proposal= function(r) rtruncnorm(n = 1, a = 0, b = Inf, mean = r, sd = 0.7)
# WOMACK PRIOR ------------------------------------------------------------

#p: number of predictors
#rho: rho parameter
womack = function(p, rho)
{
  out = rep( -Inf, p + 1 )
  names(out) = 0 : p
  out[ p + 1 ] = 0
  for( p in (p-1):0 )
  {
    j = (p+1):p
    bb = out[j+1] + lgamma(j+1) - lgamma(p+1) - lgamma(j+1-p) + log(rho)
    out[p+1] = log( sum( exp(bb - max(bb) )) ) + max(bb)
    out = out - max(out) -log(sum( exp(out - max(out) ) ) )		
  }
  exp( out + lbeta(c(0:p) + 1, p - c(0:p) + 1 ) + log(p+1) )
}

# DENSITY PROYECTED NORMAL IDENTITY COV -----------------------------------

# Evaluates the density of the proyected normal through a vector of means
# mu: mean matrix nx2
# v: direction vector nx2
# plog: T or F return log density

fdpnorm <- function(mu,v,plog=F) 
{
  U=rowSums(v*mu) 
  PU=pnorm(U,mean=0,sd=1)
  DU=dnorm(U,mean=0,sd=1)
  D=rep(NA,dim(v)[1])
  
  if(plog==T)
  {
    D=log(0.5)-log(pi)-0.5*(mu[,1]^{2}+mu[,2]^{2})+log(1+U*(DU/PU)^{-1})
    # IF RETURNING INF REPLACE 
    suma=sum(is.na(D))+sum(D==Inf)
    
    if(any(is.na(D)) || any(D==Inf))
    {
      D[is.na(D)]=log(0.5)-log(pi)-0.5*(mu[is.na(D),1]^{2}+mu[is.na(D),2]^{2})+
        log(1+U[is.na(D)]*PU[is.na(D)]/1e-16)
      D[D==Inf]=log(0.5)-log(pi)-0.5*(mu[D==Inf,1]^{2}+mu[D==Inf,2]^{2})+
        log(1+U[D==Inf]*PU[D==Inf]/1e-16)
    }
    
  } else
  {
    D=0.5/pi*exp(-0.5*(mu[,1]^{2}+mu[,2]^{2}))*(1+U*PU/DU)
    # IF RETURNING INF REPLACE 
    suma=sum(is.na(D))+sum(D==Inf)
    if(any(is.na(D)) || any(D==Inf))
    {
      D[is.na(D)]=0.5/pi*exp(-0.5*(mu[is.na(D),1]^{2}+
                                     mu[is.na(D),2]^{2}))*(1+U[is.na(D)]*PU[is.na(D)]/1e-16)
      D[D==Inf]=0.5/pi*exp(-0.5*(mu[D==Inf,1]^{2}+
                                   mu[D==Inf,2]^{2}))*(1+U[D==Inf]*PU[D==Inf]/1e-16)
    }
    
  }
  
  
  return(D)
}

# Plot Circle density ----------------------------------------------------------

CircplotEst=function(aux,datos,fyt,pr,lim)
{
  # aux: matrix with MCMC density estimations
  # fyt: true density
  # pr: control parameter bindwidth height
  # lim: x and y axis limit
  
  # points
  p = seq(0,2*pi,length.out = 100)
  # labels
  labp=c(0,pi/2,pi,3*pi/2)
  labp1=ifelse(abs(cos(labp))==1,cos(labp)-cos(labp)*0.1,cos(labp))
  labp2=ifelse(abs(sin(labp))==1,sin(labp)-sin(labp)*0.25,sin(labp))
  # width of each bin
  binwidth <- 2 * pi / 20  
  # estimated density
  fyest=apply(aux,2,mean)
  
  ci = matrix(NA,ncol=2,nrow = 100)
  for (i in 1:100) {
    ci[i,] <- quantile(aux[,i],c(.25,.975))
  }
  
  xx1=(1+fyt)*cos(p)
  xx2=(1+fyest)*cos(p)
  
  pts=(1+fyt)*sin(p)
  ptsE=(1+fyest)*sin(p)
  
  xxl1=(1+ci[,1])*cos(p)
  xxl2=(1+ci[,2])*cos(p)
  
  lim1=(1+ci[,1])*sin(p)
  lim2=(1+ci[,2])*sin(p)
  
  hist_counts <- hist(as.numeric(datos), breaks = seq(0, 2 * pi, by = binwidth), plot = FALSE)$counts
  
  plot(cos(p), sin(p), type="l", xlim=c(-lim,lim), ylim=c(-lim,lim),axes = F, xlab=" ",ylab = " ")
  text(labp1, labp2, labels=expression(0 ,frac(pi,2), pi,frac(3*pi,2)),cex=1.5)
  for (i in 1:n) {
    angle <- (i - 1) * binwidth
    start_angle <- angle - binwidth / 2
    end_angle <- angle + binwidth / 2
    a <-  c(pr * hist_counts[i]/n * cos(start_angle), pr * hist_counts[i]/n * cos(end_angle), 0)
    y <- c(pr * hist_counts[i]/n * sin(start_angle), pr * hist_counts[i]/n * sin(end_angle), 0)
    polygon(a, y, col = "snow4", border = "black")
  }
  lines(xx1,pts, col="red",lwd=2)
  lines(xx2,ptsE, col="#00028C", lwd=2,lty=2)
  # CREDIBILITY INTERVAL
  polygon(c(xxl1, rev(xxl2)),c(lim1, rev(lim2)),
          col = rgb(0,0,1,0.25), border = NA) 
}  


