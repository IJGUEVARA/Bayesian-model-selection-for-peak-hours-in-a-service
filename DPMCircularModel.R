proposal <- function(r) rtruncnorm(n = 1, a = 0, 
                                   b = Inf, mean = r, sd = 0.5)
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


# WOMACK PRIOR ------------------------------------------------------------

#K: number of predictors, poods is rho paramether
womack = function(K, podds)
{
  out = rep( -Inf, K + 1 )
  names(out) = 0 : K
  out[ K + 1 ] = 0
  for( k in (K-1):0 )
  {
    j = (k+1):K
    bb = out[j+1] + lgamma(j+1) - lgamma(k+1) - lgamma(j+1-k) + log(podds)
    out[k+1] = log( sum( exp(bb - max(bb) )) ) + max(bb)
    out = out - max(out) -log(sum( exp(out - max(out) ) ) )		
  }
  exp( out + lbeta(c(0:K) + 1, K - c(0:K) + 1 ) + log(K+1) )
}


# MCMC PROYECTED NORMAL ---------------------------------------------------
proposal= function(r) rtruncnorm(n = 1, a = 0, b = Inf, mean = r, sd = 0.7)
dpm_pnorm=function(datos,num,x,p)
{
  datos=circular(datos, units="radians", type="angles")
  # SETTING HYPERPARAMETERS -------------------------------------------------
  
  # NORMAL PRIOR FOR BETA
  tau = 0.7# VARIANCE NORMAL DISTRIBUTION
  np=dim(x)[2]
  ph=circular(p)
  n=length(datos)
  vp=cbind(cos(p),sin(p))
  v = cbind(cos(datos),sin(datos))
  invtau=1/tau # PRECISION NORMAL DISTRIBUTION
  b0_1 = rep(0,np) # PRIOR MEAN FIRST BETA
  b0_2 = rep(0,np) # PRIOR MEAN FIRST BETA
  
  delta0 = diag(invtau,np) # PRIOR PRECISION MATRIX 
  # GAMMA PRIOR CONCENTRATION PARAMETER DIRICHLET
  
  am=20 # SHAPE
  bm=0.1 # RATE
  # POSSIBLE MODELS ----------------------------------------------------------
  
  model = expand.grid(rep(list(c(0,1)), times = np-1))
  m = dim(model)[1] # TOTAL NUMBER POSSIBLE MODELS
  
  # MODEL LIST AS MATRIX OF LOGICAL VALUES
  model = as.matrix(cbind(rep(1,m),model)) 
  dimnames(model) =list(NULL, NULL) # REMOVE NAMES
  model = model==1 
  # CHOOSE RANDOM NUMBER OF CLUSTERS --------------------------------------
  
  M = 2 # CONCENTRATION PARAMETER
  Ni = rep(0,n) # SIZE OF EACH CLUSTER
  vj = rep(NA,100) # STICK BREAKING
  wj = rep(NA,100) # WEIGHTS DP
  
  # AUXILIARY VARIABLES
  u = runif(n) 
  aux =T
  
  for (i in 1:n) {
    while (aux) {
      Ni[i]=Ni[i]+1
      vj[Ni[i]]=rbeta(1,1,M)
      wj[Ni[i]]=vj[Ni[i]]*prod(1-vj[1:Ni[i]-1])
      if(sum(wj[1:Ni[i]])>1-u[i]){aux=F}
    }  
    aux=T
  }
  
  # NUMBER OF CLUSTERS
  K=max(Ni)
  # CREATE OBJECTS TO STORE RESULTS ------------------------------------------
  
  b1j=matrix(NA,ncol = np, nrow = 100) # COEFFICIENTS BETA 1 
  b2j=matrix(NA,ncol = np, nrow = 100) # COEFFICIENTS BETA 2
  pdi=matrix(NA,nrow = n,ncol=K) # CLUSTER PROBABILITY FOR EACH OBSERVATION
  di = rep(NA,K) # CLUSTER MEMBERSHIPS INDICATORS
  # STARTING VALUES MCMC ----------------------------------------------------
  covM=solve(delta0 + t(x)%*%x) # INITIAL COV MATRIX
  
  b1j[1:K,]= rmvnorm(K, b0_1, covM) # VALUES BETA 1
  b2j[1:K,]= rmvnorm(K,b0_2,covM) # VALUES BETA 2
  
  r = rep(0.5,n) # RADIUS 
  
  g1 = model[1,] # MODEL FIRST COMPONENT
  g2 = model[1,] # MODEL SECOND COMPONENT
  
  # EVALUATE OBSERVATIONS PROBABILITY
  for ( i in 1:K ) {pdi[,i] = dpnorm(datos, mu = c(x[i,]%*%b1j[i,],x[i,]%*%b2j[i,]), sigma = diag(2)) }
  
  
  # ASSIGN CLUSTER MEMBERSHIPS INDICATORS
  for (i in 1:n) {di[i]=sample(K,size = 1,prob = pdi[i,])}
  
  w = womack(np-1,1)/(3^{0:(np-1)}*choose(np-1,0:(np-1))) # PRIOR PROB MODELS 
  
  # STORE MCMC RESULTS ------------------------------------------------------
  
  samples_g1 = numeric(num) # SELECTED MODELS BETA 1
  samples_g2 = numeric(num) # SELECTED MODELS BETA 2
  #samples_b1 = vector("list", length = num) # COEFFICIENTS BETA 1 
  #samples_b2 = vector("list", length = num) # COEFFICIENTS BETA 1 
  #samples_wj = vector("list", length = num) # WEIGHTS
  
  # MCMC --------------------------------------------------------------------
  for (it in 1:num) {
    
    #cat(it, "of ",num ,"\r") # PRINT ITERATION NUMBER
    
    index = sort(unique(di)) # UNIQUE VALUES MEMBERSHIPS
    K = length(index) # NUMBER OF NON-EMPTY CLUSTERS
    pdi = matrix(NA,nrow = n,ncol=K+1) # PROBABILITY CLUSTER FOR EACH OBSERVATION
    Ni = rep(NA,K) # SIZE OF CLUSTERS
    u = rep(NA,n) # LATENT VARIABLE SLICE SAMPLER DIRICHLET PROCESS
    
    ## ADDING AND DELETING PREDICTORS ##
    # SAMPLING GAMMA MODEL MATRIX 
    ### NEW ###
    indg1 = sample(size=1,2:np)      
    G=cbind(g1, g2)
    Ga=G
    if (sum(Ga[indg1, ]) > 0) {
      Ga[indg1, ] = c(FALSE, FALSE)
    } else {
      while(1)
      {
        Ga[indg1, ] = replicate(1, sample(c(FALSE,TRUE),replace = T))
        if(sum(Ga[indg1, ]) > 0){break}
      }
    }
    
    
    ###
    
    
    # LOG LIKELIHOOD NEW PROPOSED MODEL
    total=sum(rowSums(Ga)!=0)
    likN=log(w[total])-.5*sum(Ga[,1])*log(tau)*K-.5*sum(Ga[,2])*log(tau)*K
    
    # LOG LIKELIHOOD OLD MODEL
    total=sum(rowSums(G)!=0)
    likO=log(w[total])-.5*sum(G[,1])*log(tau)*K-.5*sum(G[,2])*log(tau)*K
    ga1=Ga[,1]
    ga2=Ga[,2]
    ####
    
    
    for (i in 1:K) {
      
      # CLUSTER RADIUS
      
      rd=r[di==index[i]] 
      
      Ni[i]=length(rd) # CLUSTER SIZE
      
      vd=matrix(v[di==index[i],], ncol = 2) # CLUSTER DIRECTIONS
      
      # CLUSTER PREDICTOR VALUES 
      
      xd1=matrix(x[di==index[i],], ncol = np) # 1ST COMPONENT
      xd2=matrix(x[di==index[i],], ncol = np) # 2ND COMPONENT
      
      # M-H SAMPLING RADIUS
 
      for (t in 1:Ni[i])
      {
        # CANDIDATE
        r_star <- proposal(rd[t])
        
        # LOG PROBABILITY ACCEPT
        
        aux <- log(r_star)-log(rd[t])+.5*(rd[t]^{2}-r_star^{2})+
          (vd[t,1]*t(b1j[i,g1])%*%xd1[t,g1]+vd[t,2]*t(b2j[i,g2])%*%xd2[t,g2])*
          (r_star-rd[t])+
          log(pnorm(rd[t]))-log(pnorm(r_star))
        
        lnalpha <- min(0,aux)
        ru <- runif(1)
        
        # ACCEPT/REJECT
        if (ru <= exp(lnalpha)) 
        {
          rd[t] <- r_star
        }
      }
      
      
      r[di==index[i]]=rd # UPDATE RADIUS CLUSTER
      
      
      y1 =rd*vd[,1] # COMPUTE 1ST COMPONENT VECTOR
      y2 =rd*vd[,2] # COMPUTE 2ND COMPONENT VECTOR
      
      # DESIGN MATRIX FOR CLUSTER
      Xtx1 =t(xd1) %*% xd1 # 1ST COMPONENT
      Xtx2 =t(xd2) %*% xd2 # 2ND COMPONENT
      
      
      b1_pre =Xtx1 + diag(invtau, np) # PRECISION MATRIX BETA 1
      b2_pre =Xtx2 + diag(invtau, np) # PRECISION MATRIX BETA 2
      
      ### 1ST COMPONENT ###
      # POSTERIOR VALUES FOR NEW MODEL
      
      sg = solve(as.matrix(b1_pre[ga1,ga1])) # POSTERIOR COVARIANCE MATRIX
      bhg = sg%*%t(matrix(xd1[,ga1],nrow=Ni[i]))%*%y1 # POSTERIOR MEAN
      
      
      # POSTERIOR VALUES FOR OLD MODEL
      s1 = solve(as.matrix(b1_pre[g1,g1])) # POSTERIOR COVARIANCE MATRIX
      bh1 = s1%*%t(matrix(xd1[,g1],nrow=Ni[i]))%*%y1 # POSTERIOR MEAN
      
      # UPDATE LOG LIKELIHOOD NEW MODEL
      likN = likN+.5*log(det(sg))+.5*t(bhg)%*%b1_pre[ga1,ga1]%*%(bhg) 
      
      # UPDATE LOG LIKELIHOOD OLD MODEL
      likO = likO+.5*log(det(s1))+.5*t(bh1)%*%b1_pre[g1,g1]%*%(bh1) 
      
      ### 2ND COMPONENT ###  
      # POSTERIOR VALUES FOR NEW MODEL
      sg = solve(as.matrix(b2_pre[ga2,ga2])) # POSTERIOR COVARIANCE MATRIX
      bhg =sg%*%t(matrix(xd2[,ga2],nrow=Ni[i]))%*%y2 # POSTERIOR MEAN
      
      # POSTERIOR VALUES FOR OLD MODEL
      s2 =solve(as.matrix(b2_pre[g2,g2])) # POSTERIOR COVARIANCE MATRIX
      bh2 =s2%*%t(matrix(xd2[,g2],nrow=Ni[i]))%*%y2 # POSTERIOR MEAN
      
      # UPDATE LOG LIKELIHOOD NEW MODEL
      likN =likN+.5*log(det(sg))+.5*t(bhg)%*%b2_pre[ga2,ga2]%*%(bhg) 
      
      # UPDATE LOG LIKELIHOOD OLD MODEL
      likO =likO+.5*log(det(s2))+.5*t(bh2)%*%b2_pre[g2,g2]%*%(bh2) 
      
      
      # UPDATE BETAS
      b1j[i,] =rep(0,np)
      b2j[i,] =rep(0,np)
      
      # SAMPLE BETAS IN ACTUAL MODEL
      b1j[i,g1] = rmvnorm(1, bh1,s1)
      b2j[i,g2] = rmvnorm(1, bh2,s2)
      
      # UPDATE STICK BREAKING
      vj[i]=rbeta(1,1+Ni[i],M+sum(di>index[i]))
      wj[i]=vj[i]*prod(1-vj[1:i-1])
      
      # UPDATE LATENT VARAIBLE SLICE SAMPLER DP
      u[di==index[i]]=runif(Ni[i],min = 0,max = wj[i])
      
      # COMPUTE PROB CLUSTER
      pdi[,i]=fdpnorm(mu = cbind(x%*%b1j[i,],x%*%b2j[i,]) ,v = v, plog = T)
      
    }
    
    # UPDATE LAST CLUSTER VALUES
    
    wj[K+1]=1-sum(wj[1:K]) # WEIGHT
    
    # BETA 1
    b1j[K+1,]=rmvnorm(n = 1, b0_1, delta0) 
    b1j[K+1,!g1]=0
    
    # BETA 2
    b2j[K+1,]=rmvnorm(n = 1,b0_2, delta0) 
    b2j[K+1,!g2]=0
    
    # PROBABILITY
    pdi[,K+1]=fdpnorm(mu = cbind(x%*%b1j[K+1,],x%*%b2j[K+1,]), v = v, plog = T)
    
    # ADD NEW CLUSTER INDEX
    index[K+1]=max(index)+1 
    
    # UPDATE K  
    K=K+1 
    
    # UPDATE MEMBERSHIPS
    
    for (j in 1:n) {
      pdi[j,] =exp(pdi[j,]-max(pdi[j,]))
      pdi[j,] =(wj[1:K]>u[j])*pdi[j,]
      if(sum(pdi[j,])==0){pdi[j,]=1e-16}
      di[j]=sample(index,size = 1,prob = pdi[j,]) } 
    
    # SAVE BETA 1, BETA 2, AND DENSITY ESTIMATIONS
    
    #samples_b1[[it]][1:K,] <- b1j[1:K,]
    #samples_b2[[it]][1:K,] <- b2j[1:K,]
    #samples_wj[[it]][1:K] <- wj[1:K]
    
    
    # UPDATE MODEL 
    
    if(min(likN-likO,0)>log(runif(1))) { G = Ga } 
    g1=G[,1]
    g2=G[,2]
    
    
    # SAVE SELECTED MODEL 
    samples_g1[it] = which(apply(model,1, function(x) all.equal(x, g1)) == "TRUE")
    samples_g2[it] = which(apply(model,1, function(x) all.equal(x, g2)) == "TRUE")
    
    # UPDATE DP CONCENTRATION PARAMETER
    phi=rbeta(1,M+1,n)
    val=(am+K-1)/(n*(bm-log(phi)))
    probpi=val/(1+val)
    M=ifelse(runif(1)<=probpi,
             rgamma(1,shape = am+K,rate = bm-log(phi)),
             rgamma(1,shape = am+K-1,rate = bm-log(phi)))
  }
  # RESULTS -----------------------------------------------------------------
  nburn=num*0.5+1
  return(cbind(samples_g1[nburn:num],samples_g2[nburn:num]))
}
