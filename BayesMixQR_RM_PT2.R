###############################################
#
# Bayesian Mixtures of Quantile Regressions
#
###############################################
#
# 04-03-2016

### Load data: 'data("ReadingSkills", package = "betareg")

#library(label.switching)
library(MCMCpack)
library(LearnBayes)

# Variable names: Dependent - 'accuracy'; Response - 'dyslexia', 'iq'.

# p = quantile of interest
# s = burn in for mcmc sampler
# r = total iterations - burn in
# k = number of components
# d = number of Beta parameters


BMQR<-function(p, s, r, k ,d)
{
  
  #########################################################################
  #functions for calcs for posterior#
  #########################################################################
  
  prob.beta <- function(y, x, p, beta, sig, z, j, d)
  {
    v <- x %*% beta[j,]
    u <- (y - v)/sig[j]
    ind <- ifelse(u < 0., 1., 0.)
    g <-  - u * (p - ind)
    gnew <- z[,k]*g*(beta[j,d])^2
    return(sum(gnew))
  }
  
  prob.sigma <- function(y, x, p, beta, sig, z, j, a, b)
  {
    v <- x %*% beta[j,]
    u <- (y - v)/sig[j]
    ind <- ifelse(u < 0., 1., 0.)
    g <-  - u * (p - ind)
    gnew <- (a-z[,j]-1)*log(sig[j])-z[,j]*g-(sig[j]/b)
    return(sum(gnew))
  }
  
  prob.z <- function(y, x, p, beta, sig, phi, z, m, k)
  {
    gplus<-matrix(NA,m,k)
    #print(dim(gplus))
    for(j in 1:k)
    {
      v <- x %*% beta[j,]
      u <- (y - v)/sig[j]
      ind <- ifelse(u < 0., 1., 0.)
      g <-  - u * (p - ind)
      #print(dim(z))
      gplus[,j] <- z[j]*log((phi[j]/sig[j]*(p*(1-p))*exp(g))) 
    }
    
    gnew <- sum(gplus) + dmultinom(z,prob=phi, log = TRUE)
    return(gnew)
  }
  
  
  #########################################################################
  #Acceptance functions for posterior MH MCMC #
  #########################################################################
  
  accept1 <- function(y, x, p, denom1, thold, thnew, beta, sig, z, j, d)
  {#print("first")
    prob.beta <- function(y, x, p, beta, sig, z, j, d)
    {
      v <- x %*% beta[j,]
      u <- (y - v)/sig[j]
      ind <- ifelse(u < 0., 1., 0.)
      g <-  - u * (p - ind)
      gnew <- z[,j]*g*(beta[j,d])^2
      return(sum(gnew))
    }
    
    num1 <- prob.beta(y, x, p, beta, sig, z, j, d)
    mh1 <- (num1 - denom1)
    mh1 <- exp(mh1)
    mh1 <- min(mh1, 1.)
    q1 <- runif(1., 0., 1.)
    if(mh1 > q1) {
      th <- thnew
      den <- num1
    }
    else {
      th <- thold
      den <- denom1
    }
    list(th = th, den = den)
  }
  
  #########################################################################
  accept2 <- function(y, x, p, denom2, beta, sigold,signew, z, j, a, b)
  {#print("second")
    prob.sigma <- function(y, x, p, beta, sig, z, j, a, b)
    {
      v <- x %*% beta[j,]
      u <- (y - v)/sig[j]
      ind <- ifelse(u < 0., 1., 0.)
      g <-  - u * (p - ind)
      gnew <- (a-z[,j]-1)*log(sig[j])-z[,j]*g-(sig[j]/b)#component 1 has very small probability as only one component, so variance estimate is messed up.
      return(sum(gnew))
    }
    num2 <- prob.sigma(y, x, p, beta, signew, z, j, a, b)
    mh2 <- (num2 - denom2)
    mh2 <- exp(mh2)
    mh2 <- min(mh2, 1.)
    q2 <- runif(1., 0., 1.)
    if(mh2 > q2) {
      sig <- signew[j]
      den2 <- num2
    }
    else {
      sig <- sigold[j]
      den2 <- denom2
    }
    
    list(sig = sig, den2 = den2)
  }
  
  #########################################################################
  accept3 <- function(y, x, p, denom3, zold, znew, beta, sig, phi, m, k)
  {#print("Third")
    #print(paste("dim(zold)=",length(zold)))
    prob.z <- function(y, x, p, beta, sig, phi, z, m, k)
    {
      gplus<-matrix(NA,length(y),k)
      #print(dim(gplus))
      for(j in 1:k)
      {
        v <- x %*% beta[j,]
        u <- (y - v)/sig[j]
        ind <- ifelse(u < 0., 1., 0.)
        g <-  - u * (p - ind)
        #print(dim(z[j]*log((phi[j]/sig[j]*(p*(1-p))*exp(g)))))
        gplus[,j] <- z[j]*log((phi[j]/sig[j]*(p*(1-p))*exp(g))) 
      }
      
      gnew <- sum(gplus) + dmultinom(z,prob=phi, log = TRUE)
      return(gnew)
    }
    #denom <- prob(y, x, p, thold)
    
    num3 <- prob.z(y, x, p, beta, sig, phi, znew, m, k)
    
    print(paste("num3",num3))
    mh3 <- (num3 - denom3)
    print(paste("mh3",mh3))
    mh3 <- exp(mh3)
    print(paste("mh3",mh3))
    mh3 <- min(mh3, 1.)
    print(paste("mh3",mh3))
    q3 <- runif(1., 0., 1.)
    print(paste("q3",q3))
    if(mh3 > q3) {
      z <- znew
      den3 <- num3
    }
    else {
      z <- zold
      den3 <- denom3
    }
    list(z = z, den3 = den3)
  }
  
  
  
  
  
  ########################################################################
  # DATA###################################################################
  x <- serum.dat[, 1.:3.]
  y <- serum.dat[, 4.]
  n <- dim(serum.dat)[1.]
  #number of components
  #k=2
  #number of beta parameters
  #d=3
  #Sigma prior parameter values
  a= 1.2
  b= 1.4
  #number of iterations after burnin
  #r=5000
  #Burn in iterations
  #Burnin=1000
  
  av1 <- mean(x[, 2.])
  x[, 2.] <- x[, 2.] - av1
  x[, 3.] <- x[, 2.]^2.
  
  #beta array to store all beta parameter iterations 
  beta <- matrix(NA,k,d)
  
  b0 <- matrix(NA, r,k)
  b1 <- matrix(NA, r,k)
  b2 <- matrix(NA, r,k)
  #variance parameters for beta priors
  sig0 <- rep(1.3,k)
  sig1 <- rep(1.3,k)
  sig2 <- rep(0.1,k)
  sigma2 <- 1
  
  
  sigma <- PHI<- matrix(NA,r,k)
  z_all<-array(NA,c(n,k,r))
  
  zold<-matrix(NA,nrow=length(y),ncol=k)
  for(m in 1:n)
  {
    location <- sample(1:k,1,replace=TRUE) 
    zold[m,] <- rep(0,k)
    zold[m,location] <- 1
  }
  
  z_all[,,1]<-zold
  
  b0old <- rep(0,k)
  b1old <- rep(0,k)
  b2old <- rep(0,k)
  sigold <- rep(1,k)
  
  sig0 <- 1.3
  sig1 <- 1.3
  sig2 <- 0.1
  
  b0[1,] <- b0old
  b1[1,] <- b1old
  b2[1,] <- b2old
  beta <- cbind(b0old, b1old, b2old)
  
  phi<-c(0.5,0.5)
  PHI[1,]<-phi
  
  denom1<-denom2<-vector(mode="numeric",length=k)
  for(j in 1:k)
  {
    denom1[j] <- prob.beta(y, x, p, beta, sigold, zold, j, d)
    denom2[j] <- prob.sigma(y, x, p, beta, sigold, zold, j, a, b)
  }  
  denom3 <- prob.z(y, x, p, beta, sigold, phi, zold[m,], m, k)
  
  #PHI posterior
  
  alpha0<-vector(mode="numeric",length=k)
  
  sigma[1,]<-sigold
  
  
  
  ###########################################################  
  
  for(i in 2.:r) {
    if(i %% 100==0) {
      cat(paste0("iteration: ", i, "\n"))
    }
    
    b0new<-b1new<-b2new<-signew<-vector(mode="numeric",length=k)
    znew<-matrix(NA,nrow=length(y),ncol=k)
    
    for(j in 1:k)
    {
      
      
      b0new[j] <- rnorm(1., b0old[j], sig0)
      beta[j,] <- c(b0new[j], b1old[j], b2old[j])
      acc1 <- accept1(y, x, p, denom1, b0old, b0new, beta, sigold, zold, j, 1)
      b0old[j] <- acc1$th
      denom1[j] <- acc1$den
      
      b1new[j] <- rnorm(1., b1old[j], sig1)
      beta[j,] <- c(b0old[j], b1new[j], b2old[j])
      acc2 <- accept1(y, x, p, denom1, b1old, b1new, beta, sigold, zold, j, 2)
      b1old[j] <- acc2$th
      denom1[j] <- acc2$den
      
      b2new[j] <- rnorm(1., b2old[j], sig2)
      beta[j,] <- c(b0old[j], b1old[j], b2new[j])
      acc3 <- accept1(y, x, p, denom1, b2old, b2new, beta, sigold, zold, j, 3)
      b2old[j] <- acc3$th
      denom1[j] <- acc3$den
      
      ################
      
      signew[j] <- rnorm(1,log(sigold[j]),sigma2)
      signew[j] <- exp(signew)
      beta[j,] <- c(b0old[j], b1old[j], b2old[j])
      acc4 <- accept2(y, x, p, denom2, beta, sigold, signew, zold, j, a, b) 
      sigold[j] <- acc4$sig
      denom2 <- acc4$den2
      #
      b2[i,j] <- b2old[j]
      b1[i,j] <- b1old[j] - (2. * b2old[j] * av1)
      b0[i,j] <- b0old[j] - b1old[j] * av1 + b2old[j] * av1 * av1
      
      
      sigma[i,j] <- sigold[j]
      
      ################
      
      alpha0[j]<-1/k+sum(zold[,j])
    }  
    
    
    phi <- rdirichlet(1, alpha0)
    
    for(m in 1:n)
    {
      location <- sample(1:k,1,replace=TRUE) 
      znew[m,] <- rep(0,k)
      znew[m,location] <- 1
      
      acc5 <- accept3(y, x, p, denom3, zold[m,], znew[m,], beta, sigold, phi, m, k) #y, x, p, denom, thold, thnew, beta, sig, phi
      zold[m,] <- acc5$z
      denom3 <- acc5$den3
    } 
    ################
    
    z_all[,,i] <-zold
    
    PHI[i,]<-phi
  }
  return(list("phi"=PHI, "b0"=b0, "b1"=b1, "b2"=b2, "sigma"=sigma, "z_all"=z_all))
}

### TEST ###

paul<-BMQR(0.5, 100., 600.,2,3)

library(ggplot2)
library(grid)
  windows()
  
  ###PHI###
  
  phi.plot<-as.data.frame(paul$phi)
  names(phi.plot)<-c("component1","component2")
  phi.plot$iterations<-1:600

  plot1<-ggplot() + geom_line(data=phi.plot,aes(x=iterations,y=component1,color='lightgreen')) + theme_bw() +
    #geom_line(data=data1.f,aes(x=Iterations,Beta0_2,color='blue')) +
    geom_hline(yintercept=mean(phi.plot$component1[101:600])) + theme(legend.position="none") + labs(x="Iterations",y="phi")

  plot2<-ggplot() + geom_line(data=phi.plot,aes(x=iterations,y=component1,color='lightgreen')) + theme_bw() +
    #geom_line(data=data2.f,aes(x=Iterations,Beta1_2,color='blue')) +
    geom_hline(yintercept=mean(phi.plot$component2[101:600])) + theme(legend.position="none") + labs(x="Iterations",y="phi")


  pushViewport(viewport(layout = grid.layout(3, 1, heights = unit(c(1, 4, 4), "null"))))

  grid.text("Median Phi", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plot1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(plot2, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
  
  #####################################################################################
  
  ###B0###
  windows()
  b0.plot<-as.data.frame(paul$b0)
  names(b0.plot)<-c("component1","component2")
  b0.plot$iterations<-1:600
  
  plot1<-ggplot() + geom_line(data=b0.plot,aes(x=iterations,y=component1,color='lightgreen')) + theme_bw() +
    #geom_line(data=data1.f,aes(x=Iterations,Beta0_2,color='blue')) +
    geom_hline(yintercept=mean(b0.plot$component1[101:600])) + theme(legend.position="none") + labs(x="Iterations",y="B0")
  
  plot2<-ggplot() + geom_line(data=b0.plot,aes(x=iterations,y=component1,color='lightgreen')) + theme_bw() +
    #geom_line(data=data2.f,aes(x=Iterations,Beta1_2,color='blue')) +
    geom_hline(yintercept=mean(b0.plot$component2[101:600])) + theme(legend.position="none") + labs(x="Iterations",y="B0")
  
  
  pushViewport(viewport(layout = grid.layout(3, 1, heights = unit(c(1, 4, 4), "null"))))
  
  grid.text("Median B0", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plot1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(plot2, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
  
  #####################################################################################
  
  ###B1###
  windows()
  b1.plot<-as.data.frame(paul$b1)
  names(b1.plot)<-c("component1","component2")
  b1.plot$iterations<-1:600
  
  plot1<-ggplot() + geom_line(data=b1.plot,aes(x=iterations,y=component1,color='lightgreen')) + theme_bw() +
    #geom_line(data=data1.f,aes(x=Iterations,Beta0_2,color='blue')) +
    geom_hline(yintercept=mean(b1.plot$component1[101:600])) + theme(legend.position="none") + labs(x="Iterations",y="B1")
  
  plot2<-ggplot() + geom_line(data=b1.plot,aes(x=iterations,y=component1,color='lightgreen')) + theme_bw() +
    #geom_line(data=data2.f,aes(x=Iterations,Beta1_2,color='blue')) +
    geom_hline(yintercept=mean(b1.plot$component2[101:600])) + theme(legend.position="none") + labs(x="Iterations",y="B1")
  
  
  pushViewport(viewport(layout = grid.layout(3, 1, heights = unit(c(1, 4, 4), "null"))))
  
  grid.text("Median B1", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plot1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(plot2, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
  
  ###B2###
  windows()
  b2.plot<-as.data.frame(paul$b2)
  names(b2.plot)<-c("component1","component2")
  b2.plot$iterations<-1:600
  
  plot1<-ggplot() + geom_line(data=b2.plot,aes(x=iterations,y=component1,color='lightgreen')) + theme_bw() +
    #geom_line(data=data1.f,aes(x=Iterations,Beta0_2,color='blue')) +
    geom_hline(yintercept=mean(b2.plot$component1[101:600])) + theme(legend.position="none") + labs(x="Iterations",y="B2")
  
  plot2<-ggplot() + geom_line(data=b1.plot,aes(x=iterations,y=component1,color='lightgreen')) + theme_bw() +
    #geom_line(data=data2.f,aes(x=Iterations,Beta1_2,color='blue')) +
    geom_hline(yintercept=mean(b2.plot$component2[101:600])) + theme(legend.position="none") + labs(x="Iterations",y="B2")
  
  
  pushViewport(viewport(layout = grid.layout(3, 1, heights = unit(c(1, 4, 4), "null"))))
  
  grid.text("Median B2", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plot1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(plot2, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))

 ###################################################################################################
     
  ###Sigma###
  windows()
  sigma.plot<-as.data.frame(paul$sigma)
  names(sigma.plot)<-c("component1","component2")
  sigma.plot$iterations<-1:600
  sigma.plot<-sigma.plot[101:600,]
  
  plot1<-ggplot() + geom_line(data=sigma.plot,aes(x=iterations,y=component1,color='lightgreen')) + theme_bw() +
    #geom_line(data=data1.f,aes(x=Iterations,Beta0_2,color='blue')) +
    geom_hline(yintercept=mean(sigma.plot$component1)) + theme(legend.position="none") + labs(x="Iterations",y="Sigma")
  
  plot2<-ggplot() + geom_line(data=sigma.plot,aes(x=iterations,y=component1,color='lightgreen')) + theme_bw() +
    #geom_line(data=data2.f,aes(x=Iterations,Beta1_2,color='blue')) +
    geom_hline(yintercept=mean(sigma.plot$component2)) + theme(legend.position="none") + labs(x="Iterations",y="Sigma")
  
  
  pushViewport(viewport(layout = grid.layout(3, 1, heights = unit(c(1, 4, 4), "null"))))
  
  grid.text("Median Sigma", vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(plot1, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(plot2, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
  
  
