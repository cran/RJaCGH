####  Copyright (C) 2005-2006  Oscar Rueda Palacio and Ram� D�z-Uriarte

#### This program is free software; you can redistribute it and/or
#### modify it under the terms of the GNU General Public License
#### as published by the Free Software Foundation; either version 2
#### of the License, or (at your option) any later version.

#### This program is distributed in the hope that it will be useful,
#### but WITHOUT ANY WARRANTY; without even the implied warranty of
#### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#### GNU General Public License for more details.

#### You should have received a copy of the GNU General Public License
#### along with this program; if not, write to the Free Software
#### Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
#### USA.



.__RJACGH_DEBUG <- FALSE

#########################################################################
##Functions for Non Homogeneous Hidden Markov Models with normal
##probability emissions
#########################################################################

#########################################################################
##Function to compute the transition matrix
## q are the intercept probabilities
## beta is the vector of regression coefficients
## x is the distance to the next gene
#########################################################################
Q.NH <- function(beta, x, q=-beta) {

  k <- nrow(q)
  Q <- exp(q + beta*x)
  Q <- Q / apply(Q, 1, sum)
  Q
}

#########################################################################
## filter and likelihood of the non - homogeneous HMM
## q, beta, x are the same as before
## mu & sigma.2 are the normal distribution parameters
## stat is the initial distribution
## k is the number of hidden states
#########################################################################

normal.HMM.likelihood.NH.filter <- function(y, k, beta, x, mu, sigma.2,
                                     stat=NULL, q=-beta) {
  ## Prevents that c, denominator never gets zero
  TOL <- 1E-10
  if (is.null(stat)) stat <- rep(1/k, k)
  
  ##last x is not observed. Besides, it doesn't count on the likelihood
  x <- c(x, 0)
  
  n <- length(y)
  if (k==1) {
    loglik <- sum(dnorm(y, mu, sqrt(sigma.2), log=TRUE))
    filter.cond <- matrix(rep(1, n+1), n+1)
    filter <- matrix(rep(1, n), n)
    res <- list(loglik=loglik, filter.cond=filter.cond, filter=filter)
  }
  else {
    filter <- matrix(NA, n,k)
    filter.cond <- matrix(NA, n+1,k)
    c <- rep(NA, n)
    filter.cond[1,] <- stat
    for (i in 1:n) {
      c[i] <- sum(filter.cond[i,]*dnorm(y[i], mu, sqrt(sigma.2)))
      c[i] <- ifelse(identical(c[i], 0), TOL, c[i])
      filter[i,] <- filter.cond[i,]*dnorm(y[i], mu, sqrt(sigma.2))/c[i]
      Q <- Q.NH(q=q, beta=beta, x=x[i])
      for (j in 1:k) {
        filter.cond[i+1,j] <- sum(filter[i,]*Q[,j,drop=FALSE])
      }
    }
    res <- list(loglik=sum(log(c)), filter.cond=filter.cond, filter=filter, c=log(c))
  }
  res
}

#########################################################################
## likelihood of the non - homogeneous HMM. C version
## q, beta, x are the same as before
## mu & sigma.2 are the normal distribution parameters
## stat is the initial distribution
## k is the number of hidden states
#########################################################################

normal.HMM.likelihood.NH.C <- function(y, x, mu,
                                       sigma.2, beta, stat=NULL, q=-beta) {
  k <- length(mu)
  if (is.null(stat)) stat <- rep(1/k, k)
  ##last x is not observed. Besides, it doesn't count on the likelihood
  x <- c(x, 0)
  loglik <- .C("normalNHHMMlikelihood", y=as.double(y), k=as.integer(k),
               x=as.double(x), n=as.integer(length(y)),
               q=as.double(as.vector(q)), beta=as.double(beta),
               stat=as.double(stat), mu=as.double(mu), sigma.2=as.double(sigma.2),
               loglik=as.double(0))
  loglik$x <- loglik$x[-length(x)]
  loglik
  
}

#########################################################################
## Function to simulate hidden state sequence by backward smoothing
## The call to normal.HMM.likelihood.NH should be changed for c++ version
## Slightly modified. Doesn't sample, it chooses the state with maximal
## probability.
## now can do model averaging
## now viterbi is used instead of this one
#########################################################################
hidden.states <- function(obj, x) {
  k <- length(obj$mu)
    n <- length(obj$y)
    filter <- normal.HMM.likelihood.NH.filter(y=obj$y, k=k, q=-obj$beta,
                                       beta=obj$beta, x=x, mu=obj$mu,
                                       sigma.2=obj$sigma.2, stat=obj$stat)$filter
    states <-rep(NA, n)
    B <- array(NA, c(n, k))
    B[n,] <- filter[n,]
    states[n] <- which.max(B[n,])
    for (i in (n-1):1) {
      Q <- Q.NH(q=-obj$beta, beta=obj$beta, x=x[i])
      den <- sum(filter[i,]*Q[,states[i+1]])
      den <-ifelse(den==0, 1, den)
      B[i,] <- filter[i,]*Q[,states[i+1]]
      B[i,] <- B[i,] / den
      states[i] <- which.max(B[i,])
    }
  res <- list(states=states, prob.states=B)
  res
}

simulateRJaCGH <- function(n, x=NULL, mu, sigma.2, beta, start, q=-beta) {
  y <- rep(NA, n)
  states <- rep(NA, n)
  k <- length(mu)
  if (is.null(x)) x <- rnorm(n-1, 0, 1)
  states[1] <- start
  y[1] <- rnorm(1, mu[states[1]], sqrt(sigma.2[states[1]]))
  for (i in 2:n) {
    Q <- Q.NH(q=q, beta=beta, x=x[i-1])
    states[i] <- sample(x=1:k, size=1, prob=Q[states[i-1],])
    y[i] <- rnorm(1, mu[states[i]], sqrt(sigma.2[states[i]]))
  }
  list(states=states, y=y)
}

#########################################################################
## Plot transition probabilities
#########################################################################

plot.Q.NH <- function(x, beta, q=-beta, col=NULL,...) {
  n<- length(x)
  k <- nrow(q)
  prob <- matrix(NA, n, k)
  if (is.null(col)) {
    col <- rep(1, k)
    col[grep("[G]", rownames(beta))] <- 2
    col[grep("[L]", rownames(beta))] <- 3
  }
  for (i in 1:n) {
    prob[i,] <- diag(Q.NH(q=q, beta=beta, x=x[i]))
  }
  plot(prob[order(x),1] ~ sort(x), type="l", ylim=c(min(prob), max(prob)), col=col[1], ...)
  ## case with homogeneous probabilities
  if (min(x)==0 && max(x)==0) abline(h=prob[1,1], col=col[1])
  if (k >1) {
    for (i in 2:k) {
      lines(prob[order(x),i] ~ sort(x), col=col[i])
        ## case with homogeneous probabilities
      if (min(x)==0 && max(x)==0) abline(h=prob[1,i], col=col[i])

    }
  }
  rug(x)
}


MetropolisSweep.C <- function(y, x, k.max, Chrom, model=NULL, burnin, TOT, prob.k, pb, ps, g, ka, mu.alfa,
                              mu.beta, sigma.tau.mu, sigma.tau.sigma.2,
                              sigma.tau.beta, tau.split.mu=NULL,
                              tau.split.beta=NULL, stat, start.k,
                              RJ=TRUE, auto.label=NULL) {
  if (!is.null(auto.label) && auto.label < 0 && auto.label > 1)
    stop("'auto.label' must be NULL or a number between 0 and 1")
  n <- length(y)
  ##Size of vectors
  size.mu <- TOT * k.max*(k.max+1)/2
  size.sigma.2 <- TOT * k.max*(k.max+1)/2
  size.beta <- TOT * k.max * (k.max+1) * (2*k.max+1) / 6
  mu <- rep(0, size.mu)
  sigma.2 <- rep(0, size.sigma.2)
  beta <- rep(0, size.beta)
  probStates <- rep(0, n*(k.max^2 - k.max) / 2)
  probJointStates <- rep(0, (n-1)*(k.max^2 + k.max -2) / 2)
  loglik <- rep(0, TOT * k.max)
  if (is.null(start.k)) start.k <- 0

  ## checks
  if(k.max != length(sigma.tau.mu)) stop("k.max != length(sigma.tau.mu)")
  if(k.max != length(sigma.tau.sigma.2)) stop("k.max != length(sigma.tau.sigma.2)")
  if(k.max != length(sigma.tau.beta)) stop("k.max != length(sigma.tau.beta)")
  
  if (model=="Chrom") {
    ## Impose a restriction to the variance
    maxVar <- var(y)
    if(length(y) != (length(x) + 1))
        stop("Length of vector of distances and data vector are different in MetropolisSweep.C (Chrom)!")
    if (.__RJACGH_DEBUG) cat("\n sigma.tau.mu ", sigma.tau.mu, "\n")
    res <- .C("MetropolisSweep", y=as.double(y), x=as.double(x), genome=as.integer(1), index=as.integer(c(0, length(y))),
              kMax=as.integer(k.max), n=as.integer(length(y)), burnin=as.integer(burnin),
              TOT=as.integer(TOT), times=as.integer(rep(0, k.max)),burninTimes=as.integer(rep(0,k.max)),
              probB=as.integer(0), probD=as.integer(0),
              probS=as.integer(0),
              probC=as.integer(0), probK=as.double(prob.k),
              pb=as.double(pb), ps=as.double(ps), g=as.double(g),
              ka=as.double(ka), muAlfa=as.double(mu.alfa),
              muBeta=as.double(mu.beta), sigmaTauMu=as.double(sigma.tau.mu),
              sigmaTauSigma2=as.double(sigma.tau.sigma.2),
              sigmaTauBeta=as.double(sigma.tau.beta),
              tauSplitMu=as.double(tau.split.mu),
              tauSplitBeta=as.double(tau.split.beta),
              k=as.integer(rep(0, 2*(TOT-1)+1)), mu=as.double(mu),
              sigma2=as.double(sigma.2), beta=as.double(beta),
              stat=as.double(stat), startK=as.integer(start.k),
              RJ=as.integer(1*RJ), maxVar=as.double(maxVar),
              probStates=as.double(probStates),
              probJointStates=as.double(probJointStates),
              loglik=as.double(loglik))
  }
  else {
    index <- diff(Chrom)
    index <- which(index>0)
    index <- c(0, index, length(y))
    maxVar <- var(y)
    if(length(y) != (length(x) + 1))
        stop("Length of vector of distances and data vector are different in MetropolisSweep.C!")
    res <- .C("MetropolisSweep", y=as.double(y), x=as.double(x), genome=as.integer(length(index)-1), index=as.integer(index), 
              kMax=as.integer(k.max), n=as.integer(length(y)), burnin=as.integer(burnin),
              TOT=as.integer(TOT), times=as.integer(rep(0, k.max)),burninTimes=as.integer(rep(0,k.max)),
              probB=as.integer(0), probD=as.integer(0),
              probS=as.integer(0),
              probC=as.integer(0), probK=as.double(prob.k),
              pb=as.double(pb), ps=as.double(ps), g=as.double(g),
              ka=as.double(ka), muAlfa=as.double(mu.alfa),
              muBeta=as.double(mu.beta), sigmaTauMu=as.double(sigma.tau.mu),
              sigmaTauSigma2=as.double(sigma.tau.sigma.2),
              sigmaTauBeta=as.double(sigma.tau.beta),
              tauSplitMu=as.double(tau.split.mu),
              tauSplitBeta=as.double(tau.split.beta),
              k=as.integer(rep(0, 2*(TOT-1)+1)), mu=as.double(mu),
              sigma2=as.double(sigma.2), beta=as.double(beta),
              stat=as.double(stat), startK=as.integer(start.k),
              RJ=as.integer(1*RJ), maxVar=as.double(maxVar),
              probStates=as.double(probStates),
              probJointStates=as.double(probJointStates),
              loglik=as.double(loglik))

  }
  ##Reconstruct objects
  obj <- list()
  gc()
  indexStat <- 1
  indexStates <- 1
  indexJointStates <- 1
  indexMu <- 1
  indexBeta <- 1
  indexLoglik <- 1
  for (i in 1:k.max) {
    obj[[i]] <- list()
    obj[[i]]$stat <- res$stat[indexStat:(indexStat + i-1)]
    obj[[i]]$mu <- matrix(res$mu[indexMu: (indexMu + res$times[i]*i -1)], ncol=i, byrow=TRUE)
    obj[[i]]$sigma.2 <- matrix(res$sigma2[indexMu: (indexMu + res$times[i]*i -1)], ncol=i, byrow=TRUE)
    indexMu <-  indexMu + TOT*i

    obj[[i]]$loglik <- res$loglik[indexLoglik:(indexLoglik +
                                               res$times[i] -1)]
    indexLoglik <- indexLoglik + TOT
    obj[[i]]$beta <- array(res$beta[indexBeta: (indexBeta + res$times[i]*i*i-1)], dim=c(i, i, res$times[i]))
    indexBeta <- indexBeta + TOT*i*i
    if (burnin >0) {
      obj[[i]]$mu <- matrix(obj[[i]]$mu[-c(1:res$burninTimes[i]),], ncol=i)
      obj[[i]]$sigma.2 <- matrix(obj[[i]]$sigma.2[-c(1:res$burninTimes[i]),], ncol=i)
      obj[[i]]$beta <-
        array(obj[[i]]$beta[,,-c(1:res$burninTimes[i])], dim=c(i,i, res$times[i] - res$burninTimes[i]))
      obj[[i]]$loglik <- obj[[i]]$loglik[-c(1:res$burninTimes[i])]
    }
    indexStat <- indexStat + i
    if (i >1) {
      if (nrow(obj[[i]]$mu) > 0) {
        obj[[i]]$prob.states <- res$probStates[indexStates : (indexStates +
                                                              n*(i-1) -1)]
        obj[[i]]$prob.states <- matrix(obj[[i]]$prob.states, nrow=n, ncol=i-1)
        obj[[i]]$prob.states <- cbind(obj[[i]]$prob.states, 1 -
                                      rowSums(obj[[i]]$prob.states))
        obj[[i]]$prob.joint.states <- res$probJointStates[indexJointStates : (indexJointStates +
                                                              (n-1)*i -1)]
        obj[[i]]$prob.joint.states <- matrix(obj[[i]]$prob.joint.states, nrow=n-1, ncol=i)
      }
      else {
        obj[[i]]$prob.states <- NULL
        obj[[i]]$prob.joint.states <- NULL
      }
      indexJointStates <- indexJointStates + (n-1)*(i)
    }
    else  {
      obj[[i]]$prob.states <- matrix(rep(1, length(y)), ncol=1)
      obj[[i]]$prob.joint.states <- matrix(rep(1, length(y)-1), ncol=1)
    }
    
    indexStates <- indexStates + n*(i-1)

    
  }
  
  obj[[1]]$prob.mu <- length(unique(obj[[1]]$mu)) / length(obj[[1]]$mu)
  obj[[1]]$prob.sigma.2 <- length(unique(obj[[1]]$sigma.2)) / length(obj[[1]]$sigma.2)
  obj[[1]]$prob.beta <- length(unique(obj[[1]]$beta)) / length(obj[[1]]$beta)
  for (i in 2:k.max) {
    obj[[i]]$prob.mu <- length(unique(obj[[i]]$mu[,1])) / length(obj[[i]]$mu[,1])
    obj[[i]]$prob.sigma.2 <- length(unique(obj[[i]]$sigma.2[,1])) / length(obj[[i]]$sigma.2[,1])
    obj[[i]]$prob.beta <- length(unique(obj[[i]]$beta[2,1,])) / length(obj[[i]]$beta[2,1,])
  }
  
  obj$k <- res$k[-c(1:(2*burnin+1))]
  obj$k <- factor(obj$k , levels=1:k.max)
  ## Relabel states
  p.labels <- 0.25 ## should be a parameter??
  for (i in 1:k.max) {
    if (table(obj$k)[i] > 0) {
      obj.sum <- summary.RJaCGH(obj, k=i, point.estimator="median")
      normal.levels <- (qnorm(mean=obj.sum$mu,
                              sd=sqrt(obj.sum$sigma.2),
                              p=p.labels) < 0 &
                        qnorm(mean=obj.sum$mu, sd=sqrt(obj.sum$sigma.2),
                              p=1-p.labels) > 0)
      n.Norm <- sum(normal.levels)
      if (n.Norm <=0) {
        normal.levels <- rep(FALSE, i)
        normal.levels[which.min(abs(obj.sum$mu))] <- TRUE
        n.Norm <- sum(normal.levels)
      }
      n.Loss <- which.max(normal.levels) -1
      n.Gain <- i - max(normal.levels*(1:i))
      obj[[i]]$state.labels <- NULL
      if(n.Loss>0)
        obj[[i]]$state.labels <- c(obj[[i]]$state.labels,
                                   paste("Loss", n.Loss:1, sep="-"))
      obj[[i]]$state.labels <- c(obj[[i]]$state.labels, rep("Normal", n.Norm))
      if(n.Gain>0)
        obj[[i]]$state.labels <- c(obj[[i]]$state.labels,
                                   paste("Gain", 1:n.Gain, sep="-"))
      ## Should automatic labelling include the former steps?
      
      if (!is.null(auto.label)) {
        states <- states.RJaCGH(obj, k=i)$states
        states <- factor(states, levels=names(summary.RJaCGH(obj,
                                   k=i)$mu))
        freqs <- prop.table(table(states))
        labels <- names(freqs)
        means <- obj.sum$mu
        means.order <- order(abs(means))
        tot.norm <- freqs['Normal']
        ind.tot <- 1
        while(tot.norm <= auto.label) {
          labels[means.order][ind.tot+1] <- "Normal"
          tot.norm <- tot.norm + freqs[means.order][ind.tot+1]
          ind.tot <- ind.tot + 1
        }
        obj[[i]]$state.labels <- labels
      }
  
  ##
      
    }
  }    

  ##
  
  obj$prob.b <- res$probB
  obj$prob.d <- res$probD
  obj$prob.s <- res$probS
  obj$prob.c <- res$probC
  obj$y <- y
  obj$x <- x
  obj$x[obj$x==-1] <- NA
  rm(res)
  gc()
  attr(obj, "auto.label") <- auto.label
  obj

}


RJMCMC.NH.HMM.Metropolis <- function(y, Chrom=NULL, x=NULL, index=NULL, model=NULL, burnin=0, TOT=1000, k.max=6,
                                     stat=NULL, mu.alfa=NULL, mu.beta=NULL, ka=ka, g=g, prob.k=NULL,
                                     sigma.tau.mu, sigma.tau.sigma.2, sigma.tau.beta,
                                     tau.split.mu, tau.split.beta,
                                     start.k, RJ=RJ, auto.label=NULL) {

  if (k.max==1) {
    cat("k.max must be 2 or more\n")
    stop()
  }
  TOT <- burnin + TOT
  n <- length(y)
  k <- rep(NA, 2*TOT)
  times <- rep(2, k.max)
  res <- list()
  x.old <- x
  if (is.null(x) | isTRUE(all.equal(x.old, rep(1, length(y)-1)))) {
    x <- rep(0, length(y)-1) ## zz:??? a "rep"?
  }
  ## Scale x'2 to avoid overflow
  else if(max(x, na.rm=TRUE)!=0){
    x <- x/max(x, na.rm=TRUE)
  }

  ## convert NA's to -1
  x[is.na(x)] <- -1
  
  ##Hyperparameters
  if(is.null(mu.alfa)) mu.alfa <- median(y)
  if(is.null(mu.beta)) mu.beta <- diff(range(y))
  if(is.null(ka)) ka <- 2
  if(is.null(g)) g <- diff(range(y))^2 / 50


  ## Prior over the number of states  ####################

  if (is.null(prob.k)) {
    prob.k <- rep(1/k.max, k.max)
  }
  pb <- c(1,rep(1/2, k.max-2),0)
  ps <- c(1, rep(1/2, k.max-2),0)

  if (is.null(stat)) {
    for(i in 1:k.max)  stat <- c(stat, rep(1/i, i))
  }
if(is.null(sigma.tau.mu) | is.null(sigma.tau.sigma.2) | is.null(sigma.tau.beta)
   | is.null(tau.split.mu) | is.null(tau.split.beta)) {
  cat("Searching jump parameters...\n")
  if(length(y) != (length(x) + 1))
      stop("Length of vector of distances and data vector are different in get jump!")

  params <- get.jump(y=y, x=x, k.max=k.max, Chrom=Chrom, model=model,
                           prob.k=prob.k, pb=pb, ps=ps, g=g, 
                           ka=ka, mu.alfa=mu.alfa, mu.beta=mu.beta, stat=stat)
  if (is.null(sigma.tau.mu)) sigma.tau.mu <- params$sigma.tau.mu
  if (is.null(sigma.tau.sigma.2)) sigma.tau.sigma.2 <- params$sigma.tau.sigma.2
  if (is.null(sigma.tau.beta)) sigma.tau.beta <- params$sigma.tau.beta
  if(is.null(tau.split.mu)) tau.split.mu <- mean(params$sigma.tau.mu) ^2
  if(is.null(tau.split.beta)) tau.split.beta <- mean(params$sigma.tau.beta) ^2
}
 
  cat("Starting Reversible Jump\n")
  if(length(y) != (length(x) + 1))
      stop("Length of vector of distances and data vector are different in metropolis!")
  
  res <- MetropolisSweep.C(y=y, x=x, k.max=k.max, Chrom=Chrom,
                           model=model, burnin=burnin, TOT=TOT,
                           prob.k=prob.k, pb=pb, ps=ps, g=g, 
                           ka=ka, mu.alfa=mu.alfa, mu.beta=mu.beta,
                           sigma.tau.mu=sigma.tau.mu,
                           sigma.tau.sigma.2=sigma.tau.sigma.2,
                           sigma.tau.beta=sigma.tau.beta,
                           tau.split.mu=tau.split.mu,
                           tau.split.beta=tau.split.beta, stat=stat,
                           start.k=start.k, RJ=RJ, auto.label=auto.label)
  class(res) <- "RJaCGH"
  res
}



RJaCGH.one.array <- function(y, Chrom=NULL, Pos=NULL, Dist=NULL, model="genome", burnin=0, TOT=1000, k.max=6,
                             stat=NULL, mu.alfa=NULL, mu.beta=NULL, ka, g, prob.k=NULL,
                             sigma.tau.mu, sigma.tau.sigma.2, sigma.tau.beta,
                             tau.split.mu, tau.split.beta,
                             start.k, RJ=RJ, auto.label=NULL) {

  ## Check that Positions are absolute and not relative
  if(!is.null(Pos)) {
    if (!is.null(Chrom) && any(diff(Pos)<0)) {
      last.Pos <- 0
      for (i in 2:length(Pos)) {
        if (Chrom[i] != Chrom[i-1]) {
          last.Pos <- Pos[i-1]
        }
        Pos[i] <- Pos[i] + last.Pos
      }
    }
  }


  ## Model without Chromosomes
  if (is.null(Chrom)) {
    if (!is.null(Pos)) {
      chrom.Pos <- Pos
      chrom.Dist <- diff(chrom.Pos) -1 
    }
    else {
      chrom.Pos <- 1:length(y)
      chrom.Dist <- rep(0, length(y)-1)
    }
    if (!is.null(Dist)) chrom.Dist <- Dist
    res <- RJMCMC.NH.HMM.Metropolis(y=y, Chrom=rep(1, length(y)), x=chrom.Dist, burnin=burnin,
                                    TOT=TOT, k.max=k.max, stat=stat, mu.alfa=mu.alfa, mu.beta=mu.beta,
                                    ka=ka, g=g, prob.k=prob.k, sigma.tau.mu=sigma.tau.mu,
                                    sigma.tau.sigma.2=sigma.tau.sigma.2,
                                    sigma.tau.beta=sigma.tau.beta, model=model,
                                    tau.split.mu=tau.split.mu,  
                                    tau.split.beta=tau.split.beta,
                                    start.k=start.k, RJ=RJ, auto.label=auto.label)
      res$Pos <- chrom.Pos
      res
    }
  else {
    if (!is.numeric(Chrom)) {
      cat("Chrom must be numeric\n")
      stop()
    }
    
    ## ###################################
    ##different model for each chromosome

    if (model=="Chrom") {
      res <- list()
      chrom.names <- as.numeric(names(table(Chrom)))
      for (i in chrom.names) {
        ## x should be preprocessed here
        if (is.null(Pos)) {
          chrom.Pos <- (1:length(y))[Chrom==i]
          chrom.Dist <- rep(0, length(chrom.Pos - 1))
        }
        else {
          chrom.Pos <- Pos[Chrom==i]
          chrom.Dist <- diff(chrom.Pos) - 1
          chrom.Dist[chrom.Dist<0] <- 0
        }
        cat("Chromosome", i, "\n")
        if (!is.null(Dist)) chrom.Dist <- Dist[Chrom==i]
        res[[i]] <-
          RJMCMC.NH.HMM.Metropolis(y=y[Chrom==i], Chrom=Chrom[Chrom==i], x=chrom.Dist, burnin=burnin,
                                   TOT=TOT, k.max=k.max, stat=stat, mu.alfa=mu.alfa, mu.beta=mu.beta,
                                   ka=ka, g=g, prob.k=prob.k, sigma.tau.mu=sigma.tau.mu,
                                   sigma.tau.sigma.2=sigma.tau.sigma.2,
                                   sigma.tau.beta=sigma.tau.beta, model=model,
                                   tau.split.mu=tau.split.mu,
                                   tau.split.beta=tau.split.beta,
                                   start.k=start.k, RJ=RJ, auto.label=auto.label)
        res[[i]]$Pos <- chrom.Pos
        class(res[[i]]) <- "RJaCGH"

      }
      res$model <- model
      if (is.null(Pos)) res$Pos <- 1:length(y)
      else res$Pos <- Pos
      res$Chrom <- Chrom
      class(res) <- "RJaCGH.Chrom"
      res
    }
    
    ## ###################################
    ## Same model for each chromosome
    
    else {
      if (is.null(Pos)) {
        Pos <- 1:length(y)
        Dist <- rep(0, length(y)-1)
      }
      else {
        Dist <- diff(Pos) -1
        for (i in 1:23) 
          Dist[Chrom==i][length(Dist[Chrom==i])] <- NA
        Dist <- Dist[-length(Dist)]
      }
      ## We'll have to take out the last Dist of every Chrom
      res <- RJMCMC.NH.HMM.Metropolis(y=y, Chrom=Chrom, x=Dist, burnin=burnin, TOT=TOT,
                                      k.max=k.max, stat=stat, mu.alfa=mu.alfa, mu.beta=mu.beta,
                                      ka=ka, g=g, prob.k=prob.k, sigma.tau.mu=sigma.tau.mu,
                                      sigma.tau.sigma.2=sigma.tau.sigma.2,
                                      sigma.tau.beta=sigma.tau.beta, model=model,
                                      tau.split.mu=tau.split.mu, 
                                      tau.split.beta=tau.split.beta,
                                      start.k=start.k, RJ=RJ, auto.label=auto.label)
      res$Pos <- Pos
      res$model <- model
      res$Chrom <- Chrom
      class(res) <- "RJaCGH.genome"
      res
    }
  }
  res
}

RJaCGH <- function(y, Chrom=NULL, Pos=NULL, Dist=NULL, model="genome", burnin=10000, TOT=10000, k.max=6,
                  stat=NULL, mu.alfa=NULL, mu.beta=NULL, ka=NULL,
                   g=NULL, prob.k=NULL, jump.parameters=list(),
                   start.k=NULL, RJ=TRUE, auto.label=0.60) {
  sigma.tau.mu <- jump.parameters$sigma.tau.mu
  sigma.tau.sigma.2 <- jump.parameters$sigma.tau.sigma.2
  sigma.tau.beta <- jump.parameters$sigma.tau.beta
  tau.split.mu <- jump.parameters$tau.split.mu
  tau.split.beta <- jump.parameters$tau.split.beta
  ## Check if we have 1 array or several
  if(is.null(dim(y))) {
    res <- RJaCGH.one.array(y, Chrom=Chrom, Pos=Pos, Dist=NULL, model=model, burnin=burnin, TOT=TOT, k.max=k.max,
                            stat=stat, mu.alfa=mu.alfa, mu.beta=mu.beta,
                            ka=ka, g=g, prob.k=prob.k,
                            sigma.tau.mu=sigma.tau.mu, sigma.tau.sigma.2=sigma.tau.sigma.2,
                            sigma.tau.beta=sigma.tau.beta, tau.split.mu=tau.split.mu,
                            tau.split.beta=tau.split.beta,
                            start.k=start.k, RJ=RJ, auto.label=auto.label)
    res
  }
  else {
    res <- list()
    res$array.names <- NULL
    for (i in 1:ncol(y)) {
      cat("array", colnames(y)[i], "\n")
      res[[colnames(y)[i]]] <-
        RJaCGH.one.array(y[,i], Chrom=Chrom,
                         Pos=Pos, Dist=NULL, model=model, burnin=burnin, TOT=TOT, k.max=k.max,
                         stat=stat, mu.alfa=mu.alfa, mu.beta=mu.beta, ka=ka, g=g, prob.k=prob.k,
                         sigma.tau.mu=sigma.tau.mu, sigma.tau.sigma.2=sigma.tau.sigma.2, sigma.tau.beta=sigma.tau.beta,
                         tau.split.mu=tau.split.mu,
                         tau.split.beta=tau.split.beta,
                         start.k=start.k, RJ=RJ, auto.label=auto.label)

      res$array.names[i] <- colnames(y)[i]
    }
    class(res) <- "RJaCGH.array"
    res
  }
     
}

summary.RJaCGH <- function(object, k=NULL, point.estimator="median", ...) {
  res <- list()
  res$y <- object$y
  res$Pos <- object$Pos
  if (is.null(k)) {
    k <- as.numeric(names(which.max(table(object$k))))
  }
  res$stat <- object[[k]]$stat
  if (point.estimator=="mode") {
      dens <- apply(object[[k]]$mu, 2, density, bw="nrd0")
      res$mu <- unlist(lapply(dens, function(x) x$x[which.max(object$y)]))
      dens <- apply(object[[k]]$sigma.2, 2, density, bw="nrd0")
      res$sigma.2 <- unlist(lapply(dens, function(x) x$x[which.max(object$y)]))
      dens <- apply(x[[k]]$beta, c(1,2), density, bw="nrd0")
      res$beta <- unlist(lapply(dens, function(x) x$x[which.max(object$y)]))
      res$beta <- matrix(res$beta, k)
      diag(res$beta) <- 0
    }

  else{
      res$mu <- apply(matrix(object[[k]]$mu, ncol=k), 2, point.estimator)
      res$sigma.2 <- apply(matrix(object[[k]]$sigma.2, ncol=k), 2, point.estimator)
      res$beta <- apply(object[[k]]$beta, c(1,2), point.estimator)
    }
  if (is.null(object[[k]]$state.labels)) {
  ref <- as.numeric(names(which.max(table(apply(object[[k]]$prob.states, 1,
      which.max)))))
  names(res$mu)[ref] <- "Normal"
  names(res$sigma.2)[ref] <- "Normal"
  rownames(res$beta) <- 1:k
  colnames(res$beta) <- 1:k
  rownames(res$beta)[ref] <- "Normal"
  colnames(res$beta)[ref] <- "Normal"
  names(res$stat)[ref] <- "Normal"
  if (ref < k) {
    names(res$mu)[(ref+1):k] <- paste("Gain", 1:(k-ref), sep="-")
    names(res$sigma.2)[(ref+1):k] <- paste("Gain", 1:(k-ref), sep="-")
    rownames(res$beta)[(ref+1):k] <- paste("Gain", 1:(k-ref), sep="-")
    colnames(res$beta)[(ref+1):k] <- paste("Gain", 1:(k-ref), sep="-")
    names(res$stat)[(ref+1):k] <- paste("Gain", 1:(k-ref), sep="-")
  }
  if (ref > 1) {
    names(res$mu)[1:(ref-1)] <- paste("Loss", (ref-1):1, sep="-")
    names(res$sigma.2)[1:(ref-1)] <- paste("Loss", (ref-1):1, sep="-")
    rownames(res$beta)[1:(ref-1)] <- paste("Loss", (ref-1):1, sep="-")
    colnames(res$beta)[1:(ref-1)] <- paste("Loss", (ref-1):1, sep="-")
    names(res$stat)[1:(ref-1)] <- paste("Loss", (ref-1):1, sep="-")
  }
}
  else {
    names(res$mu) <- object[[k]]$state.labels
    names(res$sigma.2) <- object[[k]]$state.labels
    rownames(res$beta) <- object[[k]]$state.labels
    colnames(res$beta) <- object[[k]]$state.labels
    names(res$stat) <- object[[k]]$state.labels
  }
  res
}


summary.RJaCGH.Chrom <- function(object, point.estimator="median", ...) {

  res <- list()
  res$y <- object$y
  res$Pos <- object$Pos
  
  for (i in 1:max(object$Chrom)) {
    k <- as.numeric(names(which.max(table(object[[i]]$k))))
    res[[i]] <- summary(object[[i]], k=k, point.estimator=point.estimator)
  }
  res
}

summary.RJaCGH.genome <- function(object, k=NULL,
  point.estimator="median", ...) {
  res <- summary.RJaCGH(object, k, point.estimator)
  res
}

summary.RJaCGH.array <- function(object, point.estimator="median", ...) {
  res <- list()
  for (i in object$array.names) {
    res[[i]] <- summary(object[[i]], point.estimator=point.estimator)
  }
  res
}

states <- function(obj, k) {
  UseMethod("states")
}
## Now does the same with joint probs.
states.RJaCGH <- function(obj, k=NULL) {
  res <- NULL
  if (is.null(k)) {
    k <- as.numeric(names(which.max(table(obj$k))))
  }
  if (nrow(obj[[k]]$mu)==0) stop ("No observations in that HMM\n")
  res$states <- apply(obj[[k]]$prob.states, 1, which.max)
  res$states <- factor(res$states, levels=1:k)
  res$prob.states <- obj[[k]]$prob.states
  res$prob.joint.states <- obj[[k]]$prob.joint.states

  ## Region of normal, gain and loss
  if (is.null(obj[[k]]$state.labels)) {
    ref <- as.numeric(names(which.max(table(res$states))))
    colnames(res$prob.states) <- 1:k
    colnames(res$prob.states)[ref] <- "Normal"
    levels(res$states)[ref] <- "Normal"
    if (ref < k) {
      colnames(res$prob.states)[(ref+1):k] <- paste("Gain", 1:(k-ref), sep="-")
      levels(res$states)[(ref+1):k] <- paste("Gain", 1:(k-ref), sep="-")
    }
    if (ref > 1) {
      colnames(res$prob.states)[1:(ref-1)] <- paste("Loss", (ref-1):1, sep="-")
      levels(res$states)[1:(ref-1)] <- paste("Loss", (ref-1):1, sep="-")
    }
  }
  else {
    colnames(res$prob.states) <- obj[[k]]$state.labels
    levels(res$states) <- obj[[k]]$state.labels
  }
  colnames(res$prob.joint.states) <- colnames(res$prob.states)
  res <- list(states=res$states, prob.states=res$prob.states)
}

states.RJaCGH.Chrom <- function(obj, k=NULL) {
  res <- list()
  for (i in 1:max(obj$Chrom)) {
    k <- as.numeric(names(which.max(table(obj[[i]]$k))))
    res[[i]] <- states(obj[[i]], k=k)
  }
  res
}

states.RJaCGH.genome <- function(obj, k=NULL) {

 states.RJaCGH(obj, k=k)

}

states.RJaCGH.array <- function(obj, k=NULL) {
  res <- list()
  for (i in obj$array.names) {
    res[[i]] <- states(obj[[i]], k=k)
  }
  res
}

model.averaging <- function(obj) {
  UseMethod("model.averaging")
}

model.averaging.RJaCGH <-function(obj) {
  res <- list()
  n <- length(obj$y)
  probs <- prop.table(table(obj$k))
  Loss <- rep(0, n)
  Normal <- rep(0, n)
  Gain <- rep(0, n)
  joint.Loss <- rep(0, n-1)
  joint.Normal <- rep(0, n-1)
  joint.Gain <- rep(0, n-1)

  k.max <- max(as.numeric(as.character(obj$k)))
  for (i in 1:k.max) {
    if (probs[i] > 0) {
      objSummary <- states(obj, i)
      prob.states <- objSummary$prob.states
##      prob.joint.states <- objSummary$prob.joint.states
      for (j in 1:i) {
        if (length(grep("[L]", colnames(prob.states)[j]))) {
          Loss <- Loss + probs[i] * prob.states[,j]
##          joint.Loss <- joint.Loss + probs[i] * prob.joint.states[,j]
        }
        else if (length(grep("[N]", colnames(prob.states)[j]))) {
          Normal <- Normal + probs[i] * prob.states[,j]
##          joint.Normal <- joint.Normal + probs[i] * prob.joint.states[,j]
        }
        else if (length(grep("[G]", colnames(prob.states)[j]))) {
          Gain <- Gain + probs[i] * prob.states[,j]
##          joint.Gain <- joint.Gain + probs[i] * prob.joint.states[,j]
        }
      }
    }
  }
  res$prob.states <- cbind(Loss, Normal, Gain)
  ##res$prob.joint.states <- cbind(joint.Loss, joint.Normal, joint.Gain)
  colnames(res$prob.states) <- c("Loss", "Normal", "Gain")
  ##colnames(res$prob.joint.states) <- c("Loss", "Normal", "Gain")
  res$states <- apply(res$prob.states, 1, function(x) names(which.max(x)))
  res$states <- ordered(res$states, levels=c("Loss", "Normal", "Gain"))
  res <- list(states=res$states, prob.states=res$prob.states)
}

model.averaging.RJaCGH.Chrom <-function(obj) {
  res <- list()
  for (ch in 1:23) {
    res[[ch]] <- model.averaging(obj[[ch]])
  }
  res
}

model.averaging.RJaCGH.genome <- function(obj) {
  res <- list()
  n <- length(obj$y)
  probs <- prop.table(table(obj$k))
  Loss <- rep(0, n)
  Normal <- rep(0, n)
  Gain <- rep(0, n)
  joint.Loss <- rep(0, n-1)
  joint.Normal <- rep(0, n-1)
  joint.Gain <- rep(0, n-1)
  k.max <- max(as.numeric(as.character(obj$k)))  
  for (i in 1:k.max) {
    if (probs[i] > 0) {
      objSummary <- states(obj, i)
      prob.states <- objSummary$prob.states
    ##  prob.joint.states <- objSummary$prob.joint.states
      for (j in 1:i) {
        if (length(grep("[L]", colnames(prob.states)[j]))) {
          Loss <- Loss + probs[i] * prob.states[,j]
      ##    joint.Loss <- joint.Loss + probs[i] * prob.joint.states[,j]
        }
        else if (length(grep("[N]", colnames(prob.states)[j]))) {
          Normal <- Normal + probs[i] * prob.states[,j]
        ##  joint.Normal <- joint.Normal + probs[i] * prob.joint.states[,j]
        }
        else if (length(grep("[G]", colnames(prob.states)[j]))) {
          Gain <- Gain + probs[i] * prob.states[,j]
        ##  joint.Gain <- joint.Gain + probs[i] * prob.joint.states[,j]
        }
      }
    }
  }
  res$prob.states <- cbind(Loss, Normal, Gain)
  ##res$prob.joint.states <- cbind(joint.Loss, joint.Normal, joint.Gain)
  colnames(res$prob.states) <- c("Loss", "Normal", "Gain")
  ##colnames(res$prob.joint.states) <- c("Loss", "Normal", "Gain")
  res$states <- apply(res$prob.states, 1, function(x) names(which.max(x)))
  res$states <- ordered(res$states, levels=c("Loss", "Normal", "Gain"))
  res <- list(states=res$states, prob.states=res$prob.states)
}

model.averaging.RJaCGH.array <- function(obj) {
  res <- list()
  for (i in obj$array.names) {
    res[[i]] <- model.averaging(obj[[i]])
  }
  res
}


plot.RJaCGH <- function(x, k=NULL, model.averaging=TRUE, cex=1, ...)  {

  if (is.null(k)) {
    k <- as.numeric(names(which.max(table(x$k))))
  }

  if (model.averaging) {
    res <- model.averaging(x)
    ## It's needed to get the classification with k states
    ## and access the colors of the first graphics
    cols <- colnames(states(x, k)$prob.states)
  }

  else {
    res <- states(x, k)
    cols <- colnames(res$prob.states)
  }

  layout(matrix(c(1,3,5,2,4, 5),3,2), heights=c(0.25, 0.25, 0.5))
  par(cex.lab=cex*0.8, cex.main=cex*0.9, cex.axis=cex*0.8)
  par(mar=c(5,4,4,4) + 0.1)
  barplot(prop.table(table(x$k)), xlab="# of hidden states", ylab="probability",
          main="Prob. of number of hidden states", col="red", ylim=c(0,1))
  col <- rep(1, k)
  col[grep("[G]", cols)] <- 2
  col[grep("[L]", cols)] <- 3
  plot(density(x[[k]]$mu[,1], bw=sd(x[[k]]$mu[,1])), col=col[1], xlim=range(x$y),
       main="Posterior probability of mean of hidden states")
  if (k >1) for (i in 2:k) lines(density(x[[k]]$mu[,i],
                                         bw=sd(x[[k]]$mu[,i])), col=col[i])
  plot(density(x[[k]]$sigma.2[,1],
               bw=sd(x[[k]]$sigma.2[,1]), from=0), col=col[1], main="Posterior probability of variance of hidden states")
 if (k >1) for (i in 2:k) lines(density(x[[k]]$sigma.2[,i],
                                        bw=sd(x[[k]]$sigma.2[,i])), col=col[i])
  summary.obj <- summary(x, k)
  plot.Q.NH(q=-summary.obj$beta, beta=summary.obj$beta, x=x$x,
            main="Probability of permanence in the same hidden state", xlab="Distance", ylab="Prob.", col=col)

  if (model.averaging)
    main.text <- "Prediction of copy gain/loss. Bayesian Model Averaging"
  else
    main.text <- paste("Prediction of copy gain/loss. Number of hidden states: ", k)
  ylim <- range(x$y)
  margin <- diff(ylim) / 4
  ylim[1] <- ylim[1] - margin - 1
  ylim[2] <- ylim[2] + margin
  col <- rep(1, length(x$y))
  col[as.numeric(res$states) %in% grep("[G]", levels(res$states))] <- 2
  col[as.numeric(res$states) %in% grep("[L]", levels(res$states))] <- 3
  pch <- rep(16, length(x$y))
  pch[as.numeric(res$states) %in% grep("[2-9]", levels(res$states))] <- 17
  plot(x$y~x$Pos, pch=pch, col=col, ylim=ylim, ylab="Log2 ratio", xlab="Pos.Base",
       main=main.text)
  abline(h=0, lty=2)
  par(new=TRUE)
  plot(apply(res$prob.states, 1, max) ~ x$Pos, yaxt="n", xlab="", ylab="", col="blue", ylim=c(0,5), type="s")
  abline(h=c(0, 0.25, 0.5, 0.75,1), lty=2, yaxt="n")
  axis(4, summary.obj$prob.states, at=c(0, 0.5, 0.75, 1), labels=c(0, 0.5, 0.75, 1))
  mtext("Probability", 4, 2, cex=cex*0.6)
}

plot.RJaCGH.Chrom <- function(x, Chrom="genome",
                              model.averaging=TRUE, cex=1, ...)  {

  if (Chrom=="genome") {
    states <- NULL
    prob.states <- NULL
    if (model.averaging) {
      res <- model.averaging(x)
    }
    else {
      res <- states(x)
    }
    max.levels <- lapply(res, function(x) levels(x$states))
    max.levels <- names(table(unlist(max.levels)))
    states <- lapply(res, function(x) as.character(x$states))
    states <- factor(unlist(states), levels=max.levels)
    prob.states <- lapply(res, function(x) x$prob.states)
    prob.states <- do.call("rbind", prob.states)
    par(cex.lab=cex*0.8, cex.main=cex*0.9, cex.axis=cex*0.8)
    par(mar=c(5,4,4,4) + 0.1)
    main.text <- paste("Prediction of copy gain/loss.")
    if (model.averaging) main.text <- c(main.text, "\n Bayesian Model Averaging")
    else main.text <- c(main.text, paste("\n Number of hidden states: ", k))
    y <- unlist(lapply(x, function(x) x$y))
    Pos <- unlist(lapply(x, function(x) x$Pos))
    col <- rep(1, length(y))
    col[as.numeric(states) %in% grep("[G]", levels(states))] <- 2
    col[as.numeric(states) %in% grep("[L]", levels(states))] <- 3
    pch <- rep(16, length(y))
    pch[as.numeric(states) %in% grep("[2-9]", levels(states))] <- 17
    ylim <- range(y)
    margin <- diff(ylim) / 4
    ylim[1] <- ylim[1] - margin - 1
    ylim[2] <- ylim[2] + margin
    ## Start of every Chromosome
    start.Chrom <- c(Pos[1], Pos[diff(x$Chrom)!=0], Pos[length(Pos)])
    xx <- rep(start.Chrom, rep(2, length(start.Chrom)))
    yy <- rep(c(ylim, rev(ylim)), length(start.Chrom)/2)
    plot(y~Pos, pch=pch, col=col, ylim=ylim, ylab="Log2 ratio", xlab="Pos.Base",
         main=main.text, type="n")
    polygon(xx, yy, col="gray85", xpd=TRUE)
    points(y~Pos, pch=pch, col=col)
    abline(h=0, lty=2)
    text(start.Chrom[-length(start.Chrom)] + diff(start.Chrom)/2,
         rep(ylim[1], length(start.Chrom)), labels=unique(x$Chrom),
         cex=cex*0.7, pos=1)
    par(new=TRUE)
    plot(apply(prob.states, 1, max) ~ Pos, yaxt="n", xlab="", ylab="", col="blue", ylim=c(0,5), type="s")
    abline(h=c(0, 0.25, 0.5, 0.75,1), lty=2, yaxt="n")
    axis(4, prob.states, at=c(0, 0.5, 0.75, 1), labels=c(0, 0.5, 0.75, 1))
    mtext("Probability", 4, 2, cex=cex*0.6)
    
  }
  else {
    if (Chrom=="X" || Chrom=="Y") Chrom==23
    plot(x[[Chrom]], model.averaging=model.averaging, cex=cex)
  }
}



plot.RJaCGH.genome <- function(x, k=NULL,
                               model.averaging=TRUE, cex=1, ...)  {

  if (is.null(k)) {
    k <- as.numeric(names(which.max(table(x$k))))
  }

  if (model.averaging) {
    res <- model.averaging(x)
    cols <- colnames(states(x, k)$prob.states)
  }
  else {
    res <- states(x, k)
    cols <- colnames(res$prob.states)
  }

  layout(matrix(c(1,3,5,2,4, 5),3,2), heights=c(0.25, 0.25, 0.5))
  par(cex.lab=cex*0.8, cex.main=cex*0.9, cex.axis=cex*0.8)
  par(mar=c(5,4,4,4) + 0.1)
  barplot(prop.table(table(x$k)), xlab="# of hidden states", ylab="probability",
          main="Prob. of number of hidden states", col="red", ylim=c(0,1))
  if (!model.averaging) {
    col <- rep(1, k)
    col[grep("[G]", cols)] <- 2
    col[grep("[L]", cols)] <- 3
  }
  else {
    names.colors <- colnames(states(x, k)$prob.states)
    col <- rep(1, k)
    col[grep("[G]", names.colors)] <- 2
    col[grep("[L]", names.colors)] <- 3

  }
  plot(density(x[[k]]$mu[,1], bw=sd(x[[k]]$mu[,1])), col=col[1], xlim=range(x$y),
       main="Posterior probability of mean of hidden states")
  if (k >1) for (i in 2:k) lines(density(x[[k]]$mu[,i],
                                         bw=sd(x[[k]]$mu[,i])), col=col[i])
  plot(density(x[[k]]$sigma.2[,1],
               bw=sd(x[[k]]$sigma.2[,1]), from=0), col=col[1], main="Posterior probability of variance of hidden states")
  if (k >1) for (i in 2:k) lines(density(x[[k]]$sigma.2[,i],
                                         bw=sd(x[[k]]$sigma.2[,i])), col=col[i])
  summary.obj <- summary(x, k)
  plot.Q.NH(q=-summary.obj$beta, beta=summary.obj$beta, x=x$x[!is.na(x$x)],
            main="Probability of permanence in the same hidden state", xlab="Distance", ylab="Prob.", col=col)

  if (model.averaging)
    main.text <- "Prediction of copy gain/loss. Bayesian Model Averaging"
  else
    main.text <- paste("Prediction of copy gain/loss. Number of hidden states: ", k)
  ylim <- range(x$y)
  margin <- diff(ylim) / 4
  ylim[1] <- ylim[1] - margin - 1
  ylim[2] <- ylim[2] + margin
  col <- rep(1, length(x$y))
  col[as.numeric(res$states) %in% grep("[G]", levels(res$states))] <- 2
  col[as.numeric(res$states) %in% grep("[L]", levels(res$states))] <- 3
  pch <- rep(16, length(x$y))
  pch[as.numeric(res$states) %in% grep("[2-9]", levels(res$states))] <- 17
  ## Start of every Chromosome
  start.Chrom <- c(x$Pos[1], x$Pos[diff(x$Chrom)!=0], x$Pos[length(x$Pos)])
  xx <- rep(start.Chrom, rep(2, length(start.Chrom)))
  yy <- rep(c(ylim, rev(ylim)), length(start.Chrom)/2)

  plot(x$y~x$Pos, pch=pch, col=col, ylim=ylim, ylab="Log2 ratio", xlab="Pos.Base",
       main=main.text, type="n")
  polygon(xx, yy, col="gray85", xpd=TRUE)
  points(x$y~x$Pos, pch=pch, col=col)
  abline(h=0, lty=2)
  text(start.Chrom[-length(start.Chrom)] + diff(start.Chrom)/2,
       rep(ylim[1], length(start.Chrom)), labels=unique(x$Chrom),
       cex=cex*0.7, pos=1)
  par(new=TRUE)
  plot(apply(res$prob.states, 1, max) ~ x$Pos, yaxt="n", xlab="", ylab="", col="blue", ylim=c(0,5), type="s")
  abline(h=c(0, 0.25, 0.5, 0.75,1), lty=2, yaxt="n")
  axis(4, summary.obj$prob.states, at=c(0, 0.5, 0.75, 1), labels=c(0, 0.5, 0.75, 1))
  mtext("Probability", 4, 2, cex=cex*0.6)
}

## There could be a weights argument for array precision
plot.RJaCGH.array <- function(x,method="averaging", weights=NULL,
                              cex=1, ...)  {

  array.names <- x$array.names
  par(mfrow=c(1,1))
  x$array.names <- NULL
  if (is.null(weights)) weights <- rep(1/length(array.names), length(array.names))
  weights <- weights / sum(weights)
  if (method=="averaging") {
    res <- lapply(x, model.averaging)

    if(class(x[[1]])=="RJaCGH.Chrom") {
      prob.states <- list()
      y <- list()
      k <- max(as.numeric(as.character(x[[1]]$Chrom)))
        for (i in 1:k) {
          prob.states[[i]] <- matrix(0, nrow=length(x[[1]][[i]]$y), ncol=3)
          y[[i]] <- rep(0, length(x[[1]][[i]]$y))
          for (j in 1:length(array.names)) {
            prob.states[[i]] <- prob.states[[i]] + res[[j]][[i]]$prob.states * weights[j]
            y[[i]] <- y[[i]] + x[[j]][[i]]$y * weights[j]
          }
        }
      prob.states <- do.call("rbind", prob.states)
      y <- do.call("c", y)
      states <- apply(prob.states, 1, function(x) names(which.max(x)))
      states <- factor(states, levels=c("Loss", "Normal", "Gain"))
    }
    else {
      prob.states <- matrix(0, nrow=length(x[[1]]$y), ncol=3)
      y <- rep(0, length(x[[1]]$y))
      for (j in 1:length(array.names)) {
        prob.states <- prob.states + res[[j]]$prob.states * weights[j]
        y <- y + x[[j]]$y * weights[j]
      }
      states <- apply(prob.states, 1, function(x) names(which.max(x)))
      states <- factor(states, levels=c("Loss", "Normal", "Gain"))
    }

    par(cex.lab=cex*0.8, cex.main=cex*0.9, cex.axis=cex*0.8)
    par(mar=c(5,4,4,4) + 0.1)
    main.text <- "Prediction of copy gain/loss. Bayesian Model Averaging"
    Pos <- x[[1]]$Pos
    col <- rep(1, length(y))
    col[as.numeric(states) %in% grep("[G]", levels(states))] <- 2
    col[as.numeric(states) %in% grep("[L]", levels(states))] <- 3
    pch <- rep(16, length(y))
    pch[as.numeric(states) %in% grep("[2-9]", levels(states))] <- 17
    ylim <- range(y)
    margin <- diff(ylim) / 4
    ylim[1] <- ylim[1] - margin - 1
    ylim[2] <- ylim[2] + margin
    ## Start of every Chromosome
    start.Chrom <- c(Pos[1], Pos[diff(x[[1]]$Chrom)!=0], Pos[length(Pos)])
    xx <- rep(start.Chrom, rep(2, length(start.Chrom)))
    yy <- rep(c(ylim, rev(ylim)), length(start.Chrom)/2)
    plot(y~Pos, pch=pch, col=col, ylim=ylim, ylab="Mean Log Ratio of the arrays",
         xlab="Pos.Base",
         main=main.text, type="n")
    polygon(xx, yy, col="gray85", xpd=TRUE)
    points(y~Pos, pch=pch, col=col)
    abline(h=0, lty=2)
    text(start.Chrom[-length(start.Chrom)] + diff(start.Chrom)/2,
         rep(ylim[1], length(start.Chrom)),
         labels=unique(x[[1]]$Chrom),
         cex=cex*0.7, pos=1)
    par(new=TRUE)
    plot(apply(prob.states, 1, max) ~ Pos, yaxt="n", xlab="", ylab="", col="blue", ylim=c(0,5), type="s")
    abline(h=c(0, 0.25, 0.5, 0.75,1), lty=2, yaxt="n")
    axis(4, prob.states, at=c(0, 0.5, 0.75, 1), labels=c(0, 0.5, 0.75, 1))
    mtext("Probability", 4, 2, cex=cex*0.7)
  }
  ## Name of this method?
 else {
   res <- lapply(x, model.averaging)
   if(class(x[[1]])=="RJaCGH.Chrom") {
     states <- list()
     k <- max(as.numeric(as.character(x[[1]]$Chrom)))
     for (i in 1:k) {
       states[[i]] <- rep(0, length(x[[1]][[i]]$y))
       for (j in 1:length(array.names)) {
         temp.states <- rep(0, length(res[[j]][[i]]$states))
         temp.states[res[[j]][[i]]$states=="Loss"] <- -1
         temp.states[res[[j]][[i]]$states=="Gain"] <- 1
         states[[i]] <- states[[i]] + temp.states * weights[j]
       }
     }
     states <- do.call("c",states)
     states <- states * 100
   }
    else {
      states <- rep(0, length(x[[1]]$y))
      for (j in 1:length(array.names)) {
        temp.states <- rep(0, length(res[[j]]$states))
        temp.states[res[[j]]$states=="Loss"] <- -1
        temp.states[res[[j]]$states=="Gain"] <- 1
        states <- states + temp.states * weights[j]
       }
      states <- states * 100
    }
   par(cex.lab=cex*0.8, cex.main=cex*0.9, cex.axis=cex*0.8)
   par(mar=c(5,4,4,4) + 0.1)
   main.text <- "Prediction of copy gain/loss. Bayesian Model Averaging"
   Pos <- x[[1]]$Pos

   col <- rep(1, length(states))
   col[states >= 50] <- 2
   col[states <= -50] <- 3
   pch <- rep(16, length(y))
   ylim <- c(-100, 100)
   margin <- 4
   ylim[1] <- ylim[1] - margin - 1
   ylim[2] <- ylim[2] + margin
   ## Start of every Chromosome
   start.Chrom <- c(Pos[1], Pos[diff(x[[1]]$Chrom)!=0], Pos[length(Pos)])
   xx <- rep(start.Chrom, rep(2, length(start.Chrom)))
   yy <- rep(c(ylim, rev(ylim)), length(start.Chrom)/2)
   plot(states~Pos, pch=pch, col=col, ylim=ylim, ylab="Percent of copy gain/loss in all arrays",
        xlab="Pos.Base",
        main=main.text, type="n")
   polygon(xx, yy, col="gray85", xpd=TRUE)
   points(states~Pos, pch=pch, col=col)
   abline(h=0, lty=2)
   text(start.Chrom[-length(start.Chrom)] + diff(start.Chrom)/2,
         rep(ylim[1], length(start.Chrom)),
        labels=unique(x[[1]]$Chrom), cex=cex*0.7, pos=1)
 }
}

trace.plot <- function(obj, k=NULL, array=NULL, Chrom=NULL, main.text=NULL) {
  ## Something must be made with colors
  if(class(obj)=="RJaCGH.array" && is.null(array)) {
    stop("Must specify array name\n")
  }

  if(class(obj)=="RJaCGH.Chrom" && is.null(Chrom)) {
    stop("Must specify Chromosome number\n")
  }
  if(class(obj)=="RJaCGH" || class(obj) == "RJaCGH.genome") {
    if (is.null(k)) k <- as.numeric(names(which.max(table(obj$k))))
    par(mfrow=c(2, 2))
    matplot(as.numeric(as.character(obj$k)), pch=16, cex=0.2, main="Trace plot of number of states",
         xlab="iteration", ylab="Number of states", type="l")
    matplot(obj[[k]]$mu, type="l", main="Trace plot of means", xlab="iteration",
           ylab="Mean of states")
    matplot(obj[[k]]$sigma.2, type="l", main="Trace plot of variance", xlab="iteration",
         ylab="Variance of the states")
    if (k >1) {
      matplot(t(obj[[k]]$beta[1,,]), type="n", main="Trace plot of beta", xlab="iteration",
              ylab="Beta")
      for (i in 1:k)
        matplot(t(obj[[k]]$beta[i,,]), type="l", add=TRUE)
    }
    if (is.null(main.text))
      main.text <- paste("Whole genome.", k, "hidden states")
    else
      main.text <- paste(main.text, k, "hidden states")
    mtext(main.text, outer=TRUE)
  }
  if(class(obj)=="RJaCGH.Chrom") {
    main.text <- paste("Chromosome", Chrom, ".")
    diagnostic.plot(obj[[Chrom]], main.text=main.text)
  }

  if(class(obj)=="RJaCGH.array") {
    main.text <- past("Array", array, ". Chromosome", Chrom, ".")
    diagnostic.plot(obj[[array]], chrom=chrom)
  }
}

##Brooks, S P. and Gelman, A. (1998) General Methods for Monitoring
##Convergence of Iterative Simulations. Journal of Computational and
##Graphical Statistics. 7. p434-455.


gelman.brooks.plot <- function(obj, bin=1000, array=NULL, Chrom=NULL, k=NULL) {
  obj.R <- list()
  if (class(obj[[1]])=="RJaCGH.array" && is.null(array)) {
    stop("Must specify array\n")
  }
  if (class(obj[[1]])=="RJaCGH.Chrom" && is.null(Chrom)) {
    stop("Must specify chromosome\n")
  }
  
  for (i in 1:length(obj)) {
    if (!is.null(array)) obj[[i]] <- obj[[i]][[array]]
    if (!is.null(Chrom)) obj[[i]] <- obj[[i]][[Chrom]]
  }
  par(mfrow=c(2,2), oma=c(1,1,3,1))

  C <- length(obj)
  ## Graph of number of states
  T <- min(unlist(lapply(obj, function(x) length(x$k))))
  t <- floor(T / bin)
  R <- rep(0, t+1)
  
  for (i in 1:(t+1)) {

    chain.mean <- rep(0, C)
    chain.var <- rep(0, C)
    for(j in 1:C) {
      if (i>t) batch <- obj[[j]]$k
      else batch <- obj[[j]]$k[1:(bin*i)]
      chain.mean[j] <- mean(as.numeric(as.character(batch)))
      chain.var[j] <- var(as.numeric(as.character(batch)))    
    }
    n <- length(batch)
    B.k <- sum(((chain.mean-mean(chain.mean))^2) * (n)/(C-1))
    W.k <- mean(chain.var)
    R[i] <- sqrt((((n-1) / n) * W.k + B.k/n) / W.k)
  }
  obj.R$k <- R[length(R)]
  plot(seq(bin, by=bin, length=t+1), R, type="l",
       ylim=c(0.8, 1.2), main="Parameter: k", xlab="Batch size")
  abline(h=1, lty=2)
  abline(h=c(0.9, 1.1), lty=3)

  ## Choose k, the max of the first chain (for example)
  if (is.null(k)) k <- as.numeric(names(which.max(table(obj[[1]]$k))))
  for (i in 1:C)   obj[[i]] <- obj[[i]][[k]]

  
  T <- min(unlist(lapply(obj, function(x) nrow(x$mu))))
  if (T==0)
    stop ("MCMC not converged.\n Gelman-Brooks values can't be computed and graph will not be drawn")
  ## no two chains ever visited the same state
  t <- floor(T / bin)

  ## Graph of mu  
  R <- matrix(0, t+1, k)
  for (i in 1:(t+1)) {
    chain.mean <- matrix(0, C, k)
    chain.var <- matrix(0, C, k)
    ## Number of observations in every chain
    TOT <- rep(0, C)
    for (j in 1:C) {
      if (i>t) batch <- obj[[j]]$mu
      else batch <- obj[[j]]$mu[1:(bin*i), ]
      batch <- matrix(batch, ncol=k)
      TOT[j] <- nrow(batch)
      chain.mean[j,] <- apply(batch, 2, mean)
      chain.var[j,] <- apply(batch, 2, var)
    }
    mean.chain.mean <- matrix(apply(chain.mean, 2, mean), C, k, byrow=TRUE)
    B <- ((chain.mean - mean.chain.mean)^2) * matrix(TOT/(C-1), C, k)
    B <- apply(B, 2, sum)
    W <- apply(chain.var, 2, mean)
    n <- mean(TOT)
    R[i,] <- sqrt((((n-1) / n) * W + B/n) / W)

  }
  matplot(seq(bin, by=bin, length=t+1), R, type="l",
          ylim=c(0.8, 1.2), main="Parameter: mean", xlab="Batch size")
  abline(h=c(0.9, 1.1), lty=3)
  abline(h=1, lty=2)
  obj.R$mu <- R[nrow(R),]
  ## Graph of sigma2
  R <- matrix(0, t+1, k)
  for (i in 1:(t+1)) {
    chain.mean <- matrix(0, C, k)
    chain.var <- matrix(0, C, k)
    ## Number of observations in every chain
    TOT <- rep(0, C)
    for (j in 1:C) {
      if (i>t) batch <- obj[[j]]$sigma.2
      else batch <- obj[[j]]$sigma.2[1:(bin*i), ]
      batch <- matrix(batch, ncol=k)
      TOT[j] <- nrow(batch)
      chain.mean[j,] <- apply(batch, 2, mean)
      chain.var[j,] <- apply(batch, 2, var)
    }
    mean.chain.mean <- matrix(apply(chain.mean, 2, mean), C, k, byrow=TRUE)
    B <- ((chain.mean - mean.chain.mean)^2) * matrix(TOT/(C-1), C, k)
    B <- apply(B, 2, sum)
    W <- apply(chain.var, 2, mean)
    n <- mean(TOT)
    R[i,] <- sqrt((((n-1) / n) * W + B/n) / W)

  }
  matplot(seq(bin, by=bin, length=t+1), R, type="l",
          ylim=c(0.8, 1.2),main="Parameter: sigma.2", xlab="Batch size")
  abline(h=1, lty=2)
  abline(h=c(0.9, 1.1), lty=3)
  obj.R$sigma.2 <- R[nrow(R),]
  ## Graph of beta
  if (k >1) {
    R <- matrix(0, t+1, k*(k-1))  
    for (i in 1:(t+1)) {
      chain.mean <- matrix(0, C, k*(k-1))
      chain.var <- matrix(0, C, k*(k-1))

      ## Number of observations in every chain
      TOT <- rep(0, C)
      for (j in 1:C) {
        batch <-NULL
        if (i>t) {
          for (r in 1:k) {
            for (c in 1:k) {
              if (r!=c) batch <- cbind(batch,log(obj[[j]]$beta[r,c,]))
            }
          }
        }
        else {
          for (r in 1:k) {
            for (c in 1:k) {
              if (r!=c) batch <- cbind(batch,log(obj[[j]]$beta[r, c, 1:(bin*i)]))
            }
          }
        }
        TOT[j] <- nrow(batch)
        chain.mean[j,] <- apply(batch, 2, mean)
        chain.var[j,] <- apply(batch, 2, var)
      }
      mean.chain.mean <- matrix(apply(chain.mean, 2, mean), C, k*(k-1), byrow=TRUE)
      B <- ((chain.mean - mean.chain.mean)^2) * matrix(TOT/(C-1), C, k*(k-1))
      B <- apply(B, 2, sum)
      W <- apply(chain.var, 2, mean)
      n <- mean(TOT)
      R[i,] <- sqrt((((n-1) / n) * W + B/n) / W)
    }
    matplot(seq(bin, by=bin, length=t+1), R, type="l",
            ylim=c(0.8, 1.2), main="Parameter: beta", xlab="Batch size")
    abline(h=1, lty=2)
    abline(h=c(0.9, 1.1), lty=3)
  }
  obj.R$beta <- R[nrow(R),]
  mtext("Gelman-Rubin diagnostic plots", outer=TRUE)
  obj.R
}

## This function mixes several parallel chains

collapseChain <- function(obj) {
  UseMethod("collapseChain", obj[[1]])
}

collapseChain.RJaCGH <- function(obj) {
  newobj <- list()
  class(newobj) <- "RJaCGH"
  newobj$y <- NULL
  newobj$x <- NULL
  newobj$Pos <- NULL
  newobj$prob.b <- NULL
  newobj$prob.d <- NULL
  newobj$prob.s <- NULL
  newobj$prob.c <- NULL
  newobj$k <- NULL
  C <- length(obj)
  k <- max(as.numeric(as.character(levels(obj[[1]]$k))))
  for (i in 1:k) {
    newobj[[i]] <- list()
    newobj[[i]]$stat <- NULL
    newobj[[i]]$mu <- NULL
    newobj[[i]]$sigma.2 <- NULL
    newobj[[i]]$beta <- NULL
    newobj[[i]]$loglik <- NULL
    newobj[[i]]$prob.states <- matrix(0, nrow=length(obj[[1]]$y),
  ncol=i)
    newobj[[i]]$prob.joint.states <- matrix(0,
  nrow=length(obj[[1]]$y)-1, ncol=i)
    newobj[[i]]$state.labels <- NULL
  }
  newobj$prob.b <- 0
  newobj$prob.d <- 0
  newobj$prob.s <- 0
  newobj$prob.c <- 0
  for (i in 1:C) {
    newobj$k <- c(newobj$k, obj[[i]]$k)
    newobj$prob.b <- newobj$prob.b + obj[[i]]$prob.b
    newobj$prob.d <- newobj$prob.d + obj[[i]]$prob.d
    newobj$prob.s <- newobj$prob.s + obj[[i]]$prob.s
    newobj$prob.c <- newobj$prob.c + obj[[i]]$prob.c
    newobj$y <- obj[[i]]$y
    newobj$x <- obj[[i]]$x
    newobj$Pos <- obj[[i]]$Pos
    for (j in 1:k) {
      newobj[[j]]$stat <- obj[[i]][[j]]$stat
      newobj[[j]]$mu <- rbind(newobj[[j]]$mu, obj[[i]][[j]]$mu)
      newobj[[j]]$sigma.2 <- rbind(newobj[[j]]$sigma.2, obj[[i]][[j]]$sigma.2)
      newobj[[j]]$beta <- c(newobj[[j]]$beta, obj[[i]][[j]]$beta)
      newobj[[j]]$loglik <- c(newobj[[j]]$loglik,
                              obj[[i]][[j]]$loglik)

        ##
      if (!is.null(obj[[i]][[j]]$prob.states)) {
        newobj[[j]]$prob.states <- newobj[[j]]$prob.states + 
          obj[[i]][[j]]$prob.states * nrow(obj[[i]][[j]]$mu)
        newobj[[j]]$prob.joint.states <- newobj[[j]]$prob.joint.states + 
          obj[[i]][[j]]$prob.joint.states * nrow(obj[[i]][[j]]$mu)
      }
      else {
        newobj[[j]]$prob.states <- newobj[[j]]$prob.states
        newobj[[j]]$prob.joint.states <- newobj[[j]]$prob.joint.states
      }
    }
  }
  for (j in 1:k) {
    newobj[[j]]$beta <- array(newobj[[j]]$beta,
                              c(j, j, length(newobj[[j]]$beta)/(j*j)))
    newobj[[j]]$prob.mu <- length(unique(newobj[[j]]$mu[,1])) /
      length(newobj[[j]]$mu[,1])
    newobj[[j]]$prob.sigma.2 <- length(unique(newobj[[j]]$sigma.2[,1])) /
      length(newobj[[j]]$sigma.2[,1])
    if (j >1) {
      newobj[[j]]$prob.beta <- length(unique(newobj[[j]]$beta[1,2,])) /
        length(newobj[[j]]$beta[1,2,])
    }
    else newobj[[j]]$prob.beta <- NULL
    if(nrow(newobj[[j]]$mu)) {
      newobj[[j]]$prob.states <- newobj[[j]]$prob.states /
        nrow(newobj[[j]]$mu)
      newobj[[j]]$prob.joint.states <- newobj[[j]]$prob.joint.states /
        nrow(newobj[[j]]$mu)
    }
    else {
      newobj[[j]]$prob.states <- NULL
      newobj[[j]]$prob.joint.states <- NULL
    }
  }
  newobj$k <- factor(newobj$k, levels=1:k)
  ## Recompute state labels
  for (i in 1:k) {
    if (table(newobj$k)[i] > 0) {
      p.labels <- 0.25 ## should be a parameter??
      obj.sum <- summary.RJaCGH(newobj, k=i, point.estimator="median")
      normal.levels <- (qnorm(mean=obj.sum$mu,
                              sd=sqrt(obj.sum$sigma.2),
                              p=p.labels) < 0 &
                        qnorm(mean=obj.sum$mu, sd=sqrt(obj.sum$sigma.2),
                              p=1-p.labels) > 0)
      n.Norm <- sum(normal.levels)
      if (n.Norm <=0) {
        normal.levels <- rep(FALSE, i)
        normal.levels[which.min(abs(obj.sum$mu))] <- TRUE
        n.Norm <- sum(normal.levels)
      }
      n.Loss <- which.max(normal.levels) -1
      n.Gain <- i - max(normal.levels*(1:i))
      newobj[[i]]$state.labels <- NULL
      if(n.Loss>0)
        newobj[[i]]$state.labels <- c(newobj[[i]]$state.labels,
                                   paste("Loss", n.Loss:1, sep="-"))
      newobj[[i]]$state.labels <- c(newobj[[i]]$state.labels, rep("Normal", n.Norm))
      if(n.Gain>0)
        newobj[[i]]$state.labels <- c(newobj[[i]]$state.labels,
                                   paste("Gain", 1:n.Gain, sep="-"))
      ## Should automatic labelling include the former steps?
      auto.label <- attr(obj[[1]], "auto.label")
      if (!is.null(auto.label)) {
        states <- states.RJaCGH(newobj, k=i)$states
        states <- factor(states, levels=names(summary.RJaCGH(newobj,
                                   k=i)$mu))
        freqs <- prop.table(table(states))
        labels <- names(freqs)
        means <- obj.sum$mu
        means.order <- order(abs(means))
        tot.norm <- freqs['Normal']
        ind.tot <- 1
        while(tot.norm < auto.label) {
          labels[means.order][ind.tot+1] <- "Normal"
          tot.norm <- tot.norm + freqs[means.order][ind.tot+1]
          ind.tot <- ind.tot + 1
        }
        newobj[[i]]$state.labels <- labels
      }
    }
  }
  newobj
}

collapseChain.RJaCGH.genome <- function(obj) {
  newobj <- collapseChain.RJaCGH(obj)
  newobj$model <- obj[[1]]$model
  newobj$Chrom <- obj[[1]]$Chrom
  class(newobj) <- class(obj[[1]])
  newobj
}

collapseChain.RJaCGH.Chrom <- function(obj) {
  C <- length(obj)
  newobj <- list()
  for (i in 1:23) {
    obj.temp <- list()
    for (j in 1:C) {
      obj.temp[[j]] <- obj[[j]][[i]]
    }
    newobj[[i]] <- collapseChain(obj.temp)
  }
  newobj$model <- obj[[1]]$model
  newobj$Chrom <- obj[[1]]$Chrom
  newobj$Pos <- obj[[1]]$Pos
  class(newobj) <- "RJaCGH.Chrom"
  newobj
}


collapseChain.RJaCGH.array <- function(obj) {
  C <- length(obj)
  newobj <- list()
  newobj$array.names <- obj[[1]]$array.names
  for (i in newobj$array.names) {
    newobj[[i]] <- list()
    obj.temp <- list()
    for (j in 1:C) {
      obj.temp[[j]] <- obj[[j]][[i]]
    }
    newobj[[i]] <- collapseChain(obj.temp)
  }
  
  class(newobj) <- "RJaCGH.array"
  newobj
}


chainsSelect <- function(obj, nutrim = 4, trim = NULL) {
  UseMethod("chainsSelect", obj[[1]])
}

chainsSelect.RJaCGH <- function(obj, nutrim = 4, trim = NULL) {
    ## This uses too much memory: obj is passed by
    ## copy, and we make yet another partial copy at the end
    n <- length(obj)

    if((is.null(nutrim) & is.null(trim)) |
       (!is.null(nutrim) & !is.null(trim)))
        stop("Exactly one of nutrim or trim must have non-NULL values")

    if (is.null(nutrim)) {
        if(trim > 0.5)
            stop("Trim cannot be larger than 0.5")
        lo <- floor(n * trim) + 1
        hi <- n + 1 - lo
    } else {
        remove <- n - nutrim
        lo <- floor(remove/2) + 1
        hi <- n - (remove - lo + 1)
##        print(paste("lo ", lo, " hi ", hi))
    }
    meank <- unlist(lapply(obj,function(x) mean(as.numeric(as.character(x$k)))))
    keepInd <- order(meank)[lo:hi]
    print(paste("keepInd ", paste(keepInd, collapse = " ")))
    newob <- list()
    j <- 1
    for(i in keepInd) {
        newob[[j]] <- obj[[i]]
        j <- j + 1
    }
    class(newob) <- class(obj)
    return(newob)
}

chainsSelect.RJaCGH.genome <- function(obj, nutrim = 4, trim = NULL) {
    newobj <- chainsSelect.RJaCGH(obj, nutrim = nutrim, trim = trim)
    class(newobj) <- class(obj)
    newobj
}


## chainsSelect.RJaCGH.Chrom <- function(obj, nutrim = 4, trim = NULL) {
##     n <- length(obj)
##     ch <- length(obj[[1]]) - 3 ## number of chromos.
    
##     if((is.null(nutrim) & is.null(trim)) |
##        (!is.null(nutrim) & !is.null(trim)))
##         stop("Exactly one of nutrim or trim must have non-NULL values")
##     if (is.null(nutrim)) {
##         if(trim > 0.5)
##             stop("Trim cannot be larger than 0.5")
##         lo <- floor(n * trim) + 1
##         hi <- n + 1 - lo
##     } else {
##         remove <- n - nutrim
##         lo <- floor(remove/2) + 1
##         hi <- n - (remove - lo + 1)
##     }

##     k.mat <- matrix(NA, nrow = n, ncol = ch)
##     for(i in 1:n) {
##         for(j in 1:ch) {
##             k.mat[i, j] <- mean(as.numeric(as.character(obj[[i]][[j]]$k)))
##         }
##     }
##     keepInd <- apply(k.mat, 2, function(x) order(x)[lo:hi])
##     new.obj <- list()
##     for(ii in 1:nrow(keepInd)) {
##         new.obj[[ii]] <- list()
##         for(jj in 1:ch) {
##             new.obj[[ii]][[jj]] <- obj[[k.mat[ii, jj]]][[jj]]
##             class(new.obj[[ii]][[jj]]) <- class(obj[[ii]][[jj]])
##         }
##         new.obj[[ii]][[ch + 1]] <- obj[[ii]][[ch + 1]]
##         new.obj[[ii]][[ch + 2]] <- obj[[ii]][[ch + 2]]
##         new.obj[[ii]][[ch + 3]] <- obj[[ii]][[ch + 3]]
##         class(new.obj[[ii]]) <- class(obj[[ii]])
##     }
##     class(new.obj) <- class(obj)
##     return(new.obj)
## }


chainsSelect.RJaCGH.array <- function(obj, nutrim = 4, trim = NULL) {
  C <- length(obj)
  newobj <- list()
  array.names <- obj[[1]]$array.names
  for (i in array.names) {
      newobj[[i]] <- list()
      obj.temp <- list()
      for (j in 1:C) {
          obj.temp[[j]] <- obj[[j]][[i]]
      }
      newobj[[i]] <- chainsSelect(obj.temp, nutrim = nutrim, trim = trim)
      class(newobj[[i]]) <- class(obj[[j]][[i]])
  }
  class(newobj) <- class(obj)
  newobj
}




##############################
## Adapt parameters intra model
get.jump <- function(y, x, Chrom, model, k.max=6,
                     prob.k=NULL, pb=NULL, ps=NULL, g=NULL, 
                     ka=NULL, mu.alfa=NULL, mu.beta=NULL, stat)
 {
  TOL <- 0.01

  optim.mu <- rep(0, k.max)
  optim.sigma.2 <- rep(0, k.max)
  optim.beta <- rep(0, k.max)
  lim <- sd(y)
  for (k in 2:k.max) {
    ##   cat("max", sigma.tau.mu.max, "\n")
    ##   cat("mu", sigma.tau.mu, "\n")
    ##   cat("min", sigma.tau.mu.min, "\n")
    p.mu <- 0
    p.sigma.2 <- 0
    p.beta <- 0
    
    prob.lim <- c(0.2, 0.35)

    sigma.tau.mu.max <- lim / k
    sigma.tau.mu.min <- 0.0001
    sigma.tau.sigma.2.max <- lim / k
    sigma.tau.sigma.2.min <- 0.0001
    sigma.tau.beta.max <- lim
    sigma.tau.beta.min <- 0.0001

    sigma.tau.mu <- (sigma.tau.mu.max - sigma.tau.mu.min) / 2
    sigma.tau.sigma.2 <- (sigma.tau.sigma.2.max - sigma.tau.sigma.2.min) / 2
    sigma.tau.beta <- (sigma.tau.beta.max - sigma.tau.beta.min) / 2
    
    tries <- 1
    
    ## Check if min == max
  while ((!p.mu | !p.sigma.2 | !p.beta) & tries < 5) {
    fit <- MetropolisSweep.C(y=y, x=x, k.max=k.max, Chrom=Chrom,
                           model=model, burnin=0, TOT=1000,
                           prob.k=prob.k, pb=pb, ps=ps, g=g,
                           ka=ka, mu.alfa=mu.alfa, mu.beta=mu.beta,
                           sigma.tau.mu=rep(sigma.tau.mu,k.max),
                           sigma.tau.sigma.2=rep(sigma.tau.sigma.2, k.max),
                           sigma.tau.beta=rep(sigma.tau.beta, k.max),
                           stat=stat, RJ=FALSE, start.k=k)
    prob.mu <- fit[[k]]$prob.mu
    prob.sigma.2 <- fit[[k]]$prob.sigma.2
    prob.beta <- fit[[k]]$prob.beta
    ## Check mu prob. of acceptance
    if (prob.mu < prob.lim[1] & !p.mu) {
      sigma.tau.mu.max <- sigma.tau.mu
      sigma.tau.mu <- sigma.tau.mu.min + (sigma.tau.mu.max - sigma.tau.mu.min)/2
    }
    if (prob.mu > prob.lim[2] | !p.mu) {
      sigma.tau.mu.min <- sigma.tau.mu
      sigma.tau.mu <- sigma.tau.mu.min + (sigma.tau.mu.max - sigma.tau.mu.min)/2
    }
    
    if (prob.mu > prob.lim[1] & prob.mu < prob.lim[2] & !p.mu) {
      p.mu <- 1
    }
    ## Check sigma.2 prob. of acceptance
    if (prob.sigma.2 < prob.lim[1] & !p.sigma.2) {
      sigma.tau.sigma.2.max <- sigma.tau.sigma.2
      sigma.tau.sigma.2 <- sigma.tau.sigma.2.min + (sigma.tau.sigma.2.max - sigma.tau.sigma.2.min)/2
    }
    if (prob.sigma.2 > prob.lim[2] | !p.sigma.2) {
      sigma.tau.sigma.2.min <- sigma.tau.sigma.2
      sigma.tau.sigma.2 <- sigma.tau.sigma.2.min + (sigma.tau.sigma.2.max - sigma.tau.sigma.2.min)/2
    }
    
  if (prob.sigma.2 > prob.lim[1] & prob.sigma.2 < prob.lim[2] & !p.sigma.2) {
    p.sigma.2 <- 1
  }
  ## Check beta prob. of acceptance
  if (prob.beta < prob.lim[1] & !p.beta) {
    sigma.tau.beta.max <- sigma.tau.beta
    sigma.tau.beta <- sigma.tau.beta.min + (sigma.tau.beta.max - sigma.tau.beta.min)/2
  }
  if (prob.beta > prob.lim[2] | !p.beta) {
    sigma.tau.beta.min <- sigma.tau.beta
    sigma.tau.beta <- sigma.tau.beta.min + (sigma.tau.beta.max - sigma.tau.beta.min)/2
  }

  if (prob.beta > prob.lim[1] & prob.beta < prob.lim[2] & !p.beta) {
    p.beta <- 1
  }
##   cat("Iteration=", i, "\n")
##   cat("max", sigma.tau.mu.max, "\n")
##   cat("mu", sigma.tau.mu, "\n")
##   cat("min", sigma.tau.mu.min, "\n")
    tries <- tries + 1
  }
  optim.mu[k] <- sigma.tau.mu
  optim.sigma.2[k] <- sigma.tau.sigma.2
  optim.beta[k] <- sigma.tau.beta
}
  optim.mu[1] <- optim.mu[2]
  optim.sigma.2[1] <- optim.sigma.2[2]
  optim.beta[1] <- optim.beta[2]
  res <- list(sigma.tau.mu=abs(optim.mu), sigma.tau.sigma.2=abs(optim.sigma.2),
              sigma.tau.beta=abs(optim.beta))
  print(res)
  res
}


## Example for viterbi

viterbi.C <- function(y, x=NULL, Chrom=NULL, mu, sigma.2, beta, stat=NULL) {

  if (!is.null(Chrom)) {
    index <- diff(Chrom)
    index <- which(index>0)
    index <- c(0, index, length(y))
    genome <- length(index) - 1
  }
  else {
    genome <- 1
    index <- c(0, length(y))
  }
  if (is.null(x)) x <- rep(0, length(y)-1)
  k <- length(mu)
  n <- length(y)
  if (is.null(stat)) stat <- rep(1/k, k)
  states <- .C("viterbi", y=as.double(y), x=as.double(x),
               genome=as.integer(genome),
               index = as.integer(index), k =as.integer(k),
               n=as.integer(n), mu=as.double(mu),
               sigma2=as.double(sigma.2), beta=as.double(beta),
               stat=as.double(stat), states=as.integer(rep(0, n)))
  states <- states$states
  states
}

akaike <- function(logliks, param=NULL) {
  if (is.null(param)) {
    param <- 1:length(logliks)
    param <- 2 * param + param * (param-1)
  }
  AIC <- -2*logliks + 2*param
  AIC <- AIC - min(AIC)
  akaike <- exp(-0.5*AIC)
  akaike <- akaike / sum(akaike)
  akaike
}





## Examples for chainsSelect







###### A test of chainsSelect

## one array

## y <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
##        rnorm(100,0, 1)) 
## Pos <- runif(230)
## Pos <- cumsum(Pos)
## Chrom <- rep(1:23, rep(10, 23))

## jp <- list(sigma.tau.mu=rep(0.5, 4), sigma.tau.sigma.2=rep(0.3, 4),
##            sigma.tau.beta=rep(0.7, 4), tau.split.mu=0.5, tau.split.beta=0.5)

## fit.genome <- list()
## for (i in 1:8) {
##     fit.genome[[i]] <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="genome",
##                               burnin=10, TOT=1000, jump.parameters=jp, k.max = 4)
## }


## fit.chrom <- list()
## for (i in 1:8) {
##     fit.chrom[[i]] <- RJaCGH(y=y, Pos=Pos, Chrom=Chrom, model="Chrom",
##                               burnin=10, TOT=1000, jump.parameters=jp, k.max = 4)
## }





## ## many arrays

## y1 <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
##        rnorm(100,0, 1)) 
## y2 <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
##        rnorm(100,0, 1)) 
## y3 <- c(rnorm(100, 0, 1), rnorm(10, -3, 1), rnorm(20, 3, 1),
##         rnorm(100,0, 1)) 
## Y <- cbind(y1, y2, y3)


## fit.genome.arrays <- list()
## for (i in 1:8) {
##     fit.genome.arrays[[i]] <- RJaCGH(y=Y, Pos=Pos, Chrom=Chrom, model="genome",
##                               burnin=10, TOT=1000, jump.parameters=jp, k.max = 4)
## }


## fit.chrom.arrays <- list()
## for (i in 1:8) {
##     fit.chrom.arrays[[i]] <- RJaCGH(y=Y, Pos=Pos, Chrom=Chrom, model="Chrom",
##                               burnin=10, TOT=1000, jump.parameters=jp, k.max = 4)
## }


## #####

## o1 <- chainsSelect(fit.chrom.arrays)
## o2 <- chainsSelect(fit.genome.arrays)
## o3 <- chainsSelect(fit.chrom)
## o4 <- chainsSelect(fit.genome)

## MCR <- function(obj, p, method, model.averaging=TRUE, alteration="Gain",
##                 array.weights=NULL, array.freq=0.5) {
##   UseMethod("MCR")
## }

## MCR.RJaCGH <- function(obj, p, method, model.averaging=TRUE, alteration="Gain",
##                        array.weights=NULL) {
##   if (alteration !="Gain" && alteration!="Loss")
##     stop ("'alteration' must be either 'Gain' or 'Loss'")

##   n <- length(obj$y)
##   regions <- list()
##   counter <- 1
##   i <- 1
##   ## First, we relabel the joint probabilities according to the states
##   ## First, only for gains
##   avg <- model.averaging(obj)
##   marginal.probs <- avg$prob.states[, alteration]
##   joint.probs <- avg$prob.joint.states[, alteration]
##   while (i <=n) {
##     if (marginal.probs[i] >= p) {
##       regions[[counter]] <- list()
##       regions[[counter]]$start <- obj$Pos[i]
##       regions[[counter]]$end <- obj$Pos[i]
##       regions[[counter]]$prob <- marginal.probs[i]
##       while (i < n && regions[[counter]]$prob * joint.probs[i] >=p) {
##         regions[[counter]]$end <- obj$Pos[i+1]
##         regions[[counter]]$prob <- regions[[counter]]$prob * joint.probs[i]
##         i <- i + 1
##       }
##       counter <- counter + 1
##     }
##     i <- i +1
##   }
##   attr(regions, "alteration") <- alteration
##   class(regions) <- "MCR.RJaCGH"
##   regions
## }

## MCR.RJaCGH.Chrom <- function(obj, p, method, model.averaging=TRUE, alteration="Gain",
##                        array.weights=NULL) {
##   if (alteration !="Gain" && alteration!="Loss") 
##     stop ("'alteration' must be either 'Gain' or 'Loss'")

##   avg <- model.averaging(obj)
##   regions <- list()
##   for (chr in unique(obj$Chrom)) {
##     marginal.probs <- avg[[chr]]$prob.states[, alteration]
##     joint.probs <- avg[[chr]]$prob.joint.states[, alteration]
##     n <- length(obj[[chr]]$y)
##     regions[[chr]] <- list()
##     attr(regions[[chr]], "Chrom") <- chr
##     counter <- 1
##     i <- 1
##     ## First, we relabel the joint probabilities according to the states
##     ## First, only for gains
##     while (i <=n) {
##       if (marginal.probs[i] >= p) {
##         regions[[chr]][[counter]] <- list()
##         regions[[chr]][[counter]]$start <- obj$Pos[chr][i]
##         regions[[chr]][[counter]]$end <- obj$Pos[chr][i]
##         regions[[chr]][[counter]]$prob <-
##           marginal.probs[i]
##         while (i <n &&
##                regions[[chr]][[counter]]$prob * joint.probs[i] >=p) {
##           regions[[chr]][[counter]]$end <- obj$Pos[chr][i+1]
##           regions[[chr]][[counter]]$prob <-
##             regions[[chr]][[counter]]$prob * joint.probs[i]
##           i <- i + 1
##         }
##       counter <- counter + 1
##       }
##       i <- i +1
##     }
##   }
##   attr(regions, "alteration") <- alteration
##   class(regions) <- "MCR.RJaCGH.Chrom"
##   regions
## }


## MCR.RJaCGH.genome <- function(obj, p, method, model.averaging=TRUE, alteration="Gain",
##                        array.weights=NULL) {
##   if (alteration !="Gain" && alteration!="Loss") 
##     stop ("'alteration' must be either 'Gain' or 'Loss'")

##   avg <- model.averaging(obj)
##   marginal.probs <- avg$prob.states[, alteration]
##   joint.probs <- avg$prob.joint.states[, alteration]
  
##   regions <- list()
##   for (chr in unique(obj$Chrom)) {
##     n <- length(obj$y[obj$Chrom == chr])
##     regions[[chr]] <- list()
##     attr(regions[[chr]], "Chrom") <- chr
##     counter <- 1
##     i <- 1
##     ## First, we relabel the joint probabilities according to the states
##     ## First, only for gains
##     while (i <=n) {
##       if (marginal.probs[obj$Chrom==chr][i] >= p) {
##         regions[[chr]][[counter]] <- list()
##         regions[[chr]][[counter]]$start <- obj$Pos[obj$Chrom==chr][i]
##         regions[[chr]][[counter]]$end <- obj$Pos[obj$Chrom==chr][i]
##         regions[[chr]][[counter]]$prob <-
##           marginal.probs[obj$Chrom==chr][i]
##         while (i <n &&
##                regions[[chr]][[counter]]$prob * joint.probs[obj$Chrom==chr][i] >=p) {
##           regions[[chr]][[counter]]$end <- obj$Pos[obj$Chrom==chr][i+1]
##           regions[[chr]][[counter]]$prob <-
##             regions[[chr]][[counter]]$prob * joint.probs[obj$Chrom==chr][i]
##           i <- i + 1
##         }
##       counter <- counter + 1
##       }
##       i <- i +1
##     }
##   }
##   attr(regions, "alteration") <- alteration
##   class(regions) <- "MCR.RJaCGH.genome"
##   regions
## }


## ##What about joint prob of all regions (chrom/genome)?


## MCR.RJaCGH.array <- function(obj, p, method="global", model.averaging=TRUE,
##                              alteration="Gain", array.weights=NULL, array.freq=0.5) {
##   if (method=="global") {
##     regions <- MCR.RJaCGH.array.global(obj, p, model.averaging,
##                              alteration, array.weights)
##     attr(regions, "alteration") <- alteration
##   }

##   else if (method=="individual") {
##     regions <- MCR.RJaCGH.array.individual(obj, p, model.averaging,
##                              alteration, array.weights, array.freq)
##     attr(regions, "alteration") <- alteration
##   }
##   else stop("'method' must be 'global' or 'individual'\n")
##   regions
## }

## MCR.RJaCGH.array.global <- function(obj, p, model.averaging=TRUE,
##                              alteration="Gain",
##                              array.weights=NULL) {

##   if (alteration !="Gain" && alteration!="Loss")
##     stop ("'alteration' must be either 'Gain' or 'Loss'")
##   if (is.null(array.weights)) array.weights <- rep(1/length(obj$array.names),
##                                                      length(obj$array.names))
##   ar.name <- obj$array.names[[1]]
##   ## ##########################################
##   ## Model == genome
##   if (class(obj[[ar.name]]) == "RJaCGH.genome") {
##     Chromosomes <- obj[[ar.name]]$Chrom
##     avg <- model.averaging(obj)
##     marginal.probs <- lapply(avg, function(x, alteration)
##                              x$prob.states[, alteration], alteration)
##     joint.probs <- lapply(avg, function(x, alteration)
##                           x$prob.joint.states[, alteration], alteration)
##     regions <- list()
##     for (chr in unique(obj[[ar.name]]$Chrom)) {
##       ## We take the first array for reference (length, Position...)
##       n <- length(obj[[ar.name]]$y[Chromosomes==chr])
##       regions[[chr]] <- list()
##       attr(regions[[chr]], "Chrom") <- chr
##       counter <- 1
##       i <- 1
##       ## First, we relabel the joint probabilities according to the states
##       ## First, only for gains
##       while (i <=n) {
##         prob.temp <- sapply(marginal.probs, function(x)
##                             x[Chromosomes==chr][i])
##         prob.temp <- sum(prob.temp * array.weights)
##         if (prob.temp >= p) {
##           index.start <- i
##           regions[[chr]][[counter]] <- list()
##           regions[[chr]][[counter]]$start <-
##             obj[[ar.name]]$Pos[Chromosomes==chr][i]
##           regions[[chr]][[counter]]$end <-
##             obj[[ar.name]]$Pos[Chromosomes== chr][i]
##           regions[[chr]][[counter]]$prob <- prob.temp
##           prob.temp <- mapply(function(x,y,i) {
##             x[Chromosomes==chr][i] * y[Chromosomes==chr][i]
##           }, x=marginal.probs, y=joint.probs, i=i)
##           while (i < n && sum(prob.temp * array.weights) >=p) {
##             regions[[chr]][[counter]]$end <- obj[[ar.name]]$Pos[Chromosomes==chr][i+1]
##             regions[[chr]][[counter]]$prob <- sum(prob.temp * array.weights)
##             i <- i + 1
##             prob.temp <- mapply(function(x,y,i) {
##               x[Chromosomes==chr][index.start] * prod(y[Chromosomes==chr][index.start:i])
##             }, x=marginal.probs, y=joint.probs, i=i)
##           }
##           counter <- counter + 1
##         }
##         i <- i +1
##       }
##     }
##     class(regions) <- "MCR.RJaCGH.array.global"
##   }
##   ## ##########################################
##   ## Model == Chromosome
##   else if (class(obj[[ar.name]]) == "RJaCGH.Chrom") {
##     regions <- list()
##     avg <- model.averaging(obj)
##     for (chr in unique(obj[[ar.name]]$Chrom)) {
##       ## We take the first array for reference (length, Position...)
##       n <- length(obj[[ar.name]][[chr]]$y)
##       regions[[chr]] <- list()
##       attr(regions[[chr]], "Chrom") <- chr
##       counter <- 1
##       i <- 1
##       ## First, we relabel the joint probabilities according to the states
##       marginal.probs <- lapply(avg, function(x, alteration)
##                                x[[chr]]$prob.states[, alteration], alteration)
##       joint.probs <- lapply(avg, function(x, alteration)
##                             x[[chr]]$prob.joint.states[, alteration], alteration)
##       while (i <=n) {
##         prob.temp <- sapply(marginal.probs, function(x) x[i])
##         prob.temp <- sum(prob.temp * array.weights)
##         if (prob.temp >= p) {
##           index.start <- i
##           regions[[chr]][[counter]] <- list()
##           regions[[chr]][[counter]]$start <- obj[[ar.name]][[chr]]$Pos[i]
##           regions[[chr]][[counter]]$end <- obj[[ar.name]][[chr]]$Pos[i]
##           regions[[chr]][[counter]]$prob <- prob.temp
##           prob.temp <- mapply(function(x,y,i) {
##             x[i] * y[i]
##           }, x=marginal.probs, y=joint.probs, i=i)
##           while (i < n && sum(prob.temp * array.weights) >=p) {
##             regions[[chr]][[counter]]$end <- obj[[ar.name]][[chr]]$Pos[i+1]
##             regions[[chr]][[counter]]$prob <- sum(prob.temp * array.weights)
##             i <- i + 1
##             prob.temp <- mapply(function(x,y,i) {
##               x[index.start] * prod(y[index.start:i])
##             }, x=marginal.probs, y=joint.probs, i=i)
##           }
##           counter <- counter + 1
##         }
##         i <- i +1
##       }
##     }
##     class(regions) <- "MCR.RJaCGH.array.global"
##   }

##   ## ##########################################
##   ## Model == None (no Chromosome info)
##   else if (class(obj[[ar.name]]) == "RJaCGH") {
##     regions <- list()
##     avg <- model.averaging(obj)
##       ## We take the first array for reference (length, Position...)
##       n <- length(obj[[ar.name]]$y)
##       counter <- 1
##       i <- 1
##       ## First, we relabel the joint probabilities according to the states
##       marginal.probs <- lapply(avg, function(x, alteration)
##                                x$prob.states[, alteration], alteration)
##       joint.probs <- lapply(avg, function(x, alteration)
##                             x$prob.joint.states[, alteration], alteration)
##       while (i <=n) {
##         prob.temp <- sapply(marginal.probs, function(x) x[i])
##         prob.temp <- sum(prob.temp * array.weights)
##         if (prob.temp >= p) {
##           index.start <- i
##           regions[[counter]] <- list()
##           regions[[counter]]$start <- obj[[ar.name]]$Pos[i]
##           regions[[counter]]$end <- obj[[ar.name]]$Pos[i]
##           regions[[counter]]$prob <- prob.temp
##           prob.temp <- mapply(function(x,y,i) {
##             x[i] * y[i]
##           }, x=marginal.probs, y=joint.probs, i=i)
##           while (i < n && sum(prob.temp * array.weights) >=p) {
##             regions[[counter]]$end <- obj[[ar.name]]$Pos[i+1]
##             regions[[counter]]$prob <- sum(prob.temp * array.weights)
##             i <- i + 1
##             prob.temp <- mapply(function(x,y,i) {
##               x[index.start] * prod(y[index.start:i])
##             }, x=marginal.probs, y=joint.probs, i=i)
##           }
##           counter <- counter + 1
##         }
##         i <- i +1
##       }
##     class(regions) <- "MCR.RJaCGH.array.noChrom.global"
##   }
##   regions
## }

## MCR.RJaCGH.array.individual <- function(obj, p, model.averaging=TRUE,
##                              alteration="Gain",
##                              array.weights=NULL, array.freq=0.5) {
##   if (alteration !="Gain" && alteration!="Loss")
##     stop ("'alteration' must be either 'Gain' or 'Loss'")
##   if (is.null(array.weights)) array.weights <- rep(1/length(obj$array.names),
##                                                      length(obj$array.names))
##   ar.name <- obj$array.names[[1]]
##   ## ##########################################
##   ## Model == genome
##   if (class(obj[[ar.name]]) == "RJaCGH.genome") {
##     Chromosomes <- obj[[ar.name]]$Chrom
##     avg <- model.averaging(obj)
##     marginal.probs <- lapply(avg, function(x, alteration)
##                              x$prob.states[, alteration], alteration)
##     joint.probs <- lapply(avg, function(x, alteration)
##                           x$prob.joint.states[, alteration], alteration)
##     regions <- list()
##     for (chr in unique(obj[[ar.name]]$Chrom)) {
##       ## We take the first array for reference (length, Position...)
##       n <- length(obj[[ar.name]]$y[Chromosomes==chr])
##       regions[[chr]] <- list()
##       attr(regions[[chr]], "Chrom") <- chr
##       counter <- 1
##       i <- 1
##       ## First, we relabel the joint probabilities according to the states
##       ## Check that it's similar to Rouveirol et al.
##       while (i <=n) {
##         prob.temp <- sapply(marginal.probs, function(x)
##                             x[Chromosomes==chr][i])
##         if (mean(prob.temp >= p) >= array.freq) {
##           index.start <- i
##           regions[[chr]][[counter]] <- list()
##           regions[[chr]][[counter]]$start <-
##             obj[[ar.name]]$Pos[Chromosomes==chr][i]
##           regions[[chr]][[counter]]$end <-
##             obj[[ar.name]]$Pos[Chromosomes== chr][i]
##           arrays.selected <- prob.temp[prob.temp >=p]
##           regions[[chr]][[counter]]$prob <-  data.frame(
##                                                         arrays=names(arrays.selected),
##                                                         prob=arrays.selected)
          
##           prob.temp <- mapply(function(x,y,i) {
##             x[Chromosomes==chr][i] * y[Chromosomes==chr][i]
##           }, x=marginal.probs, y=joint.probs, i=i)
##           names(prob.temp) <- obj$array.names
##           ## because we've lost 'em
##           while (i < n && mean(prob.temp >=p)>= array.freq &&
##                 isTRUE(all.equal(names(prob.temp[prob.temp >p]),
##                 names(arrays.selected)))) {
##             ## Check that we keep all the arrays
##               regions[[chr]][[counter]]$end <- obj[[ar.name]]$Pos[Chromosomes==chr][i+1]
##               regions[[chr]][[counter]]$prob <- data.frame(
##                                                            arrays=names(arrays.selected),
##                                                            prob=arrays.selected)
##               i <- i + 1
##               prob.temp <- mapply(function(x,y,i) {
##                 x[Chromosomes==chr][index.start] * prod(y[Chromosomes==chr][index.start:i])
##               }, x=marginal.probs, y=joint.probs, i=i)
##               names(prob.temp) <- obj$array.names
##             }
##           counter <- counter + 1
##         }
##         i <- i +1
##       }
##     }
##     class(regions) <- "MCR.RJaCGH.array.individual"
##   }
##   ## ##########################################
##   ## Model == Chromosome
##   else if (class(obj[[ar.name]]) == "RJaCGH.Chrom") {
##     regions <- list()
##     avg <- model.averaging(obj)
##     for (chr in unique(obj[[ar.name]]$Chrom)) {
##       ## We take the first array for reference (length, Position...)
##       n <- length(obj[[ar.name]][[chr]]$y)
##       regions[[chr]] <- list()
##       attr(regions[[chr]], "Chrom") <- chr
##       counter <- 1
##       i <- 1
##       ## First, we relabel the joint probabilities according to the states
##       marginal.probs <- lapply(avg, function(x, alteration)
##                                x[[chr]]$prob.states[, alteration], alteration)
##       joint.probs <- lapply(avg, function(x, alteration)
##                             x[[chr]]$prob.joint.states[, alteration], alteration)
##       while (i <=n) {
##         prob.temp <- sapply(marginal.probs, function(x) x[i])
##         if (mean(prob.temp >= p)>=array.freq) {
##           index.start <- i
##           regions[[chr]][[counter]] <- list()
##           regions[[chr]][[counter]]$start <- obj[[ar.name]][[chr]]$Pos[i]
##           regions[[chr]][[counter]]$end <-
##             obj[[ar.name]][[chr]]$Pos[i]
##           arrays.selected <- prob.temp[prob.temp >=p]
##           regions[[chr]][[counter]]$prob <-
##             data.frame(arrays=names(arrays.selected), prob=arrays.selected)
##           prob.temp <- mapply(function(x,y,i) {
##             x[i] * y[i]
##           }, x=marginal.probs, y=joint.probs, i=i)
##           names(prob.temp) <- obj$array.names
##           while (i < n && mean(prob.temp >=p)>=array.freq &&
##                 isTRUE(all.equal(names(prob.temp[prob.temp>p]),
##                 names(arrays.selected)))) {
##             regions[[chr]][[counter]]$end <- obj[[ar.name]][[chr]]$Pos[i+1]
##             regions[[chr]][[counter]]$prob <-
##               data.frame(arrays=names(arrays.selected), prob=arrays.selected)
##             i <- i + 1
##             prob.temp <- mapply(function(x,y,i) {
##               x[index.start] * prod(y[index.start:i])
##             }, x=marginal.probs, y=joint.probs, i=i)
##             names(prob.temp) <- obj$array.names
##           }
##           counter <- counter + 1
##         }
##         i <- i +1
##       }
##     }
##     class(regions) <- "MCR.RJaCGH.array.individual"
##   }
##   ## ##########################################
##   ## what about selecting minimum array frequence
##   ## Model == none (No Chrom info)
##   else if (class(obj[[ar.name]]) == "RJaCGH") {
##     regions <- list()
##     avg <- model.averaging(obj)
##     ## We take the first array for reference (length, Position...)
##     n <- length(obj[[ar.name]]$y)
##     counter <- 1
##     i <- 1
##     ## First, we relabel the joint probabilities according to the states
##     marginal.probs <- lapply(avg, function(x, alteration)
##                              x$prob.states[, alteration], alteration)
##     joint.probs <- lapply(avg, function(x, alteration)
##                           x$prob.joint.states[, alteration], alteration)
##     while (i <=n) {
##       prob.temp <- sapply(marginal.probs, function(x) x[i])
##       if (mean(prob.temp >= p)>=array.freq) {
##         index.start <- i
##         regions[[counter]] <- list()
##         regions[[counter]]$start <- obj[[ar.name]]$Pos[i]
##         regions[[counter]]$end <-
##           obj[[ar.name]]$Pos[i]
##         arrays.selected <- prob.temp[prob.temp >=p]
##         regions[[counter]]$prob <-
##           data.frame(arrays=names(arrays.selected), prob=arrays.selected)
##         prob.temp <- mapply(function(x,y,i) {
##           x[i] * y[i]
##         }, x=marginal.probs, y=joint.probs, i=i)
##         names(prob.temp) <- obj$array.names
##         while (i < n && mean(prob.temp >=p)>=array.freq &&
##                isTRUE(all.equal(names(prob.temp[prob.temp>p]),
##                                 names(arrays.selected)))) {
##           regions[[counter]]$end <- obj[[ar.name]]$Pos[i+1]
##           regions[[counter]]$prob <-
##             data.frame(arrays=names(arrays.selected), prob=arrays.selected)
##           i <- i + 1
##           prob.temp <- mapply(function(x,y,i) {
##             x[index.start] * prod(y[index.start:i])
##           }, x=marginal.probs, y=joint.probs, i=i)
##           names(prob.temp) <- obj$array.names
##         }
##         counter <- counter + 1
##       }
##       i <- i +1
##     }
##     class(regions) <- "MCR.RJaCGH.array.noChrom.individual"
##   }
##   regions
## }


## ##What about joint prob of all regions (chrom/genome)?
## print.MCR <- function(obj) {
##   UseMethod("print.MCR")
## }

## print.MCR.RJaCGH <- function(obj) {

##   cat("Start\tEnd\tProb.", attr(obj, "alteration"), "\n")
##   res <- sapply(obj, function(x) c(x$start, x$end, x$prob)
##          )
##   res <- t(res)
##   colnames(res) <-   c("Start", "End",
##                        paste("Prob.", attr(obj, "alteration")))
##   print(res)
  
## }

## print.MCR.RJaCGH.Chrom <- function(obj) {

##   res <- sapply(obj, function(x) {
##     sapply(x, function(y) c(attr(x, "Chrom"), y$start, y$end, y$prob))
##   })
##   ## res <- do.call("cbind", res)
##   ## this does't work well
##   ## some elements are matrices and others (the empty ones) lists
##   res <- do.call("cbind", res[!sapply(res, is.list)])
##   res <- t(res)
##   colnames(res) <-   c("Chromosome", "Start", "End",
##                        paste("Prob.", attr(obj, "alteration")))
##   print(res)

## }
       


## print.MCR.RJaCGH.genome <- function(obj) {

##   res <- sapply(obj, function(x) {
##     sapply(x, function(y) c(attr(x, "Chrom"), y$start, y$end, y$prob))
##   })
##   ## res <- do.call("cbind", res)
##   ## this does't work well
##   ## some elements are matrices and others (the empty ones) lists
##   res <- do.call("cbind", res[!sapply(res, is.list)])
##   res <- t(res)
##   colnames(res) <-   c("Chromosome", "Start", "End",
##                        paste("Prob.", attr(obj, "alteration")))
##   print(res)

## }

## print.MCR.RJaCGH.array.global <- function(obj) {


##   res <- sapply(obj, function(x) {
##     sapply(x, function(y) c(attr(x, "Chrom"), y$start, y$end, y$prob))
##   })
##   ## res <- do.call("cbind", res)
##   ## this does't work well
##   ## some elements are matrices and others (the empty ones) lists
##   res <- do.call("cbind", res[!sapply(res, is.list)])
##   res <- t(res)
##   colnames(res) <-   c("Chromosome", "Start", "End",
##                        paste("Prob.", attr(obj, "alteration")))
##   print(res)
## }

## print.MCR.RJaCGH.array.individual <- function(obj) {


##   res <- sapply(obj, function(x) {
##     sapply(x, function(y) c(attr(x, "Chrom"), y$start, y$end, nrow(y$prob)))
##   })
##   ## res <- do.call("cbind", res)
##   ## this does't work well
##   ## some elements are matrices and others (the empty ones) lists
##   res <- do.call("cbind", res[!sapply(res, is.list)])
##   res <- t(res)
##   colnames(res) <-   c("Chromosome", "Start", "End",
##                        paste("#Samples with", attr(obj, "alteration")))
##   print(res)
## }

## print.MCR.RJaCGH.array.noChrom.global <- function(obj) {

##   res <- sapply(obj, function(y) c(y$start, y$end, y$prob))
##   ## res <- do.call("cbind", res)
##   ## this does't work well
##   ## some elements are matrices and others (the empty ones) lists
##   res <- t(res)
##   colnames(res) <-   c("Start", "End",
##                        paste("Prob.", attr(obj, "alteration")))
##   print(res)
## }



## print.MCR.RJaCGH.array.noChrom.individual <- function(obj) {

##   res <- sapply(obj, function(y) c(y$start, y$end, nrow(y$prob)))
##   ## res <- do.call("cbind", res)
##   ## this does't work well
##   ## some elements are matrices and others (the empty ones) lists
##   res <- t(res)
##   colnames(res) <-   c("Start", "End",
##                        paste("#Samples with", attr(obj, "alteration")))
##   print(res)
## }

## ## method='individual' does not exactly what's intented
## ## what's intended?
