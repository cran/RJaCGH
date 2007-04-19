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
  freqJointStates <- rep(0, (n-1)*((k.max*(k.max+1)*(2*k.max+1)/6)-1))
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
              freqJointStates=as.double(freqJointStates),
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
              freqJointStates=as.double(freqJointStates),
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
        obj[[i]]$freq.joint.states <- list()
        for (k in 1:i) {
          obj[[i]]$freq.joint.states[[k]] <-
            res$freqJointStates[indexJointStates: (indexJointStates +
                                               (n-1)*i -1)]
          obj[[i]]$freq.joint.states[[k]] <-
            matrix(obj[[i]]$freq.joint.states[[k]], nrow=n-1, ncol=i)
          indexJointStates <- indexJointStates + (n-1)*i
        }
      }
      else {
        obj[[i]]$prob.states <- NULL
        obj[[i]]$freq.joint.states <- NULL
        indexJointStates <- indexJointStates + (n-1)*i^2
      }
    }
    else  {
      obj[[i]]$prob.states <- matrix(rep(1, length(y)), ncol=1)
      obj[[i]]$freq.joint.states <- list()
      ## it doesn't matter to compute joint probs.
      ## so we take all ones
      obj[[i]]$freq.joint.states[[1]] <- matrix(rep(1, length(y)-1), ncol=1)
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
  p.labels <- 0.95 ## should be a parameter??
  for (i in 1:k.max) {
    if (table(obj$k)[i] > 0) {
      obj.sum <- summary.RJaCGH(obj, k=i, point.estimator="median")
      normal.levels <- (qnorm(mean=obj.sum$mu,
                              sd=sqrt(obj.sum$sigma.2),
                              p=(1-p.labels)/2) < 0 &
                        qnorm(mean=obj.sum$mu, sd=sqrt(obj.sum$sigma.2),
                              p=1-(1-p.labels)/2) > 0)
      ## bug fix to prevent a case
      ## in which normal.levels are
      ## ...TRUE FALSE TRUE...
      if (sum(normal.levels) >0) {
        normal.levels[!normal.levels][obj.sum$mu[!normal.levels] > min(obj.sum$mu[normal.levels]) &
                                      obj.sum$mu[!normal.levels] < max(obj.sum$mu[normal.levels])] <- 
                                        TRUE
      }
      ##
      ##
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
        means <- summary.RJaCGH(obj, k=i)$mu
        means.1 <- abs(means - max(means[names(means)=='Normal']))
        means.2 <- abs(means - min(means[names(means)=='Normal']))
        means <- pmin(means.1, means.2)
        means.order <- order(means)
        names(means.order) <- names(means)[order(means)]
        means.order <- means.order[names(means.order) != 'Normal']
        tot.norm <- sum(freqs[names(freqs)=='Normal'])
        ind.tot <- 1
        while(tot.norm <= auto.label) {
          labels[means.order][ind.tot] <- "Normal"
          tot.norm <- tot.norm + freqs[means.order][ind.tot]
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
  ## if all distances are the same, they should be zero
  if (is.null(x) | isTRUE(all.equal(min(x.old, na.rm=TRUE), max(x.old, na.rm=TRUE)))) {
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
  Pos.rel <- NULL
  if(!is.null(Pos)) {
    if (!is.null(Chrom) && any(diff(Pos)<0)) {
      ## save relative positions for pMCR, plots, etc.
      Pos.rel <- Pos
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
    res$Pos.rel <- Pos.rel
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
          chrom.Dist <- rep(0, length(chrom.Pos) - 1)
        }
        else {
          chrom.Pos <- Pos[Chrom==i]
          chrom.Dist <- diff(chrom.Pos) - 1
          chrom.Dist[chrom.Dist<0] <- 0
        }
        if (!is.null(Dist)) {
          chrom.Dist <- Dist[Chrom==i]
          if (length(chrom.Dist) == length(chrom.Pos)) chrom.Dist <-
            chrom.Dist[-length(chrom.Dist)]
          
        }
        cat("Chromosome", i, "\n")
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
        if (!is.null(Pos.rel))
          res[[i]]$Pos.rel <- Pos.rel[Chrom==i]
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
        chrom.Pos <- 1:length(y)
        chrom.Dist <- rep(0, length(y)-1)
      }
      else {
        chrom.Pos <- Pos
        chrom.Dist <- diff(chrom.Pos) -1
        for (i in unique(Chrom)) 
          chrom.Dist[Chrom==i][length(chrom.Dist[Chrom==i])] <- NA
        chrom.Dist <- chrom.Dist[-length(chrom.Dist)]
      }
      if (!is.null(Dist)) chrom.Dist <- Dist
      ## We'll have to take out the last Dist of every Chrom
      res <- RJMCMC.NH.HMM.Metropolis(y=y, Chrom=Chrom, x=chrom.Dist, burnin=burnin, TOT=TOT,
                                      k.max=k.max, stat=stat, mu.alfa=mu.alfa, mu.beta=mu.beta,
                                      ka=ka, g=g, prob.k=prob.k, sigma.tau.mu=sigma.tau.mu,
                                      sigma.tau.sigma.2=sigma.tau.sigma.2,
                                      sigma.tau.beta=sigma.tau.beta, model=model,
                                      tau.split.mu=tau.split.mu, 
                                      tau.split.beta=tau.split.beta,
                                      start.k=start.k, RJ=RJ, auto.label=auto.label)
      res$Pos <- chrom.Pos
      res$Pos.rel <- Pos.rel
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
    res <- RJaCGH.one.array(y, Chrom=Chrom, Pos=Pos, Dist=Dist, model=model, burnin=burnin, TOT=TOT, k.max=k.max,
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
                         Pos=Pos, Dist=Dist, model=model, burnin=burnin, TOT=TOT, k.max=k.max,
                         stat=stat, mu.alfa=mu.alfa, mu.beta=mu.beta, ka=ka, g=g, prob.k=prob.k,
                         sigma.tau.mu=sigma.tau.mu, sigma.tau.sigma.2=sigma.tau.sigma.2, sigma.tau.beta=sigma.tau.beta,
                         tau.split.mu=tau.split.mu,
                         tau.split.beta=tau.split.beta,
                         start.k=start.k, RJ=RJ, auto.label=auto.label)


    }
    res$array.names <- colnames(y)
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
      res$mu <- unlist(lapply(dens, function(x) x$x[which.max(x$y)]))
      dens <- apply(object[[k]]$sigma.2, 2, density, bw="nrd0")
      res$sigma.2 <- unlist(lapply(dens, function(x) x$x[which.max(x$y)]))
      dens <- apply(object[[k]]$beta, c(1,2), density, bw="nrd0")
      res$beta <- unlist(lapply(dens, function(x) x$x[which.max(x$y)]))
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

  res$freq.joint.states <- obj[[k]]$freq.joint.states
  res$prob.cond.states <- obj[[k]]$freq.joint.states
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
  names(res$freq.joint.states) <- colnames(res$prob.states)
  names(res$prob.cond.states) <- colnames(res$prob.states)
  for(i in 1:k) {
    res$prob.cond.states[[i]] <- matrix(rbind(apply(res$freq.joint.states[[i]], 1,
                                       function(x) {
                                         if(sum(x)>0) x / sum(x)
                                         else rep(0, k)
                                       })), ncol=k)

    colnames(res$prob.cond.states[[i]]) <- colnames(res$prob.states)
    colnames(res$freq.joint.states[[i]]) <- colnames(res$prob.states)
  }
  res <- list(states=res$states, prob.states=res$prob.states,
              freq.joint.states= res$freq.joint.states,
              prob.cond.states=res$prob.cond.states)
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
  
  k.max <- max(as.numeric(as.character(obj$k)))
  for (i in 1:k.max) {

    if (probs[i] > 0) {
      objSummary <- states(obj, i)
      prob.states <- objSummary$prob.states
      for (j in 1:i) {
        if (length(grep("[L]", colnames(prob.states)[j]))) {
          Loss <- Loss + probs[i] * prob.states[,j]
        }
        else if (length(grep("[N]", colnames(prob.states)[j]))) {
          Normal <- Normal + probs[i] * prob.states[,j]
        }
        else if (length(grep("[G]", colnames(prob.states)[j]))) {
          Gain <- Gain + probs[i] * prob.states[,j]
        }
      }
    }
  }
  
  res$prob.states <- cbind(Loss, Normal, Gain)
  colnames(res$prob.states) <- c("Loss", "Normal", "Gain")
  res$states <- apply(res$prob.states, 1, function(x) names(which.max(x)))
  res$states <- ordered(res$states, levels=c("Loss", "Normal", "Gain"))
  res <- list(states=res$states, prob.states=res$prob.states)


}

model.averaging.RJaCGH.Chrom <-function(obj) {
  res <- list()
  for (ch in unique(obj$Chrom)) {
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
  k.max <- max(as.numeric(as.character(obj$k)))  
  for (i in 1:k.max) {
    if (probs[i] > 0) {
      objSummary <- states(obj, i)
      prob.states <- objSummary$prob.states
      for (j in 1:i) {
        if (length(grep("[L]", colnames(prob.states)[j]))) {
          Loss <- Loss + probs[i] * prob.states[,j]
        }
        else if (length(grep("[N]", colnames(prob.states)[j]))) {
          Normal <- Normal + probs[i] * prob.states[,j]
        }
        else if (length(grep("[G]", colnames(prob.states)[j]))) {
          Gain <- Gain + probs[i] * prob.states[,j]
        }
      }
    }
  }
  res$prob.states <- cbind(Loss, Normal, Gain)
  colnames(res$prob.states) <- c("Loss", "Normal", "Gain")
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
                              model.averaging=TRUE, cex=1, k=NULL,...)  {

  if (Chrom=="genome") {
    states <- NULL
    prob.states <- NULL
    if (model.averaging) {
      res <- model.averaging(x)
    }
    else {
      res <- states(x,k=k)
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
    ## small hack for the case with even number of chroms
    if (length(start.Chrom) %% 2 !=0) {
      xx <- xx[-c(length(xx)-1, length(xx))]
    }
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
    plot(x[[Chrom]], model.averaging=model.averaging, cex=cex, k=k)
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
  ## small hack for the case with even number of chroms
  if (length(start.Chrom) %% 2 !=0) {
    xx <- xx[-c(length(xx)-1, length(xx))]
  }
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
    ## small hack for the case with even number of chroms
    if (length(start.Chrom) %% 2 !=0) {
      xx <- xx[-c(length(xx)-1, length(xx))]
    }
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
   pch <- rep(16, length(states))
   ylim <- c(-100, 100)
   margin <- 4
   ylim[1] <- ylim[1] - margin - 1
   ylim[2] <- ylim[2] + margin
   ## Start of every Chromosome
   start.Chrom <- c(Pos[1], Pos[diff(x[[1]]$Chrom)!=0], Pos[length(Pos)])
   xx <- rep(start.Chrom, rep(2, length(start.Chrom)))
   ## small hack for the case with even number of chroms
   if (length(start.Chrom) %% 2 !=0) {
     xx <- xx[-c(length(xx)-1, length(xx))]
   }
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
    trace.plot(obj[[Chrom]], main.text=main.text)
  }

  if(class(obj)=="RJaCGH.array") {
    main.text <- past("Array", array, ". Chromosome", Chrom, ".")
    trace.plot(obj[[array]], chrom=chrom, main.text=main.text)
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
    newobj[[i]]$freq.joint.states <- list()
    for (j in 1:i) {
      newobj[[i]]$freq.joint.states[[j]] <-
        matrix(0,nrow=length(obj[[1]]$y)-1, ncol=i)
    }
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

     
      if (!is.null(obj[[i]][[j]]$prob.states)) {
        newobj[[j]]$prob.states <- newobj[[j]]$prob.states + 
          obj[[i]][[j]]$prob.states * nrow(obj[[i]][[j]]$mu)
        for (m in 1:j) {
          newobj[[j]]$freq.joint.states[[m]] <- newobj[[j]]$freq.joint.states[[m]] + 
          obj[[i]][[j]]$freq.joint.states[[m]]
        }
      }
      else {
        newobj[[j]]$prob.states <- newobj[[j]]$prob.states
        for (m in 1:j) {
          newobj[[j]]$freq.joint.states[[m]] <-
            newobj[[j]]$freq.joint.states[[m]]
        }
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
    }
    else {
      newobj[[j]]$prob.states <- NULL
      newobj[[j]]$freq.joint.states <- NULL
    }
  }
  newobj$k <- factor(newobj$k, levels=1:k)
  ## Recompute state labels
  for (i in 1:k) {
    if (table(newobj$k)[i] > 0) {
      p.labels <- 0.95 ## should be a parameter??
      obj.sum <- summary.RJaCGH(newobj, k=i, point.estimator="median")
      normal.levels <- (qnorm(mean=obj.sum$mu,
                              sd=sqrt(obj.sum$sigma.2),
                              p=(1-p.labels)/2) < 0 &
                        qnorm(mean=obj.sum$mu, sd=sqrt(obj.sum$sigma.2),
                              p=1-(1-p.labels)/2) > 0)

      ## bug fix to prevent a case
      ## in which normal.levels are
      ## ...TRUE FALSE TRUE...
      if (sum(normal.levels) <0) {
        normal.levels[!normal.levels][obj.sum$mu[!normal.levels] > min(obj.sum$mu[normal.levels]) &
                                      obj.sum$mu[!normal.levels] < max(obj.sum$mu[normal.levels])] <- 
                                        TRUE
      }
      ##
      ##

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
      ## changed 07-04-19. Need testing
      auto.label <- attr(obj[[1]], "auto.label")
      if (!is.null(auto.label)) {
        states <- states.RJaCGH(newobj, k=i)$states
        states <- factor(states, levels=names(summary.RJaCGH(newobj,
                                   k=i)$mu))
        freqs <- prop.table(table(states))
        labels <- names(freqs)
        means <- summary.RJaCGH(newobj, k=i)$mu
        means.1 <- abs(means - max(means[names(means)=='Normal']))
        means.2 <- abs(means - min(means[names(means)=='Normal']))
        means <- pmin(means.1, means.2)
        means.order <- order(means)
        names(means.order) <- names(means)[order(means)]
        means.order <- means.order[names(means.order) != 'Normal']
        tot.norm <- sum(freqs[names(freqs)=='Normal'])
        ind.tot <- 1
        while(tot.norm < auto.label) {
          labels[means.order][ind.tot] <- "Normal"
          tot.norm <- tot.norm + freqs[means.order][ind.tot]
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
  for (i in unique(obj[[1]]$Chrom)) {
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



#######################################################################################
#######################################################################################
## Minimal Common Regions (MCR)
## For now, only BMA
#######################################################################################
## Speed up: states as an argument
## or only compute them one time.
prob.seq <- function(obj, start, end, alteration="Gain") {
  if (start > end)
    stop("'start' must be an integer < than 'end'\n")
  if (alteration!="Gain" && alteration!="Loss")
    stop("'alteration' must be either 'Loss' or 'Gain'\n")
  K <- max(as.numeric(as.character(levels(obj$k))))
  n <- length(obj$y)
  probs <- prop.table(table(obj$k))
  prob.seq <- rep(0, K)
  ## states probs
  selected.states <- list()
  inds <- list()
  ## First gene: marginal prob.
  for (j in 1:K) {
    if(probs[j] > 0) {
      selected.states[[j]] <- states(obj,k=j)
      prob.states <- selected.states[[j]]$prob.states[start,]
      inds[[j]] <- grep(paste("[", substr(alteration, 1, 1), "]"), names(prob.states))
      if (length(inds[[j]]) > 0) {
          prob.seq[j] <- sum(prob.states[inds[[j]]])
      }
    }
  }
  if (end > start) {
    for(i in start:(end-1)) {
      for (j in 1:K) {
        if(probs[j] > 0) {
          if (length(inds[[j]]) > 0) {
            ## Conditional probs.
            freqs <- selected.states[[j]]$freq.joint.states[inds[[j]]]
            freqs <- lapply(freqs, function(x) x[i,])
            freqs <- apply(do.call("rbind", freqs), 2, sum)
            if(sum(freqs>0))
              prob.cond <- freqs / sum(freqs)
            ## is not a probability measure
            else
              prob.cond <- rep(0, j)
            prob.seq[j] <- prob.seq[j] * sum(prob.cond[inds[[j]]])
          }
        }
      }
    }
  }
  prob.seq <- as.numeric(prob.seq %*% probs)
  prob.seq
}
    

## Only with model.averaging (for now)
pMCR <- function(obj, p, alteration="Gain", 
                array.weights=NULL) {
  UseMethod("pMCR")
}

pMCR.RJaCGH <- function(obj, p, alteration="Gain", 
                       array.weights=NULL) {
  if (alteration !="Gain" && alteration!="Loss")
    stop ("'alteration' must be either 'Gain' or 'Loss'")
  if (!is.null(obj$Pos.rel))
    Pos <- obj$Pos.rel
  else
    Pos <- obj$Pos
  n <- length(obj$y)
  regions <- list()
  counter <- 1
  i <- 1
  marginal.probs <- model.averaging(obj)$prob.states[,alteration]
  while (i <=n) {
    prob <- marginal.probs[i]
    if (prob >= p) {
      regions[[counter]] <- list()
      regions[[counter]]$start <- Pos[i]
      regions[[counter]]$indexStart <- i
      regions[[counter]]$indexEnd <- i
      regions[[counter]]$end <- Pos[i]
      regions[[counter]]$genes <- 1
      regions[[counter]]$prob <- prob
      if (i <n) {
        prob <- prob.seq(obj, start=regions[[counter]]$indexStart,
                         end=regions[[counter]]$indexEnd + 1, alteration=alteration)

        while (i < n && prob >=p) {
          regions[[counter]]$end <- Pos[i+1]
          regions[[counter]]$indexEnd <- i + 1
          regions[[counter]]$genes <- regions[[counter]]$genes + 1
          regions[[counter]]$prob <- prob
          if (regions[[counter]]$indexEnd < n) {
            prob <- prob.seq(obj, start=regions[[counter]]$indexStart,
                             end=regions[[counter]]$indexEnd + 1, alteration=alteration)
          }
          i <- i + 1
        }
      }
      counter <- counter + 1
    }
    i <- i + 1
  }

  attr(regions, "alteration") <- alteration
  class(regions) <- "pMCR.RJaCGH"
  regions
}

##   ## #####################################################
##   ##   ## Experimental. If two regions collapse, we try to get
##   ##   ## the partition with highest posterior probability
##   if (test) { 
##     if (length(regions) > 1) {
##       for(i in 2:length(regions)) {
##         if (regions[[i-1]]$genes > 1 &&
##             regions[[i-1]]$indexEnd +1 == regions[[i]]$indexStart) {
##           newprob2 <- marginal.probs[regions[[i-1]]$indexEnd] *
##             prod(joint.probs[regions[[i-1]]$indexEnd:(regions[[i]]$indexEnd-1)])
##           newprob1 <- marginal.probs[regions[[i-1]]$indexStart]
##           if (regions[[i-1]]$genes > 2) {
##             newprob1 <- newprob1 *
##               prod(joint.probs[regions[[i-1]]$indexStart :
##                                (regions[[i-1]]$indexEnd - 2)])
##           }
##           while ((newprob1 + newprob2) >
##               (regions[[i-1]]$prob + regions[[i]]$prob) &&
##               newprob2 >= p) { 
##             ## Take into account Positions
##             regions[[i-1]]$indexEnd <- regions[[i-1]]$indexEnd - 1
##             regions[[i-1]]$end <- obj$Pos[regions[[i-1]]$indexEnd]
##             regions[[i-1]]$genes <- regions[[i-1]]$genes - 1
##             regions[[i-1]]$prob <- newprob1
##             regions[[i]]$start <- regions[[i-1]]$end
##             regions[[i]]$indexStart <- regions[[i]]$indexStart + 1
##             regions[[i]]$genes <- regions[[i]]$genes + 1
##             regions[[i]]$prob <- newprob2
##             ## Try another join
##             if (regions[[i-1]]$genes > 1) {
##               newprob2 <- marginal.probs[regions[[i-1]]$indexEnd] *
##                 prod(joint.probs[regions[[i-1]]$indexEnd:
##                                  (regions[[i]]$indexEnd-1)])
##               newprob1 <- marginal.probs[regions[[i-1]]$indexStart]
##               if (regions[[i-1]]$genes > 2) {
##                 newprob1 <- newprob1 *
##                   prod(joint.probs[regions[[i-1]]$indexStart :
##                                    (regions[[i-1]]$indexEnd - 2)])
##               }
##             }
##           }
##         }
##       }
##     }
##   }
## ## ##################################################
##   attr(regions, "alteration") <- alteration
##   class(regions) <- "pMCR.RJaCGH"
##   regions
## }

pMCR.RJaCGH.Chrom <- function(obj, p, alteration="Gain",
                       array.weights=NULL) {
  regions <- list()
  for (chr in unique(obj$Chrom)) {
    cat("Chromosome: ", chr, "\n")
    regions[[chr]] <- pMCR.RJaCGH(obj[[chr]], p=p,
                       alteration=alteration, 
                       array.weights=array.weights)
    attr(regions[[chr]], "Chrom") <- chr
  }
  attr(regions, "alteration") <- alteration
  class(regions) <- "pMCR.RJaCGH.Chrom"
  regions
}

pMCR.RJaCGH.genome <- function(obj, p, alteration="Gain",
                       array.weights=NULL) {
  if (alteration !="Gain" && alteration!="Loss")
    stop ("'alteration' must be either 'Gain' or 'Loss'")

  if (!is.null(obj$Pos.rel))
    Pos <- obj$Pos.rel
  else
    Pos <- obj$Pos

  regions <- list()
  i.Tot <- 1
  for(chr in unique(obj$Chrom)) {
    n <- length(obj$y[obj$Chrom==chr])
    regions[[chr]] <- list()
    counter <- 1
    i <- 1
    marginal.probs <- model.averaging(obj)$prob.states[obj$Chrom==chr,alteration]
    while (i <=n) {
      prob <- marginal.probs[i]
      if (prob >= p) {
        regions[[chr]][[counter]] <- list()
        regions[[chr]][[counter]]$start <- Pos[obj$Chrom==chr][i]
        regions[[chr]][[counter]]$indexStart <- i.Tot
        regions[[chr]][[counter]]$indexEnd <- i.Tot
        regions[[chr]][[counter]]$end <- Pos[obj$Chrom==chr][i]
        regions[[chr]][[counter]]$genes <- 1
        regions[[chr]][[counter]]$prob <- prob
        if (i < n) {
          prob <- prob.seq(obj, start=regions[[chr]][[counter]]$indexStart,
                           end=regions[[chr]][[counter]]$indexEnd + 1, alteration=alteration)
          while ((i < n) && (prob >=p)) {
            regions[[chr]][[counter]]$end <- Pos[obj$Chrom==chr][i+1]
            regions[[chr]][[counter]]$indexEnd <- i.Tot + 1
            regions[[chr]][[counter]]$genes <- regions[[chr]][[counter]]$genes + 1
            regions[[chr]][[counter]]$prob <- prob
            if (i < n-1) {
              prob <- prob.seq(obj, start=regions[[chr]][[counter]]$indexStart,
                               end=regions[[chr]][[counter]]$indexEnd + 1, alteration=alteration)
            }
            i <- i + 1
            i.Tot <- i.Tot + 1
          }
        }
        counter <- counter + 1
      }
      i <- i + 1
      i.Tot <- i.Tot + 1
    }
    attr(regions[[chr]], "Chrom") <- chr
  }
  attr(regions, "alteration") <- alteration
  class(regions) <- "pMCR.RJaCGH.genome"
  regions
}
  

pMCR.RJaCGH.array <- function(obj, p, alteration="Gain", 
                       array.weights=NULL) {
  if (alteration !="Gain" && alteration!="Loss")
    stop ("'alteration' must be either 'Gain' or 'Loss'")
  if (is.null(array.weights))
    array.weights <- rep(1/length(obj$array.names),
                       length(obj$array.names))
  else
    array.weights <- array.weights / sum(array.weights)

  if (inherits(obj[[obj$array.names[1]]], "RJaCGH.Chrom")) {
    regions <- pMCR.RJaCGH.array.Chrom(obj, p, alteration,
                                       array.weights)
  }
  if (inherits(obj[[obj$array.names[1]]], "RJaCGH.genome")) {
    regions <- pMCR.RJaCGH.array.genome(obj, p, alteration,
                                        array.weights)
  }
  
  if (!is.null(obj[[obj$array.names[1]]]$Pos.rel))
    Pos <- obj[[obj$array.names[1]]]$Pos.rel
  else
    Pos <- obj[[obj$array.names[1]]]$Pos

  if (inherits(obj[[obj$array.names[1]]], "RJaCGH")) {
    n <- length(obj[[obj$array.names[1]]]$y)
    regions <- list()
    counter <- 1
    i <- 1
    prob <- rep(0, length(obj$array.names))
    marginal.probs <- list()
    count <- 1
    for(n.array in obj$array.names) {
      marginal.probs[[count]] <- model.averaging(obj[[n.array]])$prob.states[,alteration]
      count <- count + 1
    }
    while (i <=n) {
      count <- 1
      for(n.array in obj$array.names) {
        prob[count] <- marginal.probs[[count]][i]
        count <- count +1
      }
      count <- 1
      if (prob %*% array.weights >= p) {
        regions[[counter]] <- list()
        regions[[counter]]$start <- Pos[i]
        regions[[counter]]$indexStart <- i
        regions[[counter]]$indexEnd <- i
        regions[[counter]]$end <- Pos[i]
        regions[[counter]]$genes <- 1
        regions[[counter]]$prob <- prob %*% array.weights
        if (i < n) {
          for(n.array in obj$array.names) {
            prob[count] <- prob.seq(obj[[n.array]], start=regions[[counter]]$indexStart,
                                    end=regions[[counter]]$indexEnd + 1, alteration=alteration)
            count <- count + 1
          }
          count <- 1
          while ((i < n) && (prob %*%array.weights >=p)) {
            regions[[counter]]$end <- Pos[i+1]
            regions[[counter]]$indexEnd <- i + 1
            regions[[counter]]$genes <- regions[[counter]]$genes + 1
            regions[[counter]]$prob <- prob %*% array.weights
            if (regions[[counter]]$indexEnd < n) {
              for(n.array in obj$array.names) {
                prob[count] <- prob.seq(obj[[n.array]], start=regions[[counter]]$indexStart,
                                        end=regions[[counter]]$indexEnd + 1, alteration=alteration)
                count <- count + 1
              }
              count <- 1
            }
            i <- i + 1
          }
        }
        counter <- counter + 1
      }
      i <- i + 1
    }
    attr(regions, "alteration") <- alteration
    class(regions) <- "pMCR.RJaCGH.array"
  }
  regions
}

pMCR.RJaCGH.array.Chrom <- function(obj, p, alteration="Gain",
                       array.weights=NULL) {
  regions <- list()
  for (chr in unique(obj[[1]]$Chrom)) {
    cat("Chromosome: ", chr, "\n")
    if (!is.null(obj[[obj$array.names[1]]][[chr]]$Pos.rel))
      Pos <- obj[[obj$array.names[1]]][[chr]]$Pos.rel
    else
      Pos <- obj[[obj$array.names[1]]][[chr]]$Pos

    n <- length(obj[[obj$array.names[1]]][[chr]]$y)
    regions[[chr]] <- list()
    counter <- 1
    i <- 1
    prob <- rep(0, length(obj$array.names))
    marginal.probs <- list()
    count <- 1
    for(n.array in obj$array.names) {
      marginal.probs[[count]] <- model.averaging(obj[[n.array]][[chr]])$prob.states[,alteration]
      count <- count + 1
    }
    while (i <=n) {
      count <- 1
      for(n.array in obj$array.names) {
        prob[count] <- marginal.probs[[count]][i]
        count <- count +1
      }
      count <- 1
      if (prob %*% array.weights >= p) {
        regions[[chr]][[counter]] <- list()
        regions[[chr]][[counter]]$start <- Pos[i]
        regions[[chr]][[counter]]$indexStart <- i
        regions[[chr]][[counter]]$indexEnd <- i
        regions[[chr]][[counter]]$end <- Pos[i]
        regions[[chr]][[counter]]$genes <- 1
        regions[[chr]][[counter]]$prob <- prob %*% array.weights
        if (i < n) {
          for(n.array in obj$array.names) {
            prob[count] <- prob.seq(obj[[n.array]][[chr]],
                                    start=regions[[chr]][[counter]]$indexStart,
                                    end=regions[[chr]][[counter]]$indexEnd + 1,
                                    alteration=alteration)
            count <- count + 1
          }
          count <- 1
          while ((i < n) && (prob %*%array.weights >=p)) {
            regions[[chr]][[counter]]$end <- Pos[i+1]
            regions[[chr]][[counter]]$indexEnd <- i + 1
            regions[[chr]][[counter]]$genes <- regions[[chr]][[counter]]$genes + 1
            regions[[chr]][[counter]]$prob <- prob %*% array.weights
            if (regions[[chr]][[counter]]$indexEnd < n-1) {
              for(n.array in obj$array.names) {
                prob[count] <- prob.seq(obj[[n.array]][[chr]],
                                        start=regions[[chr]][[counter]]$indexStart,
                                        end=regions[[chr]][[counter]]$indexEnd + 1,
                                        alteration=alteration)
                count <- count + 1
              }
              count <- 1
            }
            i <- i + 1
          }
        }
        counter <- counter + 1
      }
      i <- i + 1
    }
    attr(regions[[chr]], "Chrom") <- chr
  }
  attr(regions, "alteration") <- alteration
  class(regions) <- "pMCR.RJaCGH.array.Chrom"
  regions

}


pMCR.RJaCGH.array.genome <- function(obj, p, alteration="Gain",
                       array.weights=NULL) {
  if (alteration !="Gain" && alteration!="Loss")
    stop ("'alteration' must be either 'Gain' or 'Loss'")
  regions <- list()
  i.Tot <- 1
  Chrom <- obj[[obj$array.names[1]]]$Chrom
  for(chr in unique(Chrom)) {
    if (!is.null(obj[[obj$array.names[1]]]$Pos.rel))
      Pos <- obj[[obj$array.names[1]]]$Pos.rel
    else
      Pos <- obj[[obj$array.names[1]]]$Pos

    n <- length(obj[[obj$array.names[1]]]$y[Chrom==chr])
    regions[[chr]] <- list()
    counter <- 1
    i <- 1
    prob <- rep(0, length(obj$array.names))
    marginal.probs <- list()
    count <- 1
    for(n.array in obj$array.names) {
      marginal.probs[[count]] <-
        model.averaging(obj[[n.array]])$prob.states[Chrom==chr,alteration]
      count <- count + 1
    }
    while (i <=n) {
      count <- 1
      for(n.array in obj$array.names) {
        prob[count] <- marginal.probs[[count]][i]
        count <- count +1
      }
      count <- 1
      if (prob %*% array.weights >= p) {
        regions[[chr]][[counter]] <- list()
        regions[[chr]][[counter]]$start <- Pos[Chrom==chr][i]
        regions[[chr]][[counter]]$indexStart <- i.Tot
        regions[[chr]][[counter]]$indexEnd <- i.Tot
        regions[[chr]][[counter]]$end <- Pos[Chrom==chr][i]
        regions[[chr]][[counter]]$genes <- 1
        regions[[chr]][[counter]]$prob <- prob %*% array.weights
        if (i < n) {
          for(n.array in obj$array.names) {
            prob[count] <- prob.seq(obj[[n.array]],
                                    start=regions[[chr]][[counter]]$indexStart,
                                    end=regions[[chr]][[counter]]$indexEnd + 1,
                                    alteration=alteration)
            count <- count + 1
          }
          count <- 1

          while ((i < n) && (prob %*%array.weights >=p)) {
            regions[[chr]][[counter]]$end <- Pos[Chrom==chr][i+1]
            regions[[chr]][[counter]]$indexEnd <- i.Tot + 1
            regions[[chr]][[counter]]$genes <- regions[[chr]][[counter]]$genes + 1
            regions[[chr]][[counter]]$prob <- prob %*% array.weights
            if (i < n-1) {
              for(n.array in obj$array.names) {
                prob[count] <- prob.seq(obj[[n.array]],
                                        start=regions[[chr]][[counter]]$indexStart,
                                        end=regions[[chr]][[counter]]$indexEnd + 1,
                                        alteration=alteration)
                count <- count + 1
              }
              count <- 1
            }
            i <- i + 1
            i.Tot <- i.Tot + 1
          }
        }
        counter <- counter + 1
      }
      i <- i + 1
      i.Tot <- i.Tot + 1
    }
    attr(regions[[chr]], "Chrom") <- chr
  }
  attr(regions, "alteration") <- alteration
  class(regions) <- "pMCR.RJaCGH.array.genome"
  regions
}



##What about joint prob of all regions (chrom/genome)?

print.pMCR.RJaCGH <- function(x,...) {

  res <- sapply(x, function(y) c(y$start, y$end, y$genes, y$prob)
         )
  res <- t(res)
  if (ncol(res)>0) {
    colnames(res) <-   c("Start", "End", "#Genes",
                         paste("Prob.", attr(x, "alteration")))
    print(res)
  }
  else {
    res <- "No common minimal regions found"
    print(paste(res, "\n", sep=""))
  }
  
}

print.pMCR.RJaCGH.Chrom <- function(x,...) {

  res <- sapply(x, function(y) {
    sapply(y, function(z) c(attr(y, "Chrom"), z$start, z$end, z$genes,
    z$prob))
  })
  if(sum(sapply(res, function(x) length(x)))>0) {
    res <- do.call("cbind", res[!sapply(res, is.list)])
    res <- t(res)
    colnames(res) <-   c("Chromosome", "Start", "End", "#Genes",
                         paste("Prob.", attr(x, "alteration")))
    print(res)
  }
  else {
    res <- "No common minimal regions found"
    print(paste(res, "\n", sep=""))
  }

}
print.pMCR.RJaCGH.genome <- function(x,...) {

  res <- sapply(x, function(y) {
    sapply(y, function(z) c(attr(y, "Chrom"), z$start, z$end, z$genes,
    z$prob))
  })
  if(sum(sapply(res, function(x) length(x)))>0) {
    res <- do.call("cbind", res[!sapply(res, is.list)])
    res <- t(res)
    colnames(res) <-   c("Chromosome", "Start", "End", "#Genes",
                         paste("Prob.", attr(x, "alteration")))
    print(res)
  }
  else {
    res <- "No common minimal regions found"
    print(paste(res, "\n", sep=""))
  }

}
       
print.pMCR.RJaCGH.array <- function(x,...) {

  res <- sapply(x, function(y) c(y$start, y$end, y$genes, y$prob)
         )
  res <- t(res)
  if (ncol(res)>0) {
    colnames(res) <-   c("Start", "End", "#Genes",
                         paste("Prob.", attr(x, "alteration")))
    print(res)
  }
  else {
    res <- "No minimal common regions found"
    print(paste(res, "\n",sep=""))
  }
  
}



print.pMCR.RJaCGH.array.Chrom <- function(x,...) {


  res <- sapply(x, function(y) {
    sapply(y, function(z) c(attr(y, "Chrom"), z$start, z$end, z$genes,
                            z$prob))
  })
  if(sum(sapply(res, function(x) length(x)))>0) {
    res <- do.call("cbind", res[!sapply(res, is.list)])
    res <- t(res)
    colnames(res) <-   c("Chromosome", "Start", "End", "#Genes",
                         paste("Prob.", attr(x, "alteration")))
    print(res)
  }
  else {
    res <- "No common minimal regions found"
    print(paste(res, "\n", sep=""))
  }
}
print.pMCR.RJaCGH.array.genome <- function(x,...) {


  res <- sapply(x, function(y) {
    sapply(y, function(z) c(attr(y, "Chrom"), z$start, z$end, z$genes,
                            z$prob))
  })
  if(sum(sapply(res, function(x) length(x)))>0) {
    res <- do.call("cbind", res[!sapply(res, is.list)])
    res <- t(res)
    colnames(res) <-   c("Chromosome", "Start", "End", "#Genes",
                         paste("Prob.", attr(x, "alteration")))
    print(res)
  }
  else {
    res <- "No common minimal regions found"
    print(paste(res, "\n", sep=""))
  }
}

genome.plot <- function(obj, col=NULL, breakpoints=NULL) {
  if (!is.null(col) & !is.null(breakpoints)) {
    if(length(col) != length(breakpoints) + 1)
      stop("length(col) must be length(breakpoints + 1\n")
    if(length(breakpoints) < 2)
      stop("length(breakpoints) must be at least 2\n")
   }

  if(inherits(obj, "RJaCGH.array")) {
    Chrom <- obj[[obj$array.name[1]]]$Chrom
    Pos <- obj[[obj$array.name[1]]]$Pos
    if (!is.null(obj[[obj$array.name[1]]]$Pos.rel))
      Pos <- obj[[obj$array.name[1]]]$Pos.rel

    probs <- model.averaging(obj)
    probs <- lapply(probs, function(x) lapply(x, function(y) y$prob.states))
    probs <- lapply(probs, function(x) do.call("rbind", x))
    average <- matrix(0, nrow=nrow(probs[[1]]), ncol=3)
    y.mean <- rep(0, nrow(probs[[1]]))
    y <- lapply(obj, function(x) lapply(x, function(y) y$y))
    y <- lapply(y, function(x) do.call("c", x))
    y$array.names <- NULL
    for (i in 1:length(probs)) {
      average <- average + probs[[i]]
      y.mean <- y.mean + y[[i]]
    }
    average <- average / length(probs)
    y.mean <- y.mean / length(probs)
    y <- y.average
    probs <- average
  }
  
  if(inherits(obj, "RJaCGH.genome") | inherits(obj, "RJaCGH.Chrom")) {
    Chrom <- obj$Chrom
    Pos <- obj$Pos
    if (!is.null(obj$Pos.rel))
      Pos <- obj$Pos.rel
    if(inherits(obj, "RJaCGH.Chrom") & !is.null(obj[[1]]$Pos.rel)) {
      Pos <- numeric(0)
      for(chr in unique(obj$Chrom))
        Pos <- c(Pos, obj[[chr]]$Pos.rel)
    }
                    
    
    probs <- model.averaging(obj)
    y <- obj$y
    if (inherits(obj, "RJaCGH.Chrom")) {
      probs <-
        do.call("rbind",
                lapply(probs, function(x) x$prob.states))
      for (chr in unique(Chrom)) {
        y <- c(y, obj[[chr]]$y)
      }
    }
    else probs <- probs$prob.states
    index <- 1*(y<0) + 3*(y>=0)
    colo <- rep(0, length(index))
    for(i in 1:length(index)) {
      colo[i] <- probs[i,index[i]]
    }
  }
  colo[index==1] <- -colo[index==1]
  colo.round <- floor(colo*10) / 10
  colo.recoded <- rep(0, length(colo.round))
  
    if (is.null(col)) {
    col <- colors()
    col <- col[c(51, 50, 86, 24, 404, 552, 555)]
  }
  if(is.null(breakpoints)) {
    breakpoints <- c(-0.9, -0.7, -0.5, 0.5, 0.7, 0.9)
  }
  MidPoint <- floor(length(col)/2)
  colo.recoded[colo.round <= -breakpoints[1]] <- col[1]
  label.legend <- paste("P.Loss >= ", -breakpoints[1], sep="")
  if (length(breakpoints) > 2) {
    for (i in 2:MidPoint) {
      colo.recoded[colo.round > breakpoints[i-1] & colo.round <=
    breakpoints[i]] <- col[i]
      label.legend <- c(label.legend, paste(breakpoints[i],
                                            " > P.Loss >= ", -breakpoints[i-1], sep=""))
    }
  }
  colo.recoded[colo.round > breakpoints[MidPoint] & colo.round <
    breakpoints[MidPoint +1]] <- col[MidPoint +1]
  label.legend <- c(label.legend, paste("P.Loss < ", -breakpoints[MidPoint],
                                        " or P.Gain >= ", breakpoints[MidPoint+1], sep=""))
  
  if (length(breakpoints) > 2) {
    for (i in (MidPoint+2):length(breakpoints)) {
      colo.recoded[colo.round >=breakpoints[i-1] & colo.round
    < breakpoints[i]] <- col[i+1]
      label.legend <- c(label.legend, paste(breakpoints[i],
                                            " > P.Gain >= ", breakpoints[i-1], sep=""))
    }
  }

  colo.recoded[colo.round >=breakpoints[length(breakpoints)]] <-
  col[length(col)]
  label.legend <- c(label.legend, paste("P.Gain >= ",
                                        breakpoints[length(breakpoints)], sep=""))
  n.chrom <- length(unique(Chrom))
  
  par(omi=c(0,0,0,0))
  par(mfrow=c(1, 2))
  xmax <- max(Pos)
  plot(0,0, type="n", xlim=c(0, xmax), ylim=c(1, ceiling(n.chrom/2)),
       axes=FALSE, ylab="Chromosome", xlab="")
  axis(side=2, at=c(1:ceiling(n.chrom/2)), labels=unique(Chrom)[1:ceiling(n.chrom/2)])

  for(i in 1:ceiling(n.chrom/2)) {

    lines(c(0, max(Pos[Chrom==i])), c(i,i))
    points(Pos[Chrom==i], i + y[Chrom==i]/ (2*max(abs(y))),
           pch=19,col=colo.recoded[Chrom==i], cex=0.5)
  }
  plot(0,0, type="n", xlim=c(0, xmax), ylim=c(1, ceiling(n.chrom/2)),
       axes=FALSE, ylab="Chromosome", xlab="")
  axis(side=2, at=1:(n.chrom - ceiling(n.chrom/2)), labels=unique(Chrom)[(ceiling(n.chrom/2)+1):n.chrom])

  for(i in (ceiling(n.chrom/2) +1):n.chrom) {
    lines(c(0, max(Pos[Chrom==i])), c(i-12,i-12))
    points(Pos[Chrom==i], i -ceiling(n.chrom/2) + y[Chrom==i]/ (2*max(abs(y))),
           pch=16,col=colo.recoded[Chrom==i], cex=0.5)
  }
  legend(max(Pos[Chrom==(floor(n.chrom/2)+3)]), 8, legend=label.legend,
         col=col, pch=19)
  
}
