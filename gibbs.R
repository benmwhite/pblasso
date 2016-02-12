#Latent variable probit implementation (laplace priors)

#gibbs sampling function for bayesian lasso
library(mnormt) # rmnorm
library(VGAM) # rinv.gaussian, rlaplace
library(miscTools) # colMedians
library(truncnorm) #rtruncnorm
library(ggplot2)

#y must be a vector of 0's and 1's, x must be a matrix with rows representing
#observations and columns representing features (e.g. genes)

pblasso <- function(y, x, max.steps = 1000, burn = 50,
                          thin = 5, r = 0.5, delta = 0.5,
                          progress = TRUE) { 
  n <- nrow(x)
  m <- ncol(x)
  n_ones <- sum(y)
  
  #centering and scaling
  x <- scale(x)
    
  XtX <- t(x) %*% x	#Time saving
  
  betaSamples <- matrix(0, max.steps, m)
  sigma2Samples <- rep(0, max.steps)
  invTau2Samples <- matrix(0, max.steps, m)
  lambdaSamples <- rep(0, max.steps)
  
  colnames(betaSamples) <- colnames(x)
  
  #initial values
  beta <- rlaplace(m)
  residuals <- drop(y - x %*% beta)
  sigma2 <- drop((t(residuals) %*% residuals) / n)
  invTau2 <- 1 / (beta * beta)
  lambda <- m * sqrt(sigma2) / sum(abs(beta))
  
  k <- 0
  while (k < max.steps) {
    k <- k + 1
    
    if ((k %% 100 == 0) & progress) {
      print(paste('Iteration:', k))
    }
    
    # sample z (latent variable)
    z <- rep(NA, n)
    for (j in 1:n) {
      if (y[j] == 1) {
        z[j] <- rtruncnorm(1, a = 0, mean = x[j, ] %*% beta)
      } else { #y_j = 0
        z[j] <- rtruncnorm(1, b = 0, mean = x[j, ] %*% beta)
      }
    } 
    xz <- t(x) %*% z
    
    # sample beta
    invD <- diag(invTau2)
    invA <- solve(XtX + invD)
    mean <- invA %*% xz
    varcov <- sigma2 * invA
    beta <- drop(rmnorm(1, mean, varcov))
    betaSamples[k,] <- beta
    
    # sample sigma2
    shape <- (n+m-1)/2
    residuals <- drop(z - x %*% beta)
    scale <- (t(residuals) %*% residuals + t(beta) %*% invD %*% beta)/2
    sigma2 <- 1/rgamma(1, shape, 1/scale)
    sigma2Samples[k] <- sigma2
    
    # sample tau2
    muPrime <- sqrt(lambda^2 * sigma2 / beta^2)
    lambdaPrime <- lambda^2
    invTau2 <- rep(0, m)
    for (i in seq(m)) {
      invTau2[i] <- rinv.gaussian(1, muPrime[i], lambdaPrime)
    }
    invTau2Samples[k, ] <- invTau2
    
    # update lambda
    shape = r + m/2
    scale = delta + sum(1/invTau2)/2
    lambda <- rgamma(1, shape, 1/scale)
    lambdaSamples[k] <- lambda
  }  
  return(betaSamples[seq(from = burn+1, to = max.steps, by = thin),])
}

#parallelized version (unix only)
mc_pblasso <- function(y, x, max.steps = 1000, burn = 50,
                       thin = 5, nchains = 2, r = 0.5, delta = 0.5,
                       tag.chain = FALSE) {
  tmp_func <- function(whatevs) {
    pblasso(y, x, max.steps, burn, thin, r, delta, progress = FALSE)
  }
  results <- mclapply(as.list(seq_len(nchains)), tmp_func, mc.cores = nchains)
  if (tag.chain == TRUE) {
    tag <- function(l) {
      for (i in seq_len(length(l))) {
        l[[i]] <- cbind("CHAIN" = i, l[[i]])
      }
      return(l)
    }
    results <- tag(results)
  }
  results <- do.call(rbind, results)
  return(results)
}