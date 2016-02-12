#diagnostic tools for simulations and model fitting
#with probit_blasso()
#source("gibbs.R")

####generating test data####
gen_testdata <- function(n = 20, p = 50, ytype = "binary", 
                         npos = 5, nneg = 5, 
                         truepos = 1, trueneg = -1) {
  #n = number of obs
  #p = number of variables
  #ytype = "binary" or "continuous"
  #npos = number of positively associated variables to generate
  #nneg = negatively associated variables
  #truepos = true beta value for positives
  #trueneg = true beta value for negatives
  
  true_betas <- c(rep(truepos, npos), rep(trueneg, nneg), 
                  rep(0, p - npos - nneg))
  
  X <- replicate(p - npos - nneg, rnorm(n))
  #genes positively associated with outcome
  pos <- rbind(replicate(npos, rnorm(floor(n/2), 2, 0.1)), 
               replicate(npos, rnorm(ceiling(n/2), -2, 0.1)))
  #genes negatively associated with outcome
  neg <- rbind(replicate(nneg, rnorm(floor(n/2), -2, 0.1)), 
               replicate(nneg, rnorm(ceiling(n/2), 2, 0.1)))
  X <- scale(cbind(pos, neg, X)) #normalizing
  Y <- drop(X %*% true_betas) + rnorm(n, sd = 0.1)
  if (ytype == "binary") {
    return(list(X = X, Y = as.numeric(Y > 0)))
  }
  if (ytype == "continuous") {
    return(list(X = X, Y = Y))
  }  
}

####Variable selection based on available posterior samples####
prob_pos <- function(vec) {
  mean(vec > 0)
}
prob_neg <- function(vec) {
  mean(vec < 0)
}
effect_select <- function(samps, cutoff) {
  pos <- apply(samps, 2, prob_pos)
  neg <- apply(samps, 2, prob_neg)
  out <- list(pos = names(which(pos > cutoff)),
              neg = names(which(neg > cutoff)))
  return(out)
}

#function for filtering columns (genes, proteins) by variance
filter_var <- function(mat, q = NULL, ntop = NULL) {
  var_vec <- apply(mat, 2, var)
  if (is.numeric(q)) {
    cutoff <- quantile(var_vec, probs = q)
  }
  if (is.numeric(ntop)) {
    cutoff <- quantile(var_vec, probs = 1 - ntop/length(var_vec))
  }
  out <- mat[ , var_vec > cutoff]
  return(out)
}

