KS_merge <- function(p){
  ks_res <- stats::ks.test(p, y = "punif", alternative = "g")
  return(ks_res$p.value)
}

Fisher_merge <- function(pval){
  #@pval a vector of p-values
  k <- length(pval)
  chi <- -2*sum(log(pval))
  return(stats::pchisq(chi, df = 2*k, lower.tail = F))
}

geometric_merge <- function(p){
  K <- length(p)
  return(min(exp(1)*(prod(p))^(1/K),1))
}

harmonic_merge <- function(p){
  K <- length(p)
  return(min((exp(1)*log(K))*(K/(sum(1/p))),1))
}

box_muller <- function(p){
  U <- stats::runif(length(p))
  Z <- sqrt(-2*log(U))* cos(2*pi*p)
  return(stats::shapiro.test(Z)$p.value)
}


merge_pvalue <- function(chain_fission, method = c("KS", "shapiro", "Fisher", "harmonic", "geometric")){
  if (method == "KS")
    pval_chain <- KS_merge(chain_fission)
  if (method == "shapiro")
    pval_chain <- box_muller(chain_fission)
  if (method == "Fisher")
    pval_chain <- Fisher_merge(chain_fission)
  if (method == "geometric")
    pval_chain <- geometric_merge(chain_fission)
  if (method == "harmonic")
    pval_chain <- harmonic_merge(chain_fission)
  return(pval_chain)
}
