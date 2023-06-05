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

merge_pvalue <- function(chain_fission, method = c("KS", "Fisher", "harmonic", "geometric")){
  if (method == "KS")
    pval_chain <- apply(chain_fission, 2, KS_merge)
  if (method == "Fisher")
    pval_chain <- apply(chain_fission, 2, Fisher_merge)
  if (method == "geometric")
    pval_chain <- apply(chain_fission, 2, geometric_merge)
  if (method == "harmonic")
    pval_chain <- apply(chain_fission, 2, harmonic_merge)
  return(pval_chain)
}
