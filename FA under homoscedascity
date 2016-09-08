################################################
## pc under specification 1 (homoscedascity) ###
################################################

# M is time series data
# r is the number of static factors to be extracted


pca1 = function(M,r) {
  if (r>ncol(M)) {
    print ("r must be less than")
    print(ncol(M))
  }
  if (r==ncol(M)) {
    R=cor(M, use="complete.obs")
    fac_scores = solve(diag(eigen(R)$values[1:r]^(1/2)))%*%t(eigen(R)$vectors[,1:r])%*%t(M)
    load = eigen(R)$vectors[,1:r]%*%diag(eigen(R)$values[1:r]^(1/2))
    
  }
  if (r<ncol(M)&r>1) {
    R=cor(M, use="complete.obs")
    load = eigen(R)$vectors[,1:r]%*%diag(eigen(R)$values[1:r]^(1/2))
    Psi_bar = (sum(diag(R))-sum(diag(eigen(R)$values[1:r])))/ncol(M)
    fac_scores = solve(diag(eigen(R)$values[1:r])+Psi_bar*diag(1,nrow=r))%*%diag(eigen(R)$values[1:r]^(1/2))%*%t(eigen(R)$vectors[,1:r])%*%t(M)
  }
  if (r==1) {
    R=cor(M, use="complete.obs")
    load = eigen(R)$vectors[,1:r]*(eigen(R)$values[1:r]^(1/2))
    Psi_bar = (sum(diag(R))-eigen(R)$values[1:r])/ncol(M)
    fac_scores = as.numeric(solve(eigen(R)$values[1:r]+Psi_bar)*(eigen(R)$values[1:r]^(1/2)))*(t(eigen(R)$vectors[,1:r]))%*%t(M)
    fac_scores = as.vector(fac_scores)
  }
  else {
    fac_scores = as.vector(fac_scores)
    fac_scores = as.matrix(fac_scores, nrow = nrow(M), ncol=r)
    fac_scores = matrix(fac_scores, nrow = nrow(M), ncol=r)
  }
  cum = sum(eigen(R)$values[1:r]/ncol(M))
  lista = list("factor_scores"=fac_scores,"loadings"=load, "Cumulative Proportion of Variance Explained"= cum)
  return(lista)
}
