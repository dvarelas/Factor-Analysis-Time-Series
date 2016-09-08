#################################################
## pc under specification 2 (heteroscedascity) ##
#################################################

# M is time series data
# r is the number of static factors to be extracted

pca2 = function(M,r) {
  if (r>ncol(M)) {
    print ("r must be less than")
    print(ncol(M))
  }
  if (r<=ncol(M)&r>1) {
    R=cor(M, use="complete.obs")
    load = eigen(R)$vectors[,1:r]%*%diag(eigen(R)$values[1:r]^(1/2))
    Psi = diag((R-load%*%t(load)))
    Psi = diag(Psi)
    fac_scores = solve(t(load)%*%solve(Psi)%*%load+diag(1,nrow=r))%*%t(load)%*%solve(Psi)%*%t(M)
    fac_scores = as.vector(fac_scores)
    fac_scores = as.matrix(fac_scores, nrow = nrow(M), ncol=r)
    fac_scores = matrix(fac_scores, nrow = nrow(M), ncol=r)
  }
  if (r==1) {
    R=cor(M, use="complete.obs")
    load = eigen(R)$vectors[,1:r]*eigen(R)$values[1:r]^(1/2)
    Psi = diag((R-load%*%t(load)))
    Psi = diag(Psi)
    fac_scores = solve(t(load)%*%solve(Psi)%*%load+diag(1,nrow=r))%*%t(load)%*%solve(Psi)%*%t(M)
    fac_scores = as.vector(fac_scores)
  }
  cum = sum(eigen(R)$values[1:r]/ncol(M))
  lista = list("factor_scores"=fac_scores,"loadings"=load, "Cumulative Proportion of Variance Explained"= cum)
  return(lista)
}
