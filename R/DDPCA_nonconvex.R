DDPCA_nonconvex <-
function(Sigma,K,max_iter_nonconvex = 15,SDD_approx = TRUE, max_iter_SDD = 20,eps = NA){
  S = Sigma
  for (i in 1:max_iter_nonconvex){
    eig_object = eigs(Sigma,K)
    if (K>1){
      D = diag(eig_object$values)
    } else {
      D = eig_object$values
    }
    V = eig_object$vectors
    L = V%*%D%*%t(V)
    if (SDD_approx) {
      A = ProjDD(Sigma - L)
      A_sym = (A + t(A))/2
    } else {
      A_sym = ProjSDD(Sigma - L, max_iter_SDD,eps)
    }
    S = Sigma - A_sym
  }
  return(list(L=L,A=A_sym))
}
