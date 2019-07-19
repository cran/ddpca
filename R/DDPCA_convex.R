DDPCA_convex <-
function(Sigma,lambda,rho = 20,max_iter_convex = 50){
  p = nrow(Sigma)
  A = matrix(0,p,p)
  E = matrix(0,p,p)
  Lambda = matrix(0,p,p)
  
  for (i in 1:max_iter_convex){
    # M-step
    svd_object = svd(Sigma - E - A - Lambda/rho)
    D = diag(svd_object$d)
    D_s = sign(D)*pmax(abs(D) - lambda/rho,0)
    L = svd_object$u%*%D_s%*%t(svd_object$v)
    
    # A-step
    temp2 = Sigma - L - E - Lambda/rho
    A = ProjSDD(temp2)
    
    # E-step
    E = rho/(rho+1)*(Sigma - A - L - Lambda/rho)
    
    # Lambda-step
    Lambda = Lambda + rho*(A + L + E - Sigma)
  }
  
  return(list(L=L,A=A))
}
