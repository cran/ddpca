DDHC <-
function(X,known_Sigma = NA,method='nonconvex',K = 1,lambda = 3,max_iter_nonconvex = 15,SDD_approx = TRUE, max_iter_SDD = 20,eps = NA,
                  rho = 20,max_iter_convex = 50,alpha=0.5,pvalcut=NA){
    n = nrow(X)
    p = ncol(X)
    X_bar = colMeans(X)
    if (is.na(known_Sigma)) {
      # need to estimate Sigma
      Sigma_1 = cov(X)/n
    } else {
      # use true Sigma
      Sigma_1 = known_Sigma/n
    }
    if (method=='nonconvex'){
      result = DDPCA_nonconvex(Sigma_1,K,max_iter_nonconvex,SDD_approx,max_iter_SDD,eps)
    } else {
      result = DDPCA_convex(Sigma_1,lambda,rho,max_iter_convex)
      K = rankMatrix(result$L)
    }
    eig_object = eigs(Sigma_1,K)
    
    if (K>1){
      D = diag(eig_object$values)
    } else {
      D = eig_object$values
    }
    V = eig_object$vectors
    b = V%*%sqrt(D)
    resi = rq(X_bar~b,tau=0.5)$residuals
    p_value = 2*(1-pnorm(abs(resi/(sqrt(diag(result$A))))))
    HC_result = HCdetection(p_value,alpha,pvalcut)
    return(HC_result)
  }