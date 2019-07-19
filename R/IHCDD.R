IHCDD <-
function(X,method='nonconvex',K = 1,lambda = 3,max_iter_nonconvex = 15,SDD_approx = TRUE, max_iter_SDD = 20,eps = NA,
           rho = 20,max_iter_convex = 50,alpha=0.5,pvalcut=NA){
    n = nrow(X)
    X_bar = colMeans(X)
    Sigma_1 = cov(X)/n
    if (method=='nonconvex'){
      result = DDPCA_nonconvex(Sigma_1,K,max_iter_nonconvex,SDD_approx,max_iter_SDD,eps)
    } else {
      result = DDPCA_convex(Sigma_1,lambda,rho,max_iter_convex)
    }
    Sigma_hat = result$L + result$A
    p_value = 2*(1 - pnorm(abs(solve(Sigma_hat,X_bar)/diag(solve(Sigma_hat)))));
    HC_result = HCdetection(p_value,alpha,pvalcut)
    return(HC_result)
  }
