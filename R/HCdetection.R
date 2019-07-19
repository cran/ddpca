HCdetection <-
function(p,alpha=0.5,pvalcut=NA){
  n = length(p)
  if (is.na(pvalcut)){pvalcut = 1/n}
  kk = seq(1,n)/(n+1)
  psort = sort(p)
  HC = sqrt(n)*(kk-psort)/sqrt(psort-psort^2)
  HC[psort<pvalcut] = -Inf
  HC[(round(alpha*n)+1):n] = -Inf
  HCT = max(HC)
  H = (HCT/sqrt(2*log(log(n)))/(1+1/sqrt(log(log(n))))>1)*1
  return(list(H=H,HCT=HCT))
}
