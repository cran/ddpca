ProjSDD <-
function(A,max_iter_SDD=20,eps=NA){
  p = nrow(A)
  G = A
  I = matrix(0,p,p)
  for (i in 1:max_iter_SDD){
    PsG = (G+t(G))/2
    newG = ProjDD(PsG-I)
    I = newG - (PsG-I)
    if (!is.na(eps) && norm(G-newG,'F')<eps){
      break
    }
    G = newG
  }
  return(PsG)
}
