ProjDD <-
function(C){
  p = nrow(C)
  A = matrix(0,p,p)
  C_d = diag(C)
  C_od = abs(C-diag(C_d))
  C_ods = rowSums(C_od)
  ind1 = (C_d>=C_ods)
  ind2 = (C_d<C_ods)&(C_d>=0)
  ind3 = (C_d>=-C_ods)&(C_d<0)&(abs(C_d)>apply(C,1,max))
  ind4 = (C_d>=-C_ods)&(C_d<0)&(!(abs(C_d)>apply(C,1,max)))
  ind5 = (C_d<(-1*C_ods))
  A[ind1,] = C[ind1,]
  A[ind3 | ind5,] = 0
  pos = which(ind2 | ind4)
  if (length(pos)==0) {return(A)}
  
  for (l in 1:length(pos)){
    j = pos[l]
    Cjc1 = C[j,]
    if (j==1){
      Cjc = Cjc1[2:p]
    } else if (j==p){
      Cjc = Cjc1[1:(p-1)]
    } else {
      Cjc = c(Cjc1[1:(j-1)],Cjc1[(j+1):p])
    }
    
    sort_object = sort(abs(Cjc),index.return=T)
    newCjc = sort_object$x
    idx = sort_object$ix
    d = cumsum(rev(newCjc)) - C_d[j]
    d = rev(d)
    d_bar = d/seq(p,2,by=-1)
    m_star_temp = which((newCjc>0) & (newCjc>=d_bar))
    m_star = m_star_temp[1]
    a = rep(0,p-1)
    ao = rep(0,p-1)
    a[1:(m_star-1)] = 0
    temp = Cjc[idx] + d_bar[m_star]*((-1)^(Cjc[idx]>0))
    a[m_star:(p-1)] = temp[m_star:(p-1)]
    ao[idx] = a
    if (j==1){
      A[j,] = c(C_d[j]+d_bar[m_star],ao[j:(p-1)])
    } else if (j==p){
      A[j,] = c(ao[1:(j-1)],C_d[j]+d_bar[m_star])
    } else {
      A[j,] = c(ao[1:(j-1)],C_d[j]+d_bar[m_star],ao[j:(p-1)])
    }
  }
  return(A)
}
