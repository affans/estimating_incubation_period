pinonparam_tr <- function(st, ed, filename, dirname){
  # revise the dataset
  non_param = matrix(0, 12, 12)
  il = 1
  mt <- matrix(0, 3, 12)
  t <- read.csv(filename, header = T)
  leave <- as.Date(as.character(t$leave), format = '%m/%d/%y')
  dura <- as.numeric(t$duration)
  sympton <- as.Date(as.character(t$sympton), format = '%m/%d/%y')
  if(is.null(st)){
    startdate = min(leave)
  }else{
    startdate = as.Date(st)
  }
  if(is.null(ed)){
    enddate = max(leave)
  }else{
    enddate = as.Date(ed)
  }
  lea_left = leave[leave %in% seq.Date(startdate, enddate, by = 'days')]
  dura_left = dura[leave %in% seq.Date(startdate, enddate, by = 'days')]
  dura_left = dura_left[dura_left >= 1]
  a = dura_left
  n = length(a)
  print(table(a))
  quar = c(0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 0.9999)
  piset = c(0, 0.05, 0.1, 0.2)
  for(pai in piset){
    obj=function(par)
    {
      te1=log(pai*(par[2]*a)^(par[1]-1)/exp(-par[2]^par[1]) + 
                (1-pai)/(gamma(x = 1/par[1])*(1-pgamma(par[2]^par[1], shape = 1/par[1], scale = 1))))-(par[2]*a)^par[1]
      te2 = log(par[1]*par[2])
      val=sum(te1)+n*te2
      return(-val)
    }
    rot=optim(c(1,1),obj)$par
    alp = rot[1]
    lam = rot[2]
    mean = 1/(lam*alp)*gamma(x = 1/alp)
    ssquar = (1/lam)*(-log(1-quar))^(1/alp)
    mt[2, ] = c(alp, lam, mean, ssquar)
    
    # nonparametric_bootstrap
    B = 1000
    matr <- matrix(0, B, 12)
    colnames(matr) <- c("alpha","lambda", "mean", paste(quar*100,"th-", "quantile", sep = ""))
    for(i in 1:B){
      att = sample(dura_left, n, replace = TRUE)
      objet=function(par)
      {
        te1=log(pai*(par[2]*att)^(par[1]-1)/exp(-par[2]^par[1]) + 
                  (1-pai)/(gamma(x = 1/par[1])*(1-pgamma(par[2]^par[1], shape = 1/par[1], scale = 1))))-(par[2]*att)^par[1]
        te2 = log(par[1]*par[2])
        val=sum(te1)+n*te2
        return(-val)
      }
      rot=optim(c(1,1),objet)$par
      alp = rot[1]
      lam = rot[2]
      mean = 1/(lam*alp)*gamma(x = 1/alp)
      ssquar = (1/lam)*(-log(1-quar))^(1/alp)
      matr[i, ] = c(alp, lam, mean, ssquar)
    }
    
    for(i in 1:12){
      st = as.numeric(matr[, i])
      st = sort(st)
      mt[1, i] = st[25]
      mt[3, i] = st[975]
    }
    row0 = c('lower', 'mean', 'upper')
    non_param[il:(il+2), ] = mt
    il = il+3
  }
  non_param = round(non_param, 2)
  colnames(non_param) <- colnames(matr)
  rownames(non_param) <- c(paste("p-0", row0, sep = ""), paste("p-0.05", row0, sep = ""), paste("p-0.1", row0, sep = ""),
                           paste("p-0.2", row0, sep = ""))
  pi_0 = vector(mode = "character", 12)
  pi_1 = vector(mode = "character", 12)
  pi_2 = vector(mode = "character", 12)
  pi_3 = vector(mode = "character", 12)
  for(j in 1:12){
    pi_0[j] = paste(non_param[2,j], "(", non_param[1, j], ", ", non_param[3, j], ")", sep = "")
    pi_1[j] = paste(non_param[5,j], "(", non_param[4, j], ", ", non_param[6, j], ")", sep = "")
    pi_2[j] = paste(non_param[8,j], "(", non_param[7, j], ", ", non_param[9, j], ")", sep = "")
    pi_3[j] = paste(non_param[11,j], "(", non_param[10, j], ", ", non_param[12, j], ")", sep = "")
  }
  param_1 = data.frame(pi_0, pi_1, pi_2, pi_3)
  colnames(param_1) = c("pi = 0", "pi = 0.05", "pi = 0.1", "pi = 0.2")
  rownames(param_1) = colnames(non_param)
  write.csv(param_1, paste(dirname, "non_param_pi.csv", sep = "/"))
}