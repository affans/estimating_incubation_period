pi_lognormal <- function(st, ed, filename, dirname){
  # revise the dataset
  non_param = matrix(0, 12, 11)
  il = 1
  mt <- matrix(0, 3, 11)
  t <- read.csv(filename, header = T)
  leave <- as.Date(as.character(t$leave))
  dura <- as.numeric(t$duration)
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
  
  a = dura_left
  n = length(a)
  x1 = unique(sort(a))
  x2 = as.integer(table(sort(a)))
  for(j in 2:length(x2)){
    x2[j-1] = x2[j-1] + ceiling(x2[j]/2)
    x2[j] = x2[j] - ceiling(x2[j]/2)
  }
  b = vector(mode = 'numeric')
  for(j in 1:length(x1)){
    b = c(b, rep(x1[j], x2[j]))
  }
  a = sample(b, length(b))
  a = a + 0.5
  quar = c(0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999)
  piset = c(0, 0.05, 0.1, 0.2)
  for(pai in piset){
    obj=function(par)
    {
      te1 <- -n*(par[1]+par[2]^2/2)
      te2 <- sum(log(1-pnorm((log(a)-par[1])/par[2])))
      return(-te1-te2)
    }
    rot=optim(c(1,1),obj)$par
    mu = rot[1]
    sigma = rot[2]
    mean = mean = exp(mu+sigma^2/2)
    ssquar = exp(sigma*qnorm(quar)+mu)
    mt[2, ] = c(mu, sigma, mean, ssquar)
    
    # nonparametric_bootstrap
    B = 1000
    matr <- matrix(0, B, 11)
    colnames(matr) <- c("mu","sigma", "mean", paste(quar*100,"th-", "quantile", sep = ""))
    for(i in 1:B){
      att = sample(dura_left, n, replace = TRUE)
      objet=function(par)
      {
        te1 <- -n*(par[1]+par[2]^2/2)
        te2 <- sum(log(1-pnorm((log(att)-par[1])/par[2])))
        return(-te1-te2)
      }
      rot=optim(c(1,1),objet)$par
      mu = rot[1]
      sigma = rot[2]
      mean = mean = exp(mu+sigma^2/2)
      ssquar = exp(sigma*qnorm(quar)+mu)
      matr[i, ] = c(mu, sigma, mean, ssquar)
    }
    
    for(i in 1:11){
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
  pi_0 = vector(mode = "character", 11)
  pi_1 = vector(mode = "character", 11)
  pi_2 = vector(mode = "character", 11)
  pi_3 = vector(mode = "character", 11)
  for(j in 1:11){
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