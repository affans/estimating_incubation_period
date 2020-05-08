plotor_lognormal <- function(st, ed, filename, dirname){
  t <- read.csv(filename, header = T)
  leave <- as.Date(as.character(t$leave))
  dura <- as.numeric(t$duration)
  sympton <- as.Date(as.character(t$sympton))
  
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
  
  print("total_number is:")
  print(length(dura_left))
  print(summary(dura_left))
  print(table(sort(dura_left)))
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
  quar = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999)
  obj=function(par)
  {
    te1 <- -n*(par[1]+par[2]^2/2)
    te2 <- sum(log(1-pnorm((log(a)-par[1])/par[2])))
    return(-te1-te2)
  }
  rot=optim(c(1,1),obj)$par
  print("mu = ")
  mu = rot[1]
  print(mu)
  print("sigma = ")
  sigma = rot[2]
  print(sigma)
  print("estimated incubation time")
  mean = exp(mu+sigma^2/2)
  print(mean)
  ssquar = exp(sigma*qnorm(quar)+mu)
  print(ssquar)
  # plotting
  xt <- seq(0.01, 20, 0.01)
  fgamma <- function(s, m, si){
    return(exp(-m-si^2/2)*(1-pnorm((log(s)-m)/si)))
  }
  hxt = fgamma(xt, mu, sigma)
  fweibull <- function(s, m, si){
    return((1/(s*si*sqrt(2*pi)))*exp(-(log(s)-m)^2/(2*si^2)))
  }
  gxt = fweibull(xt, mu, sigma)
  
  ahat = as.numeric(table(sort(a)))
  fre = ahat/sum(ahat)
  ymax = max(c(fre, hxt, gxt))
  hist_min = min(a)-0.5
  hist_max = max(a)+0.5
  png(paste(dirname, "histagram.png", sep = "/"))
  hist(a, freq = F, breaks = hist_min:hist_max, ylim = c(0, ymax), xlab ='Forward time or incubation period (days)',
       ylab = 'Density', main = NULL)
  lines(xt, hxt, lty = 1)
  lines(xt, gxt, lty = 2)
  legend('topright', c("Forward time", "Incubation period"), lty = 1:2, trace = TRUE, plot = TRUE)
  # title('Density of Forward Time and Incubation Time')
  dev.off()
}