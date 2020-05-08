t <- read.csv('final-0.csv')
leave <- as.Date(as.character(t$leave))
dura <- as.numeric(t$duration)
startdate = as.Date('2020/01/19')
enddate = as.Date('2020/1/23')
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
print(table(sort(a)))

log_logistic <- function(par){
  te1 <- n*log(sin(pi/par[2]))-n*log(par[1]*pi/par[2])
  te2 <- -sum(log(1+(a/par[1])^par[2]))
  return(-te1-te2)
}
log_normal <- function(par){
  te1 <- -n*(par[1]+par[2]^2/2)
  te2 <- sum(log(1-pnorm((log(a)-par[1])/par[2])))
  return(-te1-te2)
}
forward_gamma <- function(par){
  te=log(1- pgamma(a, shape = par[1], scale = par[2]))
  te1=-log(par[1]*par[2])
  val=sum(te)+n*te1
  return(-val)
}
# gamma
parm_gamma <- optim(c(1,1), forward_gamma)$par
print('parameter estimation of gamma')
print(parm_gamma)
alp_gamma <- parm_gamma[1]
gam_gamma <- parm_gamma[2]
f1gamma <- function(x){
  return((1-pgamma(x, shape = alp_gamma, scale = gam_gamma))/(alp_gamma*gam_gamma))
}
# log_loistic model
parm_log <- optim(c(1,1), log_logistic)$par
alp <- parm_log[1]
bet <- parm_log[2]
print('parameter estimation of log logistic')
print(parm_log)
print("estimated incubation time(of log logistic model):")
print((alp*pi/bet)/sin(pi/bet))
flogistic <- function(x){
  return(sin(pi/bet)/(alp*pi/bet*(1+(x/alp)^bet)))
}
# log_normal model
parm_nor <- optim(c(1,1), log_normal)$par
print('parameter estimation of log normal')
print(parm_nor)
mu <- parm_nor[1]
sigma <- parm_nor[2]
fnormal <- function(x){
  return((1-pnorm((log(x)-mu)/sigma))/exp(mu+sigma^2/2))
}
#weibull
weibull = function(par)
{
  te=-(par[2]*a)^par[1]
  te1=log(par[1]*par[2])-lgamma(1/par[1])
  val=sum(te)+n*te1
  return(-val)
}
parm_weibull=optim(c(1,1),weibull)$par
print('parameter estimation of weibull')
print(parm_weibull)
wei_alp = parm_weibull[1]
wei_lam = parm_weibull[2]
fweibull <- function(s){
  return(exp(-(wei_lam*s)^wei_alp)*(wei_alp*wei_lam)/gamma(x = 1/wei_alp))
}

# quartile
quar = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999)
q_logistic <- function(al, be){
  return(al*(quar/(1-quar))^(1/be))
}
print("quartile of log_logistic model:")
print(q_logistic(alp, bet))

# chi-square test
nk = 23
pj <- vector(mode = 'numeric', nk)
qj <- vector(mode = 'numeric', nk)
oj <- vector(mode = 'numeric', nk)
gj <- vector(mode = 'numeric', nk)
at <- vector(mode = 'numeric', nk)
for(i in 1:nk){
  g = i-1
  h = i
  pj[i] = integrate(flogistic, lower = g, upper = h)$value
  qj[i] = integrate(fnormal, lower = g, upper = h)$value
  oj[i] = integrate(fweibull, lower = g, upper = h)$value
  gj[i] = integrate(f1gamma, lower = g, upper = h)$value
  at[i] = length(which(a == (g+h)/2))
}
pj = pj/sum(pj)
qj = qj/sum(qj)
oj = oj/sum(oj)
gj = gj/sum(gj)

sp <- sum((n*pj-at)^2/(n*pj))
sq <- sum((n*qj-at)^2/(n*qj))
so <- sum((n*oj-at)^2/(n*oj))
sg <- sum((n*gj-at)^2/(n*gj))

print("the chisquare test(p-value) for log logistic")
print(pchisq(sp, df = nk-1, lower.tail = F))
print("the chisquare test(p-value) for log normal")
print(pchisq(sq, df = nk-1, lower.tail = F))
print("the chisquare test(p-value) for weibull")
print(pchisq(so, df = nk-1, lower.tail = F))
print("the chisquare test(p-value) for gamma")
print(pchisq(sg, df = nk-1, lower.tail = F))
rm(list=ls())
