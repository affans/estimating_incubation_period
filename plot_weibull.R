plotor <- function(st, ed, filename, dirname){
t <- read.csv(filename, header = T)
leave <- as.Date(as.character(t$leave)), format = '%m/%d/%y')
dura <- as.numeric(t$duration)
sympton <- as.Date(as.character(t$sympton)), format = '%m/%d/%y')

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
if(F){
dura_left = dura_left[dura_left <= 25]
maxt <- max(dura_left)
simseq = vector(mode = 'numeric', length = (maxt+1))
for(i in 0:maxt){
  simseq[i+1] = length(dura_left[dura_left == i])
}
aimseq = vector(mode = 'numeric')
for(i in 0:(maxt-1)){
  aimseq = c(aimseq, rep(i, ceiling(simseq[i+1]/2)+floor(simseq[i+2]/2)))
}
aimseq = c(aimseq, rep(maxt, ceiling(simseq[maxt+1]/2)))
dura_left = sample(aimseq, length(aimseq))
}
print("total_number is:")
print(length(dura_left))
print(summary(dura_left))
print(table(sort(dura_left)))
a = dura_left
n = length(a)
quar = c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999, 0.9999)
obj=function(par)
{
te=-(par[2]*a)^par[1]
te1=log(par[1]*par[2])-lgamma(1/par[1])
val=sum(te)+n*te1
return(-val)
}
rot=optim(c(1,1),obj)$par
print("alp = ")
alp = rot[1]
print(alp)
print("lam = ")
lam = rot[2]
print(lam)
print("estimated incubation time")
mean = 1/(lam*alp)*gamma(x = 1/alp)
print(mean)
ssquar = (1/lam)*(-log(1-quar))^(1/alp)
print(ssquar)

# plotting
xt <- seq(0, 20, 0.01)
fgamma <- function(s, la, al){
  return(exp(-(la*s)^al)*(la*al)/(gamma(x = 1/al)))
}
hxt = fgamma(xt, lam, alp)
fweibull <- function(s, la, al){
  return(al*la^(al)*s^(al-1)*exp(-(la*s)^(al)))
}
gxt = fweibull(xt, lam, alp)

ahat = as.numeric(table(sort(a)))
fre = ahat/sum(ahat)
ymax = max(c(fre, hxt, gxt))
hist_min = round(min(a))-0.5
hist_max = round(max(a))+0.5
png(paste(dirname, "histagram.png", sep = "/"))
hist(a, freq = F, breaks = hist_min:hist_max, ylim = c(0, ymax), xlab ='Forward time or incubation period (days)',
     ylab = 'Density', main = NULL)
lines(xt, hxt, lty = 1)
lines(xt, gxt, lty = 2)
legend('topright', c("Forward time", "Incubation period"), lty = 1:2, trace = TRUE, plot = TRUE)
# title('Density of Forward Time and Incubation Time')
dev.off()
}