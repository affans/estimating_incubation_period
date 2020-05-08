rm(list=ls())
sta_date <- as.Date('2020/1/19')
end_date <- as.Date("2020/1/23")
file_name <- "final-0.csv"
if(F){
  dir_name <- "gamma"
  if(!dir.exists(dir_name)){
    dir.create(dir_name)
  }
  source("plot_gamma.R")
  source("pi_gamma.R")
  plotor_gamma(st = sta_date, ed = end_date, filename = file_name, dirname = dir_name)
  pi_gamma(st = sta_date, ed = end_date, filename = file_name, dirname = dir_name)
}
if(F){
  dir_name <- "lognormal"
  if(!dir.exists(dir_name)){
    dir.create(dir_name)
  }
  source("plot_lognormal.R")
  source("pi_lognormal.R")
  plotor_lognormal(st = sta_date, ed = end_date, filename = file_name, dirname = dir_name)
  pi_lognormal(st = sta_date, ed = end_date, filename = file_name, dirname = dir_name)
}
if(T){
  dir_name <- "weibull"
  if(!dir.exists(dir_name)){
    dir.create(dir_name)
  }
  source("plot_weibull.R")
  #source("pi_weibull.R")
  plotor(st = sta_date, ed = end_date, filename = file_name, dirname = dir_name)
  #pinonparam(st = sta_date, ed = end_date, filename = file_name, dirname = dir_name)
}
# truncate at 1
if(F){
  dir_name <- "weibull_truncate"
  if(!dir.exists(dir_name)){
    dir.create(dir_name)
  }
  source("plot_truncate.R")
  source("pi_truncate.R")
  plotor_tr(st = sta_date, ed = end_date, filename = file_name, dirname = dir_name)
  pinonparam_tr(st = sta_date, ed = end_date, filename = file_name, dirname = dir_name)
}
if(F){
  dir_name <- "weibull-new"
  if(!dir.exists(dir_name)){
    dir.create(dir_name)
  }
  source("newplot.R")
  source("newpi.R")
  newplotor(st = sta_date, ed = end_date, filename = file_name, dirname = dir_name)
  newpi(st = sta_date, ed = end_date, filename = file_name, dirname = dir_name)
}
rm(list=ls())