# this is the code to run all 6 clinical data: pbc, veteran, lung, anz, us, melanoma clinical
.libPaths("/dskh/nobackup/yunwei")
setwd("/dskh/nobackup/yunwei/Analysis_yw/benchmark/202012summaryresult")

library(dplyr)
library(survival)
library(glmnet)
library(rms)
library(tidyverse)
library(caret)
library(pec)
library(coefplot)
library("survAUC")
library(gridExtra)
library(ggplot2)
library("survival")
library(survminer)
library(randomForestSRC)
library(ggRandomForests)
library(penalized)
library(DMwR)
library(randomForest)
library(riskRegression)
library(pROC)
library(ROCR)
library(cvTools)
library(parallel)
library(pbmcapply)
library(MTLR)
library(profmem)
library(keras)
#install_keras()
library(pseudo)
library(survivalROC)
library(survival)
#library(survcomp)#cran package moved to bioconductor
library(survAUC)
library(CoxBoost)
library(limma)
library(partykit)
library(coin)
library(compound.Cox)
library(GenAlgo)
library(survivalsvm)
library(rmatio)

source("original_cox_fun.R")
source("p_cox1_fun.R")
source("p_cox2_fun.R")
source("p_cox3_fun.R")
source("rsf1_fun.R")
source("rsf2_fun.R")
source("bw_cox1_fun.R")
source("bw_cox2_fun.R")
source("bw_cox3_fun.R")
source("mtlr_fun.R")
source("dnnsurv_fun.R")
source("getpseudoconditional.R")
source("cox_boost_fun.R")
source("ga_original_cox_fun.R")
source("ga_mtlr_fun.R")
source("ga_cox_boost_fun.R")
source("mtlr_fun2.R")
source("cox_boost_fun2.R")
source("survivalsvm_fun.R")
source("glmnet_lassocox_fun.R")
source("glmnet_ridge_fun.R")
source("glmnet_en_fun.R")

#pbc
data("pbc",package = "randomForestSRC")
dim(pbc)
colnames(pbc)
pbc1 <- pbc[,c("days", "status", "treatment", "sex", "age", "bili", "stage")]
dim(pbc1)
colnames(pbc1)[1]="time"
any(is.na(pbc1))
pbc2=na.omit(pbc1)
dim(pbc2) #312 7
table(pbc2$status) #187 125
fitform_ogl=as.formula(Surv(time, status) ~ treatment + sex + age + bili + stage)
summary(pbc2$time)
hist(pbc2$time)
#1
#original_cox_fun(1,pbc2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess)
start_time <- Sys.time()
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(pbc2$time)[2]),as.numeric(summary(pbc2$time)[5]),(as.numeric(summary(pbc2$time)[5])-as.numeric(summary(pbc2$time)[2]))/14)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, original_cox_fun,pbc2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"pbc_cox1m.rds")
saveRDS(cox1,"pbc_cox1.rds")
cox1<-readRDS("pbc_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"pbc_cox1t.rds")
#2
env <- globalenv() # Grab the global environment
env$dd <- datadist(pbc2) # Assign the datadist to it
#bw_cox1_fun(1,pbc2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox1_fun,pbc2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"pbc_bw_cox1m.rds")
saveRDS(cox1,"pbc_bw_cox1.rds")
cox1<-readRDS("pbc_bw_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"pbc_bw_cox1t.rds")
#3
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox2_fun,pbc2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"pbc_bw_cox2m.rds")
saveRDS(cox1,"pbc_bw_cox2.rds")
cox1<-readRDS("pbc_bw_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"pbc_bw_cox2t.rds")
#4
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox3_fun,pbc2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"pbc_bw_cox3m.rds")
saveRDS(cox1,"pbc_bw_cox3.rds")
cox1<-readRDS("pbc_bw_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"pbc_bw_cox3t.rds")
#5
start_time <- Sys.time()
form1=as.formula(~treatment + sex + age + bili + stage)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox1_fun,pbc2,5,form1,formula1,formula2,formula3,formula4,timess,1,0, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"pbc_p_cox1m.rds")
saveRDS(cox1,"pbc_p_cox1.rds")
cox1<-readRDS("pbc_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"pbc_p_cox1t.rds")
#6
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox2_fun,pbc2,5,form1,formula1,formula2,formula3,formula4,timess, 0,1,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"pbc_p_cox2m.rds")
saveRDS(cox1,"pbc_p_cox2.rds")
cox1<-readRDS("pbc_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"pbc_p_cox2t.rds")
#7
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox3_fun,pbc2,5,form1,formula1,formula2,formula3,formula4,timess,1,1,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"pbc_p_cox3m.rds")
saveRDS(cox1,"pbc_p_cox3.rds")
cox1<-readRDS("pbc_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"pbc_p_cox3t.rds")
#8
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,pbc2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"pbc_rsf1m.rds")
saveRDS(cox1,"pbc_rsf1.rds")
cox1<-readRDS("pbc_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"pbc_rsf1t.rds")
#9
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,pbc2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"pbc_rsf2m.rds")
saveRDS(cox1,"pbc_rsf2.rds")
cox1<-readRDS("pbc_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"pbc_rsf2t.rds")
#10
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,pbc2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"pbc_mtlrm.rds")
saveRDS(cox1,"pbc_mtlr.rds")
cox1<-readRDS("pbc_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"pbc_mtlrt.rds")
#11
# formula1=fitform_ogl
# formula2=fitform_ogl
# formula3=Surv(time,status)~1
# formula4=Surv(time,status)~1
# timess=seq(1,900,900/15)
pickedtime=c(50,100)
#t<-surv_train
#d<-cen_train
#qt=pickedtime
#getPseudoConditional(surv_train, cen_train, pickedtime)
#dnnsurv_fun(1,pbc2,2,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, dnnsurv_fun,pbc2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"pbc_dnnsurvm.rds")
saveRDS(cox1,"pbc_dnnsurv.rds")
cox1<-readRDS("pbc_dnnsurv.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"pbc_dnnsurvt.rds")

#12
start_time <- Sys.time()
time1=timess[3]
stepnumber=10
penaltynumber=100
#data=veteran
#cvK=5
#j=1
#cox_boost_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,pbc2,5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"pbc_coxboostm.rds")
saveRDS(cox1,"pbc_coxboost.rds")
cox1<-readRDS("pbc_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"pbc_coxboostt.rds")
#20
pbc2$time=as.numeric(pbc2$time)
#survivalsvm_fun(1,pbc2,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,pbc2,5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"pbc_survivalsvmm.rds")
saveRDS(cox1,"pbc_survivalsvm.rds")
cox1<-readRDS("pbc_survivalsvm.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"pbc_survivalsvmt.rds")
###############################################################




###############################################################
# veteran data
data("veteran")
colnames(veteran)
xnam <- paste(colnames(veteran)[c(1,2,5,6,7,8)], sep="")
form=as.formula(paste("Surv(time, status)~ ", paste(xnam, collapse= "+")))
fitform_ogl=form
summary(veteran$time)
hist(veteran$time)
env <- globalenv() # Grab the global environment
env$dd <- datadist(veteran) # Assign the datadist to it
#1
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(veteran$time)[2]),as.numeric(summary(veteran$time)[5]),(as.numeric(summary(veteran$time)[5])-as.numeric(summary(veteran$time)[2]))/14)

start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, original_cox_fun,veteran,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"veteran_cox1m.rds")
saveRDS(cox1,"veteran_cox1.rds")
cox1<-readRDS("veteran_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"veteran_cox1t.rds")
#2
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox1_fun,veteran,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"veteran_bw_cox1m.rds")
saveRDS(cox1,"veteran_bw_cox1.rds")
cox1<-readRDS("veteran_bw_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"veteran_bw_cox1t.rds")
#3
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox2_fun,veteran,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"veteran_bw_cox2m.rds")
saveRDS(cox1,"veteran_bw_cox2.rds")
cox1<-readRDS("veteran_bw_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"veteran_bw_cox2t.rds")
#4
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox3_fun,veteran,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"veteran_bw_cox3m.rds")
saveRDS(cox1,"veteran_bw_cox3.rds")
cox1<-readRDS("veteran_bw_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"veteran_bw_cox3t.rds")
#5
form1=as.formula(~trt + celltype + karno + diagtime + age + prior)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox1_fun,veteran,5,form1,formula1,formula2,formula3,formula4,timess,1,0,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"veteran_p_cox1m.rds")
saveRDS(cox1,"veteran_p_cox1.rds")
cox1<-readRDS("veteran_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"veteran_p_cox1t.rds")
#6
form1=as.formula(~trt + celltype + karno + diagtime + age + prior)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox2_fun,veteran,5,form1,formula1,formula2,formula3,formula4,timess,0,1, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"veteran_p_cox2m.rds")
saveRDS(cox1,"veteran_p_cox2.rds")
cox1<-readRDS("veteran_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"veteran_p_cox2t.rds")
#7
form1=as.formula(~trt + celltype + karno + diagtime + age + prior)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox3_fun,veteran,5,form1,formula1,formula2,formula3,formula4,timess,1,1,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"veteran_p_cox3m.rds")
saveRDS(cox1,"veteran_p_cox3.rds")
cox1<-readRDS("veteran_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"veteran_p_cox3t.rds")
#8
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,veteran,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"veteran_rsf1m.rds")
saveRDS(cox1,"veteran_rsf1.rds")
cox1<-readRDS("veteran_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"veteran_rsf1t.rds")
#9
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,veteran,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"veteran_rsf2m.rds")
saveRDS(cox1,"veteran_rsf2.rds")
cox1<-readRDS("veteran_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"veteran_rsf2t.rds")
#10
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,veteran,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"veteran_mtlrm.rds")
saveRDS(cox1,"veteran_mtlr.rds")
cox1<-readRDS("veteran_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"veteran_mtlrt.rds")
#11
start_time <- Sys.time()
pickedtime=c(2,3,5,7,9)
#t<-surv_train
#d<-cen_train
#qt=pickedtime
#getPseudoConditional(surv_train, cen_train, c(2,3,5,7,9))

#dnnsurv_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, dnnsurv_fun,veteran,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"veteran_dnnsurvm.rds")
saveRDS(cox1,"veteran_dnnsurv.rds")
cox1<-readRDS("veteran_dnnsurv.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"veteran_dnnsurvt.rds")
#12
start_time <- Sys.time()
time1=timess[3]
stepnumber=10
penaltynumber=100
#data=veteran
#cvK=5
#j=1
#cox_boost_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,veteran,5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"veteran_coxboostm.rds")
saveRDS(cox1,"veteran_coxboost.rds")
cox1<-readRDS("veteran_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"veteran_coxboostt.rds")

#20
veteran$time=as.numeric(veteran$time)
#survivalsvm_fun(1,veteran,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,veteran,5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"veteran_survivalsvmm.rds")
saveRDS(cox1,"veteran_survivalsvm.rds")
cox1<-readRDS("veteran_survivalsvm.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"veteran_survivalsvmt.rds")
###############################################################




###############################################################
#lung
lung=read.csv("lung.csv")
colnames(lung)
dim(lung)
head(lung)
any(is.na(lung))
library(mice)
imputation1 <- mice(lung,m=1,maxit=10,seed=500)
imputation1$method
lung1=mice::complete(imputation1,1)
xnam <- paste(colnames(lung1)[4:10], sep="")
form=as.formula(paste("Surv(time, status)~ ", paste(xnam, collapse= "+")))
fitform_ogl=form
summary(lung1$time)
hist(lung1$time)
lung1$status=as.vector(ifelse(lung1$status==1,0,1))
table(lung1$status)
env <- globalenv() # Grab the global environment
env$dd <- datadist(lung1) # Assign the datadist to it
#1
start_time <- Sys.time()
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(lung1$time)[2]),as.numeric(summary(lung1$time)[5]),(as.numeric(summary(lung1$time)[5])-as.numeric(summary(lung1$time)[2]))/14)

Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, original_cox_fun,lung1,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"lung_cox1m.rds")
saveRDS(cox1,"lung_cox1.rds")
cox1<-readRDS("lung_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lung_cox1t.rds")
#2
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox1_fun,lung1,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"lung_bw_cox1m.rds")
saveRDS(cox1,"lung_bw_cox1.rds")
cox1<-readRDS("lung_bw_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lung_bw_cox1t.rds")
#3
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox2_fun,lung1,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"lung_bw_cox2m.rds")
saveRDS(cox1,"lung_bw_cox2.rds")
cox1<-readRDS("lung_bw_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lung_bw_cox2t.rds")
#4
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox3_fun,lung1,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"lung_bw_cox3m.rds")
saveRDS(cox1,"lung_bw_cox3.rds")
cox1<-readRDS("lung_bw_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lung_bw_cox3t.rds")
#5
start_time <- Sys.time()
form1=as.formula(paste("~", paste(xnam, collapse= "+")))
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox1_fun,lung1,5,form1,formula1,formula2,formula3,formula4,timess, 1,0,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"lung_p_cox1m.rds")
saveRDS(cox1,"lung_p_cox1.rds")
cox1<-readRDS("lung_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lung_p_cox1t.rds")
#6
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox2_fun,lung1,5,form1,formula1,formula2,formula3,formula4,timess,0,1,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"lung_p_cox2m.rds")
saveRDS(cox1,"lung_p_cox2.rds")
cox1<-readRDS("lung_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lung_p_cox2t.rds")
#7
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox3_fun,lung1,5,form1,formula1,formula2,formula3,formula4,timess,1,1, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"lung_p_cox3m.rds")
saveRDS(cox1,"lung_p_cox3.rds")
cox1<-readRDS("lung_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lung_p_cox3t.rds")
#8
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,lung1,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"lung_rsf1m.rds")
saveRDS(cox1,"lung_rsf1.rds")
cox1<-readRDS("lung_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lung_rsf1t.rds")
#9
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,lung1,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"lung_rsf2m.rds")
saveRDS(cox1,"lung_rsf2.rds")
cox1<-readRDS("lung_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lung_rsf2t.rds")
#10
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,lung1,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"lung_mtlrm.rds")
saveRDS(cox1,"lung_mtlr.rds")
cox1<-readRDS("lung_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lung_mtlrt.rds")
#11
# formula1=fitform_ogl
# formula2=fitform_ogl
# formula3=Surv(time,status)~1
# formula4=Surv(time,status)~1
# timess=seq(1,900,900/15)
pickedtime=c(50,60,70,88,100)
#t<-surv_train
#d<-cen_train
#qt=pickedtime
#getPseudoConditional(surv_train, cen_train, c(2,3,5,7,9))
#dnnsurv_fun(1,lung1,2,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, dnnsurv_fun,lung1,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"lung_dnnsurvm.rds")
saveRDS(cox1,"lung_dnnsurv.rds")
cox1<-readRDS("lung_dnnsurv.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lung_dnnsurvt.rds")

#12
start_time <- Sys.time()
time1=timess[3]
stepnumber=10
penaltynumber=100
#data=veteran
#cvK=5
#j=1
#cox_boost_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,lung1,5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"lung_coxboostm.rds")
saveRDS(cox1,"lung_coxboost.rds")
cox1<-readRDS("lung_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lung_coxboostt.rds")

#20
lung1$time=as.numeric(lung1$time)
#survivalsvm_fun(1,lung1,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,lung1,5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"lung_survivalsvmm.rds")
saveRDS(cox1,"lung_survivalsvm.rds")
cox1<-readRDS("lung_survivalsvm.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lung_survivalsvmt.rds")
###############################################################




###############################################################

#anz
rawdata=readRDS("data_2017.rds")
#donor matrix
dmatrix=select(rawdata, starts_with("donor_"))
dim(dmatrix)
dfeature=c("donor_age","donor_sex","donor_height" ,"donor_weight","donor_state" ,"donor_blgroup","donor_eth","donor_diabetes" ,"donor_hypertension" ,"donor_smoker" ,"donor_death", "donor_creatinine","donor_dcd","donor_kdri")
length(dfeature)
dmatrix1=dmatrix[,dfeature]
dim(dmatrix1)

#recip matrix
rmatrix=select(rawdata, starts_with("recip_"))
dim(rmatrix)
rfeature=c("recip_age_rrtstart","recip_age", "recip_sex","recip_height" ,"recip_weight", "recip_eth_detailed" ,"recip_pra","recip_blgroup","recip_biopsy" ,"recip_creatinine", "recip_smoker","recip_lung","recip_coronary" ,"recip_pvd","recip_cvd","recip_diabetes","recip_cancer","recip_waittime","recip_id","recip_liststartdate","recip_rrtstartdate","recip_state_initial","recip_state_current","recip_epts","recip_listtime","recip_graftno","recip_state")
length(rfeature)
rmatrix1=rmatrix[,rfeature]
dim(rmatrix1)


txmatrix=select(rawdata, starts_with("tx_"))
dim(txmatrix)
colnames(txmatrix)=c("tx_id","tx_date","tx_ischemia" ,"tx_faildate" ,"tx_gstatus"  ,"tx_gperiod" ,  "tx_lstatus" , "tx_lperiod" , "tx_misa" ,    "tx_misb" ,   
                     "tx_misdr"  ,  "tx_gfailcat", "tx_organ")
txfeature=c("tx_ischemia" , "tx_gstatus", "tx_gperiod" ,"tx_lstatus" ,"tx_lperiod" ,"tx_misa","tx_misb", "tx_misdr","tx_date")
length(txfeature)
txmatrix1=txmatrix[,txfeature]
dim(txmatrix1)

# levels(txmatrix1$tx_graftfailcause) # this one need to check

#combine them and delete na, check how many rows we have
rawdata2=cbind.data.frame(dmatrix1,rmatrix1,txmatrix1)

check_na=as.data.frame(sapply(rawdata2, function(x) sum(is.na(x))))
check_na$name=rownames(check_na)
colnames(check_na)[1]="count"
check_na1=check_na%>% filter(count>0)
#check_na1

rawdata3=na.omit(rawdata2)
dim(rawdata3) # wehave 3325 rows left, good
any(is.na(rawdata3))

# check the other leftover columns
#colnames(rawdata)[!colnames(rawdata)%in% c(dfeature,rfeature,txfeature)]


rawdata3$donor_diabetes=as.character(rawdata3$donor_diabetes)
rawdata3$donor_diabetes=as.factor(ifelse(rawdata3$donor_diabetes=="No"|rawdata3$donor_diabetes=="No Diabetes","N","Y"))
levels(rawdata3$donor_diabetes)
table(rawdata3$donor_eth)
rawdata3$donor_eth=as.factor(ifelse(rawdata3$donor_eth=="Asian","Asian",ifelse(rawdata3$donor_eth=="Caucasian","Caucasian","other")))
levels(rawdata3$donor_eth)
rawdata3$donor_death=as.character(rawdata3$donor_death)
rawdata3$donor_causedeath_cva=as.factor(ifelse(rawdata3$donor_death=="Cerebral Infarct"|rawdata3$donor_death=="Spontaneous Subarachnoid Haemorrhage","Y","N"))
levels(rawdata3$donor_causedeath_cva)

rawdata3$recip_eth_detailed=as.factor(rawdata3$recip_eth_detailed)
levels(rawdata3$recip_eth_detailed)[c(10,11)]<-"Indigenous" #Oceanian - Australian, Oceanian - Australian Aboriginal
levels(rawdata3$recip_eth_detailed)[c(1:9,11:29)]<-"Non-Indigenous"
levels(rawdata3$recip_eth)

rawdata3$recip_diabetes=as.character(rawdata3$recip_diabetes)
rawdata3$recip_diabetes=as.factor(ifelse(rawdata3$recip_diabetes=="N","N","Y"))
levels(rawdata3$recip_diabetes)


rawdata3$recip_cardio=as.factor(ifelse(rawdata3$recip_coronary=="Y"|rawdata3$recip_pvd=="Y"|rawdata3$recip_cvd=="Y","Y","N"))
levels(rawdata3$recip_cardio)
rawdata3$recip_bmi=rawdata3$recip_weight/(rawdata3$recip_height/100)^2
rawdata3$donor_bmi=rawdata3$donor_weight/(rawdata3$donor_height/100)^2
rawdata3=rawdata3[!rawdata3$recip_smoker=="Unknown",]#delete this

rawdata3=mutate_if(rawdata3, is.character, as.factor)
rawdata3$recip_id=as.character(rawdata3$recip_id)

rawdata3$recip_liststart=as.factor(ifelse(rawdata3$recip_liststartdate<"2009-01-01",1,ifelse(rawdata3$recip_liststartdate>"2012-01-01",3,2)))
table(rawdata3$recip_liststart)
rawdata3$recip_waitstart=as.factor(ifelse(rawdata3$recip_rrtstartdate<"2009-01-01",1,ifelse(rawdata3$recip_rrtstartdate>"2012-01-01",3,2)))
table(rawdata3$recip_waitstart)

#colnames(rawdata3)
#str(rawdata3)

selected_var=c("donor_age","donor_sex","donor_bmi", "donor_state","donor_blgroup" ,"donor_eth","donor_diabetes" ,"donor_hypertension","donor_causedeath_cva" ,"donor_smoker" , "donor_creatinine","donor_dcd","donor_kdri","recip_age", "recip_sex","recip_bmi","recip_cardio","recip_eth_detailed" , "recip_pra" ,"recip_blgroup","recip_biopsy"  ,      
               "recip_creatinine"    , "recip_smoker" , "recip_lung"  ,"recip_diabetes"       ,"recip_cancer","recip_waittime","recip_listtime", "recip_liststart","recip_waitstart","recip_state_initial","recip_state_current","recip_epts","recip_graftno" ,"tx_ischemia" ,  "tx_misa","tx_misb", "tx_misdr","tx_gstatus", "tx_gperiod")
rawdata4=rawdata3[,selected_var]
colnames(rawdata4)[18]=c("recip_eth")
colnames(rawdata4)[c(39,40)]=c("status","time")
#colnames(rawdata4)
rawdata4[rawdata4$time==0,"time"]<-0.1
#change all factors to numerical to run certain methods
fac2=c()
for (n in colnames(rawdata4))
  if (is.factor(rawdata4[[n]])) {
    fac2=append(fac2,n)
  }
rawdata4_new=rawdata4 %>% mutate_at(fac2, funs(as.numeric))



fitform_ogl<-Surv(time,status)~donor_age+donor_sex+donor_bmi+donor_blgroup+donor_eth+donor_diabetes+donor_hypertension+donor_causedeath_cva+donor_smoker+donor_creatinine+donor_dcd+donor_kdri+recip_age+recip_sex+recip_bmi+recip_cardio+recip_eth+recip_pra+recip_blgroup+recip_biopsy+recip_creatinine+recip_smoker+recip_lung+recip_diabetes+recip_cancer+recip_waittime+recip_listtime+recip_liststart+recip_waitstart+recip_epts+recip_graftno+tx_ischemia+tx_misa+tx_misb+tx_misdr
summary(rawdata4$time)
env <- globalenv() # Grab the global environment
env$dd <- datadist(rawdata4) # Assign the datadist to it
hist(rawdata4$time)
#1
start_time <- Sys.time()
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(rawdata4$time)[2]),as.numeric(summary(rawdata4$time)[5]),(as.numeric(summary(rawdata4$time)[5])-as.numeric(summary(rawdata4$time)[2]))/14)


#original_cox_fun(1,clinical4,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, original_cox_fun,rawdata4,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"anz_cox1m.rds")
saveRDS(cox1,"anz_cox1.rds")
cox1<-readRDS("anz_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"anz_cox1t.rds")
#2
start_time <- Sys.time()
#bw_cox1_fun(26,rawdata4,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox1_fun,rawdata4,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"anz_bw_cox1m.rds")
saveRDS(cox1,"anz_bw_cox1.rds")
cox1<-readRDS("anz_bw_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"anz_bw_cox1t.rds")
#3
start_time <- Sys.time()
#bw_cox2_fun(56,rawdata4,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox2_fun,rawdata4,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"anz_bw_cox2m.rds")
saveRDS(cox1,"anz_bw_cox2.rds")
cox1<-readRDS("anz_bw_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"anz_bw_cox2t.rds")
#4
start_time <- Sys.time()
#bw_cox2_fun(56,rawdata4,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox3_fun,rawdata4,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"anz_bw_cox3m.rds")
saveRDS(cox1,"anz_bw_cox3.rds")
cox1<-readRDS("anz_bw_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"anz_bw_cox3t.rds")
#5
start_time <- Sys.time()
form1=as.formula(~donor_age+donor_sex+donor_bmi+donor_blgroup+donor_eth+donor_diabetes+donor_hypertension+donor_causedeath_cva+donor_smoker+donor_creatinine+donor_dcd+donor_kdri+recip_age+recip_sex+recip_bmi+recip_cardio+recip_eth+recip_pra+recip_blgroup+recip_biopsy+recip_creatinine+recip_smoker+recip_lung+recip_diabetes+recip_cancer+recip_waittime+recip_listtime+recip_liststart+recip_waitstart+recip_epts+recip_graftno+tx_ischemia+tx_misa+tx_misb+tx_misdr)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox1_fun,rawdata4,5,form1,formula1,formula2,formula3,formula4,timess,1,0,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"anz_p_cox1m.rds")
saveRDS(cox1,"anz_p_cox1.rds")
cox1<-readRDS("anz_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"anz_p_cox1t.rds")
#6
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox2_fun,rawdata4,5,form1,formula1,formula2,formula3,formula4,timess,0,1,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"anz_p_cox2m.rds")
saveRDS(cox1,"anz_p_cox2.rds")
cox1<-readRDS("anz_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"anz_p_cox2t.rds")
#7
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox3_fun,rawdata4,5,form1,formula1,formula2,formula3,formula4,timess,1,1,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"anz_p_cox3m.rds")
saveRDS(cox1,"anz_p_cox3.rds")
cox1<-readRDS("anz_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"anz_p_cox3t.rds")
#8
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,rawdata4,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"anz_rsf1m.rds")
saveRDS(cox1,"anz_rsf1.rds")
cox1<-readRDS("anz_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"anz_rsf1t.rds")
#9
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,rawdata4,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"anz_rsf2m.rds")
saveRDS(cox1,"anz_rsf2.rds")
cox1<-readRDS("anz_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"anz_rsf2t.rds")
#10
start_time <- Sys.time()
#mtlr_fun(1,clinical4_new,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,rawdata4,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"anz_mtlrm.rds")
saveRDS(cox1,"anz_mtlr.rds")
cox1<-readRDS("anz_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"anz_mtlrt.rds")
#11: 
start_time <- Sys.time()
pickedtime=c(100,270,350,560,1000)
#data<-clinical4_new
#cvK=5
#j=1
#t<-surv_train
#d<-cen_train
#qt=pickedtime
#getPseudoConditional(surv_train, cen_train, 45)
#getPseudoConditional(surv_train, cen_train, pickedtime)
#dnnsurv_fun(1,rawdata4_new,5,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, dnnsurv_fun,rawdata4_new,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"anz_dnnsurvm.rds")
saveRDS(cox1,"anz_dnnsurv.rds")
cox1<-readRDS("anz_dnnsurv.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"anz_dnnsurvt.rds")

#12
start_time <- Sys.time()
time1=timess[3]
stepnumber=10
penaltynumber=100
#data=veteran
#cvK=5
#j=1
#cox_boost_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,rawdata4,5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"anz_coxboostm.rds")
saveRDS(cox1,"anz_coxboost.rds")
cox1<-readRDS("anz_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"anz_coxboostt.rds")


###############################################################




###############################################################
#us
rawdata=readRDS("full_dt3.rds")
dim(rawdata)
any(is.na(rawdata))
str(rawdata)
#the response columns are tx_gperiod and tx_gstatus
#get rid of factor with 1 level
rawdata1=rawdata[,-which(colnames(rawdata)%in% c("donor_type","recip_data")) ]
dim(rawdata1)
#get rid of the date time variables: none (because this data comes form the imputation without those date objects, in the previous numero analysis, the date columns are obtained from the original data)
#sapply(rawdata1, class)
colnames(rawdata1 %>% dplyr::select_if(is.character))
colnames(rawdata1 %>% dplyr::select_if(is.factor))
#colnames(rawdata1 %>% dplyr::select_if(is.integer))
colnames(rawdata1 %>% dplyr::select_if(is.numeric)) #numeric includes those integer
table(rawdata1$tx_gstatus)
summary(rawdata1$tx_gperiod)
class(rawdata1$tx_gstatus)
rawdata1[rawdata1$tx_gperiod==0, "tx_gperiod"]=0.1
#the columns can be used for prediction, with those similar ones
predictors_name=c("recip_blg","recip_age" ,"recip_bmi","recip_state","recip_cmv_igg" ,"recip_cmv_igm","recip_cmv_status","recip_creatinine_trans","recip_kidney_diagnosis","recip_diabetes","recip_kidney_diagnosis_trans","recip_pretrans_dialysis" ,"recip_epstin_barr_virus_status" ,"recip_hightest_education" , "recip_pra_current"  ,"recip_eth" ,"recip_sex","recip_hepatitis" ,"recip_hepb_surface","recip_hepc","recip_height","recip_height_registration" ,"recip_hiv_tx","recip_age_listing", "recip_bmi_listing","recip_pra_listing", "recip_height_listing", "recip_weight_listing","recip_malignancy" ,"recip_malignancy_tx","recip_medical_condition","recip_kidney_num" ,"recip_kidney_num2","recip_svd","recip_previous_tx", "recip_previous_tx2","recip_serum_albumn",  "donor_blg" ,"donor_age","donor_antiht","donor_arginine_vasopressin","donor_infection_blood_source","donor_bmi" ,"donor_terminal_blood" , "donor_cancer_site","donor_cardiac_arrest", "donor_per_phs","donor_state","donor_clinical_infection","donor_cause_of_death" ,"donor_creatinine" ,"donor_circumstance_of_death","donor_mechanism_of_death","donor_diabetes2","donor_epstin_barr_virus_igg","donor_epstin_barr_virus_igm","donor_ecd","donor_kidney_status" ,"donor_eth","donor_cancer" ,"donor_sex","donor_hbsab","donor_hbv","donor_hepb_surface","donor_hepc","donor_height","donor_cancer2","donor_smoker","donor_cocaine","donor_diabetes3","donor_ht3","donor_other_drug","donor_cancer3","donor_kdri1","donor_biopsy","donor_baby","donor_non_heart_beating.1","donor_infection","donor_protein","donor_infection_pulmonary","donor_right_kidney_biopsy", "donor_terminal_sgotast","donor_terminal_sgotalt","donor_skin_cancer","donor_form", "donor_tatoos" ,"donor_terminal_bilirubin" , "tx_blg_match","tx_misa", "tx_misb","tx_ischemia","tx_allocation_days","tx_waittime","tx_misdr", "tx_kidney_received","tx_kidney_pump" ,"tx_kidney_received_ice","tx_graphical_type"              
                  ,"tx_procedure_type")
#paste(predictors_name, collapse=" + ")
rawdata2=rawdata1[,c(predictors_name,"tx_gperiod","tx_gstatus")]

##########################################################################
#change all factors to numerical to run certain methods
fac2=c()
for (n in colnames(rawdata2))
  if (is.factor(rawdata2[[n]])) {
    fac2=append(fac2,n)
  }
rawdata2_new=rawdata2 %>% mutate_at(fac2, funs(as.numeric))

##########################################################################

#change the time and status columns
colnames(rawdata2)[which(colnames(rawdata2)=="tx_gperiod")]="time"
colnames(rawdata2)[which(colnames(rawdata2)=="tx_gstatus")]="status"
#the overall formula
fitform_ogl<-Surv(time,status)~recip_blg + recip_age + recip_bmi + recip_state + recip_cmv_igg + recip_cmv_igm + recip_cmv_status + recip_creatinine_trans + recip_kidney_diagnosis + recip_diabetes + recip_kidney_diagnosis_trans + recip_pretrans_dialysis + recip_epstin_barr_virus_status + recip_hightest_education + recip_pra_current + recip_eth + recip_sex + recip_hepatitis + recip_hepb_surface + recip_hepc + recip_height + recip_height_registration + recip_hiv_tx + recip_age_listing + recip_bmi_listing + recip_pra_listing + recip_height_listing + recip_weight_listing + recip_malignancy + recip_malignancy_tx + recip_medical_condition + recip_kidney_num + recip_kidney_num2 + recip_svd + recip_previous_tx + recip_previous_tx2 + recip_serum_albumn + donor_blg + donor_age + donor_antiht + donor_arginine_vasopressin + donor_infection_blood_source + donor_bmi + donor_terminal_blood + donor_cancer_site + donor_cardiac_arrest + donor_per_phs + donor_state + donor_clinical_infection + donor_cause_of_death + donor_creatinine + donor_circumstance_of_death + donor_mechanism_of_death + donor_diabetes2 + donor_epstin_barr_virus_igg + donor_epstin_barr_virus_igm + donor_ecd + donor_kidney_status + donor_eth + donor_cancer + donor_sex + donor_hbsab + donor_hbv + donor_hepb_surface + donor_hepc + donor_height + donor_cancer2 + donor_smoker + donor_cocaine + donor_diabetes3 + donor_ht3 + donor_other_drug + donor_cancer3 + donor_kdri1 + donor_biopsy + donor_baby + donor_non_heart_beating.1 + donor_infection + donor_protein + donor_infection_pulmonary + donor_right_kidney_biopsy + donor_terminal_sgotast + donor_terminal_sgotalt + donor_skin_cancer + donor_form + donor_tatoos + donor_terminal_bilirubin + tx_blg_match + tx_misa + tx_misb + tx_ischemia + tx_allocation_days + tx_waittime + tx_misdr + tx_kidney_received + tx_kidney_pump + tx_kidney_received_ice + tx_graphical_type + tx_procedure_type

#fitform_ogl<-Surv(time,status)~recip_blg + recip_age + recip_bmi + recip_state + recip_cmv_igg + recip_cmv_igm + recip_cmv_status + recip_creatinine_trans + recip_kidney_diagnosis + recip_diabetes + recip_kidney_diagnosis_trans + recip_pretrans_dialysis + recip_epstin_barr_virus_status + recip_hightest_education + recip_pra_current + recip_eth + recip_sex + recip_hepatitis + recip_hepb_surface + recip_hepc + recip_height + recip_height_registration + recip_hiv_tx + recip_age_listing + recip_bmi_listing + recip_pra_listing + recip_height_listing + recip_weight_listing + recip_malignancy + recip_malignancy_tx + recip_medical_condition + recip_kidney_num + recip_kidney_num2 + recip_svd + recip_previous_tx + recip_previous_tx2 + recip_serum_albumn + donor_blg + donor_age + donor_antiht + donor_arginine_vasopressin + donor_infection_blood_source + donor_bmi + donor_terminal_blood + donor_cancer_site + donor_cardiac_arrest + donor_per_phs + donor_state + donor_clinical_infection + donor_cause_of_death + donor_creatinine + donor_circumstance_of_death + donor_mechanism_of_death + donor_diabetes2  + donor_ecd + donor_kidney_status + donor_eth + donor_cancer + donor_sex   + donor_height + donor_cancer2 + donor_smoker + donor_cocaine + donor_diabetes3 + donor_ht3 + donor_other_drug + donor_cancer3 + donor_kdri1 + donor_biopsy  + donor_non_heart_beating.1 + donor_infection + donor_protein + donor_infection_pulmonary + donor_right_kidney_biopsy + donor_terminal_sgotast + donor_terminal_sgotalt + donor_skin_cancer + donor_form + donor_tatoos + donor_terminal_bilirubin + tx_blg_match + tx_misa + tx_misb + tx_ischemia + tx_allocation_days + tx_waittime + tx_misdr + tx_kidney_received + tx_kidney_pump + tx_kidney_received_ice + tx_graphical_type + tx_procedure_type
# donor_epstin_barr_virus_igg ,donor_epstin_barr_virus_igm ,donor_hbsab,donor_hbv,donor_hepb_surface,donor_hepc  ,donor_baby
#downsampling
set.seed(202001)
rawdata3=rawdata2[sample(nrow(rawdata2),3000),]

summary(rawdata3$time)
hist(rawdata3$time)
fac2=c()
for (n in colnames(rawdata3))
  if (is.factor(rawdata3[[n]])) {
    fac2=append(fac2,n)
  }
rawdata3_new=rawdata3 %>% mutate_at(fac2, funs(as.numeric))

##########################################################################################
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(rawdata3$time)[2]),as.numeric(summary(rawdata3$time)[5]),(as.numeric(summary(rawdata3$time)[5])-as.numeric(summary(rawdata3$time)[2]))/14)
env <- globalenv() # Grab the global environment
env$dd <- datadist(rawdata3) # Assign the datadist to it
# data=rawdata3
# cvK=5
# j=1
#1
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, original_cox_fun,rawdata3,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"us_cox1m.rds")
saveRDS(cox1,"us_cox1.rds")
cox1<-readRDS("us_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"us_cox1t.rds")
#5
form1=as.formula(~recip_blg + recip_age + recip_bmi + recip_state + recip_cmv_igg + recip_cmv_igm + recip_cmv_status + recip_creatinine_trans + recip_kidney_diagnosis + recip_diabetes + recip_kidney_diagnosis_trans + recip_pretrans_dialysis + recip_epstin_barr_virus_status + recip_hightest_education + recip_pra_current + recip_eth + recip_sex + recip_hepatitis + recip_hepb_surface + recip_hepc + recip_height + recip_height_registration + recip_hiv_tx + recip_age_listing + recip_bmi_listing + recip_pra_listing + recip_height_listing + recip_weight_listing + recip_malignancy + recip_malignancy_tx + recip_medical_condition + recip_kidney_num + recip_kidney_num2 + recip_svd + recip_previous_tx + recip_previous_tx2 + recip_serum_albumn + donor_blg + donor_age + donor_antiht + donor_arginine_vasopressin + donor_infection_blood_source + donor_bmi + donor_terminal_blood + donor_cancer_site + donor_cardiac_arrest + donor_per_phs + donor_state + donor_clinical_infection + donor_cause_of_death + donor_creatinine + donor_circumstance_of_death + donor_mechanism_of_death + donor_diabetes2 + donor_epstin_barr_virus_igg + donor_epstin_barr_virus_igm + donor_ecd + donor_kidney_status + donor_eth + donor_cancer + donor_sex + donor_hbsab + donor_hbv + donor_hepb_surface + donor_hepc + donor_height + donor_cancer2 + donor_smoker + donor_cocaine + donor_diabetes3 + donor_ht3 + donor_other_drug + donor_cancer3 + donor_kdri1 + donor_biopsy + donor_baby + donor_non_heart_beating.1 + donor_infection + donor_protein + donor_infection_pulmonary + donor_right_kidney_biopsy + donor_terminal_sgotast + donor_terminal_sgotalt + donor_skin_cancer + donor_form + donor_tatoos + donor_terminal_bilirubin + tx_blg_match + tx_misa + tx_misb + tx_ischemia + tx_allocation_days + tx_waittime + tx_misdr + tx_kidney_received + tx_kidney_pump + tx_kidney_received_ice + tx_graphical_type + tx_procedure_type)
start_time <- Sys.time()
#p_cox1_fun(1,rawdata3,2,form1,formula1,formula2,formula3,formula4,timess,1,0)

Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox1_fun,rawdata3,5,form1,formula1,formula2,formula3,formula4,timess,1,0,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"us_p_cox1m.rds")
saveRDS(cox1,"us_p_cox1.rds")
cox1<-readRDS("us_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"us_p_cox1t.rds")
#6
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox2_fun,rawdata3,5, form1,formula1,formula2,formula3,formula4,timess,0,1,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"us_p_cox2m.rds")
saveRDS(cox1,"us_p_cox2.rds")
cox1<-readRDS("us_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"us_p_cox2t.rds")
#7
#p_cox3_fun(1,rawdata3,2,form1,formula1,formula2,formula3,formula4,timess,1,1)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox3_fun,rawdata3,5, form1,formula1,formula2,formula3,formula4,timess,1,1,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"us_p_cox3m.rds")
saveRDS(cox1,"us_p_cox3.rds")
cox1<-readRDS("us_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"us_p_cox3t.rds")
#8
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,rawdata3,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"us_rsf1m.rds")
saveRDS(cox1,"us_rsf1.rds")
cox1<-readRDS("us_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"us_rsf1t.rds")
#9
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,rawdata3,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"us_rsf2m.rds")
saveRDS(cox1,"us_rsf2.rds")
cox1<-readRDS("us_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"us_rsf2t.rds")
#10
start_time <- Sys.time()
#mtlr_fun(1,rawdata3,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,rawdata3,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"us_mtlrm.rds")
saveRDS(cox1,"us_mtlr.rds")
cox1<-readRDS("us_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"us_mtlrt.rds")
# #11: failed ,dont know why and fail to get the results
# start_time <- Sys.time()
# #pickedtime=c(45,50,100,150,200)
# pickedtime=c(45,100)
# #pickedtime=c(10,45)
# #data<-rawdata3_new
# #cvK=5
# #j=1
# #t<-surv_train
# #d<-cen_train
# #qt=pickedtime
# #getPseudoConditional(surv_train, cen_train, 45)
# #getPseudoConditional(surv_train, cen_train, pickedtime)
# #dnnsurv_fun(1,rawdata3_new,2,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
# #Error in names(x) <- value : 'names' attribute [26] must be the same length as the vector [1]
# Rprof(tf <- "rprof.log",memory.profiling=TRUE)
# cox1 <- pbmcapply::pbmclapply(21:40, dnnsurv_fun,rawdata3_new,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 15)
# Rprof(NULL)
# mm<-summaryRprof(tf,memory = "both")
# mm
# saveRDS(mm,"us_dnnsurvm.rds")
# saveRDS(cox1,"us_dnnsurv.rds")
# cox1<-readRDS("us_dnnsurv.rds")
# head(cox1)
# end_time <- Sys.time()
# saveRDS(end_time - start_time,"us_dnnsurvt.rds")
#12
start_time <- Sys.time()
time1=timess[3]
stepnumber=10
penaltynumber=100
#data=veteran
#cvK=5
#j=1
#cox_boost_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,rawdata3,5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"us_coxboostm.rds")
saveRDS(cox1,"us_coxboost.rds")
cox1<-readRDS("us_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"us_coxboostt.rds")

###############################################################




###############################################################
#melanoma_clinical
load("melanomaAssays.rdata")
class(melanomaAssays)
assays(melanomaAssays)

clinical1=as.data.frame(colData(melanomaAssays))
colnames(clinical1)[which(colnames(clinical1)=="Person_FUStatus")]="status"
colnames(clinical1)[which(colnames(clinical1)=="Prognosis_TimeSinceLNMet")]="time"
table(clinical1$status)
clinical1$status=as.vector(ifelse(clinical1$status=="Dead, melanoma",1,0))
table(clinical1$status) #39 60

predictors_name=c("Tumour_ID","Person_NumPrim","Person_Sex","Person_StageatBank","Age_Analysis","Tum_NumNodesInv","Tum_MetSize","Tum_Extranodal","Tum_CellType","Tum_CellSize","Tum_Necrosis","Tum_Pigment",
                  "Tum_NonTumour.","Tum_BRAFmut","Tum_NRASmut","Tum_FLT3mut","Tum_METmut","Tum_BRAFmut", "Tum_NRASmut","Tum_FLT3mut" , "Tum_METmut","Tum_PIK3CAmut", "Tum_PDGFRAmut" ,"Tum_EGFRmut","Tum_CKITmut","Prim_Site", "Prim_Site_SunExp","Prim_Stage" , "Prim_Tstage" , "Prim_NStage" , "Prim_Naevus" ,"Prim_Breslow" , "Prim_Mitos", "Prim_Clark","Prim_Histol","Prim_Regress" ,"Prim_Ulc" ,"Prim_Vasc","Prim_LymphInv","Prim_Satell","Sun_Damage_Score")

clinical1=clinical1[,colnames(clinical1)%in% c(predictors_name,"time","status")]
#check for na
check_na_fun=function(data){
  check_na=as.data.frame(sapply(data, function(x) sum(is.na(x))))
  check_na$name=rownames(check_na)
  colnames(check_na)[1]="count"
  check_na1=check_na%>% filter(count>0)
  return(check_na1)}
check_na_fun(clinical1)

#delete : Sun_Damage_Score
clinical1=clinical1[,!colnames(clinical1)%in% c("Sun_Damage_Score")]
clinical1$Tum_CellSize #factor??? difference between 1 and 1m??

#change the numericl values to numeric
clinical1$Tum_Necrosis=as.numeric(sub("%", "", clinical1$Tum_Necrosis))
clinical1$Tum_NonTumour.=as.numeric(sub("%", "", clinical1$Tum_NonTumour.))
#clinical1$Tum_CellType=as.numeric(clinical1$Tum_CellType)

#change character to factor
clinical1 <- mutate_if(clinical1, is.character, as.factor)
#str(clinical1)

#imputation: failed
#imputation1 <- mice(clinical1,m=1,maxit=10,seed=500) #Error in solve.default(xtx + diag(pen)) : system is computationally singular: reciprocal condition number = 2.16639e-16

clinical2=clinical1[complete.cases(clinical1),] #64 records left

predictors_name_new=c("Tumour_ID","Person_NumPrim","Person_Sex","Person_StageatBank","Age_Analysis","Tum_NumNodesInv","Tum_MetSize","Tum_Extranodal","Tum_CellType","Tum_CellSize","Tum_Necrosis","Tum_Pigment",
                      "Tum_NonTumour.","Tum_BRAFmut","Tum_NRASmut","Tum_FLT3mut","Tum_METmut","Tum_BRAFmut", "Tum_NRASmut","Tum_FLT3mut" , "Tum_METmut","Tum_PIK3CAmut", "Tum_PDGFRAmut" ,"Tum_EGFRmut","Tum_CKITmut","Prim_Site", "Prim_Site_SunExp","Prim_Stage" , "Prim_Tstage" , "Prim_NStage" , "Prim_Naevus")
clinical3=clinical1[,colnames(clinical1)%in% c(predictors_name_new,"time","status")]
clinical3=clinical3[complete.cases(clinical3),] #88 records left

#no need those prim columns
clinical4=clinical3[,1:23] #88*23

##########################################################################
#change all factors to numerical to run certain methods
fac2=c()
for (n in colnames(clinical4))
  if (is.factor(clinical4[[n]])) {
    fac2=append(fac2,n)
  }
clinical4_new=clinical4 %>% mutate_at(fac2, funs(as.numeric))

##########################################################################
fitform_ogl=Surv(time,status)~Person_NumPrim+Person_Sex+Person_StageatBank+Age_Analysis+Tum_NumNodesInv+Tum_MetSize+Tum_Extranodal+Tum_CellType+Tum_CellSize+Tum_Necrosis+Tum_Pigment+Tum_NonTumour.+Tum_BRAFmut+Tum_NRASmut+Tum_FLT3mut+Tum_METmut+Tum_BRAFmut+Tum_NRASmut+Tum_FLT3mut+Tum_METmut+Tum_PIK3CAmut+Tum_PDGFRAmut+Tum_EGFRmut+Tum_CKITmut
summary(clinical4$time)
hist(clinical4$time)
#1
start_time <- Sys.time()
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(clinical4$time)[2]),as.numeric(summary(clinical4$time)[5]),(as.numeric(summary(clinical4$time)[5])-as.numeric(summary(clinical4$time)[2]))/14)

#original_cox_fun(1,clinical4,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, original_cox_fun,clinical4,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,mc.cores = 5)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaclinical_cox1m.rds")
saveRDS(cox1,"melanomaclinical_cox1.rds")
cox1<-readRDS("melanomaclinical_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaclinical_cox1t.rds")
#5
start_time <- Sys.time()
form1=as.formula(~Person_NumPrim+Person_Sex+Person_StageatBank+Age_Analysis+Tum_NumNodesInv+Tum_MetSize+Tum_Extranodal+Tum_CellType+Tum_CellSize+Tum_Necrosis+Tum_Pigment+Tum_NonTumour.+Tum_BRAFmut+Tum_NRASmut+Tum_FLT3mut+Tum_METmut+Tum_BRAFmut+Tum_NRASmut+Tum_FLT3mut+Tum_METmut+Tum_PIK3CAmut+Tum_PDGFRAmut+Tum_EGFRmut+Tum_CKITmut)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox1_fun,clinical4,5, form1,formula1,formula2,formula3,formula4,timess,1,0, mc.cores = 5)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaclinical_p_cox1m.rds")
saveRDS(cox1,"melanomaclinical_p_cox1.rds")
cox1<-readRDS("melanomaclinical_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaclinical_p_cox1t.rds")
#6
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox2_fun,clinical4,5, form1,formula1,formula2,formula3,formula4,timess,0,1, mc.cores = 5)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaclinical_p_cox2m.rds")
saveRDS(cox1,"melanomaclinical_p_cox2.rds")
cox1<-readRDS("melanomaclinical_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaclinical_p_cox2t.rds")
#7
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, p_cox3_fun,clinical4,5, form1,formula1,formula2,formula3,formula4,timess,1,1, mc.cores = 5)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaclinical_p_cox3m.rds")
saveRDS(cox1,"melanomaclinical_p_cox3.rds")
cox1<-readRDS("melanomaclinical_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaclinical_p_cox3t.rds")
#8
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,clinical4,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 5)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaclinical_rsf1m.rds")
saveRDS(cox1,"melanomaclinical_rsf1.rds")
cox1<-readRDS("melanomaclinical_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaclinical_rsf1t.rds")
#9
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,clinical4,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 5)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaclinical_rsf2m.rds")
saveRDS(cox1,"melanomaclinical_rsf2.rds")
cox1<-readRDS("melanomaclinical_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaclinical_rsf2t.rds")
#10
start_time <- Sys.time()
#mtlr_fun(1,clinical4_new,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,clinical4,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 5)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaclinical_mtlrm.rds")
saveRDS(cox1,"melanomaclinical_mtlr.rds")
cox1<-readRDS("melanomaclinical_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaclinical_mtlrt.rds")
# #11: falied not at the pesudo step
# # env <- globalenv() # Grab the global environment
# # env$dd <- datadist(clinical4_new) # Assign the datadist to it
# # new<-clinical4_new[,!names(clinical4_new)%in% c("Tum_CKITmut")]
# start_time <- Sys.time()
# pickedtime=c(30,50,100,150,200)
# #data<-clinical4_new
# #cvK=5
# #j=1
# #t<-surv_train
# #d<-cen_train
# #qt=pickedtime
# #getPseudoConditional(surv_train, cen_train, 45)
# #getPseudoConditional(surv_train, cen_train, pickedtime)
# #dnnsurv_fun(10,new,5,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
# Rprof(tf <- "rprof.log",memory.profiling=TRUE)
# cox1 <- pbmcapply::pbmclapply(1:20, dnnsurv_fun,clinical4_new,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 5)
# Rprof(NULL)
# mm<-summaryRprof(tf,memory = "both")
# mm
# saveRDS(mm,"melanomaclinical_dnnsurvm.rds")
# saveRDS(cox1,"melanomaclinical_dnnsurv.rds")
# cox1<-readRDS("melanomaclinical_dnnsurv.rds")
# head(cox1)
# end_time <- Sys.time()
# saveRDS(end_time - start_time,"melanomaclinical_dnnsurvt.rds")

#12
start_time <- Sys.time()
time1=timess[3]
stepnumber=10
penaltynumber=100
#data=veteran
#cvK=5
#j=1
#cox_boost_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,clinical4,5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 5)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaclinica_coxboostm.rds")
saveRDS(cox1,"melanomaclinica_coxboost.rds")
cox1<-readRDS("melanomaclinica_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaclinica_coxboostt.rds")
#20
clinical4$time=as.numeric(clinical4$time)
#survivalsvm_fun(1,clinical4,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,clinical4,5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"clinical_survivalsvmm.rds")
saveRDS(cox1,"clinical_survivalsvm.rds")
cox1<-readRDS("clinical_survivalsvm.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"clinical_survivalsvmt.rds")






