# this is the code to run 6 published gene expression data: ng1,ng2,ng3,ng4,ng5,ng6
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

###############################################################
trail1=read.mat("sorliedataorig.mat")
dim(trail1$X)
table(trail1$delta)[1]/(dim(trail1$X)[1])
trail2=read.mat("bovelstaddata.mat")
dim(trail2$X)
table(trail2$delta)[1]/(dim(trail2$X)[1])
trail3=read.mat("Lung-Beer2002-86-x-7129.mat")
dim(trail3$X)
table(trail3$delta)[1]/(dim(trail3$X)[1])
trail4=read.mat("AML-Bullinger2004-116-x-6283.mat")
dim(trail4$X)
table(trail4$delta)[1]/(dim(trail4$X)[1])
trail5=read.mat("Veer2002-78-x-4751.mat")
dim(trail5$X)
table(trail5$delta)[1]/(dim(trail5$X)[1])
trail6=read.mat("rosenwalddataorig.mat")
dim(trail6$X)
table(trail6$delta)[1]/(dim(trail6$X)[1])


###############################################################
trail1_1=cbind.data.frame(trail1$X,trail1$delta,trail1$Y)
dim(trail1_1)
colnames(trail1_1)[(dim(trail1_1)[2]-1):dim(trail1_1)[2]]=c("status","time")
fitform_ogl=as.formula(Surv(time,status)~.)
form1=as.formula(~.)
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(trail1_1$time)[2]),as.numeric(summary(trail1_1$time)[5]),(as.numeric(summary(trail1_1$time)[5])-as.numeric(summary(trail1_1$time)[2]))/14)
hist(trail1_1$time)
trail1_1$os_class=as.vector(ifelse(trail1_1$status==1 & 
                                     trail1_1$time <50, "poor", 
                                   ifelse(trail1_1$status==0 & 
                                            trail1_1$time >50, "good", "not")))
table(trail1_1$os_class)
trail1_2=trail1_1[,-which(names(trail1_1) %in% c("os_class"))]
# data=trail1_1
# cvK=5
# j=1
# numm=(dim(trail1_1)[2]-3)
# topnumm=10

#19
#glmnet_lassocox_fun(1,trail1_2,5,formula1,formula2,formula3,formula4,timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_lassocox_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_p_cox1m.rds")
saveRDS(cox1,"ngse1_p_cox1.rds")
cox1<-readRDS("ngse1_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_p_cox1t.rds")
#20
#p_cox2_fun(11,trail1_1,5, form1,formula1,formula2,formula3,formula4,timess,0,50)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_ridge_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_p_cox2m.rds")
saveRDS(cox1,"ngse1_p_cox2.rds")
cox1<-readRDS("ngse1_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_p_cox2t.rds")
#21
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_en_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_p_cox3m.rds")
saveRDS(cox1,"ngse1_p_cox3.rds")
cox1<-readRDS("ngse1_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_p_cox3t.rds")

#8
# so strange why this is not working now
# data("veteran")
# fitform_ogl=as.formula(Surv(time,status)~.-trt)
# coxph(fitform_ogl,veteran)
# rfsrc(fitform_ogl, data = veteran,ntree = 1000,mtry = 10,tree.err=TRUE,importance = TRUE)

start_time <- Sys.time()
#rsf1_fun(1,trail1_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_rsf1m.rds")
saveRDS(cox1,"ngse1_rsf1.rds")
cox1<-readRDS("ngse1_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_rsf1t.rds")
#9
start_time <- Sys.time()
#rsf2_fun(1,trail1_1,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_rsf2m.rds")
saveRDS(cox1,"ngse1_rsf2.rds")
cox1<-readRDS("ngse1_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_rsf2t.rds")
#10
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_mtlrm.rds")
saveRDS(cox1,"ngse1_mtlr.rds")
cox1<-readRDS("ngse1_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_mtlrt.rds")
#11
start_time <- Sys.time()
pickedtime=c(5,7,12)
#t<-surv_train
#d<-cen_train
#qt=pickedtime
#getPseudoConditional(surv_train, cen_train, c(2,3,5,7,9))

#dnnsurv_fun(1,trail1_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, dnnsurv_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_dnnsurvm.rds")
saveRDS(cox1,"ngse1_dnnsurv.rds")
cox1<-readRDS("ngse1_dnnsurv.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_dnnsurvt.rds")
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
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_coxboostm.rds")
saveRDS(cox1,"ngse1_coxboost.rds")
cox1<-readRDS("ngse1_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_coxboostt.rds")

#13
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_original_cox_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_ga_cox1m.rds")
saveRDS(cox1,"ngse1_ga_cox1.rds")
cox1<-readRDS("ngse1_ga_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_ga_cox1t.rds")
#14
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_mtlr_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_ga_mtlrm.rds")
saveRDS(cox1,"ngse1_ga_mtlr.rds")
cox1<-readRDS("ngse1_ga_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_ga_mtlrt.rds")
#15
time1=timess[3]
stepnumber=10
penaltynumber=100
#ga_cox_boost_fun(1,trail1_1,5,(dim(trail1_1)[2]-3),10,20,time1,timess,stepnumber,penaltynumber)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_cox_boost_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_ga_coxboostm.rds")
saveRDS(cox1,"ngse1_ga_cocboost.rds")
cox1<-readRDS("ngse1_ga_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_ga_coxboostt.rds")

#16
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun2,trail1_1,5, (dim(trail1_1)[2]-3),10,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_limma_mtlrm.rds")
saveRDS(cox1,"ngse1_limma_mtlr.rds")
cox1<-readRDS("ngse1_limma_mtlr.rds")
colMeans(do.call(rbind, cox1))
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_limma_mtlrt.rds")
#17
time1=timess[3]
stepnumber=10
penaltynumber=100
#cox_boost_fun2(1,trail1_1,5,(dim(trail1_1)[2]-3),10,time1,timess,stepnumber,penaltynumber)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun2,trail1_1,5, (dim(trail1_1)[2]-3),10,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_limma_coxboostm.rds")
saveRDS(cox1,"ngse1_limma_coxboost.rds")
cox1<-readRDS("ngse1_limma_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_limma_coxboostt.rds")

#18
trail1_1$time=as.numeric(trail1_1$time)
#survivalsvm_fun(1,veteran,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse1_survivalsvmm.rds")
saveRDS(cox1,"ngse1_survivalsvm.rds")
cox1<-readRDS("ngse1_survivalsvm.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse1_survivalsvmt.rds")


###############################################################




###############################################################
trail1_1=cbind.data.frame(trail2$X,trail2$delta,trail2$Y)
dim(trail1_1)
colnames(trail1_1)[(dim(trail1_1)[2]-1):dim(trail1_1)[2]]=c("status","time")
fitform_ogl=as.formula(Surv(time,status)~.)
form1=as.formula(~.)
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(trail1_1$time)[2]),as.numeric(summary(trail1_1$time)[5]),(as.numeric(summary(trail1_1$time)[5])-as.numeric(summary(trail1_1$time)[2]))/14)
hist(trail1_1$time)
abline(v=c(as.numeric(summary(trail1_1$time)[2]),as.numeric(summary(trail1_1$time)[5])), col="red")
trail1_1$os_class=as.vector(ifelse(trail1_1$status==1 & 
                                     trail1_1$time <10, "poor", 
                                   ifelse(trail1_1$status==0 & 
                                            trail1_1$time >10, "good", "not")))
table(trail1_1$os_class)
trail1_2=trail1_1[,-which(names(trail1_1) %in% c("os_class"))]
# data=trail1_1
# cvK=5
# j=1
# numm=(dim(trail1_1)[2]-3)
# topnumm=10

#19
#glmnet_lassocox_fun(1,trail1_2,5,formula1,formula2,formula3,formula4,timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_lassocox_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_p_cox1m.rds")
saveRDS(cox1,"ngse2_p_cox1.rds")
cox1<-readRDS("ngse2_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_p_cox1t.rds")
#20
#p_cox2_fun(11,trail1_1,5, form1,formula1,formula2,formula3,formula4,timess,0,50)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_ridge_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_p_cox2m.rds")
saveRDS(cox1,"ngse2_p_cox2.rds")
cox1<-readRDS("ngse2_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_p_cox2t.rds")
#21
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_en_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_p_cox3m.rds")
saveRDS(cox1,"ngse2_p_cox3.rds")
cox1<-readRDS("ngse2_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_p_cox3t.rds")


#8
# so strange why this is not working now
# data("veteran")
# fitform_ogl=as.formula(Surv(time,status)~.-trt)
# coxph(fitform_ogl,veteran)
# rfsrc(fitform_ogl, data = veteran,ntree = 1000,mtry = 10,tree.err=TRUE,importance = TRUE)

start_time <- Sys.time()
#rsf1_fun(1,trail1_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_rsf1m.rds")
saveRDS(cox1,"ngse2_rsf1.rds")
cox1<-readRDS("ngse2_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_rsf1t.rds")
#9
start_time <- Sys.time()
#rsf2_fun(1,trail1_1,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_rsf2m.rds")
saveRDS(cox1,"ngse2_rsf2.rds")
cox1<-readRDS("ngse2_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_rsf2t.rds")
#10
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_mtlrm.rds")
saveRDS(cox1,"ngse2_mtlr.rds")
cox1<-readRDS("ngse2_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_mtlrt.rds")
#11
start_time <- Sys.time()
pickedtime=timess[c(1,3,5)]
#t<-surv_train
#d<-cen_train
#qt=pickedtime
#getPseudoConditional(surv_train, cen_train, c(2,3,5,7,9))

#dnnsurv_fun(1,trail1_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, dnnsurv_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_dnnsurvm.rds")
saveRDS(cox1,"ngse2_dnnsurv.rds")
cox1<-readRDS("ngse2_dnnsurv.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_dnnsurvt.rds")
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
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_coxboostm.rds")
saveRDS(cox1,"ngse2_coxboost.rds")
cox1<-readRDS("ngse2_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_coxboostt.rds")

#13
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_original_cox_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_ga_cox1m.rds")
saveRDS(cox1,"ngse2_ga_cox1.rds")
cox1<-readRDS("ngse2_ga_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_ga_cox1t.rds")
#14
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_mtlr_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_ga_mtlrm.rds")
saveRDS(cox1,"ngse2_ga_mtlr.rds")
cox1<-readRDS("ngse2_ga_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_ga_mtlrt.rds")
#15
time1=timess[3]
stepnumber=10
penaltynumber=100
#ga_cox_boost_fun(1,trail1_1,5,(dim(trail1_1)[2]-3),10,20,time1,timess,stepnumber,penaltynumber)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_cox_boost_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_ga_coxboostm.rds")
saveRDS(cox1,"ngse2_ga_cocboost.rds")
cox1<-readRDS("ngse2_ga_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_ga_coxboostt.rds")

#16

start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun2,trail1_1,5, (dim(trail1_1)[2]-3),10,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_limma_mtlrm.rds")
saveRDS(cox1,"ngse2_limma_mtlr.rds")
cox1<-readRDS("ngse2_limma_mtlr.rds")
colMeans(do.call(rbind, cox1))
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_limma_mtlrt.rds")
#17
time1=timess[3]
stepnumber=10
penaltynumber=100
#cox_boost_fun2(1,trail1_1,5,(dim(trail1_1)[2]-3),10,time1,timess,stepnumber,penaltynumber)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun2,trail1_1,5, (dim(trail1_1)[2]-3),10,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_limma_coxboostm.rds")
saveRDS(cox1,"ngse2_limma_coxboost.rds")
cox1<-readRDS("ngse2_limma_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_limma_coxboostt.rds")

#18
trail1_1$time=as.numeric(trail1_1$time)
#survivalsvm_fun(1,veteran,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse2_survivalsvmm.rds")
saveRDS(cox1,"ngse2_survivalsvm.rds")
cox1<-readRDS("ngse2_survivalsvm.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse2_survivalsvmt.rds")


###############################################################



###############################################################
trail1_1=cbind.data.frame(trail3$X,trail3$delta,trail3$Y)
dim(trail1_1)
colnames(trail1_1)[(dim(trail1_1)[2]-1):dim(trail1_1)[2]]=c("status","time")
fitform_ogl=as.formula(Surv(time,status)~.)
form1=as.formula(~.)
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(trail1_1$time)[2]),as.numeric(summary(trail1_1$time)[5]),(as.numeric(summary(trail1_1$time)[5])-as.numeric(summary(trail1_1$time)[2]))/14)
hist(trail1_1$time)
abline(v=c(as.numeric(summary(trail1_1$time)[2]),as.numeric(summary(trail1_1$time)[5])), col="red")
trail1_1$os_class=as.vector(ifelse(trail1_1$status==1 & 
                                     trail1_1$time <40, "poor", 
                                   ifelse(trail1_1$status==0 & 
                                            trail1_1$time >40, "good", "not")))
table(trail1_1$os_class)
trail1_2=trail1_1[,-which(names(trail1_1) %in% c("os_class"))]
# data=trail1_1
# cvK=5
# j=1
# numm=(dim(trail1_1)[2]-3)
# topnumm=10
#19
#glmnet_lassocox_fun(1,trail1_2,5,formula1,formula2,formula3,formula4,timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_lassocox_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse3_p_cox1m.rds")
saveRDS(cox1,"ngse3_p_cox1.rds")
cox1<-readRDS("ngse3_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse3_p_cox1t.rds")
#20
#p_cox2_fun(11,trail1_1,5, form1,formula1,formula2,formula3,formula4,timess,0,50)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_ridge_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse3_p_cox2m.rds")
saveRDS(cox1,"ngse3_p_cox2.rds")
cox1<-readRDS("ngse3_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse3_p_cox2t.rds")
#21
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_en_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse3_p_cox3m.rds")
saveRDS(cox1,"ngse3_p_cox3.rds")
cox1<-readRDS("ngse3_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse3_p_cox3t.rds")


#8
# so strange why this is not working now
# data("veteran")
# fitform_ogl=as.formula(Surv(time,status)~.-trt)
# coxph(fitform_ogl,veteran)
# rfsrc(fitform_ogl, data = veteran,ntree = 1000,mtry = 10,tree.err=TRUE,importance = TRUE)

start_time <- Sys.time()
#rsf1_fun(1,trail1_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse3_rsf1m.rds")
saveRDS(cox1,"ngse3_rsf1.rds")
cox1<-readRDS("ngse3_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse3_rsf1t.rds")
#9
start_time <- Sys.time()
#rsf2_fun(1,trail1_1,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse3_rsf2m.rds")
saveRDS(cox1,"ngse3_rsf2.rds")
cox1<-readRDS("ngse3_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse3_rsf2t.rds")
#10
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse3_mtlrm.rds")
saveRDS(cox1,"ngse3_mtlr.rds")
cox1<-readRDS("ngse3_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse3_mtlrt.rds")
#11
start_time <- Sys.time()
pickedtime=timess[c(1,3,5)]
#t<-surv_train
#d<-cen_train
#qt=pickedtime
#getPseudoConditional(surv_train, cen_train, c(2,3,5,7,9))

#dnnsurv_fun(1,trail1_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, dnnsurv_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse3_dnnsurvm.rds")
saveRDS(cox1,"ngse3_dnnsurv.rds")
cox1<-readRDS("ngse3_dnnsurv.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse3_dnnsurvt.rds")
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
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse3_coxboostm.rds")
saveRDS(cox1,"ngse3_coxboost.rds")
cox1<-readRDS("ngse3_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse3_coxboostt.rds")

#13
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_original_cox_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse3_ga_cox1m.rds")
saveRDS(cox1,"ngse3_ga_cox1.rds")
cox1<-readRDS("ngse3_ga_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse3_ga_cox1t.rds")
#14
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_mtlr_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse3_ga_mtlrm.rds")
saveRDS(cox1,"ngse3_ga_mtlr.rds")
cox1<-readRDS("ngse3_ga_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse3_ga_mtlrt.rds")
#15
time1=timess[3]
stepnumber=10
penaltynumber=100
#ga_cox_boost_fun(1,trail1_1,5,(dim(trail1_1)[2]-3),10,20,time1,timess,stepnumber,penaltynumber)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_cox_boost_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse3_ga_coxboostm.rds")
saveRDS(cox1,"ngse3_ga_cocboost.rds")
cox1<-readRDS("ngse3_ga_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse3_ga_coxboostt.rds")

#16

start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun2,trail1_1,5, (dim(trail1_1)[2]-3),10,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse3_limma_mtlrm.rds")
saveRDS(cox1,"ngse3_limma_mtlr.rds")
cox1<-readRDS("ngse3_limma_mtlr.rds")
colMeans(do.call(rbind, cox1))
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse3_limma_mtlrt.rds")
#17
time1=timess[3]
stepnumber=10
penaltynumber=100
#cox_boost_fun2(1,trail1_1,5,(dim(trail1_1)[2]-3),10,time1,timess,stepnumber,penaltynumber)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun2,trail1_1,5, (dim(trail1_1)[2]-3),10,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse3_limma_coxboostm.rds")
saveRDS(cox1,"ngse3_limma_coxboost.rds")
cox1<-readRDS("ngse3_limma_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse3_limma_coxboostt.rds")

# #18
# # Error in quadprog::solve.QP(C, -d, t(H), f, meq = meq) : 
# #   constraints are inconsistent, no solution!
# trail1_2$time=as.numeric(trail1_2$time)
# #survivalsvm_fun(1,trail1_2,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
# start_time <- Sys.time()
# Rprof(tf <- "rprof.log",memory.profiling=TRUE)
# cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
# Rprof(NULL)
# mm<-summaryRprof(tf,memory = "both")
# mm
# saveRDS(mm,"ngse3_survivalsvmm.rds")
# saveRDS(cox1,"ngse3_survivalsvm.rds")
# cox1<-readRDS("ngse3_survivalsvm.rds")
# head(cox1)
# end_time <- Sys.time()
# saveRDS(end_time - start_time,"ngse3_survivalsvmt.rds")
# 


###############################################################





###############################################################
trail1_1=cbind.data.frame(trail4$X,trail4$delta,trail4$Y)
dim(trail1_1)
colnames(trail1_1)[(dim(trail1_1)[2]-1):dim(trail1_1)[2]]=c("status","time")
fitform_ogl=as.formula(Surv(time,status)~.)
form1=as.formula(~.)
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(trail1_1$time)[2]),as.numeric(summary(trail1_1$time)[5]),(as.numeric(summary(trail1_1$time)[5])-as.numeric(summary(trail1_1$time)[2]))/14)
hist(trail1_1$time)
abline(v=c(as.numeric(summary(trail1_1$time)[2]),as.numeric(summary(trail1_1$time)[5])), col="red")

trail1_1$os_class=as.vector(ifelse(trail1_1$status==1 & 
                                     trail1_1$time <500, "poor", 
                                   ifelse(trail1_1$status==0 & 
                                            trail1_1$time >500, "good", "not")))
table(trail1_1$os_class)
trail1_2=trail1_1[,-which(names(trail1_1) %in% c("os_class"))]
# data=trail1_1
# cvK=5
# j=1
# numm=(dim(trail1_1)[2]-3)
# topnumm=10
#19
#glmnet_lassocox_fun(1,trail1_2,5,formula1,formula2,formula3,formula4,timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_lassocox_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_p_cox1m.rds")
saveRDS(cox1,"ngse4_p_cox1.rds")
cox1<-readRDS("ngse4_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_p_cox1t.rds")
#20
#p_cox2_fun(11,trail1_1,5, form1,formula1,formula2,formula3,formula4,timess,0,50)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_ridge_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_p_cox2m.rds")
saveRDS(cox1,"ngse4_p_cox2.rds")
cox1<-readRDS("ngse4_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_p_cox2t.rds")
#21
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_en_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_p_cox3m.rds")
saveRDS(cox1,"ngse4_p_cox3.rds")
cox1<-readRDS("ngse4_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_p_cox3t.rds")


#8
# so strange why this is not working now
# data("veteran")
# fitform_ogl=as.formula(Surv(time,status)~.-trt)
# coxph(fitform_ogl,veteran)
# rfsrc(fitform_ogl, data = veteran,ntree = 1000,mtry = 10,tree.err=TRUE,importance = TRUE)

start_time <- Sys.time()
#rsf1_fun(1,trail1_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_rsf1m.rds")
saveRDS(cox1,"ngse4_rsf1.rds")
cox1<-readRDS("ngse4_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_rsf1t.rds")
#9
start_time <- Sys.time()
#rsf2_fun(1,trail1_1,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_rsf2m.rds")
saveRDS(cox1,"ngse4_rsf2.rds")
cox1<-readRDS("ngse4_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_rsf2t.rds")
#10
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_mtlrm.rds")
saveRDS(cox1,"ngse4_mtlr.rds")
cox1<-readRDS("ngse4_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_mtlrt.rds")
#11
start_time <- Sys.time()
pickedtime=timess[c(1,3,5)]
#t<-surv_train
#d<-cen_train
#qt=pickedtime
#getPseudoConditional(surv_train, cen_train, c(2,3,5,7,9))

#dnnsurv_fun(1,trail1_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, dnnsurv_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_dnnsurvm.rds")
saveRDS(cox1,"ngse4_dnnsurv.rds")
cox1<-readRDS("ngse4_dnnsurv.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_dnnsurvt.rds")
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
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_coxboostm.rds")
saveRDS(cox1,"ngse4_coxboost.rds")
cox1<-readRDS("ngse4_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_coxboostt.rds")

#13
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_original_cox_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_ga_cox1m.rds")
saveRDS(cox1,"ngse4_ga_cox1.rds")
cox1<-readRDS("ngse4_ga_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_ga_cox1t.rds")
#14
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_mtlr_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_ga_mtlrm.rds")
saveRDS(cox1,"ngse4_ga_mtlr.rds")
cox1<-readRDS("ngse4_ga_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_ga_mtlrt.rds")
#15
time1=timess[3]
stepnumber=10
penaltynumber=100
#ga_cox_boost_fun(1,trail1_1,5,(dim(trail1_1)[2]-3),10,20,time1,timess,stepnumber,penaltynumber)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_cox_boost_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_ga_coxboostm.rds")
saveRDS(cox1,"ngse4_ga_cocboost.rds")
cox1<-readRDS("ngse4_ga_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_ga_coxboostt.rds")

#16

start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun2,trail1_1,5, (dim(trail1_1)[2]-3),10,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_limma_mtlrm.rds")
saveRDS(cox1,"ngse4_limma_mtlr.rds")
cox1<-readRDS("ngse4_limma_mtlr.rds")
colMeans(do.call(rbind, cox1))
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_limma_mtlrt.rds")
#17
time1=timess[3]
stepnumber=10
penaltynumber=100
#cox_boost_fun2(1,trail1_1,5,(dim(trail1_1)[2]-3),10,time1,timess,stepnumber,penaltynumber)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun2,trail1_1,5, (dim(trail1_1)[2]-3),10,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_limma_coxboostm.rds")
saveRDS(cox1,"ngse4_limma_coxboost.rds")
cox1<-readRDS("ngse4_limma_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_limma_coxboostt.rds")

#18
trail1_1$time=as.numeric(trail1_1$time)
#survivalsvm_fun(1,veteran,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse4_survivalsvmm.rds")
saveRDS(cox1,"ngse4_survivalsvm.rds")
cox1<-readRDS("ngse4_survivalsvm.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse4_survivalsvmt.rds")



###############################################################



###############################################################
trail1_1=cbind.data.frame(trail5$X,trail5$delta,trail5$Y)
dim(trail1_1)
colnames(trail1_1)[(dim(trail1_1)[2]-1):dim(trail1_1)[2]]=c("status","time")
fitform_ogl=as.formula(Surv(time,status)~.)
form1=as.formula(~.)
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(trail1_1$time)[2]),as.numeric(summary(trail1_1$time)[5]),(as.numeric(summary(trail1_1$time)[5])-as.numeric(summary(trail1_1$time)[2]))/14)
hist(trail1_1$time)
abline(v=c(as.numeric(summary(trail1_1$time)[2]),as.numeric(summary(trail1_1$time)[5])), col="red")

trail1_1$os_class=as.vector(ifelse(trail1_1$status==1 & 
                                     trail1_1$time <80, "poor", 
                                   ifelse(trail1_1$status==0 & 
                                            trail1_1$time >80, "good", "not")))
table(trail1_1$os_class)
trail1_2=trail1_1[,-which(names(trail1_1) %in% c("os_class"))]
# data=trail1_1
# cvK=5
# j=1
# numm=(dim(trail1_1)[2]-3)
# topnumm=10
#19
#glmnet_lassocox_fun(1,trail1_2,5,formula1,formula2,formula3,formula4,timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_lassocox_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_p_cox1m.rds")
saveRDS(cox1,"ngse5_p_cox1.rds")
cox1<-readRDS("ngse5_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_p_cox1t.rds")
#20
#p_cox2_fun(11,trail1_1,5, form1,formula1,formula2,formula3,formula4,timess,0,50)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_ridge_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_p_cox2m.rds")
saveRDS(cox1,"ngse5_p_cox2.rds")
cox1<-readRDS("ngse5_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_p_cox2t.rds")
#21
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_en_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_p_cox3m.rds")
saveRDS(cox1,"ngse5_p_cox3.rds")
cox1<-readRDS("ngse5_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_p_cox3t.rds")


#8
# so strange why this is not working now
# data("veteran")
# fitform_ogl=as.formula(Surv(time,status)~.-trt)
# coxph(fitform_ogl,veteran)
# rfsrc(fitform_ogl, data = veteran,ntree = 1000,mtry = 10,tree.err=TRUE,importance = TRUE)

start_time <- Sys.time()
#rsf1_fun(1,trail1_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_rsf1m.rds")
saveRDS(cox1,"ngse5_rsf1.rds")
cox1<-readRDS("ngse5_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_rsf1t.rds")
#9
start_time <- Sys.time()
#rsf2_fun(1,trail1_1,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_rsf2m.rds")
saveRDS(cox1,"ngse5_rsf2.rds")
cox1<-readRDS("ngse5_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_rsf2t.rds")
#10
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_mtlrm.rds")
saveRDS(cox1,"ngse5_mtlr.rds")
cox1<-readRDS("ngse5_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_mtlrt.rds")
#11
start_time <- Sys.time()
pickedtime=timess[c(1,3,5)]
#t<-surv_train
#d<-cen_train
#qt=pickedtime
#getPseudoConditional(surv_train, cen_train, c(2,3,5,7,9))

#dnnsurv_fun(1,trail1_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, dnnsurv_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_dnnsurvm.rds")
saveRDS(cox1,"ngse5_dnnsurv.rds")
cox1<-readRDS("ngse5_dnnsurv.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_dnnsurvt.rds")
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
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_coxboostm.rds")
saveRDS(cox1,"ngse5_coxboost.rds")
cox1<-readRDS("ngse5_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_coxboostt.rds")

#13
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_original_cox_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_ga_cox1m.rds")
saveRDS(cox1,"ngse5_ga_cox1.rds")
cox1<-readRDS("ngse5_ga_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_ga_cox1t.rds")
#14
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_mtlr_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_ga_mtlrm.rds")
saveRDS(cox1,"ngse5_ga_mtlr.rds")
cox1<-readRDS("ngse5_ga_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_ga_mtlrt.rds")
#15
time1=timess[3]
stepnumber=10
penaltynumber=100
#ga_cox_boost_fun(1,trail1_1,5,(dim(trail1_1)[2]-3),10,20,time1,timess,stepnumber,penaltynumber)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_cox_boost_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_ga_coxboostm.rds")
saveRDS(cox1,"ngse5_ga_cocboost.rds")
cox1<-readRDS("ngse5_ga_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_ga_coxboostt.rds")

#16

start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun2,trail1_1,5, (dim(trail1_1)[2]-3),10,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_limma_mtlrm.rds")
saveRDS(cox1,"ngse5_limma_mtlr.rds")
cox1<-readRDS("ngse5_limma_mtlr.rds")
colMeans(do.call(rbind, cox1))
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_limma_mtlrt.rds")
#17
time1=timess[3]
stepnumber=10
penaltynumber=100
#cox_boost_fun2(1,trail1_1,5,(dim(trail1_1)[2]-3),10,time1,timess,stepnumber,penaltynumber)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun2,trail1_1,5, (dim(trail1_1)[2]-3),10,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_limma_coxboostm.rds")
saveRDS(cox1,"ngse5_limma_coxboost.rds")
cox1<-readRDS("ngse5_limma_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_limma_coxboostt.rds")

#18
trail1_1$time=as.numeric(trail1_1$time)
#survivalsvm_fun(1,veteran,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse5_survivalsvmm.rds")
saveRDS(cox1,"ngse5_survivalsvm.rds")
cox1<-readRDS("ngse5_survivalsvm.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse5_survivalsvmt.rds")




###############################################################




###############################################################
trail1_1=cbind.data.frame(trail6$X,trail6$delta,trail6$Y)
dim(trail1_1)
colnames(trail1_1)[(dim(trail1_1)[2]-1):dim(trail1_1)[2]]=c("status","time")
fitform_ogl=as.formula(Surv(time,status)~.)
form1=as.formula(~.)
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(trail1_1$time)[2]),as.numeric(summary(trail1_1$time)[5]),(as.numeric(summary(trail1_1$time)[5])-as.numeric(summary(trail1_1$time)[2]))/14)
hist(trail1_1$time)
trail1_1$os_class=as.vector(ifelse(trail1_1$status==1 & 
                                     trail1_1$time <5, "poor", 
                                   ifelse(trail1_1$status==0 & 
                                            trail1_1$time >5, "good", "not")))
table(trail1_1$os_class)
trail1_2=trail1_1[,-which(names(trail1_1) %in% c("os_class"))]
# data=trail1_1
# cvK=5
# j=1
# numm=(dim(trail1_1)[2]-3)
# topnumm=10
#19
#glmnet_lassocox_fun(1,trail1_2,5,formula1,formula2,formula3,formula4,timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_lassocox_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_p_cox1m.rds")
saveRDS(cox1,"ngse6_p_cox1.rds")
cox1<-readRDS("ngse6_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_p_cox1t.rds")
#20
#p_cox2_fun(11,trail1_1,5, form1,formula1,formula2,formula3,formula4,timess,0,50)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_ridge_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_p_cox2m.rds")
saveRDS(cox1,"ngse6_p_cox2.rds")
cox1<-readRDS("ngse6_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_p_cox2t.rds")
#21
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_en_fun,trail1_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_p_cox3m.rds")
saveRDS(cox1,"ngse6_p_cox3.rds")
cox1<-readRDS("ngse6_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_p_cox3t.rds")


#8
# so strange why this is not working now
# data("veteran")
# fitform_ogl=as.formula(Surv(time,status)~.-trt)
# coxph(fitform_ogl,veteran)
# rfsrc(fitform_ogl, data = veteran,ntree = 1000,mtry = 10,tree.err=TRUE,importance = TRUE)

start_time <- Sys.time()
#rsf1_fun(1,trail1_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_rsf1m.rds")
saveRDS(cox1,"ngse6_rsf1.rds")
cox1<-readRDS("ngse6_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_rsf1t.rds")
#9
start_time <- Sys.time()
#rsf2_fun(1,trail1_1,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_rsf2m.rds")
saveRDS(cox1,"ngse6_rsf2.rds")
cox1<-readRDS("ngse6_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_rsf2t.rds")
#10
start_time <- Sys.time()
#bw_cox3_fun(1,veteran,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_mtlrm.rds")
saveRDS(cox1,"ngse6_mtlr.rds")
cox1<-readRDS("ngse6_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_mtlrt.rds")
#11
start_time <- Sys.time()
pickedtime=timess[c(1,3,5)]
#t<-surv_train
#d<-cen_train
#qt=pickedtime
#getPseudoConditional(surv_train, cen_train, c(2,3,5,7,9))

#dnnsurv_fun(1,trail1_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, dnnsurv_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_dnnsurvm.rds")
saveRDS(cox1,"ngse6_dnnsurv.rds")
cox1<-readRDS("ngse6_dnnsurv.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_dnnsurvt.rds")
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
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_coxboostm.rds")
saveRDS(cox1,"ngse6_coxboost.rds")
cox1<-readRDS("ngse6_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_coxboostt.rds")

#13
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_original_cox_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_ga_cox1m.rds")
saveRDS(cox1,"ngse6_ga_cox1.rds")
cox1<-readRDS("ngse6_ga_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_ga_cox1t.rds")
#14
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_mtlr_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_ga_mtlrm.rds")
saveRDS(cox1,"ngse6_ga_mtlr.rds")
cox1<-readRDS("ngse6_ga_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_ga_mtlrt.rds")
#15
time1=timess[3]
stepnumber=10
penaltynumber=100
#ga_cox_boost_fun(1,trail1_1,5,(dim(trail1_1)[2]-3),10,20,time1,timess,stepnumber,penaltynumber)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_cox_boost_fun,trail1_1,5, (dim(trail1_1)[2]-3),10,20,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_ga_coxboostm.rds")
saveRDS(cox1,"ngse6_ga_cocboost.rds")
cox1<-readRDS("ngse6_ga_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_ga_coxboostt.rds")

#16

start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun2,trail1_1,5, (dim(trail1_1)[2]-3),10,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_limma_mtlrm.rds")
saveRDS(cox1,"ngse6_limma_mtlr.rds")
cox1<-readRDS("ngse6_limma_mtlr.rds")
colMeans(do.call(rbind, cox1))
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_limma_mtlrt.rds")
#17
time1=timess[3]
stepnumber=10
penaltynumber=100
#cox_boost_fun2(1,trail1_1,5,(dim(trail1_1)[2]-3),10,time1,timess,stepnumber,penaltynumber)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun2,trail1_1,5, (dim(trail1_1)[2]-3),10,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_limma_coxboostm.rds")
saveRDS(cox1,"ngse6_limma_coxboost.rds")
cox1<-readRDS("ngse6_limma_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_limma_coxboostt.rds")

#18
trail1_1$time=as.numeric(trail1_1$time)
#survivalsvm_fun(1,veteran,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,trail1_2,5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"ngse6_survivalsvmm.rds")
saveRDS(cox1,"ngse6_survivalsvm.rds")
cox1<-readRDS("ngse6_survivalsvm.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"ngse6_survivalsvmt.rds")
















