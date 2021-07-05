# this shows how to use SurvBenchmark 
# this is an example to run those methods and get corresponding evaluation results 
# here, we use the published gene expression dataset: Ovarian data
# updated by Yunwei Zhang, 20210705

.libPaths("/dskh/nobackup/yunwei") #setting R library path, you need to change it correspondingly
setwd("/dskh/nobackup/yunwei/Analysis_yw/benchmark/202012summaryresult") #setting R working directory, you need to change it correspondingly

#load all the packages first
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
library(pseudo)
library(survivalROC)
library(survival)
library(survAUC)
library(CoxBoost)
library(limma)
library(partykit)
library(coin)
library(compound.Cox)
library(GenAlgo)
library(survivalsvm)
library(rmatio)

# import all those source files from the folder <functions>
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

# get the Ovarian gene expression data 
library(curatedOvarianData)
data("GSE49997_eset")
expmatrix2=exprs(GSE49997_eset)
dim(expmatrix2)

expmatrix2_1=t(expmatrix2)
cancerdt2=cbind.data.frame(expmatrix2_1,GSE49997_eset$vital_status,GSE49997_eset$days_to_death)
colnames(cancerdt2)[16049:16050]=c("status","time")
cancerdt2$status=as.vector(ifelse(cancerdt2$status=="living",0,1))
table(cancerdt2$status)
summary(cancerdt2$time)
cancerdt2_1=cancerdt2[complete.cases(cancerdt2),]
dim(cancerdt2_1)

#check the survival time distribution
#hist(cancerdt2_1$time)
cancerdt2_1$os_class=as.vector(ifelse(cancerdt2_1$status==1 & 
                                        cancerdt2_1$time <2*365, "poor", 
                                      ifelse(cancerdt2_1$status==0 & 
                                               cancerdt2_1$time >2*365, "good", "not")))


fitform_ogl=Surv(time,status)~.
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
form1=as.formula(~.)
timess=seq(as.numeric(summary(cancerdt2_1$time)[2]),as.numeric(summary(cancerdt2_1$time)[5]),(as.numeric(summary(cancerdt2_1$time)[5])-as.numeric(summary(cancerdt2_1$time)[2]))/14)

# lasso_cox model
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_lassocox_fun,cancerdt2_1[,-dim(cancerdt2_1)[2]],5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_p_cox1m.rds")
saveRDS(cox1,"gse1_p_cox1.rds")
cox1<-readRDS("gse1_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_p_cox1t.rds")
# ridge_cox model
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_ridge_fun,cancerdt2_1[,-dim(cancerdt2_1)[2]],5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_p_cox2m.rds")
saveRDS(cox1,"gse1_p_cox2.rds")
cox1<-readRDS("gse1_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_p_cox2t.rds")
# elastic_net_cox method
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_en_fun,cancerdt2_1[,-dim(cancerdt2_1)[2]],5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_p_cox3m.rds")
saveRDS(cox1,"gse1_p_cox3.rds")
cox1<-readRDS("gse1_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_p_cox3t.rds")
# random survival forest
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,cancerdt2_1[,-dim(cancerdt2_1)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_rsf1m.rds")
saveRDS(cox1,"gse1_rsf1.rds")
cox1<-readRDS("gse1_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_rsf1t.rds")
# random survival forest
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,cancerdt2_1[,-dim(cancerdt2_1)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_rsf2m.rds")
saveRDS(cox1,"gse1_rsf2.rds")
cox1<-readRDS("gse1_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_rsf2t.rds")
# mtlr
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,cancerdt2_1[,-dim(cancerdt2_1)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_mtlrm.rds")
saveRDS(cox1,"gse1_mtlr.rds")
cox1<-readRDS("gse1_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_mtlrt.rds")
# CoxBoost
start_time <- Sys.time()
time1=timess[3]
stepnumber=10
penaltynumber=100
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,cancerdt2_1[,-dim(cancerdt2_1)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_coxboostm.rds")
saveRDS(cox1,"gse1_coxboost.rds")
cox1<-readRDS("gse1_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_coxboostt.rds")

# mtlr (GA) 
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_mtlr_fun,cancerdt2_1,5, 16047,5,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_ga_mtlrm.rds")
saveRDS(cox1,"gse1_ga_mtlr.rds")
cox1<-readRDS("gse1_ga_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_ga_mtlrt.rds")
# CoxBoost (GA)
start_time <- Sys.time()
time1=timess[3]
stepnumber=10
penaltynumber=100
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_cox_boost_fun,cancerdt2_1,5, 16047,5,20,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_ga_coxboostm.rds")
saveRDS(cox1,"gse1_ga_coxboost.rds")
cox1<-readRDS("gse1_ga_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_coxboostt.rds")

# Cox (GA)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_original_cox_fun,cancerdt2_1,5, 16047,10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_ga_cox1m.rds")
saveRDS(cox1,"gse1_ga_cox1.rds")
cox1<-readRDS("gse1_ga_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_ga_cox1t.rds")

# SurvivalSVM
cancerdt2_1$time=as.numeric(cancerdt2_1$time)
#survivalsvm_fun(1,cancerdt2_1,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,cancerdt2_1[,-dim(cancerdt2_1)[2]],5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_survivalsvmm.rds")
saveRDS(cox1,"gse1_survivalsvm.rds")
cox1<-readRDS("gse1_survivalsvm.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_survivalsvmt.rds")


# mtlr (DE)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun2,cancerdt2_1,5, 16047,1000,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_limma_mtlrm.rds")
saveRDS(cox1,"gse1_limma_mtlr.rds")
cox1<-readRDS("gse1_limma_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_limma_mtlrt.rds")
# CoxBoost (DE)
start_time <- Sys.time()
time1=timess[3]
stepnumber=10
penaltynumber=100
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun2,cancerdt2_1,5, 16047,1000,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_limma_coxboostm.rds")
saveRDS(cox1,"gse1_limma_coxboost.rds")
cox1<-readRDS("gse1_limma_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_limma_coxboostt.rds")

# deepsurv
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, deepsurv_fun,cancerdt2_1[,-dim(cancerdt2_1)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_deepsurvm.rds")
saveRDS(cox1,"gse1_deepsurv.rds")
cox1<-readRDS("gse1_deepsurv.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_deepsurvt.rds")

# deephit
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, deephit_fun,cancerdt2_1[,-dim(cancerdt2_1)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"gse1_deephitm.rds")
saveRDS(cox1,"gse1_deephit.rds")
cox1<-readRDS("gse1_deephit.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"gse1_deephitt.rds")



