# this is the code to run melanoma omics data: melanoma itraq, melanoma swath, melanoma nano (this is without extra fs)
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

#melanoma itraq
load("melanomaAssays.rdata")
class(melanomaAssays)
assays(melanomaAssays)
clinical1=as.data.frame(colData(melanomaAssays))

protein1=as.data.frame(assays(melanomaAssays)[["iTRAQprotein"]])
dim(protein1) #1898 41
#check for na
check_na_fun=function(data){
  check_na=as.data.frame(sapply(data, function(x) sum(is.na(x))))
  check_na$name=rownames(check_na)
  colnames(check_na)[1]="count"
  check_na1=check_na%>% filter(count>0)
  return(check_na1)}
check_na_fun(protein1)
#get rid of missing values
protein2=na.omit(protein1)
dim(protein2) #308 41
#add surv time and status
protein3=t(protein2)
protein3=as.data.frame(protein3)
dim(protein3) # 41 308
any(is.na(protein3))
protein3$name=rownames(protein3)
protein3$name=as.numeric(protein3$name)
protein3_1=merge(protein3,clinical1,by.x="name",by.y="Tumour_ID")
colnames(protein3_1)
protein3_2=protein3_1[,c(1:309,313,316)]
protein3_2$os_class=as.vector(ifelse(protein3_2$Person_FUStatus=="Dead, melanoma" & protein3_2$Prognosis_TimeSinceLNMet <365, "poor", ifelse(protein3_2$Person_FUStatus %in% c("Alive NSR") & protein3_2$Prognosis_TimeSinceLNMet >4*365, "good", "not")))
table(protein3_2$os_class)
colnames(protein3_2)[c(310,311)]=c("status","time")
protein3_2$status=as.vector(ifelse(protein3_2$status=="Dead, melanoma",1,0))
table(protein3_2$status)
colnames(protein3_2)=gsub("\\-", "_", colnames(protein3_2))
protein3_2$os_class=as.vector(ifelse(protein3_2$status==1 & 
                                       protein3_2$time <365*3, "poor", 
                                     ifelse(protein3_2$status==0 & 
                                              protein3_2$time >365*3, "good", "not")))


#fitform_ogl=as.formula(paste("Surv(time, status)~ ", paste(colnames(protein3_2)[2:308], collapse= "+")))
summary(protein3_2$time)
fitform_ogl=as.formula(Surv(time, status)~.)
form1=as.formula(~.)
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1


#############################################################################################
timess=seq(as.numeric(summary(protein3_2$time)[2]),as.numeric(summary(protein3_2$time)[5]),(as.numeric(summary(protein3_2$time)[5])-as.numeric(summary(protein3_2$time)[2]))/14)

#data<-protein3_2
#cvK=5
#j=1
#r=1
#original_cox_fun(1,protein3_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)

#5
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_lassocox_fun,protein3_2[,-dim(protein3_2)[2]],5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaitraq_p_cox1m.rds")
saveRDS(cox1,"melanomaitraq_p_cox1.rds")
cox1<-readRDS("melanomaitraq_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaitraq_p_cox1t.rds")
#6
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_ridge_fun,protein3_2[,-dim(protein3_2)[2]],5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaitraq_p_cox2m.rds")
saveRDS(cox1,"melanomaitraq_p_cox2.rds")
cox1<-readRDS("melanomaitraq_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaitraq_p_cox2t.rds")
#7
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_en_fun,protein3_2[,-dim(protein3_2)[2]],5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaitraq_p_cox3m.rds")
saveRDS(cox1,"melanomaitraq_p_cox3.rds")
cox1<-readRDS("melanomaitraq_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaitraq_p_cox3t.rds")

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
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun2,protein3_2,5,309,20,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaitraq_limma_coxboostm.rds")
saveRDS(cox1,"melanomaitraq_limma_coxboost.rds")
cox1<-readRDS("melanomaitraq_limma_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaitraq_limma_coxboostt.rds")


#10
start_time <- Sys.time()
#mtlr_fun(1,clinical4_new,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun2,protein3_2,5,309,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaitraq_limma_mtlrm.rds")
saveRDS(cox1,"melanomaitraq_limma_mtlr.rds")
cox1<-readRDS("melanomaitraq_limma_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaitraq_limma_mtlrt.rds")

#15
#ga_original_cox_fun(1,protein3_2,5,309,10,20,timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_original_cox_fun,protein3_2,5, 309,10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"itraq_ga_cox1m.rds")
saveRDS(cox1,"itraq_ga_cox1.rds")
cox1<-readRDS("itraq_ga_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"itraq_ga_cox1t.rds")
#10
start_time <- Sys.time()
#mtlr_fun(1,clinical4_new,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_mtlr_fun,protein3_2,5,309,5,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaitraq_ga_mtlrm.rds")
saveRDS(cox1,"melanomaitraq_ga_mtlr.rds")
cox1<-readRDS("melanomaitraq_ga_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaitraq_ga_mtlrt.rds")

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
cox1 <- pbmcapply::pbmclapply(1:20, ga_cox_boost_fun,protein3_2,5,309,5,20,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaitraq_ga_coxboostm.rds")
saveRDS(cox1,"melanomaitraq_ga_coxboost.rds")
cox1<-readRDS("melanomaitraq_ga_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaitraq_ga_coxboostt.rds")


#20
protein3_2$time=as.numeric(protein3_2$time)
#survivalsvm_fun(1,protein3_2,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,protein3_2[,-dim(protein3_2)[2]],5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"itraq_survivalsvmm.rds")
saveRDS(cox1,"itraq_survivalsvm.rds")
cox1<-readRDS("itraq_survivalsvm.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"itraq_survivalsvmt.rds")

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
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,protein3_2[,-dim(protein3_2)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaitraq_coxboostm.rds")
saveRDS(cox1,"melanomaitraq_coxboost.rds")
cox1<-readRDS("melanomaitraq_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaitraq_coxboostt.rds")

#8
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,protein3_2[,-dim(protein3_2)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaitraq_rsf1m.rds")
saveRDS(cox1,"melanomaitraq_rsf1.rds")
cox1<-readRDS("melanomaitraq_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaitraq_rsf1t.rds")
#9
#rsf2_fun(1,protein3_2[,-dim(protein3_2)[2]],5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,protein3_2[,-dim(protein3_2)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaitraq_rsf2m.rds")
saveRDS(cox1,"melanomaitraq_rsf2.rds")
cox1<-readRDS("melanomaitraq_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaitraq_rsf2t.rds")
#10
start_time <- Sys.time()
#mtlr_fun(1,clinical4_new,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,protein3_2[,-dim(protein3_2)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaitraq_mtlrm.rds")
saveRDS(cox1,"melanomaitraq_mtlr.rds")
cox1<-readRDS("melanomaitraq_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaitraq_mtlrt.rds")
#11
start_time <- Sys.time()
pickedtime=c(30,45)
#data<-protein3_2
#cvK=5
#j=1
#t<-surv_train
#d<-cen_train
#qt=pickedtime
#getPseudoConditional(surv_train, cen_train, 30)
#getPseudoConditional(surv_train, cen_train, c(30,45))
#dnnsurv_fun(1,protein3_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, dnnsurv_fun,protein3_2[,-dim(protein3_2)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaitraq_dnnsurvm.rds")
saveRDS(cox1,"melanomaitraq_dnnsurv.rds")
cox1<-readRDS("melanomaitraq_dnnsurv.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaitraq_dnnsurvt.rds")

###############################################################




###############################################################
#swath
protein2=as.data.frame(assays(melanomaAssays)[["SWATHprotein"]])
dim(protein2) #2307 70
#check for na
check_na_fun=function(data){
  check_na=as.data.frame(sapply(data, function(x) sum(is.na(x))))
  check_na$name=rownames(check_na)
  colnames(check_na)[1]="count"
  check_na1=check_na%>% filter(count>0)
  return(check_na1)}
check_na_fun(protein2)
#get rid of missing values
protein2_1=na.omit(protein2)
dim(protein2_1) #1626 70
#add surv time and status
protein2_2=t(protein2_1)
protein2_2=as.data.frame(protein2_2)
dim(protein2_2) # 70 1626
any(is.na(protein2_2))
protein2_2$name=rownames(protein2_2)
protein2_2$name=as.numeric(protein2_2$name)
protein2_3=merge(protein2_2,clinical1,by.x="name",by.y="Tumour_ID")
colnames(protein2_3)[1626:1635]
protein2_4=protein2_3[,c(1:1627,1631,1634)]
colnames(protein2_4)[c(1628,1629)]=c("status","time")
protein2_4$status=as.vector(ifelse(protein2_4$status=="Dead, melanoma",1,0))
table(protein2_4$status)
#colnames(protein2_4)=gsub("\\-", "_", colnames(protein2_4))
#have to change the name because still, some begin with numbers
xnam <- paste(colnames(protein2_4)[2:1627], sep="")
xnam<-data.frame(xnam)
namedata<-xnam%>% separate(xnam, into = c("pre", "post"))
cols=c("pre","post")
namedata[,c(1,2)]=namedata[,c(2,1)]
namedata$x <- apply( namedata[ , cols ] , 1 , paste , collapse = "_" )
head(namedata)
fitform_ogl=as.formula(Surv(time, status)~.)
form1=as.formula(~.)
colnames(protein2_4)[2:1627]=namedata$x
summary(protein2_4$time)
hist(protein2_4$time)
protein2_4=protein2_4[,2:1629]
protein2_4$os_class=as.vector(ifelse(protein2_4$status==1 & 
                                       protein2_4$time <365, "poor", 
                                     ifelse(protein2_4$status==0 & 
                                              protein2_4$time >365*4, "good", "not")))
table(protein2_4$os_class)

formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(protein2_4$time)[2]),as.numeric(summary(protein2_4$time)[5]),(as.numeric(summary(protein2_4$time)[5])-as.numeric(summary(protein2_4$time)[2]))/14)


#5
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_lassocox_fun,protein2_4[,-dim(protein2_4)[2]],5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaswath_p_cox1m.rds")
saveRDS(cox1,"melanomaswath_p_cox1.rds")
cox1<-readRDS("melanomaswath_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaswath_p_cox1t.rds")
#6
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_ridge_fun,protein2_4[,-dim(protein2_4)[2]],5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaswath_p_cox2m.rds")
saveRDS(cox1,"melanomaswath_p_cox2.rds")
cox1<-readRDS("melanomaswath_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaswath_p_cox2t.rds")
#7
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_en_fun,protein2_4[,-dim(protein2_4)[2]],5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaswath_p_cox3m.rds")
saveRDS(cox1,"melanomaswath_p_cox3.rds")
cox1<-readRDS("melanomaswath_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaswath_p_cox3t.rds")

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
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun2,protein2_4,5, 1626,100,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaswath_limma_coxboostm.rds")
saveRDS(cox1,"melanomaswath_limma_coxboost.rds")
cox1<-readRDS("melanomaswath_limma_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaswath_limma_coxboostt.rds")

#10
start_time <- Sys.time()
#mtlr_fun(1,protein2_4,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun2,protein2_4,5, 1626,100,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaswath_limma_mtlrm.rds")
saveRDS(cox1,"melanomaswath_limma_mtlr.rds")
cox1<-readRDS("melanomaswath_limma_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaswath_limma_mtlrt.rds")

#2
#ga_original_cox_fun(1,protein2_4,5,1626,10,20,timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_original_cox_fun,protein2_4,5, 1626,10,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaswath_ga_cox1m.rds")
saveRDS(cox1,"melanomaswath_ga_cox1.rds")
cox1<-readRDS("melanomaswath_ga_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaswath_ga_cox1t.rds")
#10
start_time <- Sys.time()
#mtlr_fun(1,clinical4_new,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, ga_mtlr_fun,protein2_4,5,1626,5,20,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaswath_ga_mtlrm.rds")
saveRDS(cox1,"melanomaswath_ga_mtlr.rds")
cox1<-readRDS("melanomaswath_ga_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaswath_ga_mtlrt.rds")

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
cox1 <- pbmcapply::pbmclapply(1:20, ga_cox_boost_fun,protein2_4,5,1626,5,20,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaswath_ga_coxboostm.rds")
saveRDS(cox1,"melanomaswath_ga_coxboost.rds")
cox1<-readRDS("melanomaswath_ga_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaswath_ga_coxboostt.rds")

#4
protein2_4$time=as.numeric(protein2_4$time)
#survivalsvm_fun(1,protein2_4,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,protein2_4[,-dim(protein2_4)[2]],5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaswath_survivalsvmm.rds")
saveRDS(cox1,"melanomaswath_survivalsvm.rds")
cox1<-readRDS("melanomaswath_survivalsvm.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaswath_survivalsvmt.rds")

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
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,protein2_4[,-dim(protein2_4)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaswath_coxboostm.rds")
saveRDS(cox1,"melanomaswath_coxboost.rds")
cox1<-readRDS("melanomaswath_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaswath_coxboostt.rds")


#8
start_time <- Sys.time()
#rsf1_fun(1,protein2_4,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,protein2_4[,-dim(protein2_4)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaswath_rsf1m.rds")
saveRDS(cox1,"melanomaswath_rsf1.rds")
cox1<-readRDS("melanomaswath_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"lmelanomaswath_rsf1t.rds")
#9
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,protein2_4[,-dim(protein2_4)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaswath_rsf2m.rds")
saveRDS(cox1,"melanomaswath_rsf2.rds")
cox1<-readRDS("melanomaswath_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaswath_rsf2t.rds")
#10
start_time <- Sys.time()
#mtlr_fun(1,protein2_4,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,protein2_4[,-dim(protein2_4)[2]],5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomaswath_mtlrm.rds")
saveRDS(cox1,"melanomaswath_mtlr.rds")
cox1<-readRDS("melanomaswath_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomaswath_mtlrt.rds")

###############################################################




###############################################################
#melanoma nanostring
# load("melanomaAssays.RData")
# class(melanomaAssays)
# assays(melanomaAssays)
# 
# clinical1=as.data.frame(colData(melanomaAssays))
nano1=as.data.frame(assays(melanomaAssays)[["NanoString"]])
dim(nano1) #204 45
#check for na
any(is.na(nano1))
#add surv time and status
nano2=t(nano1)
nano2=as.data.frame(nano2)
dim(nano2) 
any(is.na(nano2))
nano2$name=rownames(nano2)
nano2$name=as.numeric(nano2$name)
nano2_1=merge(nano2,clinical1,by.x="name",by.y="Tumour_ID")
colnames(nano2_1)
nano2_2=nano2_1[,c(1:205,209,212)]
colnames(nano2_2)[c(206,207)]=c("status","time")
nano2_2$status=as.vector(ifelse(nano2_2$status=="Dead, melanoma",1,0))
table(nano2_2$status) #19 26
colnames(nano2_2)=gsub("\\-", "_", colnames(nano2_2))
colnames(nano2_2)=gsub("\\(", "_", colnames(nano2_2))
colnames(nano2_2)=gsub("\\)", "_", colnames(nano2_2))
fitform_ogl=as.formula(paste("Surv(time, status)~ ", paste(colnames(nano2_2)[2:205], collapse= "+")))
summary(nano2_2$time)

# nano2_2$os_class=as.vector(ifelse(nano2_2$status==1 & 
#                                     nano2_2$time <365*3, "poor", 
#                                   ifelse(nano2_2$status==0 & 
#                                            nano2_2$time >365*3, "good", "not")))
# table(nano2_2$os_class)
formula1=fitform_ogl
formula2=fitform_ogl
formula3=Surv(time,status)~1
formula4=Surv(time,status)~1
timess=seq(as.numeric(summary(nano2_2$time)[2]),as.numeric(summary(nano2_2$time)[5]),(as.numeric(summary(nano2_2$time)[5])-as.numeric(summary(nano2_2$time)[2]))/14)

#2
bw_cox1_fun(1,nano2_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess)
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox1_fun,nano2_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomanano_bw_cox1m.rds")
saveRDS(cox1,"melanomanano_bw_cox1.rds")
cox1<-readRDS("melanomanano_bw_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomanano_bw_cox1t.rds")
#3
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox2_fun,nano2_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomanano_bw_cox2m.rds")
saveRDS(cox1,"melanomanano_bw_cox2.rds")
cox1<-readRDS("melanomanano_bw_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomanano_bw_cox2t.rds")
#4
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, bw_cox3_fun,nano2_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomanano_bw_cox3m.rds")
saveRDS(cox1,"melanomanano_bw_cox3.rds")
cox1<-readRDS("melanomanano_bw_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomanano_bw_cox3t.rds")
#5
#data<-nano2_2
#cvK=5
#j=1
#r=1
#lambda_1=1
#labmda_2=0
#p_cox1_fun(1,nano2_2,5,form1,formula1,formula2,formula3,formula4,timess,1,0)
start_time <- Sys.time()
form1=as.formula(paste("~ ", paste(colnames(nano2_2)[1:28], collapse= "+")))
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_lassocox_fun,nano2_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomanano_p_cox1m.rds")
saveRDS(cox1,"melanomanano_p_cox1.rds")
cox1<-readRDS("melanomanano_p_cox1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomanano_p_cox1t.rds")
#6
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_ridge_fun,nano2_2,5,formula1,formula2,formula3,formula4,timess,mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomanano_p_cox2m.rds")
saveRDS(cox1,"melanomanano_p_cox2.rds")
cox1<-readRDS("mmelanomanano_p_cox2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomanano_p_cox2t.rds")
#7
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, glmnet_en_fun,nano2_2,5,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomanano_p_cox3m.rds")
saveRDS(cox1,"melanomanano_p_cox3.rds")
cox1<-readRDS("melanomanano_p_cox3.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomanano_p_cox3t.rds")
#8
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf1_fun,nano2_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomanano_rsf1m.rds")
saveRDS(cox1,"melanomanano_rsf1.rds")
cox1<-readRDS("melanomanano_rsf1.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomanano_rsf1t.rds")
#9
start_time <- Sys.time()
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, rsf2_fun,nano2_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomanano_rsf2m.rds")
saveRDS(cox1,"melanomanano_rsf2.rds")
cox1<-readRDS("melanomanano_rsf2.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomanano_rsf2t.rds")
#10
start_time <- Sys.time()
#mtlr_fun(1,clinical4_new,5,fitform_ogl,formula1,formula2,formula3,formula4,timess)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, mtlr_fun,nano2_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomanano_mtlrm.rds")
saveRDS(cox1,"melanomanano_mtlr.rds")
cox1<-readRDS("melanomanano_mtlr.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomanano_mtlrt.rds")
#11: 
start_time <- Sys.time()
pickedtime=c(530,600)
#data<-nano2_2
#cvK=5
#j=1
#t<-surv_train
#d<-cen_train
#qt=pickedtime
#getPseudoConditional(surv_train, cen_train, c(30,100))
#getPseudoConditional(surv_train, cen_train, pickedtime)
#dnnsurv_fun(1,nano2_2,5,fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime)
Rprof(tf <- "rprof.log",memory.profiling=TRUE)
cox1 <- pbmcapply::pbmclapply(1:20, dnnsurv_fun,nano2_2,5, fitform_ogl,formula1,formula2,formula3,formula4,timess,pickedtime, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomanano_dnnsurvm.rds")
saveRDS(cox1,"melanomanano_dnnsurv.rds")
cox1<-readRDS("melanomanano_dnnsurv.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomanano_dnnsurvt.rds")

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
cox1 <- pbmcapply::pbmclapply(1:20, cox_boost_fun,nano2_2,5, fitform_ogl,formula1,formula2,formula3,formula4,time1,timess,stepnumber,penaltynumber, mc.cores = 15)
Rprof(NULL)
mm<-summaryRprof(tf,memory = "both")
mm
saveRDS(mm,"melanomanano_coxboostm.rds")
saveRDS(cox1,"melanomanano_coxboost.rds")
cox1<-readRDS("melanomanano_coxboost.rds")
head(cox1)
end_time <- Sys.time()
saveRDS(end_time - start_time,"melanomanano_coxboostt.rds")


# #4
# nano2_2$time=as.numeric(nano2_2$time)
# #survivalsvm_fun(1,nano2_2,5, fitform_ogl,formula1,formula2,formula3, formula4, timess)
# # Error in quadprog::solve.QP(C, -d, t(H), f, meq = meq) : 
# # constraints are inconsistent, no solution! 
# start_time <- Sys.time()
# Rprof(tf <- "rprof.log",memory.profiling=TRUE)
# cox1 <- pbmcapply::pbmclapply(1:20, survivalsvm_fun,nano2_2,5, fitform_ogl,formula1,formula2,formula3, formula4, timess, mc.cores = 15)
# Rprof(NULL)
# mm<-summaryRprof(tf,memory = "both")
# mm
# saveRDS(mm,"nano_survivalsvmm.rds")
# saveRDS(cox1,"nano_survivalsvm.rds")
# cox1<-readRDS("nano_survivalsvm.rds")
# head(cox1)
# end_time <- Sys.time()
# saveRDS(end_time - start_time,"nano_survivalsvmt.rds")








