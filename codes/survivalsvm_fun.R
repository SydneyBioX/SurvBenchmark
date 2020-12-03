# data("veteran")
# library(survivalsvm)
# survsvm <- survivalsvm(data = veteran, time.variable.name = "diagtime",
#                            status.variable.name = "status",
#                            type = "hybrid", gamma.mu = c(0.1,0.1),
#                            opt.meth = "quadprog", diff.meth = "makediff3",
#                            kernel = "lin_kernel",
#                            sgf.sv = 5, sigf = 7, maxiter = 20,
#                            margin = 0.05, bound = 10)
# survsvm$var.names
# print(survsvm)
# predd=predict(survsvm,veteran) #they claim it should be the risk rank, but what is a risk rank? it seems to be the larger the value, the lower the risk
# predd$predicted #we have negative values here, but in their github, it is treated as the probability to get cindex
# rcorr.cens(predd$predicted,with(veteran,Surv(diagtime,status)))
# lp<- -predd$predicted
# lpnew <- -predd$predicted
# Surv.rsp <- Surv(veteran$diagtime, veteran$status)
# Surv.rsp.new <- Surv(veteran$diagtime, veteran$status)    
# AUC.uno(Surv.rsp, Surv.rsp.new, lpnew, c(1,11)) #according to the definition the auc is based on comparison,so maybe we can calculate auc?
# 


survivalsvm_fun=function(r,data,cvK,fitform_ogl,formula1,formula2,formula3, formula4, timess){
  set.seed(r)
  print(r)
  hc_acc5 = c()
  bc_acc5=c()
  unoc_acc5=c()
  ghc_acc5=c()
  bs1=c()
  bs2=c()
  bs3=c()
  bs4=c()
  bs5=c()
  bs6=c()
  auc1=auc2=auc3=auc4=auc5=auc6=auc7=auc8=auc9=auc10=auc11=auc12=auc13=auc14=auc15=auc=c()
  cvSets = cvFolds(nrow(data), cvK)  # permute all the data, into 5 folds
  for (j in 1:cvK) {
    print(j)
    test_id = cvSets$subsets[cvSets$which == j]
    test = data[test_id, ]
    train = data[-test_id, ]
    survsvm <- survivalsvm(data = train, time.variable.name = "time",
                               status.variable.name = "status",
                               type = "hybrid", gamma.mu = c(0.1,0.1),
                               opt.meth = "quadprog", diff.meth = "makediff3",
                               kernel = "lin_kernel",
                               sgf.sv = 5, sigf = 7, maxiter = 20,
                               margin = 0.05, bound = 10)
    predd=predict(survsvm,test)$predicted
    #predd<-predictSurvProb(original_cox1,newdata=test,times=seq(365,365*15,365))
    #harrel cindex
    harrelC1 <- rcorr.cens(predd,with(test,Surv(time,status)))
    hc_acc5[j]<-harrelC1["C Index"]
    #begg cindex
    lp<- -predict(survsvm,train)$predicted
    lpnew <- -predd
    Surv.rsp <- Surv(train$time, train$status)
    Surv.rsp.new <- Surv(test$time, test$status) 
    # beggcfun<- possibly(function(lp){rr=BeggC(Surv.rsp, Surv.rsp.new,lp, lpnew)
    # return(rr)}
    # ,otherwise = NA)
    bc_acc5[j] <- NA
    #uno cindex
    unoc_acc5[j]<-UnoC(Surv.rsp, Surv.rsp.new, lpnew)
    #gh cindex
    #ghc_acc5[j]<-GHCI(lpnew)
    ghc_acc5[j]<-NA
    #br
    briers1 <- predErr(Surv.rsp, Surv.rsp.new, lp, lpnew,times=test$time, type = "brier", int.type = "unweighted")$error
    bs1[j]<-sum(na.omit(briers1))
    briers2<-predErr(Surv.rsp, Surv.rsp.new, lp, lpnew,times=test$time, type = "brier", int.type = "weighted")$error
    bs2[j]<-sum(na.omit(briers2))
    ibsfun1=possibly(function(modell){
      briers3 <- pec(list("cox1"=modell),data=test,formula=formula1,cens.model="cox")
      return(crps(briers3)[2])
    },otherwise = NA)
    bs3[j]<-ibsfun1(original_cox1)
    ibsfun2=possibly(function(modell){
      briers4 <- pec(list("cox1"=modell),data=test,formula=formula2,cens.model="marginal")
      return(crps(briers4)[2])
    },otherwise = NA)
    bs4[j]<-ibsfun2(original_cox1)
    ibsfun3=possibly(function(modell){
      briers5 <- pec(list("cox1"=modell),data=test,formula=formula3,cens.model="cox")
      return(crps(briers5)[2])
    },otherwise = NA)
    bs5[j]<-ibsfun3(original_cox1)
    ibsfun4=possibly(function(modell){
      briers6 <- pec(list("cox1"=modell),data=test,formula=formula4,cens.model="marginal")
      return(crps(briers6)[2])
    },otherwise = NA)
    bs6[j]<-ibsfun4(original_cox1)
    #time-dependent auc
    times <- timess
    AUC_CD <- AUC.uno(Surv.rsp, Surv.rsp.new, lpnew, times)
    auc1[j]=AUC_CD$auc[1]
    auc2[j]=AUC_CD$auc[2]
    auc3[j]=AUC_CD$auc[3]
    auc4[j]=AUC_CD$auc[4]
    auc5[j]=AUC_CD$auc[5]
    auc6[j]=AUC_CD$auc[6]
    auc7[j]=AUC_CD$auc[7]
    auc8[j]=AUC_CD$auc[8]
    auc9[j]=AUC_CD$auc[9]
    auc10[j]=AUC_CD$auc[10]
    auc11[j]=AUC_CD$auc[11]
    auc12[j]=AUC_CD$auc[12]
    auc13[j]=AUC_CD$auc[13]
    auc14[j]=AUC_CD$auc[14]
    auc15[j]=AUC_CD$auc[15]
    auc[j]=AUC_CD$iauc
  }
  want=cbind.data.frame(hc_acc5,bc_acc5,unoc_acc5,ghc_acc5,bs1,bs2,bs3,bs4,bs5,bs6,auc1,auc2,auc3,auc4,auc5,auc6,auc7,auc8,auc9,auc10,auc11,auc12,auc13,auc14,auc15,auc)
  return(want)}
