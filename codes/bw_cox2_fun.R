bw_cox2_fun=function(r,data,cvK,fitform_ogl,formula1,formula2,formula3, formula4, timess){
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
    
    original_cox1=cph(fitform_ogl,x=TRUE,y=TRUE,data=train)
    z<-fastbw(original_cox1, rule="p",sls=.05)
    if(isEmpty(z$names.kept)){coxx=coxph(Surv(time,status)~1,data = train)} else{
    coxx=coxph(as.formula(paste("Surv(time,status)~",paste(z$names.kept, collapse = "+"),sep = "")),data = train)}
    predd<-predict(object=coxx,newdata=test,type="risk")
    #predd<-predictSurvProb(original_cox1,newdata=test,times=seq(365,365*15,365))
    #harrel cindex
    harrelC1 <- rcorr.cens(-predd,with(test,Surv(time,status)))
    hc_acc5[j]<-harrelC1["C Index"]
    #begg cindex
    lp<- predict(coxx)
    lpnew <- predict(coxx, newdata=test)
    Surv.rsp <- Surv(train$time, train$status)
    Surv.rsp.new <- Surv(test$time, test$status)              
    bc_acc5[j] <- BeggC(Surv.rsp, Surv.rsp.new,lp, lpnew)
    #uno cindex
    unoc_acc5[j]<-UnoC(Surv.rsp, Surv.rsp.new, lpnew)
    #gh cindex
    ghc_acc5[j]<-GHCI(lpnew)
    #br
    briers1 <- predErr(Surv.rsp, Surv.rsp.new, lp, lpnew,times=test$time, type = "brier", int.type = "unweighted")$error
    bs1[j]<-sum(na.omit(briers1))
    briers2<-predErr(Surv.rsp, Surv.rsp.new, lp, lpnew,times=test$time, type = "brier", int.type = "weighted")$error
    bs2[j]<-sum(na.omit(briers2))
    ibsfun1=possibly(function(modell){
      briers3 <- pec(list("cox1"=modell),data=test,formula=formula1,cens.model="cox")
      return(crps(briers3)[2])
    },otherwise = NA)
    
    bs3[j]<-ibsfun1(coxx)
    ibsfun2=possibly(function(modell){
      briers4 <- pec(list("cox1"=modell),data=test,formula=formula2,cens.model="marginal")
      return(crps(briers4)[2])
    },otherwise = NA)
    bs4[j]<-ibsfun2(coxx)
    ibsfun3=possibly(function(modell){
      briers5 <- pec(list("cox1"=modell),data=test,formula=formula3,cens.model="cox")
      return(crps(briers5)[2])
    },otherwise = NA)
    bs5[j]<-ibsfun3(coxx)
    ibsfun4=possibly(function(modell){
      briers6 <- pec(list("cox1"=modell),data=test,formula=formula4,cens.model="marginal")
      return(crps(briers6)[2])
    },otherwise = NA)
    bs6[j]<-ibsfun4(coxx)
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