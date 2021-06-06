rsf1_fun=function(r,data,cvK,fitform_ogl,formula1,formula2,formula3, formula4, timess){
  set.seed(r)
  print(r)
  cvSets = cvFolds(nrow(data), cvK)  # permute all the data, into 5 folds
  bicfun=possibly(function(j){
    test_id = cvSets$subsets[cvSets$which == j]
    test = data[test_id, ]
    train = data[-test_id, ]
    set.seed(123)
    out.rsf.1 <- rfsrc(fitform_ogl, 
                       data = train,
                       ntree = 1000,
                       mtry = 10,
                       tree.err=TRUE,
                       importance = TRUE)
    pred.test.fin = predict( out.rsf.1, 
                             newdata = test, 
                             importance = "none" ,type="response",block.size = 1)
    pred.test.fin2 = predict( out.rsf.1, 
                              newdata = train, 
                              importance = "none" ,type="response",block.size = 1)
    pred_te=pred.test.fin$predicted
    pred_tr=pred.test.fin2$predicted
    
    #predd<-predictSurvProb(original_cox1,newdata=test,times=seq(365,365*15,365))
    #harrel cindex
    harrelC1 <- rcorr.cens(-pred_te,with(test,Surv(time,status)))
    hc<-harrelC1["C Index"]
    #begg cindex
    lp<- pred_tr
    lpnew <- pred_te
    Surv.rsp <- Surv(train$time, train$status)
    Surv.rsp.new <- Surv(test$time, test$status)              
    bc <- BeggC(Surv.rsp, Surv.rsp.new,lp, lpnew)
    #uno cindex
    unoc<-UnoC(Surv.rsp, Surv.rsp.new, lpnew)
    #gh cindex
    ghc<-GHCI(lpnew)
    #br
    briers1 <- predErr(Surv.rsp, Surv.rsp.new, lp, lpnew,times=test$time, type = "brier", int.type = "unweighted")$error
    br1<-sum(na.omit(briers1))
    briers2<-predErr(Surv.rsp, Surv.rsp.new, lp, lpnew,times=test$time, type = "brier", int.type = "weighted")$error
    br2<-sum(na.omit(briers2))
    ibsfun1=possibly(function(modell){
      briers3 <- pec(list("cox1"=modell),data=test,formula=formula1,cens.model="cox")
      return(crps(briers3)[2])
    },otherwise = NA)
    #briers3 <- pec(list("cox1"=original_cox1),data=test,formula=Surv(tx_gperiod,tx_gstatus)~recip_sex+recip_eth+recip_age+recip_height+recip_weight+recip_smoker+recip_lung+recip_coronary+recip_pvd+recip_cvd+recip_diabetes+recip_waittime+donor_age+donor_sex+donor_height+donor_weight+donor_causedeath_cva+donor_dcd+donor_diabetes+donor_ht+donor_smoker+donor_creatinine+tx_ischaemia+tx_misa+tx_misb+tx_misdr,cens.model="cox")
    #bs3[j]<-crps(briers3)[2]
    br3<-ibsfun1(out.rsf.1)
    ibsfun2=possibly(function(modell){
      briers4 <- pec(list("cox1"=modell),data=test,formula=formula2,cens.model="marginal")
      return(crps(briers4)[2])
    },otherwise = NA)
    br4<-ibsfun2(out.rsf.1)
    ibsfun3=possibly(function(modell){
      briers5 <- pec(list("cox1"=modell),data=test,formula=formula3,cens.model="cox")
      return(crps(briers5)[2])
    },otherwise = NA)
    br5<-ibsfun3(out.rsf.1)
    ibsfun4=possibly(function(modell){
      briers6 <- pec(list("cox1"=modell),data=test,formula=formula4,cens.model="marginal")
      return(crps(briers6)[2])
    },otherwise = NA)
    br6<-ibsfun4(out.rsf.1)
    #time-dependent auc
    times <- timess
    AUC_CD <- AUC.uno(Surv.rsp, Surv.rsp.new, lpnew, times)
    a1=AUC_CD$auc[1]
    a2=AUC_CD$auc[2]
    a3=AUC_CD$auc[3]
    a4=AUC_CD$auc[4]
    a5=AUC_CD$auc[5]
    a6=AUC_CD$auc[6]
    a7=AUC_CD$auc[7]
    a8=AUC_CD$auc[8]
    a9=AUC_CD$auc[9]
    a10=AUC_CD$auc[10]
    a11=AUC_CD$auc[11]
    a12=AUC_CD$auc[12]
    a13=AUC_CD$auc[13]
    a14=AUC_CD$auc[14]
    a15=AUC_CD$auc[15]
    a=AUC_CD$iauc
    return(c(hc,bc,unoc,ghc,br1,br2,br3,br4,br5,br6,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a))},otherwise=NA)
  
  cv5_result=rbind.data.frame(bicfun(1),bicfun(2),bicfun(3),bicfun(4),bicfun(5))   
  colnames(cv5_result)=c("hc","bc","unoc","ghc","br1","br2","br3","br4","br5","br6","a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15","a")
  
  # want=cbind.data.frame(hc_acc5,bc_acc5,unoc_acc5,ghc_acc5,bs1,bs2,bs3,bs4,bs5,bs6,auc1,auc2,auc3,auc4,auc5,auc6,auc7,auc8,auc9,auc10,auc11,auc12,auc13,auc14,auc15,auc)
  return(cv5_result)}