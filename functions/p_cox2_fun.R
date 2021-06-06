p_cox2_fun=function(r,data,cvK,form1,formula1,formula2,formula3, formula4, timess,lambda_1,lambda_2){
  set.seed(r)
  print(r)
  cvSets = cvFolds(nrow(data), cvK)  # permute all the data, into 5 folds
  bicfun=possibly(function(j){
    test_id = cvSets$subsets[cvSets$which == j]
    test = data[test_id, ]
    train = data[-test_id, ]
    
    fit=penalized(Surv(time,status),penalized = form1,lambda1 = lambda_1,lambda2 = lambda_2,data = train)
    predd=predict(fit,penalized=form1,data=test)
    predd1=as.matrix(predd)
    
    lpp<- predict(fit,penalized=form1,data=train)
    lpp1=as.matrix(lpp)
    
    #cindex
    predd2=predd1[,round(dim(predd1)[2]/2,digits = 0)]
    harrelC1 <- rcorr.cens(predd2,Surv(test$time,test$status))
    hc_1<-harrelC1["C Index"]
    lpp2=-lpp1[,round(dim(lpp1)[2]/2,digits = 0)]
    lpnew <- -predd2
    Surv.rsp <- Surv(train$time, train$status)
    Surv.rsp.new <- Surv(test$time, test$status) 
    bc_1 <- BeggC(Surv.rsp, Surv.rsp.new,lpp2, lpnew)
    unoc_1<-UnoC(Surv.rsp, Surv.rsp.new, lpnew)
    ghc_1<-GHCI(lpnew)
    
    
    #br
    lp=lpp2
    briers1 <- predErr(Surv.rsp, Surv.rsp.new, lp, lpnew,times=test$time, type = "brier", int.type = "unweighted")$error
    br1<-sum(na.omit(briers1))
    briers2<-predErr(Surv.rsp, Surv.rsp.new, lp, lpnew,times=test$time, type = "brier", int.type = "weighted")$error
    br2<-sum(na.omit(briers2))
    ibsfun1=possibly(function(modell){
      briers3 <- pec(list("cox1"=modell),data=test,formula=fitform_ogl1,cens.model="cox")
      return(crps(briers3)[2])
    },otherwise = NA)
    br3<-NA
    ibsfun2=possibly(function(modell){
      briers4 <- pec(list("cox1"=modell),data=test,formula=fitform_ogl2,cens.model="marginal")
      return(crps(briers4)[2])
    },otherwise = NA)
    br4<-NA
    ibsfun3=possibly(function(modell){
      briers5 <- pec(list("cox1"=modell),data=test,formula=fitform_ogl3,cens.model="cox")
      return(crps(briers5)[2])
    },otherwise = NA)
    br5<-NA
    ibsfun4=possibly(function(modell){
      briers6 <- pec(list("cox1"=modell),data=test,formula=fitform_ogl4,cens.model="marginal")
      return(crps(briers6)[2])
    },otherwise = NA)
    br6<-NA
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
    return(c(hc_1,bc_1,unoc_1,ghc_1,br1,br2,br3,br4,br5,br6,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a))},otherwise=NA)
  
  cv5_result=rbind.data.frame(bicfun(1),bicfun(2),bicfun(3),bicfun(4),bicfun(5))   
  colnames(cv5_result)=c("hc_1","bc_1","unoc_1","ghc_1","br1","br2","br3","br4","br5","br6","a1","a2","a3","a4","a5","a6","a7","a8","a9","a10","a11","a12","a13","a14","a15","a")
  
  
  return(cv5_result)}