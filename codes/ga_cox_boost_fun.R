ga_cox_boost_fun=function(r,data,cvK,numm,topnumm, generation_num,time1,timess,stepnumber,penaltynumber){
  set.seed(r)
  print(r)
  cvSets = cvFolds(nrow(data), cvK)  # permute all the data, into 5 folds
  bicfun=possibly(function(j){
    test_id = cvSets$subsets[cvSets$which == j]
    test = data[test_id, ]
    train = data[-test_id, ]
    train2=train[!train$os_class=="not",]
    test2=test[!test$os_class=="not",]
    
    ##ga feature selection
    Data=t(train2[,1:numm]) #gene by patient
    myContext <- list(dataset=Data, gps=train2$os_class)
    n.individuals <- round(numm*0.8,digits = 0) #searching into 80% of genes
    n.features <- topnumm #how many genes we would like to select
    y <- matrix(0, n.individuals, n.features)
    for (i in 1:n.individuals) {
      
      y[i,] <- sample(1:nrow(Data), n.features)
    }
    mahaFitness <- function(arow, context) {
      maha(t(context$dataset[arow,]), context$gps, method='var')
    }
    my.ga <- GenAlg(y, mahaFitness, selectionMutate, myContext, 0.001, 0.75) #there might be sigularity issues
    for (i in 1:generation_num) {
      my.ga <- newGeneration(my.ga)
    }
    #summary(my.ga)
    selectedname <- rownames(Data[my.ga@best.individual,])
    #print(selectedname)
    
    
    
    train=train[,colnames(train)%in%c(selectedname,"status","time")]
    test=test[,colnames(test)%in%c(selectedname,"status","time")]
    
    
    #fitform_ogl=as.formula(paste("Surv(time, status)~ ", paste(colnames(train)[1:(dim(train)[2]-2)], collapse= "+")))
    fitform_ogl=Surv(time,status)~.
    form1=as.formula(~.)
    formula1=fitform_ogl
    formula2=fitform_ogl
    formula3=Surv(time,status)~1
    formula4=Surv(time,status)~1
    
    #form1=as.formula(paste("~ ",paste(colnames(train)[1:(dim(train)[2]-2)], collapse= "+")))
    
    tr_predictormatrix=train[,-which(colnames(train)%in% c("time","status"))]
    tr_predictormatrix=data.matrix(tr_predictormatrix)
    te_predictormatrix=test[,-which(colnames(test)%in% c("time","status"))]
    te_predictormatrix=data.matrix(te_predictormatrix)
    coxboost=CoxBoost(time=train$time,status=train$status,x=tr_predictormatrix,stepno=stepnumber,penalty=penaltynumber)
    lpnew=predict(coxboost,type="risk",times=time1,newdata=te_predictormatrix)
    #it is the survival probability #the predicted probability of not yet having had the event at the time points given in times
    lpnew=-lpnew #change to risk
    
    lp<- predict(coxboost,type="risk",times=time1,newdata=tr_predictormatrix)
    lp=-lp
    #cindex
    harrelC1 <- rcorr.cens(-lpnew,Surv(test$time,test$status))
    hc_1<-harrelC1["C Index"]
    Surv.rsp <- Surv(train$time, train$status)
    Surv.rsp.new <- Surv(test$time, test$status) 
    bc_1 <- BeggC(Surv.rsp, Surv.rsp.new,lp, lpnew)
    unoc_1<-UnoC(Surv.rsp, Surv.rsp.new, lpnew)
    ghc_1<-GHCI(lpnew)
    
    
    #br
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