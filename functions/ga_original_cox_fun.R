ga_original_cox_fun=function(r,data,cvK,numm,topnumm,generation_num, timess){
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
    
    train2=train[!train$os_class=="not",]
    test2=test[!test$os_class=="not",]
    
    #ohhh, i think the problem is their gene data is column:individual
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
    
    original_cox1=coxph(fitform_ogl,x=TRUE,data=train)
    predd<-predict(object=original_cox1,newdata=test,type="risk")
    #predd<-predictSurvProb(original_cox1,newdata=test,times=seq(365,365*15,365))
    #harrel cindex
    harrelC1 <- rcorr.cens(-predd,with(test,Surv(time,status)))
    hc_acc5[j]<-harrelC1["C Index"]
    #begg cindex
    lp<- predict(original_cox1)
    lpnew <- predict(original_cox1, newdata=test)
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
