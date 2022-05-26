This file contains the tables for the datasets and methods we benchmarked.

**Table1.** Datasets summary

| **Datasets summary** |
| **Dataset (name usedinthispaper)** | **Number ofobservations** | **No. ofvariables** | **Type ofdata** | **Censoringrate(rounded to 4decimalplaces)**|**Reference** |
| --- | --- | --- | --- | --- |--- |
| Melanoma\_itraq | 41 | 642 | Omics | 0.4146 | Wang,K.Y.X. et al. Cross-Platform Omics Prediction procedure: agame changer for implementing precision medicine in patientswithstage-IIImelanoma.bioRxiv2020.12.09.415927;doi:[https://doi.org/10.1101/2020.12.09.415927](https://doi.org/10.1101/2020.12.09.415927) |
| Melanoma\_nano | 45 | 206 | Omics | 0.4222 | Wang,K.Y.X. et al. Cross-Platform Omics Prediction procedure: agame changer for implementing precision medicine in patientswithstage-IIImelanoma.bioRxiv2020.12.09.415927;doi:[https://doi.org/10.1101/2020.12.09.415927](https://doi.org/10.1101/2020.12.09.415927) |
| Ovarian\_2 | 58 | 19818 | Omics | 0.3793 | Ganzfried,B.F.etal.(2013)curatedOvarianData:clinicallyannotateddatafortheovariancancertranscriptome.Database,2013. |
| GE\_5 | 78 | 4753 | Omics | 0.5641 | van&#39;tVeer,L.J.etal.(2002)Geneexpressionprofilingpredictsclinical outcomeofbreast cancer.Nature,415,530–536. |
| GE\_3 | 86 | 6288 | Omics | 0.7209 | Bullinger,L.etal.(2004)UseofGene-ExpressionProfilingtoIdentifyPrognostic Subclasses in Adult Acute Myeloid Leukemia. NewEnglandJournalofMedicine,350, 1605–1616. |
| Melanoma\_clinical | 88 | 16 | Clinical | 0.3939 | Wang,K.Y.X. et al. Cross-Platform Omics Prediction procedure: agame changer for implementing precision medicine in patientswithstage-IIImelanoma.bioRxiv2020.12.09.415927;doi:[https://doi.org/10.1101/2020.12.09.415927.](https://doi.org/10.1101/2020.12.09.415927) |
| GE\_1 | 115 | 551 | Omics | 0.6670 | Sorlie,T. et al. (2003) Repeated observation of breast tumor subtypesin independent gene expression data sets. Proc. Natl. Acad. Sci.U.S. A., 100, 8418–8423. |
| GE-\_4 | 116 | 6285 | Omics | 0.5641 | van de Vijver,M.J. et al. (2002) A gene-expression signature as apredictorofsurvivalinbreastcancer.N.Engl.J.Med.,347,1999–2009. |
| Veteran | 137 | 8 | Clinical | 0.0657 | Kalbfleisch,J.D.andPrentice,R.L.(2002)TheStatisticalAnalysisofFailureTimeData.WileySeriesinProbabilityandStatistics. |
| Ovarian\_1 | 194 | 16050 | Omics | 0.7062 | Ganzfried,B.F.etal.(2013)curatedOvarianData:clinicallyannotateddatafortheovariancancertranscriptome.Database,2013. |
| Lung | 228 | 9 | Clinical | 0.2763 | Loprinzi,C.L.etal.(1994)Prospectiveevaluationofprognosticvariables from patient-completed questionnaires. North CentralCancerTreatment Group.J. Clin.Oncol., 12,601–607. |
| GE\_6 | 240 | 7401 | Omics | 0.4250 | Van Houwelingen,H.C. (2004) The Elements of Statistical Learning,Data Mining, Inference, and Prediction. Trevor Hastie, RobertTibshirani and Jerome Friedman, Springer, New York, 2001. No.of pages: xvi 533. ISBN 0-387-95284-5. Statistics in Medicine,23, 528–529. |
| GE\_2 | 295 | 4921 | Omics | 0.7322 | Beer,D.G.etal.(2002)Gene-expressionprofilespredictsurvivalofpatientswithlungadenocarcinoma.Nat.Med.,8,816–824. |
| PBC | 312 | 7 | Clinical | 0.5994 | Fleming,T.R.andHarrington,D.P.(2005)CountingProcessesandSurvivalAnalysis.WileySeriesinProbabilityandStatistics. |
| UNOS\_Kidney | 3000 | 101 | Clinical | 0.7350 | OPTNdata (https://optn.transplant.hrsa.gov/) |
| ANZ | 3323 | 40 | Clinical | 0.8739 | ANZDATA (https://www.anzdata.org.au/) |



**Table2.** Summary of methods used in this study

| **Methodname** | **Methodnameinthispaper** | **Rfunctionname** | **Rpackagename** | **Parameters (default)** |
| --- | --- | --- | --- | --- |
| Cox | Cox | coxph | survival | NA |
| Cox withbackwardeliminationusingAIC | Cox\_bw\_AIC | cph,fastbw | rms | rule=&quot;aic&quot;,sls=.05,k.aic=2 |
| Coxwithbackwardeliminationusingpvalue | Cox\_bw\_p | cph,fastbw | rms | rule=&quot;p&quot;,sls=.05 |
| Cox withbackwardeliminationusingBIC | Cox\_bw\_BIC | cph,fastbw | rms | rule=&quot;aic&quot;,sls=.05,k.aic=log(as.numeric(table(train$status)[2])) |
| Lassocox (for clinical datasets) | Lasso\_Cox | penalized | penalized | Lambda1=1,lambda2=0 |
| Ridgecox (for clinical datasets) | Ridge\_Cox | penalized | penalized | Lambda1=0,lambda2=1 |
| Elasticnetcox (for clinical datasets) | EN\_Cox | penalized | penalized | Lambda1=1,lambda2=1 |
| Lassocox (for omics datasets) | Lasso\_Cox | glmnet | glmnet | alpha=1,nfolds=5,type.measure=&quot;C&quot; |
| Ridgecox (for omics datasets) | Ridge\_Cox | glmnet | glmnet | alpha=0,nfolds=5,type.measure=&quot;C&quot; |
| Elasticnetcox (for omics datasets) | EN\_Cox | glmnet | glmnet | alpha=0.5,nfolds=5,type.measure=&quot;C&quot; |
| Randomsurvivalforest | RSF | rfsrc | RandomSurvivalForest | Default:ntree=1000,mtry=10 |
| Multitasklogisticregressionmethod | MTLR | mtlr | MTLR | C1=1 |
|DNNSurv (Deeplearningsurvivalmodel) |DNNSurv | multiplefunctionsasinGithubcodes | DNNSurv|Default: no parameter arguments to be changed by users |
| Boostingcoxmodel | CoxBoost | coxboost | CoxBoost | stepnumber=10, penalty number=100 |
| Coxmodelwithgeneticalgorithmasfeatureselectionmethod | Cox (GA) | GenAlg | GenAlgo | n.features=10(foromics),n.features=4(forclinical),generation\_num=20 |
| Multitasklogistic regressionmodelwithgenetic algorithmas featureselectionmethod | MTLR(GA) | GenAlg | GenAlgo | n.features=10 (foromics),n.features=4 (forclinical),generation\_num=20 |
| Boostingcoxmodelwithgeneticalgorithmasfeatureselectionmethod | CoxBoost (GA) | GenAlg | GenAlgo | n.features=10(foromics),n.features=4(forclinical),generation\_num=20 |
| Multitasklogisticregressionmodelwithrankingbasedmethodasfeatureselectionmeth | MTLR(DE) | lmFit,eBayes | limma | n.features=10(foromics),n.features=4(forclinical) |
| Boostingcox modelwithrankingbasedmethodas feature selectionmethod | CoxBoost (DE) | lmFit,eBayes | limma | n.features=10(foromics),n.features=4(forclinical) |
| Survival supportvectormachine | SurvivalSVM | survivalsvm | survivalsvm | Default: sgf.sv = 5, sigf = 7, maxiter = 20, margin = 0.05, bound = 10, eig.tol = 1e-06, conv.tol = 1e-07, posd.tol = 1e-08 |
|DeepSurv(Deeplearningsurvival model) |DeepSurv |deepsurv |survivalmodels | Default:frac=0.3,activation=&quot;relu&quot;,num\_nodes=c(4L,8L,4L,2L),dropout=0.1,early\_stopping=TRUE,epochs=100L,batch\_size=32L |
|DeepHit(Deeplearningsurvival model) |DeepHit |deephit |survivalmodels | Default:frac=0.3,activation=&quot;relu&quot;,num\_nodes=c(4L,8L,4L,2L),dropout=0.1,early\_stopping=TRUE,epochs=100L,batch\_size=32L |

