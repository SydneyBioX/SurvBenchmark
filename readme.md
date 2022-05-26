# SurvBenchmark: comprehensive benchmarking study of survival analysis methods using both omics data and clinical data


This is the work for SurvBenchmark (202205 updated) and the associated paper can be found: doi (to be added). 

# Introduction 
We develop a benchmarking design, SurvBenchmark, that evaluates a diverse collection of survival models for both clinical and omics datasets. SurvBenchmark not only focuses on classical approaches such as the Cox model, but it also evaluates state-of-art machine learning survival models. 
There are 16 datasets (https://github.com/SydneyBioX/SurvBenchmark/blob/main/tables/table1.docx) 
**Table1.** Datasets summary

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

and 20 survival methods (https://github.com/SydneyBioX/SurvBenchmark/blob/main/tables/table2.docx) 
benchmarked in this study. 

##############################################################################

# Files in this repo 

In this repo, all the high resolution figures related to the paper can be found under folder <figures>.

The folder <functions> contains functions to run all methods.

The folder <datasets> contains all datasets benchmarked in our paper.

The folder <figures_data> contains all figure data used to generate the figures in our paper. 

The github_example.R file gives an example to get the results using methods in <functions> on the Ovarian dataset. 

For the datasets we used, please check this Table1 in our paper, this is also under <tables> table1.

For the survival methods we benchmarked, please check Table2 in our paper, this is also under <tables> table2.

The R package is available at https://github.com/yunwezhang/SurvBenchmark_package, on-going work will be updated continuously. 


###############################################################################

# Installation 

```r
library(devtools)
devtools::install_github("SydneyBioX/SurvBenchmark_package")
library(SurvBenchmark)
```


## Requirements


You may need to install the following dependencies first:   

```r
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
library(survcomp)
library(survAUC)
library(CoxBoost)
library(limma)
library(partykit)
library(coin)
library(compound.Cox)
library(GenAlgo)
library(survivalsvm)
library(rmatio)
library(survivalmodels)
library(reticulate)
```


## Visualise the results

The comparison of survival models can be visualized using heatmap as the below example.  


<img src="figures/github_figure.png" align="center" width="400"/>



# Reference

paper (to be added)


# License

```r
Copyright [2022] [Yunwei Zhang]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
```



   
