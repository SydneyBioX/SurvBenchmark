# SurvBenchmark: comprehensive benchmarking study of survival analysis methods using both omics data and clinical data


This is the work for SurvBenchmark (202205 updated) and the associated paper can be found: doi (to be added). 

# Introduction 
We develop a benchmarking design, SurvBenchmark, that evaluates a diverse collection of survival models for both clinical and omics datasets. SurvBenchmark not only focuses on classical approaches such as the Cox model, but it also evaluates state-of-art machine learning survival models. 
There are 16 datasets (https://github.com/SydneyBioX/SurvBenchmark/blob/main/tables/table1.docx) and 20 survival methods (https://github.com/SydneyBioX/SurvBenchmark/blob/main/tables/table2.docx) benchmarked in this study. 

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



   
