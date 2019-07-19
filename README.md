# CEPIM

CEPIM: a comparison and evaluation platform for integration methods in cancer subtyping.

## How to install

>Download [CEPIM_1.2.1.tar.gz](https://github.com/GaoLabXDU/CEPIM/releases/download/1.2.1/CEPIM_1.2.1.tar.gz). 

Install R package locally from R studio.

We have prepared init.R to help install the packages that CEPIM depends on.

>Download [init.R](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/init.R).

Source init.R in R studio.
```R
source(init.R)
```

## How to use

Just enter the data and some parameters to get the report.
```R
CEPIM(datalist, datatype = c("gaussian", "gaussian", "gaussian"),
    functionList = list('PINS', 'LRA', 'SNF', 'iCluster', 'PFA'), kMax=5)
```
**Please make sure that MATLAB has been installed on your computer if you want to run PFA.**   
For more details, please see [Supplementary Materials](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/CEPIM_Supplementary_materials.rar) or the help in this R package.

## Reports
Here, we present four reports generated by CEPIM for four scenarios. 

>Download [CEPIM Reports for all scenarios and Supplementary Materials](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/CEPIM_Supplementary_materials.rar).

The reports generated by CEPIM consist of **two parts**.   
***The first part*** shows the **comparison of subtyping results across different methods and different numbers of subtypes**, including *time consumption, Cox p-value, NMI, ARI, silhouette coefficient*. Then, we present an *overall samples similarity heatmap* representing a robust prediction for sample pairwise similarities. 

### Time Comsumption
>![TimeConsumption](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/TimeConsumption.png)

### Cox P-value
>![CoxPvalue](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/CoxPvalue.png)

### Silhouette Coefficient
>![SI](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/SI.png)

There are two comparison strategies for the NMI and ARI depending on the **availability of true labels of patients**. If you upload an empirical pre-determination of subtypes for samples by experts or clinicians based on clinical phenotypes, images or experience, CEPIM will **take it as gold standard to compare**. If you don’t have any subtype information, CEPIM will calculate the NMI and ARI **between every two methods from k=2 to the k-max**.

### NMI and ARI
- True label available

#### NMI
>![NMI](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/NMI.png)

#### ARI
>![ARI](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/ARI.png)

- True label not available
>![Comparison](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/Comparison.png)

### Overall samples similarity heatmap 
>![SamplesSimilarity](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/SamplesSimilarity.png)

***The second part*** shows the **performance of each method**, including *summary information of all metrics (when the true labels are available in Scenario 2 and 4, see [Supplementary Materials](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/CEPIM_Supplementary_materials.rar) for more details), KM survival curves, and patient similarity heatmaps* for different subtype numbers of each method.

### Summary Information
>![metrics](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/metrics.png)

### Kaplan-Meier Survival Curves
>![KM](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/KM.png)

### Patient Similarity HeatMap
>![heatmap](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/heatmap.png)

