# CEPIM

CEPIM: a comparison and evaluation platform for integration methods in cancer subtyping.

## How to install

Download [CEPIM_1.0.0.tar.gz](https://github.com/GaoLabXDU/CEPIM/releases/download/v1.0.0/CEPIM_1.0.0.tar.gz). 

Install R package locally from R studio.

## How to use

Just enter the data and some parameters to get the report.
```R
CEPIM(datalist, datatype = c("gaussian", "gaussian", "gaussian"),
    functionList = list('PINS', 'LRA', 'SNF', 'iCluster', 'PFA'),
    evalList = list('Heatmap','KM', 'SI', 'ARI', 'RI', 'NMI'), kMax=5)
```
For more details, please see this Guideline ***(Guideline file)*** or the help in this R package.

## Reports
Here, we present four reports generated by CEPIM for four scenarios. 

>[S1.html](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/S1.html)
[S2.html](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/S1.html)
[S3.html](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/S1.html)
[S4.html](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/S1.html)

***Please download the report to your local computer, and modify the extension name to ".html" before open the file.***

The reports generated by CEPIM include **two parts**.   
***The first part*** shows the **comparison of subtyping results across different methods and different number of subtypes**, including *time consumption, Cox p-value, NMI, ARI, silhouette coefficient*. Then, we present an *overall samples similarity heatmap* representing a robust prediction for sample pairwise similarities. 

There are two comparison strategies for the NMI and ARI depending on the **availability of true labels of patients**. If you upload an empirical pre-determination of subtypes for samples by experts or clinicians based on clinical phenotypes, images or experience, CEPIM will **take it as gold standard to compare**. If you don’t have any subtype information, CEPIM will calculate the NMI and ARI **between every two methods from k=2 to the k-max**.

Time Comsumption
>![TimeConsumption](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/TimeConsumption.png)

Cox P-value
>![CoxPvalue](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/CoxPvalue.png)

Silhouette Coefficient
>![SI](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/SI.png)

NMI and ARI
- True label available
>![NMI](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/NMI.png)
***An ARI Figure***
- True label not available
>![Comparison](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/Comparison.png)

Overall samples similarity heatmap 
>![SamplesSimilarity](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/SamplesSimilarity.png)

***The second part*** shows the **performance of each method**, including *summary information of all metrics (only in scenario 2 and 4, see Guideline for more details), KM survival curves, and patient similarity heatmaps* for different subtype numbers of each method.

Summary Information
>![metrics](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/metrics.png)

Kaplan-Meier Survival Curves
>![KM](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/KM.png)

Patient Similarity HeatMap
>![heatmap](https://github.com/GaoLabXDU/CEPIM/raw/master/documents/pic/heatmap.png)

