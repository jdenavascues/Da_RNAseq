---
title: "BatchQC Report"
date: "2023-03-21"
output: 
  html_vignette:
    toc: true
    toc_depth: 2
    template: batchQC.html
    self_contained: no
    lib_dir: libs
---


Summary
=======
## Confounding
### Number of samples in each Batch and Condition

-------------------------------------------
        &nbsp;           Batch a   Batch b 
----------------------- --------- ---------
 **Condition Control**      3         3    

 **Condition DaDaOE**       0         3    

  **Condition DaKD**        3         0    

  **Condition DaOE**        2         0    

  **Condition ScOE**        0         3    
-------------------------------------------

### Measures of confounding between Batch and Condition

----------------------------------------------------------------------
            &nbsp;                Standardized Pearson     Cramer's V 
                                 Correlation Coefficient              
------------------------------- ------------------------- ------------
  **Confounding Coefficients             0.8859              0.8036   
 (0=no confounding, 1=complete                                        
        confounding)**                                                
----------------------------------------------------------------------

## Variation Analysis
### Variation explained by Batch and Condition
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />


----------------------------------------------------------
   &nbsp;      Full (Condition+Batch)   Condition   Batch 
------------- ------------------------ ----------- -------
  **Min.**             0.108              0.085       0   

 **1st Qu.**           63.36              63.3      4.024 

 **Median**            81.56              81.55     13.65 

  **Mean**             74.74              74.72     15.69 

 **3rd Qu.**           91.68              91.68     24.09 

  **Max.**             99.98              99.98     60.31 
----------------------------------------------------------

## P-value Analysis
### Distribution of Batch and Condition Effect p-values Across Genes

---------------------------------------------------------------------------------------------
         &nbsp;            Min.     1st Qu.     Median     Mean     3rd Qu.   Max.   Ps<0.05 
------------------------ -------- ----------- ---------- --------- --------- ------ ---------
   **Batch P-values**     0.1697    0.9169      0.9545    0.9408    0.9797     1        0    

 **Condition P-values**     0      2.416e-05   0.001446   0.07725   0.04215    1     0.7626  
---------------------------------------------------------------------------------------------

<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />


Differential Expression
=======================
## Expression Plot
Boxplots for all values for each of the samples and are colored by batch membership.

<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

## LIMMA

--------------------------------------------------------------------------------------------------------------------
  &nbsp;    Condition: DaDaOE (logFC)   Condition: DaKD (logFC)   Condition: DaOE (logFC)   Condition: ScOE (logFC) 
---------- --------------------------- ------------------------- ------------------------- -------------------------
 **6710**             91.35                     -30.25                    -26.17                     7479           

 **7194**             5.069                     -34.05                    -62.12                     4722           

 **1991**             557.3                     -281.8                     944.4                     17996          

 **5129**             23.33                     -13.25                     130.4                     1599           

 **3771**             1572                       -1203                     1042                      17409          

 **2117**            -3.686                     -19.05                    -11.24                     1222           

 **462**             -21.06                     -175.7                    -185.6                     7064           

 **813**              57.37                     -12.22                     196.6                     5075           

 **6013**            -7.343                      4.333                      164                      1549           

 **8637**            -41.86                     -90.75                    -59.15                     10430          
--------------------------------------------------------------------------------------------------------------------

Table: Table continues below

 
----------------------------------------------------
  &nbsp;    AveExpr     F      P.Value    adj.P.Val 
---------- --------- ------- ----------- -----------
 **6710**    1556     10594   1.088e-20   1.086e-16 

 **7194**    1068     4948    9.277e-19   4.629e-15 

 **1991**    3538     3932    3.552e-18   1.182e-14 

 **5129**    370.9    2238    9.528e-17   2.377e-13 

 **3771**    6489     1691    4.876e-16   9.732e-13 

 **2117**    282.1    1632    5.998e-16   9.977e-13 

 **462**     1412     1542    8.36e-16    1.192e-12 

 **813**     967.4    1504    9.654e-16   1.204e-12 

 **6013**    318.8    1415    1.382e-15   1.533e-12 

 **8637**    2040     1380    1.599e-15   1.596e-12 
----------------------------------------------------


Median Correlations
===================
This plot helps identify outlying samples.
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />


Heatmaps
========
## Heatmap
This is a heatmap of the given data matrix showing the batch effects and variations with different conditions.
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

## Sample Correlations
This is a heatmap of the correlation between samples.
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />


Circular Dendrogram
===================
This is a Circular Dendrogram of the given data matrix colored by batch to show the batch effects.
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />


PCA: Principal Component Analysis
=================================
## PCA
This is a plot of the top two principal components colored by batch to show the batch effects.
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-20-1.png" style="display: block; margin: auto;" />

## Explained Variation

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  &nbsp;    Proportion of Variance (%)   Cumulative Proportion of   Percent Variation Explained by   Percent Variation Explained by   Condition Significance   Percent Variation Explained by   Batch Significance (p-value) 
                                               Variance (%)           Either Condition or Batch                Condition                    (p-value)                      Batch                                             
---------- ---------------------------- -------------------------- -------------------------------- -------------------------------- ------------------------ -------------------------------- ------------------------------
 **PC1**               40.2                        40.2                           99                               99                           0                           26.2                           0.9809            

 **PC2**              15.55                       55.74                          98.9                             98.9                          0                           19.6                           0.9863            

 **PC3**              13.55                       69.29                          95.3                             95.3                          0                           13.3                           0.9239            

 **PC4**              6.589                       75.88                          95.6                             95.6                          0                           4.4                            0.9511            

 **PC5**              4.193                       80.08                          3.1                              3.1                         0.9852                        0.1                            0.9999            

 **PC6**              3.816                       83.89                          2.5                              2.5                         0.9917                        0.2                            0.9722            

 **PC7**              2.909                        86.8                          0.8                              0.8                         0.9989                         0                             0.9927            

 **PC8**              2.511                       89.31                          1.1                               1                          0.9986                        0.2                            0.9519            

 **PC9**              2.247                       91.56                          0.3                              0.3                         0.9998                         0                             0.9972            

 **PC10**             1.921                       93.48                          0.8                              0.8                         0.999                          0                             0.9627            

 **PC11**             1.623                        95.1                          1.7                              1.7                         0.9956                        0.1                            0.9691            

 **PC12**             1.502                        96.6                          0.2                              0.2                         0.9999                         0                             0.9899            

 **PC13**             1.328                       97.93                          0.5                              0.5                         0.9995                         0                             0.9875            

 **PC14**             1.107                       99.04                          0.2                              0.2                           1                            0                             0.9869            

 **PC15**             0.9443                      99.98                          0.1                              0.1                           1                            0                             0.986             

 **PC16**            0.01639                       100                           99.9                              0                            0                           35.8                             0               

 **PC17**            3.77e-29                      100                           26.7                             25.9                        0.4688                        1.3                            0.7288            
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


Shape
=====
This is a heatmap plot showing the variation of gene expression mean, variance, skewness and kurtosis between samples grouped by batch to see the batch effects variation
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-23-1.png" style="display: block; margin: auto;" />

```
## Note: Sample-wise p-value is calculated for the variation across samples on the measure across genes. Gene-wise p-value is calculated for the variation of each gene between batches on the measure across each batch. If the data is quantum normalized, then the Sample-wise measure across genes is same for all samples and Gene-wise p-value is a good measure.
```


Combat Plots
============
This is a plot showing whether parametric or non-parameteric prior is appropriate for this data. It also shows the Kolmogorov-Smirnov test comparing the parametric and non-parameteric prior distribution.

```
## Found 2 batches
## Adjusting for 4 covariate(s) or covariate level(s)
## Standardizing Data across genes
## Fitting L/S model and finding priors
```

<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-25-1.png" style="display: block; margin: auto;" /><img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-25-2.png" style="display: block; margin: auto;" /><img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-25-3.png" style="display: block; margin: auto;" /><img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/combat_batchqc_report_files/figure-html/unnamed-chunk-25-4.png" style="display: block; margin: auto;" />

```
## Batch mean distribution across genes: Normal vs Empirical distribution
## Two-sided Kolmogorov-Smirnov test
## Selected Batch: 1
## Statistic D = 0.0438
## p-value = 0
## 
## 
## Batch Variance distribution across genes: Inverse Gamma vs Empirical distribution
## Two-sided Kolmogorov-Smirnov test
## Selected Batch: 1
## Statistic D = 0.1219
## p-value = 0Note: The non-parametric version of ComBat takes much longer time to run and we recommend it only when the shape of the non-parametric curve widely differs such as a bimodal or highly skewed distribution. Otherwise, the difference in batch adjustment is very negligible and parametric version is recommended even if p-value of KS test above is significant.
```


SVA
===
## Summary

```
## Number of Surrogate Variables found in the given data: 1
```
