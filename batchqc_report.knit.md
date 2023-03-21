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
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />


----------------------------------------------------------
   &nbsp;      Full (Condition+Batch)   Condition   Batch 
------------- ------------------------ ----------- -------
  **Min.**             0.934              0.786       0   

 **1st Qu.**            72.2              71.55     11.19 

 **Median**            85.44              85.02     26.79 

  **Mean**             79.53              79.06     27.21 

 **3rd Qu.**            92.4              92.19     42.36 

  **Max.**             99.92              99.92     67.31 
----------------------------------------------------------

## P-value Analysis
### Distribution of Batch and Condition Effect p-values Across Genes

-------------------------------------------------------------------------------------------------
         &nbsp;             Min.      1st Qu.     Median     Mean     3rd Qu.    Max.    Ps<0.05 
------------------------ ---------- ----------- ---------- --------- --------- -------- ---------
   **Batch P-values**      0.2207     0.5178      0.6917    0.6795    0.8512    0.9999      0    

 **Condition P-values**   2.22e-16   4.333e-05   0.001321   0.05563   0.02435   0.9986   0.8088  
-------------------------------------------------------------------------------------------------

<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />


Differential Expression
=======================
## Expression Plot
Boxplots for all values for each of the samples and are colored by batch membership.

<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />

## LIMMA

-------------------------------------------------------------------------------------------------
     &nbsp;        Condition: DaDaOE (logFC)   Condition: DaKD (logFC)   Condition: DaOE (logFC) 
----------------- --------------------------- ------------------------- -------------------------
 **FBgn0029708**             169.3                     -109.9                    -17.88          

 **FBgn0266521**             41.15                     -137.7                      -44           

 **FBgn0037697**            -27.44                     -151.5                    -160.1          

 **FBgn0028516**             12.03                     -28.43                     49.35          

 **FBgn0085382**             35.17                     -39.39                    -7.476          

 **FBgn0000591**             389.6                     -161.5                     957.8          

 **FBgn0000575**             63344                      -1358                     42641          

 **FBgn0003326**             53.8                      -33.96                     105.3          

 **FBgn0039883**             203.1                     -318.4                    -119.8          

 **FBgn0015558**             26.24                      -32.9                    -8.971          
-------------------------------------------------------------------------------------------------

Table: Table continues below

 
------------------------------------------------------------------------------------
     &nbsp;        Condition: ScOE (logFC)   AveExpr    F      P.Value    adj.P.Val 
----------------- ------------------------- --------- ------ ----------- -----------
 **FBgn0029708**            5872              1196     2918   2.214e-17   1.148e-13 

 **FBgn0266521**            7845              1521     2888   2.349e-17   1.148e-13 

 **FBgn0037697**            5488              1101     2686   3.585e-17   1.168e-13 

 **FBgn0028516**            2192              418.7    2002   1.983e-16   4.846e-13 

 **FBgn0085382**            1123               239     1690   5.318e-16   1.04e-12  

 **FBgn0000591**            12944             2624     1596   7.398e-16   1.205e-12 

 **FBgn0000575**            1350              18229    1384   1.696e-15   2.368e-12 

 **FBgn0003326**            1459              316.3    1272   2.777e-15   3.394e-12 

 **FBgn0039883**            10820             2311     1238   3.239e-15   3.519e-12 

 **FBgn0015558**            1206              246.4    1115   5.959e-15   5.597e-12 
------------------------------------------------------------------------------------


Median Correlations
===================
This plot helps identify outlying samples.
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />


Heatmaps
========
## Heatmap
This is a heatmap of the given data matrix showing the batch effects and variations with different conditions.
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

## Sample Correlations
This is a heatmap of the correlation between samples.
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />


Circular Dendrogram
===================
This is a Circular Dendrogram of the given data matrix colored by batch to show the batch effects.
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />


PCA: Principal Component Analysis
=================================
## PCA
This is a plot of the top two principal components colored by batch to show the batch effects.
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-20-1.png" style="display: block; margin: auto;" />

## Explained Variation

-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  &nbsp;    Proportion of Variance (%)   Cumulative Proportion of   Percent Variation Explained by   Percent Variation Explained by   Condition Significance   Percent Variation Explained by   Batch Significance (p-value) 
                                               Variance (%)           Either Condition or Batch                Condition                    (p-value)                      Batch                                             
---------- ---------------------------- -------------------------- -------------------------------- -------------------------------- ------------------------ -------------------------------- ------------------------------
 **PC1**               50.3                        50.3                           98                              97.7                          0                           49.7                           0.2648            

 **PC2**              15.33                       65.62                          96.2                             95.8                          0                           7.1                            0.3513            

 **PC3**              10.51                       76.14                          95.7                             95.7                          0                           8.9                            0.9169            

 **PC4**               5.48                       81.62                          80.6                             79.3                       0.00077                         3                             0.4144            

 **PC5**              3.392                       85.01                          13.3                             12.2                        0.7986                        0.4                            0.7179            

 **PC6**                3                         88.01                          6.8                              6.2                         0.9328                         0                             0.7944            

 **PC7**              2.372                       90.38                          1.9                              1.2                         0.9939                         0                             0.783             

 **PC8**              2.013                       92.39                          4.1                              1.2                         0.9751                        0.1                            0.5774            

 **PC9**              1.746                       94.14                          3.3                              2.9                         0.9824                         0                             0.8397            

 **PC10**             1.544                       95.68                           10                               5                          0.8703                        0.1                            0.4505            

 **PC11**             1.239                       96.92                          0.8                              0.6                         0.9991                        0.1                            0.886             

 **PC12**             0.9847                      97.91                          1.8                              1.4                         0.9947                         0                             0.832             

 **PC13**             0.8356                      98.74                           1                               0.5                         0.9983                         0                             0.8245            

 **PC14**             0.7031                      99.45                          0.1                              0.1                           1                            0                             0.9853            

 **PC15**             0.5538                       100                           0.1                              0.1                           1                            0                             0.9992            

 **PC16**           6.649e-27                      100                            86                              0.1                        0.00067                        28.1                           1e-05             

 **PC17**           1.926e-29                      100                           83.2                             79.3                       0.01252                        50.9                           0.1356            
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


Shape
=====
This is a heatmap plot showing the variation of gene expression mean, variance, skewness and kurtosis between samples grouped by batch to see the batch effects variation
<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-23-1.png" style="display: block; margin: auto;" />

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

<img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-25-1.png" style="display: block; margin: auto;" /><img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-25-2.png" style="display: block; margin: auto;" /><img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-25-3.png" style="display: block; margin: auto;" /><img src="/Users/JQ/Documents/_CODE REPOS/GitHub/Da_RNAseq/svar_batchqc_report_files/figure-html/unnamed-chunk-25-4.png" style="display: block; margin: auto;" />

```
## Batch mean distribution across genes: Normal vs Empirical distribution
## Two-sided Kolmogorov-Smirnov test
## Selected Batch: 1
## Statistic D = 0.03163
## p-value = 6.423e-09
## 
## 
## Batch Variance distribution across genes: Inverse Gamma vs Empirical distribution
## Two-sided Kolmogorov-Smirnov test
## Selected Batch: 1
## Statistic D = 0.1865
## p-value = 0Note: The non-parametric version of ComBat takes much longer time to run and we recommend it only when the shape of the non-parametric curve widely differs such as a bimodal or highly skewed distribution. Otherwise, the difference in batch adjustment is very negligible and parametric version is recommended even if p-value of KS test above is significant.
```


SVA
===
## Summary

```
## Number of Surrogate Variables found in the given data: 1
```
