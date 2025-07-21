# crossRanger

**crossRanger: A Flexible Random Forest Toolkit for Microbiome Meta-Analysis**


## Introduction

crossRanger is an R package that extends the ranger engine (a fast implementation of random forests [ranger](https://github.com/imbs-hl/ranger) ) to enable classification or regression on a single target variable across multiple microbiome datasets. It supports:

* Building and validating models within or across studies.

* Transferring pre-trained models to new datasets.

* Streamlining meta-analyses of microbiome–phenotype relationships.

**Key Features**

1. Unified Random Forest Framework

* Supports classification (e.g., disease status) and regression (e.g., age) tasks.

* Handles high-dimensional microbiome data (e.g., taxa, pathways) efficiently.

2. Cross-Dataset Analysis

* Train a model on one dataset and predict the same target in others.

* Compare performance metrics (e.g., AUROC, RMSE) across studies.

3. Stratified Analysis

Account for covariates like sex or body site via stratified modeling.

4. Performance Evaluation

* Classification: AUROC, AUPRC, accuracy.

* Regression: MAE, MSE, R², MAPE.

5. Downstream Tools

* Identify microbial associations (e.g., Wilcoxon tests on CLR-transformed abundances).

* Generate publication-ready visualizations.

## Authors ##
Shi Huang, UC San Diego 

## Installation ##

The development version is maintained on GitHub and can be downloaded as follows:
``` r 
## install.packages('devtools') # if devtools not installed
devtools::install_github('shihuang047/crossRanger')
```

## Examples ##
* Compositional microbiome data simulation using multinomial distribution
``` r
set.seed(123)
df <- data.frame(rbind(t(rmultinom(7, 75, c(.201,.5,.02,.18,.099))),
            t(rmultinom(8, 75, c(.201,.4,.12,.18,.099))),
            t(rmultinom(15, 75, c(.011,.3,.22,.18,.289))),
            t(rmultinom(15, 75, c(.091,.2,.32,.18,.209))),
            t(rmultinom(15, 75, c(.001,.1,.42,.18,.299)))))
metadata<-data.frame(f_s=factor(c(rep("A", 15), rep("B", 15), rep("C", 15), rep("D", 15))),
                     f_c=factor(c(rep("C", 7), rep("H", 8), rep("C", 7), rep("H", 8),
                                  rep("C", 7), rep("H", 8), rep("C", 7), rep("H", 8))),
                     f_d=factor(rep(c(rep("a", 5), rep("b", 5), rep("c", 5)), 4)))
```
* RF modeling on one responsive variable
``` r
rf.out.of.bag(x=df, y=metadata$f_s, imp_pvalues=FALSE)
```
* RF modeling on one responsive variable stratified by another covariate
``` r
res_list<-rf_clf.by_datasets(df, metadata, s_category='f_s', c_category='f_c', positive_class="C")
summary(res_list)
```
*  Apply RF models to other datasets
``` r 
cross_rf<-rf_clf.cross_appl(rf_model_list=res_list$rf_model_list, 
                            x_list=res_list$x_list, 
                            y_list=res_list$y_list, positive_class="C")
cross_rf$perf_summ
# the output
Train_data Test_data   Validation_type  Accuracy       AUC       Kappa Sensitivity Specificity Pos_Pred_Value Neg_Pred_Value
1           A         A   self_validation 0.9333333 1.0000000  0.86725664   1.0000000       0.875      0.8750000      1.0000000
2           A         B cross_application 0.4666667 0.3303571  0.00000000   1.0000000       0.000      0.4666667            NaN
3           A         C cross_application 0.4666667 0.6964286  0.00000000   1.0000000       0.000      0.4666667            NaN
4           A         D cross_application 0.4666667 0.4285714  0.00000000   1.0000000       0.000      0.4666667            NaN
5           B         B   self_validation 0.5333333 0.4464286  0.05405405   0.4285714       0.625      0.5000000      0.5555556
6           B         A cross_application 0.4666667 0.4017857  0.00000000   1.0000000       0.000      0.4666667            NaN
7           B         C cross_application 0.4666667 0.3482143  0.00000000   1.0000000       0.000      0.4666667            NaN
8           B         D cross_application 0.4666667 0.4553571 -0.03448276   0.7142857       0.250      0.4545455      0.5000000
9           C         C   self_validation 0.5333333 0.3571429  0.05405405   0.4285714       0.625      0.5000000      0.5555556
10          C         A cross_application 0.4666667 0.4642857 -0.11111111   0.1428571       0.750      0.3333333      0.5000000
11          C         B cross_application 0.4666667 0.6071429  0.00000000   1.0000000       0.000      0.4666667            NaN
12          C         D cross_application 0.4000000 0.5267857 -0.13445378   0.8571429       0.000      0.4285714      0.0000000
13          D         D   self_validation 0.4000000 0.3214286 -0.21621622   0.2857143       0.500      0.3333333      0.4444444
14          D         A cross_application 0.5333333 0.4464286  0.10256410   0.8571429       0.250      0.5000000      0.6666667
15          D         B cross_application 0.6000000 0.5267857  0.19642857   0.5714286       0.625      0.5714286      0.6250000
16          D         C cross_application 0.4000000 0.3482143 -0.26168224   0.0000000       0.750      0.0000000      0.4615385
   Precision    Recall        F1 Prevalence Detection_Rate Detection_Prevalence Balanced_Accuracy
1  0.8750000 1.0000000 0.9333333  0.4666667     0.46666667            0.5333333         0.9375000
2  0.4666667 1.0000000 0.6363636  0.4666667     0.46666667            1.0000000         0.5000000
3  0.4666667 1.0000000 0.6363636  0.4666667     0.46666667            1.0000000         0.5000000
4  0.4666667 1.0000000 0.6363636  0.4666667     0.46666667            1.0000000         0.5000000
5  0.5000000 0.4285714 0.4615385  0.4666667     0.20000000            0.4000000         0.5267857
6  0.4666667 1.0000000 0.6363636  0.4666667     0.46666667            1.0000000         0.5000000
7  0.4666667 1.0000000 0.6363636  0.4666667     0.46666667            1.0000000         0.5000000
8  0.4545455 0.7142857 0.5555556  0.4666667     0.33333333            0.7333333         0.4821429
9  0.5000000 0.4285714 0.4615385  0.4666667     0.20000000            0.4000000         0.5267857
10 0.3333333 0.1428571 0.2000000  0.4666667     0.06666667            0.2000000         0.4464286
11 0.4666667 1.0000000 0.6363636  0.4666667     0.46666667            1.0000000         0.5000000
12 0.4285714 0.8571429 0.5714286  0.4666667     0.40000000            0.9333333         0.4285714
13 0.3333333 0.2857143 0.3076923  0.4666667     0.13333333            0.4000000         0.3928571
14 0.5000000 0.8571429 0.6315789  0.4666667     0.40000000            0.8000000         0.5535714
15 0.5714286 0.5714286 0.5714286  0.4666667     0.26666667            0.4666667         0.5982143
16 0.0000000 0.0000000       NaN  0.4666667     0.00000000            0.1333333         0.3750000
```
![cross_appl_MAE](cross_appl_MAE_plot.png)
## References ##
* Breiman, L. (2001). Random forests. Mach Learn, 45:5-32. https://doi.org/10.1023/A:1010933404324.
* Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. https://doi.org/10.18637/jss.v077.i01.

## License ##
All source code freely availale under [GPL-3 License](https://www.gnu.org/licenses/gpl-3.0.en.html). 

## Documentation ##
Each exported function in the package has been documented and we have also written an introductory vignette that is accessible by calling 
``` r
vignette('crossRanger--intro', package='crossRanger')
```

## Bugs/Feature requests ##
I appreciate bug reports and feature requests. Please post to the github issue tracker [here](https://github.com/shihuang047/crossRanger/issues). 

## Acknowledgements

 This work is supported by IBM Research AI through the AI Horizons Network. For
 more information visit the [IBM AI Horizons Network website](https://www.research.ibm.com/artificial-intelligence/horizons-network/).

