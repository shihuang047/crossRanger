# crossRanger
*The cross-application and visualization of random forests (ranger) models of compositional microbiome data*

***
## Introduction
crossRanger = Cross-application of Ranger models.

[ranger](https://github.com/imbs-hl/ranger) is a fast implementation of random forests (Breiman 2001) or recursive partitioning, particularly suited for high dimensional data. Classification, regression, and survival forests are supported. 

This R package provides functions for the Random Forests (ranger) modeling of multiple microbiome datasets. Specifically, this package allows cross-application and comparison of Random Forest models where the models will be cross-applied, the model performance in both training and testing will be analysis and important features of multiple models will be selected respectively for cross-datasets comparisons and visualization. 

## Authors ##
Shi Huang, UC San Diego 

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

## Installation ##

The development version is maintained on GitHub and can be downloaded as follows:
``` r 
## install.packages('devtools') # if devtools not installed
devtools::install_github('shihuang047/crossRanger')
```


## Bugs/Feature requests ##
I appreciate bug reports and feature requests. Please post to the github issue tracker [here](https://github.com/shihuang047/crossRanger/issues). 
