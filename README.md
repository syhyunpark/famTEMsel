## famTEMsel (R package version 0.1.0)
Functional Additive models for Treatment Effect-Modifier Selection 

An implementation of a *functional additive regression model* which is uniquely modified and constrained to model nonlinear interaction effects between a categorical treatment variable and a potentially large number of functional/scalar pretreatment covariates on their effects on a scalar-valued outcome. The model generalizes functional additive models by incorporating treatment-specific components into additive effect components, however, a structural constraint is imposed on the treatment-specific components, to give a class of orthogonal main and interaction effect additive models. If primary interest is in interactions, one can avoid estimating main effects, obviating the need to specify their form and thereby avoiding the issue of model misspecification. Utilizing this orthogonality, one can effectively conduct treatment effect-modifier variable selection. The selected covariates can be used to make individualized treatment recommendations. We refer to Park, Petkova, Tarpey, and Ogden (2020) <doi:10.1016/j.jspi.2019.05.008> and Park, Petkova, Tarpey, and Ogden (2020) "Constrained functional additive models for estimating interactions between a treatment and functional covariates" (pre-print) for detail of the method. The main function of this package is famTEMsel(). 


#### Description

* famTEMsel - `famTEMsel` main function
* cv.famTEMsel - `famTEMsel` cross-validation function for tuning parameter selection 
* predict_famTEMsel - `famTEMsel` prediction function
* make_ITR_famTEMsel - make individualized treatment recommendations (ITRs) based on a `famTEMsel` object
* plot_famTEMsel -  plot component functions from a `famTEMsel` object 


#### To run: 
To install an R package, start by installing the "devtools" package (from CRAN). On R, type: 
```
install.packages("devtools")  # install the devtools package from CRAN
library(devtools)
```

To install the "famTEMsel" package from github, type: 
```
devtools::install_github("syhyunpark/famTEMsel")  # install the famTEMsel package from github 
library(famTEMsel)   # load the samTEMsel package to R 
```

To see some of the example codes appearing in the "help" menu, type:  
```
?famTEMsel   
?cv.famTEMsel
```

