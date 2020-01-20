# MDR
The MDR package provides functions to test departure from Loewe additivity for binary toxicological or pharmacological mixtures with the Model Deviation Ratio (MDR) method.
It provides tools for computing the model deviation ratio (MDR) statistic and assessing the significance 
of the observed departure from Loewe additivity. 

It also provides tools for the analysis of single compound dose-response experiment data. 
It is a companion package to the article “Dose-Response Analysis of Toxicological and Pharmacological Mixtures with the Model Deviation Ratio Method: Problems and
Solutions” by Macacu, A., S. Lasagni, E. Carnesecchi, J.L Dorne, and G. Guillot (2020). 

Windows users should first install Rtools: Rtools  https://cran.r-project.org/bin/windows/Rtools/

For installation of the MDR package under R: 

`if( ! 'devtools' %in% installed.packages() )
  { install.packages('devtools') }`
`devtools::install_github('gilles-guillot/MDR', build_vignettes = TRUE)`

Load the package with: 
`library(MDR)`

And explore its feature with:
`vignette('MDR')`
