## Some R commandes to help building the package

options(usethis.protocol = "https")

library(devtools)
library(roxygen2)
library(usethis)

setwd("/home/gilles/Dropbox/iPRI/projects/ToxPharm/MDR_package/")
getwd()

## create Rstudio project and packahe structure
create('MDR')

## make folder a git project
use_git()

## connect local project ot github repo
use_github(auth_token = "37e8f5fd3b4bf0f1a407a694602c1f47b2122b43")


document()

