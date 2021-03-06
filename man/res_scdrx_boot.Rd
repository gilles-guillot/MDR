% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/res_scdrx_boot.R
\docType{data}
\name{res_scdrx_boot}
\alias{res_scdrx_boot}
\title{Object for illustrating documentation, cf vignette}
\format{a list as produced by function `Fit_Single_Hill`}
\source{
simulated by example below
}
\usage{
res_scdrx_boot
}
\description{
Object for illustrating documentation, cf vignette
}
\examples{
 \dontrun{
 subs = simdat$y==0
res_scdrx_boot = Fit_Single_Hill(x=simdat$x[subs],
                            y=simdat$z[subs], ## z contains the response
                            n=simdat$n[subs],
                            SCDR_model='Hill2',
                            bootstrap = TRUE,B=1000)
 }
}
\keyword{datasets}
