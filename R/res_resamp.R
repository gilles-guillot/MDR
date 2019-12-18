#' Object for illustrating documentation, cf vignette
#' @format a list as produced by function `Resample_MDR`
#' @source simulated by example below
#' @examples
#'  \dontrun{
#'  res_resamp = Resample_MDR(x=simdat$x,y=simdat$y,
#'  z=simdat$z,n=simdat$n,
#'  SCDR_model='Hill2',
#'  par_Hill=NULL,
#'  design='single-ray',
#'  lower = c(0.0001,0.0001,0.0001,0.0001) ,
#'  upper = c(10,10,1,1) ,
#'  n_resamp=1000,
#'  distribution='binom')
#'  }
"res_resamp"
