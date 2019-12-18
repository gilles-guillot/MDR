#' Object for illustrating documentation, cf vignette
#' @format a list as produced by function `Fit_Single_Hill`
#' @source simulated by example below
#' @examples
#'  \dontrun{
#'  subs = simdat$y==0
#' res_scdrx_boot = Fit_Single_Hill(x=simdat$x[subs],
#'                             y=simdat$z[subs], ## z contains the response
#'                             n=simdat$n[subs],
#'                             SCDR_model='Hill2',
#'                             bootstrap = TRUE,B=1000)
#'  }
"res_scdrx_boot"
