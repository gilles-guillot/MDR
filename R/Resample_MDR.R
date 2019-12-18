#' @title Perform Monte Carlo simulation of MDR statistic
#' @description Perform Monte Carlo simulation of MDR statistic
#' under generative model and experimental designed known/assumed
#' for the data
#' @details Input data should include design (x,y,n) AND
#' EITHER vector of response mortality values z
#' OR estimate of dose repose curves parameters par_Hill
#' @param x sequence of doses first compound
#' @param y sequence of doses second compound
#' @param z sequence of responses
#' @param n number of binomial samples at various dose values
#' @param SCDR_model family of the single compound dose-response curve
#' @param par_Hill parameters of the Hill functions.
#' With par_Hill = (a1,b1,a2,b2)
#' a1: a in 1/(1+(a/x)**b) for first compound
#' b1: b in  1/(1+(a/x)**b) for first compound
#' a2: a in  1/(1+(a/x)**b) for second compound
#' b2: b in 1/(1+(a/x)**b) for second compound
#' @param design  any of the below
#' 'single-parallel': single ray parallel to an axis (dose of one of the two compounds is fixed), single compound experiment data not used,
#' 'single-ray': single ray design (ratio of doses of two compounds is fixed), single compound experiment data not used
#' 'general' (TODO): fully arbitrary design, all data used including  single compound experiment data
#' @param lower lower limits for the search of parameter
#' @param upper upper limits for the search of parameter
#' @param n_resamp number of dataset resampled
#' @param distribution distribution of response conditional on Hill function value
#' @param verbose will print the percentage of computation carried out if TRUE
#' @return a list with an element named MDR_resamp
#' @export
Resample_MDR <- function(x,y,z=NULL,n,
                         SCDR_model='Hill2',
                         par_Hill=NULL,
                         design,
                         lower = c(0.0001,0.0001,0.0001,0.0001) ,
                         upper = c(10,10,1,1) ,
                         n_resamp, distribution='binom',verbose=TRUE)
{
  if( is.null(z) & is.null(par_Hill) )
    stop("Please provide at z or par_Hill")
  ## estimating parameters of single-compond DR curves
  if(is.null(par_Hill))
  {
    subs <- y==0
    res_x = Fit_Single_Hill(x=x[subs],y=z[subs],n=n[subs],
                            SCDR_model=SCDR_model,
                            lower = lower,
                            upper = upper ,
                            bootstrap=FALSE, B=10,
                            select=FALSE)
    subs <- x==0
    res_y = Fit_Single_Hill(x=y[subs],y=z[subs],n=n[subs],
                            SCDR_model=SCDR_model,
                            lower  = lower ,
                            upper = upper ,
                            bootstrap=FALSE, B=10,
                            select=FALSE)
    par_Hill = c(res_x$par_hat$Hill2,res_y$par_hat$Hill2)
  }
  ## Loop on simulated datasets
  MDR_resamp = rep(NA,n_resamp)
  for(isamp in 1:n_resamp)
  {
    if(verbose) { Sys.sleep(0.07); cat("\r","# Percentage resampling iterations:",
                         signif(100*isamp/n_resamp,digits=3),"%") }
    simdat = Sim_Mixture(x=x,y=y,n=n,
                         SCDR_model=SCDR_model,
                         par_Hill=par_Hill,
                         interaction_model='Loewe',
                         par_int=NULL,
                         distribution=distribution,par_dist=1)
    res_MDR = Compute_MDR(x=simdat$x,y=simdat$y,
                          z=simdat$z,n=simdat$n,
                          SCDR_model=SCDR_model,
                          design=design,
                          lower = lower ,
                          upper = upper ,
                          do.plot=FALSE )
    MDR_resamp[isamp] = res_MDR$MDR
  }
  return(list(MDR_resamp=MDR_resamp))
}
