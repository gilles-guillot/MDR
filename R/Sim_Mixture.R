## Simulate toxicological/pharmacological mixture experiment data
#' @title Simulate toxicological/pharmacological mixture experiment data
#' @description Simulate toxicological/pharmacological mixture experiment data
#'
#' @param x sequence of doses first compound
#' @param y sequence of doses second compound
#' @param n number of binomial samples at various dose values
#' @param par_Hill parameters of the Hill function
#' @param interaction_model 'Loewe','Bliss','Hewlett','Voelund' or 'JonkerSA'
#' @param par_int parameter for interaction model
#' @param distribution distribution of response conditional on Hill function value
#' @param par_dist param of the previous distribution
#'
#' @return list with input parameters and simulated responses
#'
#' @importFrom stats rlnorm uniroot
Sim_Mixture <- function(x,y,n,par_Hill,interaction_model,par_int=NULL,distribution,par_dist=1)
{
  if(!(distribution %in% c('binom','lnorm'))) stop('distribution not recognized')
  if(!(interaction_model %in% c('Loewe','Bliss','Hewlett','Voelund','JonkerSA'))) stop('interaction model not recognized')

  a1 = par_Hill[1]
  b1 = par_Hill[2]
  a2 = par_Hill[3]
  b2 = par_Hill[4]
  h = rep(NA,length(x)) # h[i] mean response value for combination (x[i] , y[i])
  for(i in 1:length(x)) # compute bivariate DR function at combination (x[i] , y[i])
  {
    if(x[i]==0 & y[i]==0){ h[i] = 0 }
    if(x[i]!=0 & y[i]==0){ h[i] = Hill(x=x[i],a=a1,b=b1,c=1,d=0) }
    if(x[i]==0 & y[i]!=0){ h[i] = Hill(x=y[i],a=a2,b=b2,c=1,d=0) }
    if(x[i]!=0 & y[i]!=0)
    {
      if(interaction_model == 'Bliss')
      { h[i] = DR_Bliss_Hill2(d1=x[i],d2=y[i],a1=a1,b1=b1,a2=a2,b2=b2) }
      if(interaction_model == 'Loewe')
        { h[i] = DR_Voelund_Hill2(d1=x[i],d2=y[i],a1=a1,b1=b1,a2=a2,b2=b2,lambda1=1,lambda2=1) }
      if(interaction_model == 'Hewlett')
        { h[i] = DR_Voelund_Hill2(d1=x[i],d2=y[i],a1=a1,b1=b1,a2=a2,b2=b2,lambda1=par_int[1],lambda2=par_int[1]) }
      if(interaction_model == 'Voelund')
        { h[i] = DR_Voelund_Hill2(d1=x[i],d2=y[i],a1=a1,b1=b1,a2=a2,b2=b2,
                                                                   lambda1=par_int[1],lambda2=par_int[2]) }
      if(interaction_model == 'JonkerSA')
        { h[i] = DR_JonkerSA_Hill2(d1=x[i],d2=y[i],a1=a1,b1=b1,a2=a2,b2=b2,alpha=par_int[1]) }
    }

  }
  if(distribution == 'binom')  { z = rbinom(n= length(x), size= n, prob= h) }
  if(distribution == 'lnorm')
    {
    sdlog = par_dist
    z = rlnorm(n= length(x), meanlog = h, sdlog= sdlog)
    }

  return(list(x=x,y=y,n=n,h=h,z=z,par_Hill,interaction_model,par_int,distribution,par_dist))

}
