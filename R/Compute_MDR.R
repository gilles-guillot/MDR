Compute_MDR <- function(x,y,z,n,models=c('Hill2'),design,
                lower = c(0.0001,0.0001,0.0001,0.0001) ,
                upper = c(10,10,1,1) ,
                do.plot=FALSE # plot fit of DR curve for the two-compound data
                )
{
  #' @title Compute MDR
  #' @description computing Model Deviation Ratio for data obtained under three types of designs
  #' @param x sequence of doses first compound
  #' @param y sequence of doses second compound
  #' @param z sequence of responses
  #' @param n number of binomial samples at various dose values
  #' @param models family of the single compound dose-response curve
  #' @param design any of the three below
  #' 'single-parallel': single ray parallel to an axis (dose of one of the two compounds is fixed), single compound experiment data not used,
  #' 'single-ray': single ray design (ratio of doses of two compounds is fixed), single compound experiment data not used
  #' TODO 'general': fully arbitrary design, all data used including  single compound experiment data
  #' Each of the design above has to contain single-compound experiment data
  #' @param do.plot logical
  #' @param lower lower limits for the search of parameter
  #' @param upper upper limits for the search of parameter
  #'
  #' @return estimates of doses or concentrations eliciting 50\% max effect and of MDR
  #'
  #' @importFrom graphics  abline lines plot

  if(! (design %in% c('single-parallel','single-ray','general'))) stop('design not recognized')

  if( (design == 'single-parallel') | (design == 'single-ray') )
  {
    ## compound 1 alone
    subs <- y==0
    res = Fit_Single_Hill(x=x[subs],y=z[subs],n=n[subs],models=c('Hill2'),
                          lower = lower,
                          upper = upper ,
                          bootstrap=FALSE, B=10,
                          select=FALSE)
    X50_hat = res$ed50_hat$Hill2
    ## compound 2 alone
    subs <- x==0
    res = Fit_Single_Hill(x=y[subs],y=z[subs],n=n[subs],models=c('Hill2'),
                         lower  = lower ,
                          upper = upper ,
                          bootstrap=FALSE, B=10,
                          select=FALSE)
    Y50_hat = res$ed50_hat$Hill2
    if(design == 'single-parallel')
    {
      ## Find out which of x and y varies when the other remains constant
      subs <- (x>0 & y>0)
      xtmp <- x[subs]
      ytmp <- y[subs]
      which.var <- which.max(c(length(unique(xtmp)),length(unique(ytmp))))
      ## find constant dose
      if(which.var == 1){ d_fixed = unique(ytmp) } # x varies
      if(which.var == 2){ d_fixed = unique(xtmp) } # y varies

      ## single compound parameter estimated plugged in Loewe equation
      if(which.var == 1){  d50_DA =  X50_hat * (1- d_fixed/Y50_hat) }
      if(which.var == 2){  d50_DA =  Y50_hat * (1- d_fixed/X50_hat) }

      ## remove single-compound data
      ## we may have removed a data point (x,0) or (0,y)
      if(which.var == 1)
      {
        subs <- (y == unique(ytmp))
        d <- x[subs]
        ntmp <- n[subs]
      }
      if(which.var == 2)
      {
        subs <- (x == unique(xtmp))
        d <- y[subs]
        ntmp <- n[subs]
      }
      res = Fit_Single_Hill(x=d,y=z[subs],n=n[subs],models=c('Hill2'),
                           lower  = lower ,
                            upper = upper ,
                            bootstrap=FALSE, B=10,
                            select=FALSE)
      d50_obs = res$ed50_hat$Hill2
      if(do.plot)
      {
        plot(d,z[subs]/n[subs])
        dd = seq(min(d),max(d),l=1000)
        lines(dd,1/(1+(res$par_hat$Hill2[1]/dd)**res$par_hat$Hill2[2]))
        abline(h=0.5)
        abline(v=d50_DA)
        abline(v=d50_obs)
      }
    }
    if(design == 'single-ray')
    {
      ## remove single-compound data
      subs <- (x>0 & y>0) | (x==0 & y==0)
      xtmp <- x[subs]
      ytmp <- y[subs]
      ztmp <- z[subs]
      ntmp <- n[subs]
      ## find constant ratio a = x/y
      ss = !is.na(xtmp/ytmp)
      a = mean(unique(xtmp[ss]/ytmp[ss]))
      ## d50_DA: dose x in combination with y=x/a eliciting 50% effect under DA
      d50_DA = 1/(1/X50_hat + 1/(a*Y50_hat))
      res = Fit_Single_Hill(x=xtmp,y=ztmp,n=ntmp,models=c('Hill2'),
                           lower  = lower ,
                            upper = upper ,
                            bootstrap=FALSE, B=10,
                            select=FALSE)
      d50_obs = res$ed50_hat$Hill2
      if(do.plot)
      {
        plot(xtmp,ztmp/ntmp)
        dd = seq(min(xtmp),max(xtmp),l=1000)
        lines(dd,1/(1+(res$par_hat$Hill2[1]/dd)**res$par_hat$Hill2[2]))
        abline(h=0.5)
        abline(v=d50_DA)
        abline(v=d50_obs)
      }
    }
    MDR =  d50_DA / d50_obs

  }
  ## res = list(X50_hat=X50_hat, Y50_hat=Y50_hat, d50_DA=d50_DA, d50_obs = d50_obs, MDR = MDR)
  res <- c(X50_hat, Y50_hat, d50_DA, d50_obs, MDR)
  names(res) <-  c('X50_hat', 'Y50_hat', 'd50_DA', 'd50_obs', 'MDR')

  return(res)
}

