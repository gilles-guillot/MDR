Fit_Single_Hill <- function(x,y,n,models=c('Hill2','Hill3','Hill4'),
                            lower = c(0.0001,0.0001,0.0001,0.0001) ,
                            upper = c(10,10,1,1) ,
                            bootstrap=FALSE, B=100,
                            select=FALSE,
                            eps = 1e-15)
{
  #' @title  Estimate parameters in   Hill model families
  #' @description
  #' Estimate parameters in single-compound dose-reponse models under some of the Hill model families (2,3,4 parameters)
  #'
  #' Parameters
  #' @param x sequence of dose input
  #' @param y sequence of responses
  #' @param n number of binomial samples at various dose values
  #' @param models one of 2-4 parameter Hill function
  #' @param lower lower limits for the search of parameter
  #' @param upper upper limits for the search of parameter
  #' @param bootstrap logical resampling for estimating parameter variability
  #' @param B number of bootstrap resampling
  #' @param select TRUE/FALSE perform model selection
  #' @param eps used to replace 0's by small numerical values
  #'
  #' @return
  #' Estimate of Hill model(s) parameters and bootstrap sample for this estimate
  #' Perform model selection by cross-validation
  #' Return estimate of 50% effect doses and bootstrap sample for this estimate
  #' NB:  50% effect dose may not exist under Hill4 for some combinations (c,d)
  #'
  #' Detail
  #' Model assumed: y ~ B(n,f(x))
  #' Hill models:
  #' Hill 4: f(x) = d + (c-d) / ( 1+ (a/x)**b)
  #' Hill 3: d=0
  #' Hill 2: d=0 ; c=1
  #' ML estimation, optimization of log likelihood with Generalized Simulated Annealing (GenSA::GenSA)
  #'
  #' 0 doses bring numerical issues
  #' 0's in x and y replaced by min(x) * 1e-9 / min(y) * 1e-9
  #'
  #' Model selection with leave-"some"-out likelihood cross-validation
  #' "some" = all obs receiving a certain dose level
  #'
  #' @author Gilles Guillot \email{gilles.b.guillot@gmail.com}
  #'
  #' @importFrom stats rbinom
  #' @importFrom GenSA GenSA

  llcv = par_hat = ed50_hat = par_boot = ed50_boot = list(Hill2=NA,Hill3=NA,Hill4=NA)
  ## negative log-likelihoods
  nll2 <-function(par,x,y,n)
  {
    a = par[1];   b = par[2]
    p = 1/(1+ (a/x)**b)
    -sum(y*log(p) + (n-y)*log(1-p))
  }
  ##
  nll3 <-function(par,x,y,n)
  {
    a = par[1] ;  b = par[2] ;   c = par[3]
    p = c/( 1 + (a/x)**b)
    -sum(y*log(p) + (n-y)*log(1-p))
  }
  ##
  nll4 <-function(par,x,y,n)
  {
    a = par[1];   b = par[2] ;   c = par[3] ;   d = par[4]
    p = d + (c-d)/( 1 + (a/x)**b)
    -sum(y*log(p) + (n-y)*log(1-p))
  }
  ## END negative log-likelihoods

  ## setting 0 values to eps to avoid numerical issues
  x[x==0] = min(x[x>0]) * eps
  y[y==0] = min(y[y>0]) * eps


  ## Estimation
  if('Hill2' %in% models){
    ## Estimation Hill2
    # print('Estimation Hill2')
    res = GenSA(c(1,1),fn=nll2,lower=lower[1:2],upper=upper[1:2],x=x,y=y,n=n)
    par_hat$Hill2 = res$par[1:2]
    ed50_hat$Hill2 = res$par[1]
    if(bootstrap){
      print('Bootstrap Hill2')
      par2_boot = matrix(nrow =B,ncol=2)
      ed50_2_boot = rep(NA,B)
      for(iboot in 1:B){
        print(iboot)
        y_boot = rbinom(n=length(x),size=n,prob=y/n)
        res = GenSA(c(1,1),fn=nll2,lower=lower[1:2],upper=upper[1:2],
                    x=x,y=y_boot,n=n)
        par2_boot[iboot,] = res$par[1:2]
        ed50_2_boot[iboot] = res$par[1]
      }
      par_boot$Hill2 = par2_boot
      ed50_boot$Hill2 = ed50_2_boot
    }
  }
  if('Hill3' %in% models){
    ## Estimation Hill3
    print('Estimation Hill3')
    res = GenSA(c(1,1,1),fn=nll3,lower=lower[1:3],upper=upper[1:3],x=x,y=y,n=n)
    par_hat$Hill3 = res$par[1:3]
    ed50_hat$Hill3 = res$par[1]
    if(bootstrap){
      print('Bootstrap Hill3')
      par3_boot = matrix(nrow =B,ncol=3)
      ed50_3_boot = rep(NA,B)
      for(iboot in 1:B){
        print(iboot)
        y_boot = rbinom(n=length(x),size=n,prob=y/n)
        res = GenSA(c(1,1,1),fn=nll3,lower=lower[1:3],upper=upper[1:3],
                    x=x,y=y_boot,n=n)
        par3_boot[iboot,] = res$par[1:3]
        ed50_3_boot[iboot] = res$par[1]
      }
      par_boot$Hill3 = par3_boot
      ed50_boot$Hill3 = ed50_3_boot
    }
  }
  if('Hill4' %in% models){
    ## Estimation Hill4
    print('Estimation Hill4')
    res = GenSA(c(1,1,1,0),fn=nll4,lower=lower[1:4],upper=upper[1:4],x=x,y=y,n=n)
    par_hat$Hill4 = res$par[1:4]
    ed50_hat$Hill4 = res$par[1] * ( (res$par[3] - 2*res$par[4]) / res$par[3] )**(1/res$par[2])
    if(bootstrap){
      print('Bootstrap Hill4')
      par4_boot = matrix(nrow =B,ncol=4)
      ed50_4_boot = rep(NA,B)
      for(iboot in 1:B){
        print(iboot)
        y_boot = rbinom(n=length(x),size=n,prob=y/n)
        res = GenSA(c(1,1,1,0),fn=nll3,lower=lower[1:4],upper=upper[1:4],
                    x=x,y=y_boot,n=n)
        par4_boot[iboot,] = res$par[1:4]
        ed50_4_boot[iboot] = res$par[1] * ( (res$par[3] - 2*res$par[4]) / res$par[3] )**(1/res$par[2])
      }
      par_boot$Hill4 = par4_boot
      ed50_boot$Hill4 = ed50_4_boot
    }
  }
  ## Model selection
  if(select){
    if('Hill2' %in% models){
      print('CV Hill2')
      llcv2 = 0
      for(k in 1:length(x)){
        res = GenSA(c(1,1),fn=nll2,lower=lower[1:2],upper=upper[1:2],
                    x=x[-k],y=y[-k],n=n[-k])
        llcv2 = llcv2 - nll2(res$par[1:2],x[k],y[k],n[k])
        # print(llcv2)
      }
      llcv$Hill2=llcv2
    }
    if('Hill3' %in% models){
      print('CV Hill3')
      llcv3 = 0
      for(k in 1:length(x)){
        res = GenSA(c(1,1,1),fn=nll3,lower=lower[1:3],upper=upper[1:3],
                    x=x[-k],y=y[-k],n=n[-k])
        llcv3 = llcv3 - nll3(res$par[1:3],x[k],y[k],n[k])
      }
      llcv$Hill3=llcv3
    }
    if('Hill4' %in% models){
      print('CV Hill4')
      llcv4 = 0
      for(k in 1:length(x)){
        res = GenSA(c(1,1,1,0),fn=nll4,lower=lower[1:4],upper=upper[1:4],
                    x=x[-k],y=y[-k],n=n[-k])
        llcv4 = llcv4 - nll4(res$par[1:4],x[k],y[k],n[k])
        # print(llcv4)
      }
      llcv$Hill4=llcv4
    }
  } ## END model selection
  res = list( par_hat = par_hat,  ed50_hat = ed50_hat ,
              par_boot = par_boot, ed50_boot = ed50_boot ,
              llcv = llcv)
  if(select){ res = append(res,
                           list(best.model = c('Hill2','Hill3','Hill4')[which.max(llcv)]))}
  return(res)
  }## END function Fit_Single_Hill



