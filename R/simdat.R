#' Documentation for simdat
#' @format a list as produced by function `Sim_Mixture`
#' @source simulated by example below
#' @examples
#'  \dontrun{
#'  a1 = 5 ; b1 = .3
#'  a2 = 7 ; b2 = .7
#
#'  d =  c(0,.4,.8,1.5,3,6,12,20)
#'  w = 1
#'  xy = cbind(w*d,(1-w)*d)
#'  w = 0
#'  xy = rbind(xy,cbind(w*d,(1-w)*d))
#'  w = .7
#'  xy = rbind(xy,cbind(w*d,(1-w)*d))
#
#'  x = xy[,1]
#'  y = xy[,2]
#'  n = rep(36,length(x))
#'
#'  simdat = Sim_Mixture(x=x,y=y,n=n,
#'                        par_Hill=c(a1,b1,a2,b2),
#'                        interaction_model='Loewe',
#'                        par_int=NULL,
#'                        distribution='binom')
#'  plot(simdat$x,simdat$y,
#'  xlab='dose compound 1',ylab='dose compound 2',
#'  sub='Dose combinations')
#'                        }
"simdat"
