#' @title  Hill dose response function
#' @description computation of 2-4 parameter Hill dose response function
#'
#' @param x dose (vector)
#' @param a a in d + (c-d)/(1+(a/x)**b)  ; 50\% max effect
#' @param b b in  d + (c-d)/(1+(a/x)**b)
#' @param c c in d + (c-d)/(1+(a/x)**b)
#' @param d d in d + (c-d)/(1+(a/x)**b)
#'
#' @return response
#' @export
Hill <- function(x,a,b,c=1,d=0)
{  d + (c-d)/(1+(a/x)**b)  }


#' @title reciprocal of 2-paramter Hill function
#' @description computes doses elicting effects y when dose reponse function is Hill2
#'
#' @param y response (numeric vector)
#' @param a a a in d + (c-d)/(1+(a/x)**b)  ; 50\% max effect
#' @param b b b in  d + (c-d)/(1+(a/x)**b)
#'
#' @return doses elicting effects y
Hill2_inv <- function(y,a,b)
{  a*(y/(1-y))**(1/b)  }
