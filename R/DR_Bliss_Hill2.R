#' @title Bivariate dose-response function under Bliss independance model
#' @description Compute bivariate dose-response function when each mixture compoent follows a 2-parameter Hill dose-response function and the Bliis independance  model holds.
#' @param d1 dose of first compound: a numeric vector
#' @param d2 dose of second compound: a numeric vector
#' @param a1 parameter a in Hill-2 DR fuction 1/(1+(a/x)**b) for first compound
#' @param b1 parameter b in Hill-2 DR fuction 1/(1+(a/x)**b) for first compound
#' @param a2 parameter a in Hill-2 DR fuction 1/(1+(a/x)**b) for second compound
#' @param b2 parameter b Hill-2 DR fuction in 1/(1+(a/x)**b) for second compound#'
#' @return vector of numerical responses
#' @export
DR_Bliss_Hill2 <- function(d1,d2,a1,b1,a2,b2)
{
  f = Hill(x=d1,a=a1,b=b1,c=1,d=0)
  g = Hill(x=d2,a=a2,b=b2,c=1,d=0)
  return(f+g-f*g)
}
