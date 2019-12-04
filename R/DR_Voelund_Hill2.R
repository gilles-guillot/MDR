#' @title Bivariate dose-response function under Voelund model
#' @description Compute bivariate dose-response function when each mixture compoent follows a 2-parameter Hill dose-response function and the Voelund model holds.
#' @references Eq. 6  p. 388 in Soerensen et al. Environ. Ecol. Stat. 2007
#'
#' @param d1 dose of first compound: a single numeric entry
#' @param d2 dose of second compound: a single numeric entry
#' @param a1 parameter a in Hill-2 DR fuction 1/(1+(a/x)**b) for first compound
#' @param b1 parameter b in Hill-2 DR fuction 1/(1+(a/x)**b) for first compound
#' @param a2 parameter a in Hill-2 DR fuction 1/(1+(a/x)**b) for second compound
#' @param b2 parameter b Hill-2 DR fuction in 1/(1+(a/x)**b) for second compound
#' @param lambda1 parameter Voelund model
#' @param lambda2 parameter Voelund model
#'
#' @return response: a single numerical value
#' @export
DR_Voelund_Hill2 <- function(d1,d2,a1,b1,a2,b2,lambda1,lambda2)
{
  phi = function(E,d1,d2,a1,b1,a2,b2,lambda1,lambda2)
  { ( d1/Hill2_inv(E,a1,b1) )**(1/lambda1) + ( d2/Hill2_inv(E,a2,b2) )**(1/lambda2)  - 1  }
  uniroot(f=phi,interval = c(0,1),
          d1,d2,a1,b1,a2,b2,lambda1,lambda2)$root
}


