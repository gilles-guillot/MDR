#' @title Bivariate dose-response function under Jonker S/A model
#' @description Compute bivariate dose-response function when each mixture compoent follows a 2-parameter Hill dose-response function and Jonker S/A model holds.
#' @references Jonker model S/A (Eq. 7 p. 2703 Jonker et al. Env Tox & Chem 24(10) 2005)
#' @param d1 dose of first compound: a numeric vector
#' @param d2 dose of second compound: a numeric vector
#' @param a1 parameter a in Hill-2 DR fuction 1/(1+(a/x)**b) for first compound
#' @param b1 parameter b in Hill-2 DR fuction 1/(1+(a/x)**b) for first compound
#' @param a2 parameter a in Hill-2 DR fuction 1/(1+(a/x)**b) for second compound
#' @param b2 parameter b Hill-2 DR fuction in 1/(1+(a/x)**b) for second compound
#' @param alpha parameter quantifying departure from Loewe additivity in Jonker S/A model
#' @return vector of numerical responses

DR_JonkerSA_Hill2 <- function(d1,d2,a1,b1,a2,b2,alpha)
{
  phi = function(E,d1,d2,a1,b1,a2,b2,alpha) {
    EC1 = Hill2_inv(0.5,a1,b1)
    EC2 = Hill2_inv(0.5,a2,b2)
    TU1 = d1/EC1
    TU2 = d2/EC2
    STU = TU1 + TU2
    z1 = TU1/STU
    z2 = TU2/STU
    G = alpha*z1*z2
    return( d1/Hill2_inv(E,a1,b1) + d2/Hill2_inv(E,a2,b2) - exp(G)  )
  }
  uniroot(f=phi,interval = c(0,1),
          d1,d2,a1,b1,a2,b2,alpha)$root
}
