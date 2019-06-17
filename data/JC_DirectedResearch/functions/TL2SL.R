#' @title TL2SL
#' @description Convert Total Length (cm) to Standard Length (cm)
#' 
#' @param TL A fish's total length, in centimeters
#' @param a,b The a, and b parameters. See details for further specifications.
#' 
#' @details a must always be greater than 1, since it is the slope in a conversion with the form TL = b + a*SL. When solving for SL, we obtain SL = (TL - b)/a
#'
#' @return SL a fish's standard length, in cm
#'
#' @export

TL2SL <- function(TL, aTL, bTL){
  SL <- (TL - bTL)/aTL
  
  return(SL)
}