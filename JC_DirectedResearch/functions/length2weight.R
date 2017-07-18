#' @title length2weight
#' @description Convert Total Length (cm) to Standard Length (cm)
#' 
#' @param count Fish abundance
#' @param length Fish length, in cm
#' @param a,b the allometric parameters
#' @param length_in the type of length (SL or TL) for which a and b can convert
#' @param la,lb the a and b parameters to convert TL to SL using the reverse equation where SL = (TL - lb) / la
#' 
#' @return w 

length2weight <- function(count, length, a, b, length_in, la, lb){
  w <- ifelse(length_in == "SL",
              count * TL2W(TL2SL(length, la, lb),a, b),
              count * TL2W(length, a, b))
  return(w)
}