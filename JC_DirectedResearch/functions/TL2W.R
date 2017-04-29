#' @title TL2W
#' @description Convert Total Length (cm) to Weight (gr)
#' 
#' @param TL A fish's total length, in centimeters
#' @param a,b The a, and b parameters, specificaly for converting from cm to gr!
#'
#' @return W a fish's weight, in grams
#'
#' @export

TL2W <- function(TL, a, b){
  W = a*TL^b
  
  return(W)
}