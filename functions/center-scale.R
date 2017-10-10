center_scale <- function(x){
  if(all(is.finite(x[is.na(x) == F])) & min(x, na.rm = T) != 0){

    y <- (x - mean(x, na.rm = T)) / (2 * sd(x, na.rm = T))
  } else {
    y <- x
  }

  return(y)

}