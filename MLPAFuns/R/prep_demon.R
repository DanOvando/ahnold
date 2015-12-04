#' Prepare objects required for Laplaces Demon
#'
#' \code{prep_demon} takes data and options and
#' returns a list with components required by Laplaces Demon
#' @param demondat data to be used
#' @param pos_vars possible variables to include in the model
#' @param  scale_numerics scale T or F to center and scale variables


prep_demon <- function(demondat, pos_vars,scale_numerics = T,constant = T)
{

  dep_vars <- demondat[,pos_vars] #pull out things you need

  var_types <- sapply(dep_vars,class)

  factors <- pos_vars[which(var_types == 'character' | var_types == 'factor') ]

  numerics <- pos_vars[which(var_types == 'numeric') ]

  if (scale_numerics == T)
  {
    for (j in 1:length(numerics)){
      dep_vars[,numerics[j]] <- CenterScale(as.matrix(dep_vars[,numerics[j]]))
    }
  }

  for (f in 1:length(factors))
  {
    dep_vars <- spread_factor(dep_vars,var = factors[f])
  }

  if (constant == T){
    dep_vars$constant <- 1
  }

return(dep_vars)
}