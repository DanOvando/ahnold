
spread_factor <- function(dat,var,omit = 'first')
{

  thing <- eval(parse(text = paste('dat$',var, sep = '')))

  if (class(thing) == 'character'){ thing <- as.factor(thing)}

  num_levels <- length(levels(thing))

  names <- paste(var,levels(thing), sep = '.')

  names <- names[2:length(names)]

  thingframe <- matrix(0,nrow = dim(dat)[1], ncol = length(names))

  colnames(thingframe) <- names
  for (i in 1:dim(dat)[1])
  {

    which_thing <- grep(as.character(thing[i]),names)

    thingframe[i,which_thing] <- 1

  }

  dat <- dat[, colnames(dat) %in% var == F]

  dat <- cbind(dat,thingframe)

  return(dat)
}
