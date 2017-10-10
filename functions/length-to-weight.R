length_to_weight <-
  function(mean_length,
           min_length,
           max_length,
           count,
           weight_a,
           weight_b,
           length_units = 'cm',
           weight_units = 'g',
           length_for_weight_units = 'mm',
           length_type_for_weight,
           tl_sl_a,
           tl_sl_b,
           tl_sl_type,
           tl_sl_formula) {
    #
    # generate_lengths <- function(count,mean_length, min_length, max_length){

    if (is.na(count) | count == 0){

      outweight <-  0
    }  else {

      if (is.na(min_length) |
          is.na(max_length)) {
        #generate distribution of lengths

        lengths <-  rep(mean_length, count)

      } else{
        # lengths <-  pmax(min_length,pmin(max_length,rpois(count, lambda = mean_length)))
        lengths <- runif(count, min = min_length, max = max_length)
      }

      if (length_type_for_weight == 'SL') {
        if (tl_sl_type  == 'TYPICAL') {
          weight_lengths <-  lengths * tl_sl_a + tl_sl_b
        } else{
          weight_lengths <- (lengths - tl_sl_b) / tl_sl_a

        }

      } else {
        weight_lengths <-  lengths
      }

      if (length_units == 'cm' & length_for_weight_units == 'mm') {
        weight_lengths <- weight_lengths * 10
      }

      weight <-  weight_a * weight_lengths ^ weight_b

      if (weight_units == 'kg') {
        weight <- weight * 1000
      }
      outweight = sum(weight)
    }
    return(outweight)
  } #close function
