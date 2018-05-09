#' assign tuned variables
#'
#' tunes a control variable to achieve a target depletion level
#'
#' @param control_variable
#' @param target_depletion
#' @param fish
#' @param fleet
#' @param sim_years
#' @param burn_years
#' @param num_patches
#'
#' @return ss for depletion
#' @export
#'
#' @examples
tune_fleet <- function( f_v_m,
                               fish,
                               fleet,
                               sim_years = 50,
                               burn_years = 10,
                               num_patches = 1){


  # tuned <-
  #   nlminb(
  #     fish$m * .8 / fleet$q,
  #     spasm::tune_for_depletion,
  #     lower = 0,
  #     target_depletion = target_depletion,
  #     fish = fish,
  #     fleet = fleet,
  #     sim_years = sim_years,
  #     burn_years = burn_years,
  #     num_patches = num_patches
  #   )

  if (fleet$fleet_model == "constant-catch"){

    # fleet$target_catch <-   tuned$par
    message("nope")

  } else if (fleet$fleet_model == "constant-effort"){

    fleet$initial_effort <-   (fish$m * f_v_m) / fleet$q

  } else if (fleet$fleet_model == "open-access"){

    message("Not yet you don't")
  }

  return(fleet)

}




