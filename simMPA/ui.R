



# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#


shinyUI(fluidPage(
  # Application title
  titlePanel("simMPA"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      sliderInput(
        "linf",
        "Linf:",
        min = 10,
        max = 200,
        value = 100
      ),
      sliderInput(
        "m",
        "Natural Mortality",
        min = 0.1,
        max = .9,
        value = 0.2
      ),
      sliderInput(
        "steepness",
        "Steepness",
        min = 0.201,
        max = .99,
        value = 0.6
      ),
      sliderInput(
        "adult_movement",
        "Adult Movement",
        min = 0,
        max = 100,
        value = 2
      ),
      sliderInput(
        "larval_movement",
        "Larval Movement:",
        min = 0,
        max = 100,
        value = 10
      ),
      selectInput(
        inputId = 'lhi_type',
        label = 'Prince LHI Type',
        choices = c(1,2,3),
        selected = 1
      ),
      selectInput(
        inputId = 'density_dependence_form',
        label = 'Density Dependence',
        choices = c(
          'Density dependence occurs independently in each spatial area' = 1,
          "Density dependence occurs at the regional level and recruits are distributed evenly across areas, or larvae are distributed evenly followed by local density dependence based on the local number of larvae" = 2,
          "Density dependence occurs within spatial areas, but recruits are spread evenly across areas
          " = 3,
          "The pooled larvae from all areas are distributed evenly to each area, and then density dependence occurs based on the number of spawners in each area" = 4,
          "Recruitment is independent in each area, but a fraction of the recruits in each area drift to the adjacent areas before settling" = 5
        ),
        selected = 1
      ),
      sliderInput(
        "mpa_size",
        "% of Fishery in MPA",
        min = 0,
        max = 1,
        value = 0.25
      ),
      sliderInput(
        "year_mpa",
        "Year MPA Implemented",
        min = 1,
        max = 100,
        value = 10
      ),
      sliderInput(
        "sim_years",
        "Simulation Years",
        min = 10,
        max = 100,
        value = 25
      ),
      sliderInput(
        "initial_effort",
        "Initial Effort",
        min = 0,
        max = 1000,
        value = 100
      ),
      sliderInput(
        "target_catch",
        "Target Catch",
        min = 0,
        max = 10000,
        value = 100
      ),
      sliderInput(
        "length_sel",
        "Selectivity at Length",
        min = 0,
        max = 200,
        value = c(1, 2)
      ),
      sliderInput(
        "price",
        "Price per kg",
        min = 0,
        max = 100,
        value = 1
      ),
      sliderInput(
        "theta",
        "Open Access Rate",
        min = 0,
        max = 1,
        value = .1
      ),
      sliderInput(
        "q",
        "Catchability",
        min = 0,
        max = .1,
        value = 1e-3
      ),
      selectInput(
        "fleet_model",
        "Fleet Model:",
        choices = c("constant-catch", "constant-effort", "open-access"),
        selected = "constant-effort"
      ),
      selectInput(
        "effort_allocation",
        "Effort Allocation:",
        choices = c("simple", "gravity"),
        selected = "gravity"
      )
    ),
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(type = 'tabs',
      tabPanel("Fishery Metric",plotOutput("comp_plot")),
      tabPanel("Population Structure", plotOutput("age_dist_plot"),
               plotOutput("patch_ssb_plot"), plotOutput("mpa_age_dist_plot")),
      tabPanel("Movement", plotOutput("adult_movement_plot"),
               plotOutput("larval_movement_plot"))
      ) #close tabsetPanel
      ) #close mainPanel
    ) # close sidebarLayout
) # close fluidpage
) #close app
