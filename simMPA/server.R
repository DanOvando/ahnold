
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
library(shiny)
library(spasm)
library(tidyverse)
library(hrbrthemes)
library(scales)
library(viridis)
set.seed(123)

shiny_theme <-  theme_ipsum(base_size = 20,
                            axis_title_size = 16,
                            strip_text_size = 16
)
theme_set(shiny_theme)


shinyServer(function(input, output) {

num_patches <-  50
mpa_experiment <- reactive({

fish <- create_fish(max_age = input$max_age,
                    age_mature = input$age_mat,
                    steepness = input$steepness,
                    adult_movement = input$adult_movement,
                    larval_movement = input$larval_movement,
                    vbk = input$vbk,
                    density_dependence_form = input$density_dependence_form,
                    scientific_name = 'fakeish fishis',
                    query_fishbase = F,
                    linf = 100,
                    t0 = -.2,
                    weight_a = 1e-4,
                    weight_b = 3)

fleet <- create_fleet(length_50_sel = input$length_sel[1],
                      length_95_sel = input$length_sel[2],
                      mpa_reaction = "concentrate",
                      price = input$price,
                      theta = input$theta,
                      q = input$q,
                      fleet_model = input$fleet_model,
                      effort_allocation = input$effort_allocation,
                      initial_effort = input$initial_effort,
                      target_catch = input$target_catch,
                      fish = fish)

experiment <- mpa_counterfactual(fish = fish, fleet = fleet, year_mpa = input$year_mpa,
                       mpa_size = input$mpa_size,
                       sim_years = input$sim_years,
                       num_patches = num_patches,
                       burn_years = 50)

comparison <- experiment$outcomes

min_year <- min(comparison$year)

mpa_size <-  max(comparison$percent_mpa)

year_mpa <- input$year_mpa

comp_plot <- comparison %>%
  filter(year > (min_year + year_mpa - 3) & year < (min_year + year_mpa + 10)) %>%
  gather("metric","value", ssb, catch, profits, effort) %>%
  ggplot(aes(year,value)) +
  geom_line(aes(color = experiment)) +
  geom_point(size = 2,shape = 21,aes(fill = experiment)) +
  geom_vline(aes(xintercept = year_mpa + min_year - 1)) +
  facet_wrap(~metric, scales = 'free_y') +
  labs(caption = paste0(100*mpa_size,'% MPA'))


adult_move_grid <-
  expand.grid(source = 1:num_patches,
              sink = 1:num_patches) %>%
  mutate(
    distance = source - sink,
    prob = 1 / ((2 * pi) ^ (1 / 2) * fish$adult_movement) * exp(-(distance) ^
                                                                  2 / (2 * fish$adult_movement ^ 2))
  ) %>%
  group_by(source) %>%
  mutate(prob_move = prob / sum(prob))

adult_movement_plot <- ggplot(adult_move_grid, aes(source, sink, fill = prob_move)) +
  geom_tile() +
  scale_fill_viridis(name = "% Move", labels = percent) +
  labs(title = "Adults")

larval_move_grid <-
  expand.grid(source = 1:num_patches,
              sink = 1:num_patches) %>%
  mutate(
    distance = source - sink,
    prob = 1 / ((2 * pi) ^ (1 / 2) * fish$larval_movement) * exp(-(distance) ^
                                                                  2 / (2 * fish$larval_movement ^ 2))
  ) %>%
  group_by(source) %>%
  mutate(prob_move = prob / sum(prob))

larval_movement_plot <- ggplot(larval_move_grid, aes(source, sink, fill = prob_move)) +
  geom_tile() +
  scale_fill_viridis(name = "% Move", labels = percent) +
  labs(title = "Larvae")


age_dist_plot <- experiment$wi_mpa %>%
  filter(year == max(year)) %>%
  ggplot(aes(age,numbers)) +
  geom_col()

patch_ssb_plot <- experiment$wi_mpa %>%
  filter(year == max(year)) %>%
  group_by(patch) %>%
  summarise(ssb = sum(ssb), mpa = unique(mpa)) %>%
  ggplot(aes(patch,ssb, fill = mpa)) +
  geom_col()


mpa_age_dist_plot <- experiment$wi_mpa %>%
  filter(year == max(year)) %>%
  group_by(mpa,age) %>%
  summarise(numbers = sum(numbers)) %>%
  ggplot(aes(x = age,y =numbers, fill = mpa)) +
  geom_col(position = 'dodge')


return(list(comparison = comparison, comp_plot = comp_plot,
            adult_movement_plot = adult_movement_plot,
            larval_movement_plot = larval_movement_plot,
            age_dist_plot = age_dist_plot,
            patch_ssb_plot = patch_ssb_plot,
            mpa_age_dist_plot = mpa_age_dist_plot))

})

mpa_experiment_debounced <- mpa_experiment %>% debounce(millis = 1000)


  output$comp_plot <- renderPlot({
    mpa_experiment_debounced()$comp_plot
  })
  output$adult_movement_plot <- renderPlot({
    mpa_experiment_debounced()$adult_movement_plot
  })
  output$larval_movement_plot <- renderPlot({
    mpa_experiment_debounced()$larval_movement_plot
  })
  output$age_dist_plot <- renderPlot({
    mpa_experiment_debounced()$age_dist_plot
  })
  output$patch_ssb_plot <- renderPlot({
    mpa_experiment_debounced()$patch_ssb_plot
  })
  output$mpa_age_dist_plot <- renderPlot({
    mpa_experiment_debounced()$mpa_age_dist_plot
  })

})
