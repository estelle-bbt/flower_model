#' Targets plan
#' 

## Attach required packages ----

library(targets)
library(tarchetypes)
library(ggplot2)

tar_option_set(
  packages = c("dplyr","tidyr","ggplot2","forcats")  # load dplyr in each environement
)

tar_source()

## Load Project R Functions ----

source(here::here("R", "functions.R"))

## Analyses pipeline ----

list(
  
  ## Manage data ----
  
  tar_target(proliks,install_proliks()),
  
  tar_target(data_obs,load_data("data/pollinator_observations.txt")),
  
  tar_target(data_paternities,load_data("data/paternities.txt")),
  
  tar_target(seeds_info,format_seeds_info(data_obs = data_obs, data_paternities = data_paternities)),
  
  tar_target(flowers_list,format_flowers_list(data_obs = data_obs)),
  
  tar_target(bumblebees_list,format_bumblebees_list (data_obs = data_obs)),
  
  tar_target(lik_and_ci,get_lik_and_ci(seeds_info = seeds_info,
                                       flowers_list = flowers_list,
                                       bumblebees_list = bumblebees_list)),
  
  ## Quarto ----
  
  tarchetypes::tar_quarto(index, "index.qmd")
  
)



