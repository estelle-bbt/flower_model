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
  
  
  ## Quarto ----
  
  tarchetypes::tar_quarto(index, "index.qmd")
  
)



