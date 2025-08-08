#' Generic function to just read dataset in .txt form
#'
#' @description 
#' This function allow to directly load the .txt dataset previously obtained.
#' Note that the sampling of genotypes, the script is available on dryad. We here 
#' directly provide our sampling dataset in a purpose of reproducible results.
#'
#' @param file a character of length 1. The path to the .txt file.
#'
#' @return A `table` containing data. 
#' 
#' @import dplyr
#' 
#' @export

load_data <- function(file_path){
  
  data <- read.table(file_path,head=T) 
  
  return(data)
}

#' Install proliks package
#'
#' @description 
#' This function install proliks (likelihood-based inference) 
#' contact F. Rousset for more informations
#'
#' @param none
#'
#' @return none
#' 
#' @import dplyr
#' 
#' @export

install_proliks <- function(){
  
  if (!requireNamespace("proliks", quietly = TRUE)) {
    install.packages("proliks_0.1.2.tar.gz", type = "source", repos = NULL)
  }
  
}

