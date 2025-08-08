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

#' Format seeds info
#'
#' @description 
#' This function format seed data:
#' seeds_data should be a data frame with three variables for each seed: 
#' the flower_id, the actual_father, and the mother_id 
#' each individual being given the same index when acting as father and as mother.
#' It also remove flowers that has not been recorded as visited.
#'
#' @param table data_obs and data_paternities
#'
#' @return seeds_info
#' 
#' @import dplyr
#' 
#' @export

format_seeds_info <- function(data_obs = data_obs, data_paternities = data_paternities){
  
  true_visit <- data_obs$id_flow_full
  
  seeds_info <- data_paternities %>%
    rename(actual_father=candidate_id) %>%
    mutate(mother=known_id,
           # remove the letters (H for hampe, R for receptive) to match observation data format
           flower_id = stringr::str_replace(id_flow, "_H([0-9]+)_R([0-9]+)", "_\\1_\\2"),
           # remove the last letter in case of double fruit to match observation data
           flower_id = sub("[A-Za-z]$", "", flower_id)) %>%
    select(flower_id,mother,actual_father) %>%
    filter(flower_id %in% true_visit)
  
  return(seeds_info)
}

#' Format flower visits
#'
#' @description 
#' flowers_list should be a list indexed by flower_id, 
#' whose each element is a small data frame with the following variables:
#' the identities bbb_ids of the bumblebee in each successive visit, 
#' and for each bumblebee its bbb_prev_visits (“vb(j)”)
#'
#' @param table data_obs
#'
#' @return list flowers_list
#' 
#' @import dplyr
#' 
#' @export

format_flowers_list <- function(data_obs = data_obs){
  
  flowers <- data_obs %>%
    mutate(bbb_ids=paste0(session,"_",people)) %>%
    group_by(bbb_ids) %>%
    mutate(bbb_prev_visits=row_number()) %>%
    ungroup() %>% 
    dplyr::select(id_flow_full,bbb_ids,bbb_prev_visits)
  
  # change dataframe into list
  flowers_list <- split(flowers,f=flowers$id_flow_full)
  flowers_list <- lapply(flowers_list, as.data.frame)
  
  return(flowers_list)
}

#' Format bumblebee visits
#'
#' @description 
#' 
#' bumblebees_list should be a list indexed by bbb_id, whose each element is a data frame
#' describing the visit history by the focal bumblebee in terms of variables 
#' visited_plant and flower_age, ordered as the visits 
#' (for easy selection of the bumblebee’s bbb_prev_visits previous visits).
#'
#' @param table data_obs
#'
#' @return list bbbs
#' 
#' @import dplyr
#' 
#' @export

format_bumblebees_list <- function(data_obs = data_obs){
  
  bbbs <- data_obs %>%
    arrange(session,time) %>% # to calculate flower age (visits needed to be sequentially arrange, independently of the pollinators)
    pivot_wider(values_from=id_flow_full,names_from=id_flow_full,names_prefix="newflo_",values_fn = list(id_flow_full = ~ 1), values_fill = list(id_flow_full = 0)) %>%
    mutate(across(starts_with("newflo"),~ifelse(duplicated(cumsum(.)),NA,cumsum(.)))) %>%
    rowwise() %>%
    mutate(flower_age=sum(across(starts_with("newflo")),na.rm=T)-1) %>%
    select(!starts_with("newflo")) %>%
    mutate(bbb_ids=paste0(session,"_",people)) %>%
    rename(visited_plant=ID_full) %>%
    arrange(session,people,time) %>% # re-arrange by focal pollinator for easy selection of previous visits
    dplyr::select(bbb_ids,visited_plant,flower_age)
  
  # dataframe into list
  bbbs_list <- split(bbbs,f=bbbs$bbb_ids)
  bbbs_list <- lapply(bbbs_list, as.data.frame)
  
  return(bbbs_list)
}

#' Function for log likelihood calculation
#'
#' @description 
#' 
#' log likelihood calculation for given flower 'flower_id' 
#' on given mother 'mother' plant:
#'
#' @param 
#'
#' @return 
#' 
#' @import dplyr
#' 
#' @export

calc_logl <- function(S, alpha, beta, gamma, seed_info,
                      flower_id=seed_info$flower_id,
                      actual_father=seed_info$actual_father,
                      mother=seed_info$mother,
                      n_males=10L,
                      flowers_list = flowers_list,
                      bumblebees_list = bumblebees_list) {
  focal_flower_info <- flowers_list[[flower_id]]
  nb_visits_focal_flower <- nrow(focal_flower_info)
  bbb_ids <- focal_flower_info$bbb_ids
  actual_father_rel_contrib <- numeric(nb_visits_focal_flower)
  format_n_males <- formatC(seq_len(n_males),width=2,flag="0")
  maleNames <- paste0(substr(mother,1,5),".",format_n_males)
  for (j in seq(nb_visits_focal_flower)) {
    focal_bbb <- bbb_ids[j]
    bbb_prev_visits <- focal_flower_info$bbb_prev_visits[j]
    focal_bbb_info_all <- bumblebees_list[[focal_bbb]]
    focal_bbb_info <- focal_bbb_info_all[seq_len(bbb_prev_visits),,drop=FALSE]
    # unique(focal_bbb_info_all[seq_len(bbb_prev_visits),"visited_plant",drop=FALSE])
    male_contribs <- numeric(n_males)
    names(male_contribs) <- maleNames
    for (male in maleNames){
      # male <- paste0(substr(mother,1,5),".",formatC(m,width=2,flag="0")) # to grep the correct ID according to session
      Xi <-  focal_bbb_info$visited_plant==male
      whichXi <- which(Xi)
      # exponential:
      bbb_load <- alpha^whichXi
      male_contribs[male] <- (1-(male==mother)*(1-S)) * sum(bbb_load * beta^focal_bbb_info$flower_age[whichXi])
    }
    relative_male_contribs <- male_contribs/sum(male_contribs,na.rm=T)
    actual_father_rel_contrib[j] = relative_male_contribs[actual_father]
  }
  # exponential:
  # visit_weights <- gamma^seq_len(nb_visits_focal_flower)
  # saturating:
  visit_weights <- gamma^(nb_visits_focal_flower-seq_len(nb_visits_focal_flower))
  
  # if male==mother the male contrib is S So S>0 necessarily 
  # logLik of actual father for given seed:
  return(log(sum(visit_weights*actual_father_rel_contrib,na.rm=T)) -log(sum(visit_weights,na.rm=T)))
}

#' Function to get likelihood and needed for proliks
#'
#' @description 
#' 
#' will be used in proliks latter to get confidence intervals
#'
#' @param 
#'
#' @return 
#' 
#' @import dplyr
#' 
#' @export

logL <- function(params, seeds_info, flowers_list, bumblebees_list) {
  S      <- exp(unlist(params["logS"]))
  alpha  <- exp(unlist(params["logalpha"]))
  beta   <- exp(unlist(params["logbeta"]))
  gamma  <- exp(unlist(params["loggamma"]))
  
  nb_seeds <- nrow(seeds_info)
  logL_vals <- numeric(nb_seeds)
  
  for (seed_id in seq_len(nb_seeds)) {
    logL_vals[seed_id] <- calc_logl(
      S = S,
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      seed_info = seeds_info[seed_id, ],
      flowers_list = flowers_list,
      bumblebees_list = bumblebees_list
    )
    
    # in case the male was not visited before, relative contribution = 0 so log = infinite
    # These impossible observatiosn are removed from the fitted data =>
    # we replace by 0 (to not be accounted for in the sum)
    if (is.infinite(logL_vals[seed_id])) {
      logL_vals[seed_id] <- 0
    }
  }
  
  sum(logL_vals)
}

#' Function to get likelihood and confidence intervals
#'
#' @description 
#' 
#' the package proliks is needed 
#' to do: add a function to enter parameters values
#'
#' @param 
#'
#' @return 
#' 
#' @import dplyr
#' 
#' @export

get_lik_and_ci <- function(seeds_info = seeds_info,
                           flowers_list = flowers_list,
                           bumblebees_list = bumblebees_list) {
  
  params = list("logS" = log(0.5),
                "logalpha" = log(0.5),
                "logbeta" = log(0.5),
                "loggamma" = log(0.5))
  
  lower = list("logS" = log(0),
               "logalpha" = log(0),
               "logbeta" = log(0),
               "loggamma" = log(0))
  
  upper = list("logS" = log(10),
               "logalpha" = log(10),
               "logbeta" = log(1),
               "loggamma" = log(10))
  
  # Pass a function, not a number
  Likobj <- proliks::as_Lik(
    # add "par" to get a function - otherwise we directly get the result
    logLfn = function(par) logL(par, seeds_info, flowers_list, bumblebees_list),
    lower = pmax(unlist(lower), c(-10,-10,-10,-30)),
    upper = pmin(unlist(upper), 10)
  )
  
  log_ci <- proliks::doMLfit(Likobj, init = c(-5, 0, -4, -10))
  output <- exp(summary(log_ci, verbose = FALSE)$intervals)
  
  return(output)
}

# get_lik_and_ci(seeds_info = seeds_info)
