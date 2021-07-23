#' Get crude birth rates for MX states
#'
#' \code{get_population_birth_mx} returns a data.frame with all of the
#' processed information on crude birth rate.
#'
#' @param reload Flag (default is FALSE) of whether to redownload and process 
#' the data file.
#'
#' @export
get_population_birth_mx <- function(reload = FALSE) {
  if (reload == TRUE | !exists("df_pop_birth_mx")) {
    download_population_birth_mx()
    process_population_birth_mx()
    #load(file = "./data/df_pop_birth_us.rda")
    load(file = GLOBAL_MX_POPULATION_BIRTH_FILE)
    
  }
  return(df_pop_birth_mx)
}

#' Get crude birth rates
#'
#' \code{get_population_birth} wrapper function that allows user to
#' subset data from various sources.
#' If only county is specified (other parameters = ""), 
#' it provides national level data.
#' If county is also specified (and the dataset contains counties),
#' then only those counties within the state(s) are returned.
#' 
#' @param country Country of desired data.
#' @param state State of desired data.
#' @param county County of desired data.
#' @return 
#' A data.frame with all of the processed information on crude birth rate.
#' @export
get_population_birth <- function(country = "", 
                                 state = "",
                                 county = "") {
  
  country <- stringi::stri_trans_general(country, "latin-ascii")
  state <-   stringi::stri_trans_general(state, "latin-ascii")
  county <-  stringi::stri_trans_general(county, "latin-ascii")
  
  match_country <- get_countrycode(country)
  df_return <- get_population_birth_mx()
    
  match_state <- get_fipscodes_mx(stringi::stri_trans_general(state, "latin-ascii"))
  # filter our return for the right country, state(s)
  df_return <- df_return %>%
    dplyr::mutate(s_code = get_fipscodes_mx(stringi::stri_trans_general(state, "latin-ascii"))) %>%
    dplyr::filter(get_countrycode(stringi::stri_trans_general(country, "latin-ascii")) %in% match_country &
                    s_code == match_state)
      
  return(df_return)
}