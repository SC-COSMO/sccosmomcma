#' Get population and age structure for MX subnational regions
#'
#' \code{get_population_ages_mx} returns a data.frame with all of the 
#' processed information on population age structure. 
#' 
#' @param reload Flag (default is FALSE) of whether to redownload and process 
#' the data file
#'
#' @export
get_population_ages_mx <- function(reload = FALSE) {
  if (reload == TRUE | !exists("df_pop_age_mx")) {
    download_population_ages_mx()
    process_population_ages_mx()
    load(file = "./data/df_pop_age_mx.rda")
  }
  return(df_pop_age_mx)
}

#' Get population and age structure 
#'
#' \code{get_population_ages} wrapper function that allows user to
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
#' A data.frame with all of the processed information on population age structure.
#'
#' @export
get_population_ages <- function(country = "", 
                                state = "",
                                county = "") {
  
  country <- stringi::stri_trans_general(country, "latin-ascii")
  state <-   stringi::stri_trans_general(state, "latin-ascii")
  county <-  stringi::stri_trans_general(county, "latin-ascii")
  
  match_country <- get_countrycode(country)

  df_return <- get_population_ages_mx()
  
  match_state <- get_fipscodes_mx(state)
  
  # filter our return for the right country, state(s)
  df_return <-  df_return %>%
    dplyr::mutate(country_code = get_countrycode(stringi::stri_trans_general(country, "latin-ascii"))) %>%
    dplyr::filter(country_code == match_country) %>%
    dplyr::mutate(s_code = get_fipscodes_mx(state = stringi::stri_trans_general(state, "latin-ascii"))) %>%
    dplyr::filter(s_code == match_state)
  
  return(df_return)
}