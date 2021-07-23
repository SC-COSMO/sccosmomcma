#' Get population density for Mexico subnational regions
#'
#' \code{get_densities_mx} returns a data.frame with all of the 
#' processed information on weighted population density. 
#' 
#' @param reload Flag (default is FALSE) of whether to 
#' redownload and process the data file.
#'
#' @export
get_densities_mx <- function(reload = FALSE) {
  if (reload == TRUE | !exists("df_densities_mx")) {
    download_densities_mx()
    process_densities_mx()
    load(file = "./data/df_densities_mx.rda")
  }
  return(df_densities_mx)
}

#' Get population density 
#'
#' \code{get_densities} wrapper function that allows user to subset data from 
#' various sources. 
#' If only county is specified (other parameters = ""), it provides
#' national level data.
#' Currently for US, if state(s) is specified, then it returns
#' lowest geographical level (state or county) that it has
#' within those specified state(s). 
#' If county is also specified (and the dataset contains counties), then only 
#' those counties within the state(s) are returned.
#' If county is also specified (and the dataset contains counties), then only 
#' those counties within the state(s) are returned.
#' 
#' @param country Country of desired data.
#' @param state State of desired data.
#' @param county County of desired data.
#' @param reload Flag (default is FALSE) of whether to redownload and process 
#' the data file.
#' @return 
#' A data.frame with population density of desired country/state/county.
#'
#' @export
get_densities <- function(country = "", 
                          state = "",
                          county = "") {
  
  country <- stringi::stri_trans_general(country, "latin-ascii")
  state <-   stringi::stri_trans_general(state, "latin-ascii")
  county <-  stringi::stri_trans_general(county, "latin-ascii")
  
  df_return <- get_densities_mx()
  
  match_country <- get_countrycode(country)
  match_state <- get_fipscodes_mx(stringi::stri_trans_general(state, "latin-ascii"))
  
  # filter our return for the right country, state(s)
  df_return <- df_return %>%
    dplyr::mutate(s_code = get_fipscodes_mx(stringi::stri_trans_general(state, "latin-ascii"))) %>%
    dplyr::filter(get_countrycode(stringi::stri_trans_general(country, "latin-ascii")) %in% match_country &
             s_code == match_state)
  
  return(df_return)
}

#' Get a three letter code for a country's name
#'
#' \code{get_countrycode} returns a string with a 3 letter code
#' for a country's name to help with matching between datasets. 
#'
#' @export
get_countrycode <- function(country) {
  return(suppressWarnings(countrycode::countrycode(country, origin = 'country.name', destination = 'iso3c')))
}

#' Get FIPS codes for Mexican states
#'
#' \code{get_fipscodes_mx} returns a string with the FIPS code 
#' for Mexican states.
#' 
#' @param state State of desired data.
#'
#' @export
get_fipscodes_mx <- function(state){
  tmp <- c()
  for(i in state){
    tmp <- c(tmp, df_MEX_codes$s_code[df_MEX_codes$state==i])
  }
  return(tmp)
}