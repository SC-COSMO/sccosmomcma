#' Get lifetables for MX states
#'
#' \code{get_lifetables_mx} returns a dataframe with all of the 
#' processed information on lifetables. 
#' 
#' @param reload a flag (default is FALSE) of whether to 
#' redownload and process the data file
#'
#' @export
get_lifetables_mx <- function(reload = FALSE) {
  if (reload == TRUE | !exists("df_lifetables_mx")) {
    download_lifetables_mx()
    process_lifetables_mx()
    load(file = "./data/df_lifetables_mx.rda")
  }
  return(df_lifetables_mx)
}

#' Get Lifetables 
#'
#' \code{get_lifetables} wrapper function that allows user to
#' subset data from various sources.
#' If only county is specified (other parameters = ""), 
#' it provides national level data.
#' If county is also specified (and the dataset contains counties),
#' then only those counties within the state(s) are returned
#' 
#' @param reload a flag (default is FALSE) of whether to 
#' redownload and process the data file
#'
#' @export
get_lifetables <- function(country = "", 
                           state = "",
                           county = "") {
  
  country <- stringi::stri_trans_general(country, "latin-ascii")
  state <-   stringi::stri_trans_general(state, "latin-ascii")
  county <-  stringi::stri_trans_general(county, "latin-ascii")
  
  
  match_country <- get_countrycode(country)
  df_return <- get_lifetables_mx() 
    
  match_state <- state
  # filter our return for the right country, state(s)
  df_return <- df_return %>%
    dplyr::filter(get_countrycode(stringi::stri_trans_general(country, "latin-ascii")) %in% match_country &
             stringi::stri_trans_general(state, "latin-ascii") == match_state)
      
  return(df_return)
}
