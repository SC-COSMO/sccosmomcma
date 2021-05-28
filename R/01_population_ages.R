#' Get population and age structure for MX subnational regions
#'
#' \code{get_population_ages_mx} returns a dataframe with all of the 
#' processed information on population age structure. 
#' 
#' @param reload a flag (default is FALSE) of whether to 
#' redownload and process the data file
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
