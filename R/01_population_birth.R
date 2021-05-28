#' Get crude birth rates for MX states
#'
#' \code{get_population_birth_mx} returns a dataframe with all of the
#' processed information on crude birth rate.
#'
#' @param reload a flag (default is FALSE) of whether to
#' redownload and process the data file
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