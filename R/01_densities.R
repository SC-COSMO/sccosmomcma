#' Get population density for MX subnational regions
#'
#' \code{get_densities_mx} returns a dataframe with all of the 
#' processed information on weighted population density. 
#' 
#' @param reload a flag (default is FALSE) of whether to 
#' redownload and process the data file
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