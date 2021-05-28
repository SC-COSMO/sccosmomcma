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