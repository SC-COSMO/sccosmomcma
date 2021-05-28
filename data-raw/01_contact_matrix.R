library(sccosmoData)
library(stringi)
library(usmap)
library(countrycode)

### Define country and state
choose_country <- "Mexico"
choose_state <- "MCMA"

### Get contact matrices from sccosmoData
l_contact_matrices <- sccosmoData::get_contact_matrix(country = choose_country, 
                                                      state = choose_state, 
                                                      density = 141341)
### Create .rda object of setting-specific contact matrices for an exemplary 
### population and stor it in 'data' folder
usethis::use_data(l_contact_matrices, overwrite = TRUE)

#' Get contact matrices df
#'
#' \code{get_contact_matrices_world} returns a dataframe with all of the 
#' processed contact matrices for all countries. 
#' 
#' @param reload a flag (default is FALSE) of whether to 
#' redownload and process the data file
#'
#' @export
get_contact_matrices_world <- function(reload = FALSE) {
  if (reload == TRUE | !exists("df_contacts_src_world")) {
    download_contact_matrices_world()
    process_contact_matrices_world()
    load(file = "./data/df_contacts_src_world")
  } 
  return(df_contacts_src_world)
}

#' Get contact matrix main function
#'
#' \code{get_contact_matrix} wrapper function that allows user to
#' subset data from various sources.
#' If only county is specified (other parameters = ""), 
#' it provides national level data
#' Currently for US, if state(s) is specified, then it returns
#' lowest geographical level (state or county) that it has
#' within those specified state(s)
#' If county is also specified (and the dataset contains counties),
#' then only those counties within the state(s) are returned
#' 
#' @param v_init_age_grps vector that specifies the age bins to aggregate and return
#' @param country country of desired data
#' @param state state of desired data
#' @param county county of desired data
#' @param src_country country to use for contact matrix structure 
#'                    if country is not in dataframe
#' @param density population density of country/state/county 
#' @param v_custom_prop_pop vector that specifies the proportion of total population for each age group.
#' @param recache TRUE/FALSE. If true, will run function and replace existing cache. If false, will check for cache
#' and use if cached matrix exists.

#' Must be same length as v_init_age_grps+1 and be listed in the same order as v_init_age_grps.
#'
#' @export
get_contact_matrix <- function(v_init_age_grps = c(0, 5, 15, 25, 45, 55, 65, 70),
                               country = "", 
                               state = "",
                               county = "",
                               src_country = "",
                               density = "",
                               v_custom_prop_pop = "",
                               recache = FALSE) {
  
  n_init_age_grps <- length(v_init_age_grps)
  
  # Clean up input strings
  match_country <- stri_trans_general(country, "latin-ascii")
  match_state <-   stri_trans_general(state, "latin-ascii")
  match_county <-  stri_trans_general(county, "latin-ascii")
  src_country <-  stri_trans_general(src_country, "latin-ascii")
  
  if (recache == F) {
    if (is_in_cache(GLOBAL_CONTACT_MATRICES_CACHE, match_country, match_state, match_county, density = density)) {
      return (get_from_cache(GLOBAL_CONTACT_MATRICES_CACHE, match_country, match_state, match_county, density = density))
    }
  }
  
  ##############################
  #### Error checking for inputs
  ##############################
  if (length(match_country) > 1 | length(match_state) > 1 | length(match_county) > 1 | length(density) > 1 | length(src_country) > 1) {
    stop("This function processes only one input at a time. If running for multiple locations, consider using apply or a loop.")
  }
  
  if (!is.na(match("", match_country))) {
    stop("At minimum, must supply country")
  }
  
  if (is.na(match("", match_country)) & is.na(match("", src_country))) {
    # Check whether the requested country exists in the source country list
    if ((match_country %in% unique(sccosmoData::df_contacts_src_world$country))) {
      warning("Supplied both src_country and country. Using src_country ('", paste0(src_country),"') for contact matrix.")
    }
  }
  
  if (is.na(match("", match_country)) & !is.na(match("", src_country))) {
    src_country <- match_country
  }  
  
  # Check whether the requested country exists in the source country list
  if (!(src_country %in% unique(df_contacts_src_world$country))) {
    stop("'", src_country, "' is not in contact matrix country list. Need to supply an alternate src_country.")
  }
  
  ##############################
  #### Prepare Contact Matrix
  ##############################
  
  ## Pull contact matrix for home (does not need to be scaled)
  df_contacts_home <- wrangle_contact_matrix(setting = "home",
                                             src_country = src_country,
                                             v_init_age_grps = v_init_age_grps)
  
  # Get populations and ages of target population
  if (is.na(match("", v_custom_prop_pop))) {
    if (length(v_custom_prop_pop) != n_init_age_grps) {
      stop("Length of v_custom_prop_pop is ", length(v_custom_prop_pop), ", but must equal length of v_init_age_grps (", n_init_age_grps,").")
    }
    if (sum(v_custom_prop_pop) != 1) {
      stop("v_custom_prop_pop must sum to 1. Currently sums to: ", sum(v_custom_prop_pop))
    }
    df_pop_ages <- data.frame(AgeGrpContact = cut(x = v_init_age_grps,
                                                  breaks = c(v_init_age_grps, Inf), 
                                                  include.lowest = TRUE, right = FALSE),
                              prop_pop = v_custom_prop_pop)
    message("Custom age structure applied: ")
    print(df_pop_ages)
  } else {
    message(paste0("Pulling age structure for: ", match_country, " ", match_state, " ", match_county))
    df_pop_ages <- aggregate_pop_ages(match_country = match_country,
                                      match_state = match_state,
                                      match_county = match_county,
                                      ages_to_cut = v_init_age_grps)
    
    ## Supplied state and county, only national age_pop
    if (is.na(match("", match_state)) & df_pop_ages$state[1] == "" & is.na(match("", match_county)) & df_pop_ages$county[1] == "") {
      warning("State-level and county-level populations are not available from sccosmoData. Using country-level populations instead.")
    } 
    ## Supplied state and county, only state age_pop
    if (is.na(match("", match_state)) & df_pop_ages$state[1] != "" & is.na(match("", match_county)) & df_pop_ages$county[1] == "") {
      warning("County-level populations is not available from sccosmoData. Using state-level populations instead.")
    } 
    ## Supplied state, only national age_pop
    if (is.na(match("", match_state)) & df_pop_ages$state[1] == "" & !is.na(match("", match_county))) {
      warning("State-level populations is not available from sccosmoData. Using country-level populations instead.")
    } 
    
    df_pop_ages <- df_pop_ages %>%
      ungroup() %>%
      dplyr::select(AgeGrpContact, prop_pop)
  }
  
  # Get density of target population
  if (is.na(match("", density))) {
    density_goal <- density
  } else {
    message(paste0("Pulling density for: ", match_country, " ", match_state, " ", match_county))
    density_goal <- suppressWarnings(get_densities(country = match_country,
                                                   state = match_state,
                                                   county = match_county))
    
    ## Supplied state and county, only national density
    if (is.na(match("", match_state)) & density_goal$state == "" & is.na(match("", match_county)) & density_goal$county == "") {
      warning("State-level and county-level densities are not available from sccosmoData. Using country-level density instead.")
    } 
    ## Supplied state and county, only state density
    if (is.na(match("", match_state)) & density_goal$state != "" & is.na(match("", match_county)) & density_goal$county == "") {
      warning("County-level density is not available from sccosmoData. Using state-level density instead.")
    } 
    ## Supplied state, only national density
    if (is.na(match("", match_state)) & density_goal$state == "" & !is.na(match("", match_county))) {
      warning("State-level density is not available from sccosmoData. Using country-level density instead.")
    } 
    density_goal <- pull(density_goal, density)
  }
  
  # Compute ratio between predicted contacts based on density of goal population and source population
  density_scale_factor <- (df_contacts_home$non_home_intercept_coeff[1] +
                             (log(density_goal))*
                             df_contacts_home$non_home_dens_coeff[1]) /
    (df_contacts_home$non_home_intercept_coeff[1] +
       (log(df_contacts_home$density[1]))*
       df_contacts_home$non_home_dens_coeff[1])
  message(paste0("Density Scale Factor: ", density_scale_factor))
  
  
  tmp_contacts <- list()
  for (loc in c("work", "school", "other_locations")) {
    # Pull location-specific contact_matrix
    message(paste0("Pulling ", loc, " contacts for source country: ", src_country))
    tmp_contacts[[loc]] <- wrangle_contact_matrix(setting = loc,
                                                  src_country = src_country,
                                                  v_init_age_grps = v_init_age_grps) %>%
      left_join(df_pop_ages) %>%
      # Adjust for density and population age structure
      mutate(contacts = standardized_contacts * density_scale_factor * prop_pop) %>%
      dplyr::select(country, AgeGrp, AgeGrpContact, contacts) %>%
      # Clean up to return expected format for sccosmo
      mutate(age_start_contact = as.character(AgeGrpContact),
             age_start_contact = gsub("\\[", "", age_start_contact),
             age_start_contact = gsub(",.*", "", age_start_contact),
             age_start_contact = as.numeric(age_start_contact)) %>%
      dplyr::select(-AgeGrpContact) %>%
      pivot_wider(id_cols = c(country, AgeGrp), names_from = age_start_contact, values_from = contacts,
                  names_prefix = "cd-") %>%
      dplyr::select(-c(country, AgeGrp)) %>%
      as.matrix()
  }
  
  tmp_contacts[["home"]] <- df_contacts_home %>%
    left_join(df_pop_ages) %>%
    # Adjust only for population age structure
    mutate(contacts = standardized_contacts * prop_pop) %>%
    dplyr::select(country, AgeGrp, AgeGrpContact, contacts) %>%
    # Clean up to return expected format for sccosmo
    mutate(age_start_contact = as.character(AgeGrpContact),
           age_start_contact = gsub("\\[", "", age_start_contact),
           age_start_contact = gsub(",.*", "", age_start_contact),
           age_start_contact = as.numeric(age_start_contact)) %>%
    dplyr::select(-AgeGrpContact) %>%
    pivot_wider(id_cols = c(country, AgeGrp), names_from = age_start_contact, values_from = contacts,
                names_prefix = "cd-") %>%
    dplyr::select(-c(country, AgeGrp)) %>%
    as.matrix()
  
  tmp_contacts[["all_locations"]] <- Reduce('+', tmp_contacts)
  
  # Create average contacts. Note that this is not population-weighted.     
  avg_contact <- mean(rowSums(tmp_contacts[["all_locations"]]))
  
  l_contact_matrix <- list(m_contact = tmp_contacts[["all_locations"]],
                           n_avg_ct = avg_contact,
                           m_contact_work = tmp_contacts[["work"]],
                           m_contact_school = tmp_contacts[["school"]],
                           m_contact_other_locations = tmp_contacts[["other_locations"]],
                           m_contact_home = tmp_contacts[["home"]])
  
  add_to_cache(GLOBAL_CONTACT_MATRICES_CACHE, country, state, county, l_contact_matrix, density = density)
  
  return(l_contact_matrix)
  
}

#' Get wrangle contact matrix
#'
#' \code{wrangle_contact_matrix} returns a dataframe with all of the 
#' processed information on population age structure. 
#' 
#' @param setting character that specifies which contact matrix to return. Options: all_locations, home, non_home, school, work, other_locations.
#' @param src_country character that specifies which country contact matrix to return
#' @param v_init_age_grps vector that specifies the age bins to aggregate and return
#'
#' @export
wrangle_contact_matrix <- function(setting = "all_locations", src_country = src_country, v_init_age_grps = v_init_age_grps) {
  
  ## Load contact matrices
  df_contacts_src_country <- get_contact_matrices_world() %>%
    ## Filter to country and setting
    filter(country == src_country & location == setting) %>%
    ## Expand age to one-year bins
    rowwise() %>%
    mutate(age = list(seq(age_start, age_end, 1))) %>%
    unnest(age) %>%
    ## Expand contact age to one-year bins
    mutate(age_start_contact = (contact_grpnum*5)-5,
           age_end_contact = (contact_grpnum*5)-1,
           age_end_contact = replace(age_end_contact, age_end_contact == 79, 99)) %>%
    rowwise() %>%
    mutate(age_contact = list(seq(age_start_contact, age_end_contact, 1))) %>%
    unnest(age_contact) %>%
    ## Create age groups of interest (individual)
    mutate(AgeGrp = cut(age, breaks = c(v_init_age_grps, Inf), 
                        include.lowest = TRUE, right = FALSE)) %>%
    dplyr::select(location, country, density, standardized_contacts, non_home_dens_coeff, non_home_intercept_coeff, 
                  # residual, 
                  age_contact, AgeGrp) %>%
    group_by(location, country, age_contact, AgeGrp) %>%
    summarize_all(mean) %>%
    ## Create age groups of interest (contact)
    mutate(AgeGrpContact = cut(age_contact, breaks = c(v_init_age_grps, Inf), 
                               include.lowest = TRUE, right = FALSE)) %>%
    ungroup() %>%
    dplyr::select(location, country, density, standardized_contacts, non_home_dens_coeff, non_home_intercept_coeff, 
                  # residual, 
                  AgeGrpContact, AgeGrp) %>%
    group_by(location, country, AgeGrp, AgeGrpContact) %>%
    summarize_all(mean) %>%
    ungroup()
  return(df_contacts_src_country)
}

#' Aggregate population ages into age groups of interest
#'
#' \code{aggregate_pop_ages} returns a dataframe with all of the 
#' processed information on population age structure. 
#' 
#' @param match_country character that specifies which country population age structure to return
#' @param match_state character that specifies which state population age structure to return
#' @param match_county character that specifies which county population age structure to return
#' @param ages_to_cut vector that specifies the age bins to aggregate and return
#'
#' @export
aggregate_pop_ages <- function(match_country = "",
                               match_state = "",
                               match_county = "",
                               ages_to_cut = v_init_age_grps) {
  
  df_pop_ages <- suppressWarnings(get_population_ages(country = match_country, state = match_state, county = match_county)) %>%
    mutate(AgeGrpContact = cut(age, breaks = c(ages_to_cut, Inf), include.lowest = TRUE, right = FALSE)) %>%
    dplyr::select(country, state, county, AgeGrpContact, age_pop) %>%
    group_by(country, state, county, AgeGrpContact) %>%
    summarize_all(sum) %>%
    mutate(prop_pop = age_pop/sum(age_pop)) %>%
    dplyr::select(-age_pop)
  
  return(df_pop_ages)
}