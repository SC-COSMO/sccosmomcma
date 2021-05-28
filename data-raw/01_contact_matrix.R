library(sccosmoData)
library(stringi)
library(usmap)
library(countrycode)

### Define country and state
choose_country <- "Mexico"
choose_state <- "Mexico City"

### Get contact matrices from sccosmoData
l_contact_matrices <- sccosmoData::get_contact_matrix(country = choose_country, 
                                                      state = choose_state, 
                                                      density = 141341)
### Create .rda object of setting-specific contact matrices for an exemplary 
### population and stor it in 'data' folder
usethis::use_data(l_contact_matrices, overwrite = TRUE)
