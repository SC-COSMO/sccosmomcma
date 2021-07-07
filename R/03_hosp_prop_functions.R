#' Generate hospitalization proportions.
#'
#' \code{get_MCMA_prop_hosp} wrangles individual COVID19 data for MCMA and calculates
#' hospitalization proportions.
#' 
#' @param n_date_ini Character string. Initial date to compute proportions.
#' @param n_date_last Character string. Last date to compute proportions.
#' @param save_data Logical. Saves data in a csv file if TRUE.
#' @return 
#' A data.frame with cumulative mortality risk expected by 30 and 60 days.
#' @export
get_MCMA_prop_hosp <- function(n_date_ini  = "2020-02-24",
                               n_date_last = "2020-12-07",
                               save_data   = TRUE){
  

# Load data ---------------------------------------------------------------

  df_ssa_covid_raw <- fread("data-raw/201213COVID19MEXICO.csv", header=TRUE)
  load("data-raw/df_pop_county.Rdata")
  load("data/df_hosp_MCMA.RData")
  load("data-raw/df_pop_ZMVM.RData")
  
  n_time_stamp <- unique(df_ssa_covid_raw$FECHA_ACTUALIZACION)
  

# Wrangle data  -----------------------------------------------------------
  
  df_ssa_covid <- df_ssa_covid_raw %>%
    select(ENTIDAD_RES, MUNICIPIO_RES, FECHA_INGRESO, FECHA_SINTOMAS, 
           FECHA_DEF, EDAD, SEXO, TIPO_PACIENTE, INTUBADO, UCI, CLASIFICACION_FINAL,
           RESULTADO_LAB) %>%
    mutate(state = case_when(ENTIDAD_RES==1 ~ "Aguascalientes",
                             ENTIDAD_RES==2 ~ "Baja California",
                             ENTIDAD_RES==3 ~ "Baja California Sur",
                             ENTIDAD_RES==4 ~ "Campeche",
                             ENTIDAD_RES==5 ~ "Coahuila",
                             ENTIDAD_RES==6 ~ "Colima",
                             ENTIDAD_RES==7 ~ "Chiapas",
                             ENTIDAD_RES==8 ~ "Chihuahua",
                             ENTIDAD_RES==9 ~ "Mexico City",
                             ENTIDAD_RES==10 ~ "Durango",
                             ENTIDAD_RES==11 ~ "Guanajuato",
                             ENTIDAD_RES==12 ~ "Guerrero",
                             ENTIDAD_RES==13 ~ "Hidalgo",
                             ENTIDAD_RES==14 ~ "Jalisco",
                             ENTIDAD_RES==15 ~ "State of Mexico",
                             ENTIDAD_RES==16 ~ "Michoacan",
                             ENTIDAD_RES==17 ~ "Morelos",
                             ENTIDAD_RES==18 ~ "Nayarit",
                             ENTIDAD_RES==19 ~ "Nuevo Leon",
                             ENTIDAD_RES==20 ~ "Oaxaca",
                             ENTIDAD_RES==21 ~ "Puebla",
                             ENTIDAD_RES==22 ~ "Queretaro",
                             ENTIDAD_RES==23 ~ "Quintana Roo",
                             ENTIDAD_RES==24 ~ "San Luis Potosi",
                             ENTIDAD_RES==25 ~ "Sinaloa",
                             ENTIDAD_RES==26 ~ "Sonora",
                             ENTIDAD_RES==27 ~ "Tabasco",
                             ENTIDAD_RES==28 ~ "Tamaulipas",
                             ENTIDAD_RES==29 ~ "Tlaxcala",
                             ENTIDAD_RES==30 ~ "Veracruz",
                             ENTIDAD_RES==31 ~ "Yucatan",
                             ENTIDAD_RES==32 ~ "Zacatecas")) %>%
    mutate(entidad = case_when(ENTIDAD_RES==1 ~ "Aguascalientes",
                               ENTIDAD_RES==2 ~ "Baja California",
                               ENTIDAD_RES==3 ~ "Baja California Sur",
                               ENTIDAD_RES==4 ~ "Campeche",
                               ENTIDAD_RES==5 ~ "Coahuila",
                               ENTIDAD_RES==6 ~ "Colima",
                               ENTIDAD_RES==7 ~ "Chiapas",
                               ENTIDAD_RES==8 ~ "Chihuahua",
                               ENTIDAD_RES==9 ~ "Ciudad de México",
                               ENTIDAD_RES==10 ~ "Durango",
                               ENTIDAD_RES==11 ~ "Guanajuato",
                               ENTIDAD_RES==12 ~ "Guerrero",
                               ENTIDAD_RES==13 ~ "Hidalgo",
                               ENTIDAD_RES==14 ~ "Jalisco",
                               ENTIDAD_RES==15 ~ "Estado de México",
                               ENTIDAD_RES==16 ~ "Michoacán",
                               ENTIDAD_RES==17 ~ "Morelos",
                               ENTIDAD_RES==18 ~ "Nayarit",
                               ENTIDAD_RES==19 ~ "Nuevo León",
                               ENTIDAD_RES==20 ~ "Oaxaca",
                               ENTIDAD_RES==21 ~ "Puebla",
                               ENTIDAD_RES==22 ~ "Querétaro",
                               ENTIDAD_RES==23 ~ "Quintana Roo",
                               ENTIDAD_RES==24 ~ "San Luis Potosí",
                               ENTIDAD_RES==25 ~ "Sinaloa",
                               ENTIDAD_RES==26 ~ "Sonora",
                               ENTIDAD_RES==27 ~ "Tabasco",
                               ENTIDAD_RES==28 ~ "Tamaulipas",
                               ENTIDAD_RES==29 ~ "Tlaxcala",
                               ENTIDAD_RES==30 ~ "Veracruz",
                               ENTIDAD_RES==31 ~ "Yucatán",
                               ENTIDAD_RES==32 ~ "Zacatecas"))
  
  # Create a variable = 1 if COVID-19 is positive
  df_ssa_covid$covid <- ifelse(df_ssa_covid$CLASIFICACION_FINAL %in% c(1,2,3), 1, 0)
  # df_ssa_covid$covid <- ifelse(df_ssa_covid$RESULTADO_LAB == 1, 1, 0)
  
  # Subset the data for COVID-19 cases only
  df_ssa_covid <- subset(df_ssa_covid, covid==1)
  
  # Format date variables as.Date() and create date_dx and date_sx variables
  df_ssa_covid$date_dx <- as.Date(df_ssa_covid$FECHA_INGRESO, format = "%Y-%m-%d")
  df_ssa_covid$date_sx <- as.Date(df_ssa_covid$FECHA_SINTOMAS, format = "%Y-%m-%d")
  df_ssa_covid$date_dead <- as.Date(df_ssa_covid$FECHA_DEF, format = "%Y-%m-%d")
  
  # Create variable hosp_ind = 1 if patient was hospitalized
  # Create variable icu_ind = 1 if patient required ICU
  # Create variable vent_ind = 1 if patient was in ventilator
  # Create variable dead_ind = 1 if patient died
  df_ssa_covid <- df_ssa_covid %>%
    mutate(hosp_ind = ifelse(TIPO_PACIENTE == 2, 1, 0)) %>%
    mutate(icu_ind  = ifelse(UCI == 1, 1, 0)) %>%
    mutate(vent_ind = ifelse(INTUBADO == 1, 1, 0)) %>%
    mutate(novent_ind = ifelse(INTUBADO != 1, 1, 0)) %>%
    mutate(dead_ind = ifelse(is.na(date_dead), 0, 1)) %>%
    mutate(time_to_dead = date_dead-date_dx)
  
  ####  Paste counties and population info
  df_pop_county <- df_pop_county %>%
    mutate(county_id = as.numeric(county_id))
  
  df_ssa_covid <- df_ssa_covid %>% 
    mutate(county_id = formatC(MUNICIPIO_RES, width = 3, flag = "0"), 
           state_id = formatC(ENTIDAD_RES, width = 2, flag = "0"), 
           county_id = paste0(state_id, county_id), 
           county_id= as.numeric(county_id)) %>%
    rename(age = EDAD) %>%
    left_join(df_pop_county, by = c("county_id" = "county_id")) %>%
    rename(entidad = "entidad.x")
  
  # Vector of ZMVM  county ids
  v_zmvm_id <- c("9002", "9003", "9004", "9005", "9006",
                 "9007", "9008", "9009", "9010", "9011",
                 "9012", "9013", "9014", "9015", "9016",
                 "9017", "13069",
                 "15002", "15009", "15010", "15011",
                 "15013", "15015", "15016", "15017", 
                 "15020", "15022", "15023", "15024", 
                 "15025", "15028", "15029", "15030", 
                 "15031", "15033", "15034", "15035", 
                 "15036", "15037", "15038", "15039", 
                 "15044", "15046", "15050", "15053",
                 "15057", "15058", "15059", "15060",
                 "15061", "15065", "15068", "15069",
                 "15070", "15075", "15081", "15083",
                 "15084", "15089", "15091", "15092",
                 "15093", "15094", "15095", "15096",
                 "15099", "15100", "15103", "15104",
                 "15108", "15109", "15112", "15120",
                 "15121", "15122", "15125")
  
  # Create dummy for each county in ZMVM and get a subset of the
  # data for MCMA only
  df_ssa_covid_MCMA <- df_ssa_covid %>%
    mutate(zmvm = ifelse(county_id %in% v_zmvm_id, 1, 0)) %>%
    filter(zmvm == 1) %>%
    mutate(entidad = case_when(zmvm == 1 ~ "ZMVM"), 
           state = case_when(zmvm == 1 ~ "MCMA"))
  
  # Add ZMVM population
  df_pop_ZMVM <- df_pop_ZMVM %>%
    select(entidad, population) 
  
  df_ssa_covid_MCMA <- df_ssa_covid_MCMA %>%
    left_join(df_pop_ZMVM, by = "entidad") %>%
    rename(population = "population.y") 
  

# Age groups --------------------------------------------------------------

  # # Vector of level names for each age group based on SC-COSMO model
  # v_init_age_grps <- c(0, 5, 15, 25, 45, 55, 65, 70)
  # v_names_ages <- ordered(c(paste(v_init_age_grps[-length(v_init_age_grps)], 
  #                                 (v_init_age_grps[-1]-1), sep = "-"), 
  #                           paste0(v_init_age_grps[length(v_init_age_grps)], "+")),
  #                         c(paste(v_init_age_grps[-length(v_init_age_grps)], 
  #                                 (v_init_age_grps[-1]-1), sep = "-"), 
  #                           paste0(v_init_age_grps[length(v_init_age_grps)], "+")))
  # v_names_age_groups <- paste(v_names_ages, "years")  
  # 
  # df_ssa_covid_MCMA <- df_ssa_covid_MCMA %>%
  #   mutate(age_groups = cut(age, c(v_init_age_grps, Inf), 
  #                         include.lowest = TRUE, right = FALSE))
  # levels(df_ssa_covid_MCMA$age_groups) <- v_names_age_groups
  # 
  # df_ssa_covid_MCMA_idx_gps <- df_ssa_covid_MCMA %>%
  #   filter(date_dx >= n_date_ini & date_dx <= n_date_last) %>%
  #   group_by(date_dx, age_groups) %>%
  #   summarise(idx  = n(),
  #             dead = sum(dead_ind),
  #             hosp = sum(hosp_ind),
  #             vent = sum(vent_ind),
  #             novent = sum(novent_ind),
  #             icu  = sum(icu_ind)) %>%
  #   rename(date = date_dx) %>%
  #   ungroup() %>%
  #   left_join(df_pop_size, by = c("age_groups", "date"))
  # 
  # df_ssa_covid_MCMA_grouped <- df_ssa_covid_MCMA %>%
  #   filter(date_dx >= n_date_ini & date_dx <= n_date_last) %>%
  #   group_by(age_groups) %>%
  #   summarise(idx  = n(),
  #             dead = sum(dead_ind),
  #             hosp = sum(hosp_ind),
  #             vent = sum(vent_ind),
  #             novent = sum(novent_ind),
  #             icu  = sum(icu_ind)) %>%
  #   mutate(p_dead   = dead/idx,
  #          p_hosp   = hosp/idx,
  #          p_vent   = vent/idx,
  #          p_novent = novent/idx,
  #          p_icu    = icu/idx) %>%
  #   mutate(rr_death8 = p_dead/p_dead[8],
  #          rr_hosp8 = p_hosp/p_hosp[8],
  #          rr_vent8 = p_vent/p_vent[8])
  # 
  # df_pop_size_grouped <- df_pop_size %>%
  #   group_by(age_groups) %>%
  #   summarise(mean_age_pop = mean(age_pop),
  #             mean_tot_pop = mean(tot_pop)) %>%
  #   mutate(prop_age = mean_age_pop/mean_tot_pop)
  # 
  # df_ssa_covid_MCMA_grouped <- df_ssa_covid_MCMA_grouped %>%
  #   left_join(df_pop_size_grouped, by ="age_groups")

# Generate proportions ----------------------------------------------------

  # Drop observations with time-to-dead < 0
  # df_ssa_covid_MCMA <- df_ssa_covid_MCMA %>%
  # filter(time_to_dead >= 0 | is.na(time_to_dead))
  
  df_ssa_covid_MCMA_idx <- df_ssa_covid_MCMA %>%
    filter(date_dx >= n_date_ini & date_dx <= n_date_last) %>%
    group_by(date_dx) %>%
    summarise(idx  = n(),
              dead = sum(dead_ind),
              hosp = sum(hosp_ind),
              vent = sum(vent_ind),
              novent = sum(novent_ind),
              icu  = sum(icu_ind)) %>%
    complete(date_dx = seq.Date(as.Date(n_date_ini), as.Date(n_date_last), by="day"), # Create a sequence of dates
             fill = list(idx  = 0,
                         dead = 0,
                         hosp = 0,
                         vent = 0,
                         novent = 0,
                         icu  = 0)) %>% # Fill the dates without cases with zero
    ungroup() %>%
    rename(date = date_dx)
  
  df_ssa_covid_MCMA_3wdead <- df_ssa_covid_MCMA %>%
    filter(date_dx >= n_date_ini & date_dx <= n_date_last & 
             !is.na(date_dead)) %>%
    mutate(dead_ind_3w = ifelse(time_to_dead <= 21, 1, 0)) %>%
    group_by(date_dx) %>%
    summarise(dead_3w = sum(dead_ind_3w)) %>%
    complete(date_dx = seq.Date(as.Date(n_date_ini), as.Date(n_date_last), by="day"), # Create a sequence of dates
             fill = list(dead_3w  = 0)) %>% # Fill the dates without cases with zero
    ungroup() %>%
    rename(date = date_dx)
  
  df_ssa_covid_MCMA_dead_hosp <- df_ssa_covid_MCMA %>%
    filter(date_dx >= n_date_ini & date_dx <= n_date_last &  
             !is.na(date_dead) & hosp_ind == 1) %>%
    group_by(date_dx) %>%
    summarise(dead_hosp = sum(dead_ind)) %>%
    complete(date_dx = seq.Date(as.Date(n_date_ini), as.Date(n_date_last), by="day"), # Create a sequence of dates
             fill = list(dead_hosp  = 0)) %>% # Fill the dates without cases with zero
    ungroup() %>%
    rename(date = date_dx)
  
  df_ssa_covid_MCMA_dead_hosp_vent <- df_ssa_covid_MCMA %>%
    filter(date_dx >= n_date_ini & date_dx <= n_date_last &  
             !is.na(date_dead) & hosp_ind == 1 & vent_ind == 1) %>%
    group_by(date_dx) %>%
    summarise(dead_hosp_vent = sum(dead_ind)) %>%
    complete(date_dx = seq.Date(as.Date(n_date_ini), as.Date(n_date_last), by="day"), # Create a sequence of dates
             fill = list(dead_hosp_vent  = 0)) %>% # Fill the dates without cases with zero
    ungroup() %>%
    rename(date = date_dx)
  
  df_ssa_covid_MCMA_dead_hosp_novent <- df_ssa_covid_MCMA %>%
    filter(date_dx >= n_date_ini & date_dx <= n_date_last &  
             !is.na(date_dead) & hosp_ind == 1 & vent_ind == 0) %>%
    group_by(date_dx) %>%
    summarise(dead_hosp_novent = sum(dead_ind)) %>%
    complete(date_dx = seq.Date(as.Date(n_date_ini), as.Date(n_date_last), by="day"), # Create a sequence of dates
             fill = list(dead_hosp_novent  = 0)) %>% # Fill the dates without cases with zero
    ungroup() %>%
    rename(date = date_dx)
  
  df_ssa_covid_MCMA_prop <- left_join(df_ssa_covid_MCMA_idx, 
                                      df_ssa_covid_MCMA_3wdead,
                                      by = "date") %>%
    left_join(df_ssa_covid_MCMA_dead_hosp,
              by = "date") %>%
    left_join(df_ssa_covid_MCMA_dead_hosp_vent,
              by = "date") %>%
    left_join(df_ssa_covid_MCMA_dead_hosp_novent,
              by = "date") %>%
    mutate(dead             = ifelse(is.na(dead), 0 , dead),
           dead_3w          = ifelse(is.na(dead_3w), 0 , dead_3w),
           dead_hosp        = ifelse(is.na(dead_hosp), 0 , dead_hosp),
           dead_hosp_vent   = ifelse(is.na(dead_hosp_vent), 0 , dead_hosp_vent),
           dead_hosp_novent = ifelse(is.na(dead_hosp_novent), 0 , dead_hosp_novent),
           p_hosp           = ifelse(idx != 0, hosp/idx, 0),
           p_vent_hosp      = ifelse(hosp != 0, vent/hosp, 0),
           p_novent_hosp    = ifelse(hosp != 0, novent/hosp, 0),
           p_icu_hosp       = ifelse(hosp != 0, icu/hosp, 0),
           p_vent           = ifelse(idx != 0, vent/idx, 0),
           p_icu            = ifelse(idx != 0, icu/idx, 0),
           cfr              = ifelse(idx != 0, dead/idx, 0),
           p_dead_hosp      = ifelse(idx != 0, dead_hosp/idx, 0))
  
  # Order data.frame
  df_ssa_covid_MCMA_prop <- df_ssa_covid_MCMA_prop[order(df_ssa_covid_MCMA_prop$date),]
  
  df_ssa_covid_MCMA_prop <- df_ssa_covid_MCMA_prop %>%
    mutate(date0 = as.numeric(date - date[1]))
  
  # ADIP's data
  df_adip_covid_MCMA_prop <- df_hosp_MCMA %>%
    mutate(date = as.Date(dates)) %>%
    left_join(df_ssa_covid_MCMA_idx,
              by = "date") %>%
    filter(date >= n_date_ini & date <= n_date_last) %>%
    mutate(p_hosp = ifelse(idx != 0, hosp_tot/idx, 0),
           p_vent_hosp = bed_icu_tot/hosp_tot,
           p_novent_hosp = bed_noicu_tot/hosp_tot,
           date0 = as.numeric(date - date[1])) %>%
    select(date, date0, idx, dead, hosp, hosp_tot, vent, novent, icu, p_hosp, p_vent_hosp, p_novent_hosp)
  
  l_hosp_prop <- list(ssa  = df_ssa_covid_MCMA_prop,
                      adip = df_adip_covid_MCMA_prop)
  
  if(save_data){
    save(l_hosp_prop, file = paste0("data/l_hosp_prop_",n_time_stamp,".RData"))
  }
  
  return(l_hosp_prop)
  
}

#' Plot hospitalization proportions
#'
#' \code{plot_hosp_prop} plot observed and estimated hospitalization proportions
#' for MCMA. 
#' 
#' @param df_hosp_prop Data.frame. Data of observed and estimated hospitalization
#' proportions.
#' @param save_plot Logical. Saves the plot if TRUE.
#' @return 
#' A ggplot object.
plot_hosp_prop <- function(df_hosp_prop,
                           save_plot = FALSE){
  
  df_hosp_prop$type <- as.factor(df_hosp_prop$type)
  abbrev_state <- unique(df_hosp_prop$state)
  Outcome <- unique(df_hosp_prop$Outcome)
  
  # Name of plot
  if(Outcome == "Total hospitalizations"){
    plot_name <- "TotHospProp_"
  }else if(Outcome == "Hospitalizations with ventilator"){
    plot_name <- "ICUHospProp_"
  }else if(Outcome == "Hospitalizations without ventilator"){
    plot_name <- "NoICUHospProp_"
  }
  
  # Colors
  v_colors_names <- levels(df_hosp_prop$type)
  v_colors <- c("red", "black")
  names(v_colors) <- v_colors_names
  
  # Fill
  v_colors_fill <- c("red", NA)
  names(v_colors_fill) <- v_colors_names
  
  # Linetype
  v_linetype <- c("solid", 'blank')
  names(v_linetype) <- v_colors_names
  
  # Shape 
  v_shape <- c(NA, 8)
  names(v_shape) <- v_colors_names
  
  plot_temp <- ggplot(subset(df_hosp_prop, type == "Estimated"),
         aes(x = date, y = prop, 
             ymin = lb, ymax = ub,
             color = type, fill = type, shape = type)) +
    geom_line(size = 1.2) +
    geom_point(data = subset(df_hosp_prop, type == "Observed"),
               size = 2) +
    geom_ribbon(aes(x = date,
                    y = prop,
                    ymin = lb,
                    ymax = ub,
                    fill = type),
                alpha = 0.2) +
    scale_shape_manual(values = v_shape) +
    scale_color_manual(values = v_colors) +
    scale_fill_manual(values = v_colors_fill) +
    guides(fill = guide_legend(title = ""),
           shape = guide_legend(title = ""),
           color = guide_legend(title = "",
                                ncol = 2,
                                override.aes = list(
                                  shape    = v_shape,
                                  linetype = v_linetype,
                                  color    = v_colors,
                                  fill     = v_colors_fill))) +
    scale_x_date("",
                 date_labels = "%m/%d", #Change the format of the date in the x axis
                 breaks = number_ticks(14)) +
    ylab("Proportion") +
    theme(text = element_text(size = 13),
          plot.caption = element_text(hjust = 0,
                                      size = 12),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     colour = "black"),
          axis.text.y = element_text(colour = "black"),
          panel.background = element_rect(fill = "white",
                                          colour = "gray",
                                          size = 0.15,
                                          linetype = "solid"),
          panel.border = element_rect(colour = "black",
                                      fill = NA,
                                      size = 0.7),
          legend.position = "bottom",
          legend.text = element_text(size = 10),
          legend.justification = "center",
          legend.direction = "horizontal",
          legend.title = element_blank(),
          legend.key = element_rect(fill   = "transparent",
                                    colour = "transparent")) +
    labs(title = Outcome)
  
  if(save_plot){
    ggsave(paste0("figs/",plot_name,"observed_vs_estimated_",abbrev_state,".pdf"),
           width = 10, height = 7)
  }
  
  return(plot_temp)
  
}

