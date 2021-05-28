#' Generate household possibilities
#' 
#' \code{gen_hh_n} generates household transmission (i.e., 
#' force of infection).
#' 
#' @param n_hhsize Household size
#' @param v_names_states_hh vector with state names of household model
#' @return 
#' A data.frame with possible combinations of household members within each 
#' epidemic compartment.
#' @export
gen_hh_n <- function(n_hhsize, v_names_states_hh){
  n_states <- length(v_names_states_hh)
  df_possibilities <- data.frame(id = 1, V1 = 0:n_hhsize)
  for(i in 1:(n_states - 1)){ # i <- 2
    df_possibilities <- left_join(df_possibilities, data.frame(id = 1, temp = 0:n_hhsize), by = "id")
    colnames(df_possibilities)[i+2] <- paste0("V", i+1)
  }
  df_possibilities <-  df_possibilities %>%
    dplyr::select(-id)
  v_names_cols <- colnames(df_possibilities)
  df_possibilities$new_n_hhsize <- rowSums(df_possibilities)
  df_possibilities <- df_possibilities %>% 
    filter(new_n_hhsize == n_hhsize) %>%
    dplyr::select(-new_n_hhsize) %>% 
    group_by_all() %>%
    summarise_all(list(mean)) %>%
    ungroup() %>%
    arrange_all(desc)
  colnames(df_possibilities) <- v_names_states_hh
  return(df_possibilities)
}

#' Generate household matrices
#' 
#' \code{gen_household_matrices} generates transition matrices for within 
#' household epidemics.
#' @param n_hhsize Household size
#' @param n_hh_mod Number of household model states
#' @param v_names_states_hh Vector with names of household MC SEIR model
#' @param v_names_exp_hh Vector with names of exposed states
#' @param v_names_inf_hh Vector with names of infectious states
#' @param v_index_hh_sus Indexing column vector of susceptibles
#' @param v_index_hh_exp Indexing column vector of exposed
#' @param v_index_hh_inf Indexing column vector of infectious
#' @param v_index_hh_rec Indexing column vector of recovered
#' @param m_possibilities Matrix with possible combinations of household 
#' @param r_sigma  Daily rate of progression of exposed individuals 
#' (Latent period)
#' @param r_gamma Daily rate of recovery of infectious individuals 
#' (Infectiousness period).
#' @return 
#' A list with three transition matrices: 1) community transmission matrix; 2) 
#' household transmission matric; and 3) hosehold recovered matrix
#' @export
gen_household_matrices_mc_seir <- function(n_hhsize, 
                                           n_hh_mod,
                                           v_names_states_hh,
                                           v_names_exp_hh,
                                           v_names_inf_hh,
                                           v_names_rec_hh,
                                           v_index_hh_sus,
                                           v_index_hh_exp,
                                           v_index_hh_inf,
                                           v_index_hh_rec,
                                           df_possibilities,
                                           r_sigma,
                                           r_gamma
                                           ){
  m_possibilities <- as.matrix(df_possibilities)
  #### Naming vectors ####
  v_hh_names <- as.matrix(unite(df_possibilities, col = "names", sep = ""))
  # paste("HH",m_possibilities[, 1:n_states], m_possibilities[,2], sep = "")
  ### Names of household members by class names
  v_names_hh  <- paste("HH", v_hh_names, sep = "")
  n_hh_mod    <- length(v_names_hh)
  ### Derivative names of household members by class
  v_names_dhh <- paste("dHH", v_hh_names, sep = "")
  
  n_exp_states <- length(v_names_exp_hh)
  n_inf_states <- length(v_names_inf_hh)
  
  v_names_exp_inf_hh <- c(v_names_exp_hh, v_names_inf_hh)
  v_names_inf_rec_hh <- c(v_names_inf_hh, v_names_rec_hh)
  
  #### Indexing vectors ####
  ## Column index for susceptibles and infectious
  v_index_sus_inf_hh <- c(v_index_hh_sus, v_index_hh_inf)
  ## Column index for exposed and infectious
  v_index_exp_inf_hh <- c(v_index_hh_exp, v_index_hh_inf)
  ## Column index for infectious and recovered
  v_index_inf_rec_hh <- c(v_index_hh_inf, v_index_hh_rec)
  ## Column index for susceptibles, exposed and infectious
  v_index_sus_exp_inf_hh <- c(v_index_hh_sus, v_index_hh_exp, v_index_hh_inf)
  ## Indices with active susceptibles
  v_index_keep_sus    <- (m_possibilities[, v_index_hh_sus] > 0)
  ## Indices with active exposed
  v_index_keep_exp       <- rowSums((m_possibilities[, v_index_hh_exp, drop=FALSE] > 0)) > 0
  v_index_keep_exp_first <- (m_possibilities[, v_index_hh_exp[1]] > 0)
  v_index_keep_exp_last  <- (m_possibilities[, v_index_hh_exp[n_exp_states]] > 0)
  m_index_keep_exp       <- as.matrix((m_possibilities[, v_index_hh_exp] > 0), 
                                      ncol = n_exp_states) # matrix of indices to keep exposed by class
  ## Indices with active infectious
  v_index_keep_inf       <- rowSums((m_possibilities[, v_index_hh_inf, drop=FALSE] > 0)) > 0 
  v_index_keep_inf_first <- (m_possibilities[, v_index_hh_inf[1]] > 0)
  v_index_keep_inf_last  <- (m_possibilities[, v_index_hh_inf[n_inf_states]] > 0)
  m_index_keep_inf       <- as.matrix((m_possibilities[, v_index_hh_inf] > 0), 
                                      ncol = n_inf_states) # matrix of indices to keep exposed by class
  ## Indices with active exposed and infectious
  m_index_keep_exp_inf  <- (m_possibilities[, v_index_exp_inf_hh] > 0) # matrix of indices to keep exposed by classv_index_exp_inf_hh
  ## Indices with active infectious and recovered
  m_index_keep_inf_rec  <- (m_possibilities[, v_index_inf_rec_hh] > 0) # matrix of indices to keep exposed by classv_index_inf_rec_hh
  ## Indices with active exposed for HH transmission
  v_index_keep_exp_hh <- ((rowSums(m_possibilities[, v_index_hh_exp, drop=FALSE]) > 0) & 
                            (rowSums(m_possibilities[, v_index_hh_inf, drop=FALSE]) > 0))
  v_index_keep_exp_hh_first <- ((m_possibilities[, v_index_hh_exp[1]] > 0) & 
                                  (rowSums(m_possibilities[, v_index_hh_inf, drop=FALSE]) > 0))
  
  ## Indices with active infectious for HH transmission
  v_index_keep_inf_hh <- (m_possibilities[, v_index_hh_inf] > 1)
  ## Indices with active recovered
  v_index_keep_rec    <- (m_possibilities[, v_index_hh_rec] > 0)
  ## Indices with active susceptibles and infectious
  v_index_keep_tau <- ((m_possibilities[, v_index_hh_sus] > 0) & 
                         (rowSums(m_possibilities[, v_index_hh_inf, drop=FALSE]) > 0))
  
  #### Community transmission matrix ####
  m_comm_trans <- matrix(0, 
                         nrow = n_hh_mod, ncol = n_hh_mod, 
                         dimnames = list(v_names_hh, v_names_hh))
  diag(m_comm_trans)[v_index_keep_sus] <- -1* # r_beta*
    m_possibilities[v_index_keep_sus, v_index_hh_sus]
  if(is.null(dim(m_comm_trans[v_index_keep_sus, v_index_keep_exp_first]))){ 
    m_comm_trans[v_index_keep_sus, v_index_keep_exp_first] <- 1* #r_beta*
      m_possibilities[v_index_keep_sus, v_index_hh_sus]  
  } else {
    diag(m_comm_trans[v_index_keep_sus, v_index_keep_exp_first]) <- 1* #r_beta*
      m_possibilities[v_index_keep_sus, v_index_hh_sus]  
  }
  
  
  #### Household matrices ####
  ### Household transmission matrix
  m_hh_trans <- matrix(0, 
                       nrow = n_hh_mod, ncol = n_hh_mod, 
                       dimnames = list(v_names_hh, v_names_hh))
  diag(m_hh_trans)[v_index_keep_tau] <- -1* # -tau_avg*
    m_possibilities[v_index_keep_tau, v_index_hh_sus] * 
    rowSums(m_possibilities[v_index_keep_tau, v_index_hh_inf, drop=FALSE])
  if(is.null(dim(m_hh_trans[v_index_keep_tau, v_index_keep_exp_hh_first]))){ 
    # If resulting object is a scalar, don't use 'diag'
    m_hh_trans[v_index_keep_tau, v_index_keep_exp_hh_first] <- 1* # tau_avg*
      m_possibilities[v_index_keep_tau, v_index_hh_sus] * 
      rowSums(m_possibilities[v_index_keep_tau, v_index_hh_inf, drop=FALSE])  
  } else{ # If resulting object is a matrix, use 'diag'
    diag(m_hh_trans[v_index_keep_tau, v_index_keep_exp_hh_first]) <- 1* # tau_avg*
      m_possibilities[v_index_keep_tau, v_index_hh_sus] * 
      rowSums(m_possibilities[v_index_keep_tau, v_index_hh_inf, drop=FALSE])  
  }
  
  ### Household progression matrix
  a_hh_prog_out <- array(0, 
                         dim = c(n_hh_mod, n_hh_mod, n_exp_states),
                         dimnames = list(v_names_hh, v_names_hh, v_names_exp_hh))
  m_hh_prog_out_flat <- array(0, 
                              dim = c(n_hh_mod, n_hh_mod),
                              dimnames = list(v_names_hh, v_names_hh))
  for(i in 1:n_exp_states){ # i <- 1
    diag(a_hh_prog_out[, , i])[m_index_keep_exp[, i]] <- diag(a_hh_prog_out[, , i])[m_index_keep_exp[, i]] - 
      (r_sigma*n_exp_states)*m_possibilities[m_index_keep_exp[, i], v_index_hh_exp[i]]
    m_hh_prog_out_flat <- m_hh_prog_out_flat + a_hh_prog_out[, , i]
  }
  
  m_hh_prog_in_flat <- compute_inflows(n_hhsize = n_hhsize,
                                       n_hh_mod = n_hh_mod,
                                       v_names_hh = v_names_hh,
                                       v_names_states = v_names_states_hh,
                                       v_names_source = v_names_exp_hh,
                                       v_names_destiny = v_names_inf_hh,
                                       v_names_source_destiny = v_names_exp_inf_hh[-length(v_names_exp_inf_hh)], #### CHECK,
                                       v_index_hh_source = v_index_hh_exp,
                                       v_index_hh_destiny = v_index_hh_inf,
                                       v_index_keep_source = v_index_keep_exp,
                                       v_index_souce_destiny_hh = v_index_exp_inf_hh,
                                       m_possibilities = m_possibilities,
                                       r_flow = r_sigma,
                                       n_source_states = n_exp_states,
                                       m_hh_out = m_hh_prog_out_flat)
  
  m_hh_prog <- m_hh_prog_out_flat + m_hh_prog_in_flat
  
  ### Household recovered matrix
  a_hh_recov_out <- array(0, 
                          dim = c(n_hh_mod, n_hh_mod, n_inf_states),
                          dimnames = list(v_names_hh, v_names_hh, v_names_inf_hh))
  
  m_hh_recov_out_flat <- array(0, 
                               dim = c(n_hh_mod, n_hh_mod),
                               dimnames = list(v_names_hh, v_names_hh))
  
  for(i in 1:n_inf_states){ # i <- 1
    diag(a_hh_recov_out[, , i])[m_index_keep_inf[, i]] <- diag(a_hh_recov_out[, , i])[m_index_keep_inf[, i]] - 
      (r_gamma*n_inf_states)*m_possibilities[m_index_keep_inf[, i], v_index_hh_inf[i]]
    
    m_hh_recov_out_flat <- m_hh_recov_out_flat + a_hh_recov_out[, , i]
  }
  
  m_hh_recov_in_flat <- compute_inflows(n_hhsize = n_hhsize,
                                        n_hh_mod = n_hh_mod ,
                                        v_names_states = v_names_states_hh,
                                        v_names_hh = v_names_hh,
                                        v_names_source = v_names_inf_hh,
                                        v_names_destiny = v_names_rec_hh,
                                        v_names_source_destiny = v_names_inf_rec_hh,
                                        v_index_hh_source = v_index_hh_inf,
                                        v_index_hh_destiny = v_index_hh_rec,
                                        v_index_keep_source = v_index_keep_inf,
                                        v_index_souce_destiny_hh = v_index_inf_rec_hh,
                                        m_possibilities = m_possibilities,
                                        r_flow = r_gamma,
                                        n_source_states = n_inf_states,
                                        m_hh_out = m_hh_recov_out_flat)
  
  m_hh_recov <- m_hh_recov_out_flat + m_hh_recov_in_flat
  
  
  ### Household waning matrix
  m_hh_waning <- matrix(0, 
                        nrow = n_hh_mod, ncol = n_hh_mod, 
                        dimnames = list(v_names_hh, v_names_hh))
  diag(m_hh_waning)[v_index_keep_rec] <- -1* #-r_omega*
    m_possibilities[v_index_keep_rec, v_index_hh_rec]
  if(is.null(dim(m_hh_waning[v_index_keep_rec, v_index_keep_sus]))){
    # If resulting object is a scalar, don't use 'diag'
    m_hh_waning[v_index_keep_rec, v_index_keep_sus] <- 1* #r_omega*
      m_possibilities[v_index_keep_rec, v_index_hh_rec]
  } else {
    # If resulting object is a matrix, use 'diag'
    diag(m_hh_waning[v_index_keep_rec, v_index_keep_sus]) <- 1* #r_omega*
      m_possibilities[v_index_keep_rec, v_index_hh_rec]
  }
  
  # rowSums(m_comm_trans)
  # rowSums(m_hh_trans)
  # rowSums(m_hh_prog)
  # rowSums(m_hh_recov)
  # rowSums(m_hh_waning)
  
  
  return(list(m_comm_trans = m_comm_trans,
              m_hh_trans   = m_hh_trans,
              m_hh_prog    = m_hh_prog,
              m_hh_recov   = m_hh_recov,
              m_hh_waning  = m_hh_waning))
}

#' Compute inflows considering competing risks using convolution of binomials

#' \code{compute_inflows} computes inflows 
#' @param n_hhsize Household size
#' @param n_hh_mod Number of household model states
#' @param v_names_hh Vector with all names of household MC SEIR model
#' @param v_names_states Vector with names of household MC SEIR model
#' @param v_names_source Vector with names of source states
#' @param v_names_destiny Vector with names of destiny states
#' @param v_names_source_destiny Vector with names of source and destiny states
#' @param v_index_hh_source Vector with indexes of source states
#' @param v_index_hh_destiny Vector with indexes of destiny states
#' @param v_index_keep_source Vector with indexes of destiny states to keep
#' @param v_index_souce_destiny_hh Vector with indexes of source and destiny 
#' states of the household model
#' @param m_possibilities Matrix with possible combinations of household 
#' @param r_flow Daily rate of flow
#' @param n_source_states Number of source states
#' @param m_hh_out Matrix with outflows from the source states, which are the 
#' flows that we are going to balance with the inflows we will compute inside
#' the function
#' @return 
#' A matrix with inflow rates to the destiny states
#' @export
compute_inflows <- function(n_hhsize,
                            n_hh_mod,
                            v_names_hh,
                            v_names_states,
                            v_names_source,
                            v_names_destiny,
                            v_names_source_destiny,
                            v_index_hh_source,
                            v_index_hh_destiny,
                            v_index_keep_source,
                            v_index_souce_destiny_hh,
                            m_possibilities,
                            r_flow,
                            n_source_states,
                            m_hh_out){
  m_hh_in  <- array(0, 
                    dim = c(n_hh_mod, n_hh_mod),
                    dimnames = list(v_names_hh, v_names_hh))
  
  v_index_keep_source_hh_destiny_first_or <- (rowSums(m_possibilities[, v_index_hh_source[-1], drop=FALSE] > 0) | 
                                                (rowSums(m_possibilities[, v_index_hh_destiny[1], drop=FALSE]) > 0))
  
  v_index_hh_signature <- which(!(v_names_states %in% v_names_source_destiny))
  
  rownames(m_possibilities) <- v_names_hh
  df_possibilities <- as.data.frame(m_possibilities)
  
  m_flow_source  <- m_possibilities[v_index_keep_source, ]
  m_flow_destiny <- m_possibilities[v_index_keep_source_hh_destiny_first_or, ]
  m_flow_destiny <- rbind(m_flow_destiny, 
                          m_flow_source[rownames(m_flow_source)[(!rownames(m_flow_source) %in% rownames(m_flow_destiny))], ]) # these rows represent the states where theres is no progression. It could have been progression but there wasn't for this particular case. i.e., the fraction of people that doesn't move
  
  v_names_HH_flow_source <-  paste("HH", 
                                   as.matrix(tidyr::unite(df_possibilities[v_index_keep_source, ], 
                                                          col = "names", sep = "")), 
                                   sep = "")
  v_names_HH_flow_destiny <-  paste("HH", 
                                    as.matrix(tidyr::unite(df_possibilities[v_index_keep_source_hh_destiny_first_or, ], 
                                                           col = "names", sep = ""))
                                    , sep = "")
  v_names_HH_flow_destiny <- c(v_names_HH_flow_destiny, 
                               v_names_HH_flow_source[!(v_names_HH_flow_source %in% v_names_HH_flow_destiny)])
  
  m_p_flow <- matrix(0, 
                     nrow = length(v_names_HH_flow_source), 
                     ncol = length(v_names_HH_flow_destiny), 
                     dimnames = list(v_names_HH_flow_source, 
                                     v_names_HH_flow_destiny))
  # View(m_p_flow)
  
  p_flow <- 1-exp(-r_flow*n_source_states*1)
  
  for(i in 1:nrow(m_flow_source)){ # i = 41
    # print(rownames(m_flow_source)[i])
    # We want to determine whether the row is binomial or a convolution of binomials
    # and how many columns we are going to need and which columns we need.
    curr_row      <- m_flow_source[i, ]
    name_curr_row <- rownames(m_flow_source)[i]
    # For all the columns that are not first element in the source through the 
    # first element in the destiny, we want to fix their values and total them 
    # up. The household size minus that total is the number of people that can 
    # move.
    signature <- curr_row[v_index_hh_signature] # c(1, (v_index_hh_destiny[1]+1))
    n_movers  <- n_hhsize - sum(signature)
    
    # We want to scan across all elements of the source if ony one of those 
    # numbers is greater than 0, then we are binomial; otherwise, we are 
    # convolution.
    binomial_flag <- FALSE
    if(n_movers == 1){
      binomial_flag <- TRUE
    } else{
      binomial_flag <- sum(curr_row[v_names_source] == n_movers) == 1
    }
    
    # If we are binomial, then the number of states that we need to find is 
    # number of elements in source + 1, which is the diagonal state. 
    if(binomial_flag){
      index_source  <- which(names(curr_row) == names(which(curr_row[v_names_source]>0)))
      index_destiny <- index_source + 1
      # v_index_not_movers <- names(curr_row)[!(names(curr_row) %in% c(names(signature), 
      #                                                                names(curr_row[index_source]), 
      #                                                                names(curr_row[index_destiny])))]
      m_destination_cols <- matrix(0, 
                                   nrow = n_movers + 1,
                                   ncol = length(curr_row), 
                                   dimnames = list(1:(n_movers+1), names(curr_row)))
      t_temp <- t(m_destination_cols)
      t_temp[names(signature), ] <- signature
      m_destination_cols         <- t(t_temp)
      for(mover in 0:n_movers){
        m_destination_cols[mover + 1, index_source]  <- n_movers - mover
        m_destination_cols[mover + 1, index_destiny] <- mover
      }
      v_names_hh_destiny <- paste("HH", 
                                  as.matrix(tidyr::unite(as.data.frame(m_destination_cols),
                                                         col = "names", sep = "")), 
                                  sep = "")
      
      m_p_flow[name_curr_row, v_names_hh_destiny] <- dbinom(0:n_movers, 
                                                            n_movers, 
                                                            prob = p_flow)
    } else{
      # Loop over `m_flow_destiny` columns, gather a list of them that we have to
      # test using the backwards algorithm.
      # Rules:
      #  1. Reject if 
      #     a. first Source column greater than curr_row first source column,
      #     b. destiny signature columns are not equal. to source signature columns
      # For the shorter list, we will use the backwards algorithm
      m_keep_destiny <- cbind(m_flow_destiny, Keep = 0)
      for(j in 1:nrow(m_keep_destiny)){ # j <- 7
        curr_row_destiny <- m_keep_destiny[j, ]
        flag1 <- curr_row_destiny[v_names_source[1]] <= curr_row[v_names_source[1]] 
        flag2 <- sum(curr_row_destiny[names(signature)] == curr_row[names(signature)]) == length(names(signature))
        m_keep_destiny[j, "Keep"] <- flag1 & flag2
      }
      m_keep_destiny_limited <- m_keep_destiny[m_keep_destiny[, "Keep"]==1, ]
      # v_names_source_destiny <- c(v_names_source, v_names_destiny[1])
      for(j in 1:nrow(m_keep_destiny_limited)){ # j <- 2
        v_curr_delta <- vector(mode = "numeric", length = (length(v_names_source_destiny)-1))
        curr_delta <- 0
        for(k in n_source_states:1){# k = 3
          curr_col_source  <- curr_row[v_names_source_destiny[k+1]]
          curr_col_destiny <- m_keep_destiny_limited[j, v_names_source_destiny[k+1]] 
          curr_delta <- curr_col_destiny + curr_delta - curr_col_source
          v_curr_delta[k] <- curr_delta
          if(curr_delta < 0){
            m_keep_destiny_limited[j, "Keep"] <- FALSE
            break
          }
          if(curr_delta > curr_row[v_names_source_destiny[k]]){
            m_keep_destiny_limited[j, "Keep"] <- FALSE
            break
          }
        }
        
        if(m_keep_destiny_limited[j, "Keep"]==TRUE){
          binom_conv <- 1
          for(k in 1:n_source_states){ # k = 2
            # print(curr_row[k])
            # print(v_curr_delta[k])
            binom_conv <- binom_conv*dbinom(x = v_curr_delta[k], 
                                            size = curr_row[v_names_source[k]], 
                                            prob = p_flow)
            # print(binom_conv)
          }
          m_p_flow[name_curr_row, rownames(m_keep_destiny_limited)[j]] <- binom_conv
        }
      }
      
    }
    # t(t(m_keep_destiny[m_keep_destiny[, "Keep"]==1, ]) - c(curr_row, 0))
  }
  
  m_r_flow <- -log(1-m_p_flow)
  diag(m_r_flow[v_names_HH_flow_source, v_names_HH_flow_source]) <- 0
  m_r_prop_flow <- m_r_flow/rowSums(m_r_flow)
  
  m_hh_in[v_names_HH_flow_source, v_names_HH_flow_destiny] <- -1*(diag(m_hh_out[v_names_HH_flow_source, v_names_HH_flow_source]) * m_r_prop_flow)
  
  return(m_hh_in)
}