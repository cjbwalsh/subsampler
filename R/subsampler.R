#' Simulate one or more sub-samples from a 100-cell subsampler
#'
#' @param ss The number cells that make up the subsample (equal to the %
#'     subsample size)
#' @param N The total number of specimens in the sample
#' @param n_iter the number of repeated subsamples (with replacement) to be
#'     simulated
#' @param set_random_seed If NULL, R's random seed generator is used.
#'     Alternatively, a random seed can be specified to achieve a reproducible
#'     result. The function resets to R's generator after running.
#' @return a vector of length n_iter of counts in the subsample
#' @details the function simulates the process of evenly spreading the sample
#'     across the 100 cells (allocating each of the N specimens to a cell drawn
#'     randomly from a uniform distribution between 1 and 100), and then
#'     randomly selecting ss cells to draw the subsample from.
#' @examples
#'      box_ss(10,200,10)

box_ss <- function(ss = 10, N = 200, n_iter = 1e5, set_random_seed = NULL){
  set.seed(set_random_seed)
  for(i in 1:n_iter){
    # Distribute the N randomly across the 100 cells of the subsampler
    N_by_cell <- sample(1:100,size = N, replace = TRUE)
    # Number of specimens in each cell with >0 specimens
    mixi <- aggregate(N_by_cell, by = list(cell = N_by_cell), FUN = length)
    # sample ss of the 100 cells in the subsampler
    ssi <- sample(1:100, size = ss, replace = FALSE)
    # sum the total number of specimens in the sampled cells
    if(i == 1){
      count_ss <- sum(mixi$x[mixi$cell %in% ssi])
    }else{
      count_ss <- c(count_ss,sum(mixi$x[mixi$cell %in% ssi]))
    }
  }
  set.seed(NULL)
  count_ss
}

#' Calculate Gamma distribution parameters for a sub-sample from a 100-cell subsampler
#'
#' @param count_in_subsample The number of specimens collected in a subsample
#' @param ss_proportion The size of the subsample as a proportion of the whole sample
#' @return a list of 4 values
#'    \code{d_var} variance of the gamma distribution
#'    \code{mo} mode/mean, calculated as count_in_subsample/ss_proportion
#'    \code{mo_adj} adjusted mode/mean used for determining rate and shape
#'    \code{rate} the rate (beta) parameter of the distribution
#'    \code{shape} the shape (alpha) parameter of the distribution
#' @details the mo estimate is based on the general relationship between mo and
#'    the combination of count_in_subsample and ss_proportion, derived by
#'    simulation. Shape and rate are estimated by moment method from d_var and
#'    mo_adj. See manuscript for more details.
#' @examples
#'      dist_0_10 <- gamma_params_from_s_c(0.1,0)
#'      plot(0:40, dgamma(0:40, dist_0_10$shape, dist_0_10$rate),
#'           type = 'b',xlab = "T", ylab = "Density", las = 1)
#'      dist_10_25 <- gamma_params_from_s_c(0.25,10)
#'      plot(10:100, dgamma(10:100, dist_10_25$shape, dist_10_25$rate),
#'           type = 'l',xlab = "T", ylab = "Density", las = 1)

gamma_params_from_s_c <- function(ss_proportion,count_in_subsample){
  if(ss_proportion < 0.1)  warning('Function not validated for values of s < 0.1.')
  if(ss_proportion < 1)
  d_var <- T_var(ss_proportion, count_in_subsample)
  mo <- count_in_subsample/ss_proportion
  if(count_in_subsample > 0 | ss_proportion == 1){
    if(ss_proportion == 1){
      d_var <- 0
      mo <- count_in_subsample
      rate <- 100
      shape <- ifelse(count_in_subsample == 0, 1, rate * count_in_subsample)
      # This ensures that 100% of probability density is within +/- 1% of mode
       }else{
        rate <- (sqrt(4*d_var + mo^2) + mo)/(2*d_var)
        shape <- d_var * rate^2}
     }else{
    #if count = 0 and s < 1, shape = 1, and pgamma(1,1,rate) == s
    shape <- 1
    rate <- rate_for_c_zero(ss_proportion)
    d_var <- NA; mo <- NA;
  }
  list(d_var = d_var, mo = mo,
       rate = rate, shape = shape)
}


#' Calculate rate parameter of Gamma distribution for a sub-sample of proportion s with count = 0
#'
#' @param s The subsample proportion
#' @return The rate parameter of the gamma distribution
#' @details Internal function for ss_count_gamma_params()
rate_for_c_zero <- function(s){
  if(s == 1){
    # boundary condition of full sample -
    # see ss_count_gamma_params
    rate <- 100
  }else{
    # Search for rate value where pgamma(1,1,rate) = s
    # Narrower searches after exploratory full search
    upper_rate_adj <- ifelse(s <= 0.4, 1.3,
                             ifelse(s <= 0.5, 1.4,
                                    ifelse(s <= 0.7,1.75,
                                           ifelse(s <= 0.8,2.1,
                                                  ifelse(s <= 0.9, 2.6, 5)))))
    #First coarse search
    rate_seq <- seq(s*1.005,s*upper_rate_adj, length = 1e2) #rate for s=0.01 is 0.01005
    search_range <- s*upper_rate_adj - s*1.005
    p_seq <- rep(NA,1e5)
    for(j in 1:length(rate_seq)){
      p_seq[j] <- pgamma(1,1,rate_seq[j])
    }
    #Refined search
    best_coarse <- rate_seq[which.min(abs(p_seq - s))]
    rate_seq <- seq(best_coarse - search_range/200,
                    best_coarse + search_range/200, length = 1e5)
    p_seq <- rep(NA,1e5)
    for(j in 1:length(rate_seq)){
      p_seq[j] <- pgamma(1,1,rate_seq[j])
    }

    rate <- rate_seq[which.min(abs(p_seq - s))]
  }
  return(rate)
}

#' A generalised approximation of variance of the distribution of T|(s, c)
#' @param c The number of specimens counted in a subsample
#' @param s The subsample proportion
#' @return The variance of the gamma distribution
#' @details Internal function for ss_count_gamma_params(). Coefficients used to
#'   calculate var_slope and var_intercept were determined as described by Walsh
#'   (2022) Measurement-error models improve prediction from subsampled data
#'   (See code chunk c_var_params, and Fig S1-2 in Appendix S1.

T_var <- function(ss_proportion,count_in_subsample){
  if(ss_proportion < 0.1)  warning('T_var function not validated for values of s < 0.1')
  if(ss_proportion > 0.97)  warning('T_var function not validated for values of s > 0.97')
    out_var <- (count_in_subsample + 1)*(0.3139303 + 0.1246125*(boot::logit(ss_proportion) + 5))^-11
  out_var
}

