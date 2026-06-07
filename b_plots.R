library(tidyverse)
source("helper.R")

set.seed(123)

# =====================================================
# estimates of b(q) derived from the simulated cascade
# =====================================================

lambda <- exp(seq(log(1e-6), 0, length.out = 100))
orders <- seq(0.25, 6, 0.25)

N_r <- 100
N_c <- 100

b_values_1 <- imap_dfr(
  gamma2_vals, 
  function(gamma2, nm){
    casc <- lognorm_casc_2d(gamma2 = gamma2, ncol_out = N_c, nrow_out = N_r)
    moments <- matrix(NA, nrow = length(lambda), ncol = length(orders))
  
    for(i in seq_along(lambda)){
      l <- lambda[i]
      k_c <- max(1, floor(l * N_c))
      k_r <- max(1, floor(l * N_r))
      
      xi <- as.vector(sapply(
        1:(N_r - k_r + 1),
        function(j){
          sapply(
            1:(N_c - k_c + 1),
            function(k){
              sum(casc[j:(j+k_r-1), k:(k+k_c-1)])
            }
          )
        }
      )) 
      
      moments[i, ] <- sapply(orders, function(q) mean(xi^q))
    }
  
  xi_moments <- data.frame(groups = lambda, moments)
  colnames(xi_moments)[-1] <- paste0("M", orders)
  
  b_tbl <- compute_b(xi_moments)
  b_tbl |>
    mutate(
      gamma2_hat = coef(lm(b ~ 0 + offset(q) + I(q - q^2), data = b_tbl))[1],
      sheet = factor(nm, levels = names(titles))
    )
})

curve_data_1 <- b_values_1 |>
  distinct(sheet, gamma2_hat) |>
  expand_grid(q = seq(0.25, 6, length.out = 200)) |>
  mutate(b_fit = q * (1 + gamma2_hat) - q^2 * gamma2_hat)

plot_b(data = b_values_1, curve_data = curve_data_1)

# =====================================================
# estimates of b(q) derived from the actual data
# =====================================================

b_values_2 <- map(
  names(titles), 
  function(sheet){
    TL_data <- readxl::read_excel("TL_data.xlsx", sheet = sheet)
    moments <- compute_moments(TL_data, orders = orders)
    b_tbl <- compute_b(moments)
    b_tbl |>
      mutate(
        gamma2_hat = coef(lm(b ~ 0 + offset(q) + I(q - q^2), data = b_tbl))[1],
        sheet = factor(sheet, levels = names(titles))
      )
  }
) |> bind_rows()

curve_data_2 <- b_values_2 |>
  distinct(sheet, gamma2_hat) |>
  expand_grid(q = seq(0.25, 6, length.out = 200)) |>
  mutate(b_fit = q * (1 + gamma2_hat) - q^2 * gamma2_hat)

plot_b(data = b_values_2, curve_data = curve_data_2)
