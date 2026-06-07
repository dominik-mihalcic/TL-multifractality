# =====================================================
# moment computation
# =====================================================
compute_moments <- function(data, orders){
  data |> 
    group_by(groups) |>
    summarise(
      across(
        .cols = counts,
        .fns = set_names(
          map(orders, ~ function(x) mean(abs(x)^.x)),
          paste0("M", orders)
        ),
        .names = "{.fn}"
      ),
      .groups = "drop"
    )
}

# =====================================================
# b(q) estimation
# =====================================================
compute_b <- function(data){
  data |> 
    pivot_longer(
      cols = -c(groups, M1),
      values_to = "value",
      names_to = "moment"
    ) |> 
    group_by(moment) |> 
    summarise(
      b = coef(lm(log10(value) ~ log10(M1)))[2],
      .groups = "drop"
    ) |> 
    mutate(q = as.numeric(substring(moment, 2))) |> 
    bind_rows(tibble(moment = "M1", b = 1, q = 1)) 
}

# =====================================================
# 1D cascade
# =====================================================

lognorm_casc_1d <- function(gamma2, init_val = 1, n_out, tmax = 1){
  
  rmin <- tmax/n_out # minimum resolution 
  t_pos <- seq(0, tmax, length.out = n_out) # points on t-axis
  
  # compound Poisson random measure
  lambda <- (1/rmin - 1)*(tmax + 1) 
  N <- rpois(1, lambda) 
  
  # N random centers
  ti <- runif(N, -0.5, tmax + 0.5) 
  
  # N random radii
  umax <- 1/rmin - 1
  ui <- runif(N, 0, umax)
  ri <- (1 + ui)^(-1)                         
  
  # lognormal cascade
  wi <- rlnorm(N, -gamma2/2, sqrt(gamma2)) 
  casc <- rep(init_val, n_out)
  
  for(i in seq_along(wi)) {
    ind <- which(abs(t_pos - ti[i]) <= ri[i]/2)
    casc[ind] <- casc[ind]*wi[i]
  }
  
  return(casc)
}

# =====================================================
# 2D cascade
# =====================================================

lognorm_casc_2d <- function(gamma2, ncol_out, nrow_out, init_val = 1){
  
  # size of rectangle 
  xmax <- 1     
  ymax <- 1
  
  # points in the xy plane
  x_pos <- seq(0, xmax, length.out = ncol_out)
  y_pos <- seq(0, ymax, length.out = nrow_out)
  X <- matrix(
    rep(x_pos, each = nrow_out), 
    nrow = nrow_out,
    ncol = ncol_out
  )
  Y <- matrix(
    rep(y_pos, times = ncol_out),
    nrow = nrow_out,
    ncol = ncol_out
  )
  
  # minimum resolutions
  rmin_x <- xmax / ncol_out
  rmin_y <- ymax / nrow_out
  
  # compound Poisson random measure
  lambda <- 2/pi*((1/rmin_x)*(1/rmin_y) - 1)*(xmax + 1)*(ymax + 1)
  N <- rpois(1, lambda)
  
  # N random centers 
  xi <- runif(N, -0.5, xmax + 0.5)
  yi <- runif(N, -0.5, ymax + 0.5)
  
  # N random radii
  umax <- 1 / (rmin_x * rmin_y) - 1
  ui <- runif(N, 0, umax) 
  ri <- (1 / (rmin_x * rmin_y) - ui)^(-0.5)
  
  # lognormal cascade 
  wi <- rlnorm(N, -gamma2/2, sqrt(gamma2))
  casc <- matrix(init_val, nrow = nrow_out, ncol = ncol_out)
  
  for(i in seq_along(wi)){
    ind <- which((X - xi[i])^2 + (Y - yi[i])^2 <= (ri[i]^2)/4)
    casc[ind] <- casc[ind]*wi[i]
  }
  
  return(casc)
}

# =====================================================
# visual settings
# =====================================================

titles <- c(
  "Beall" = "Colorado potato beetle",
  "Kleczkowski" = "Tobacco mosaic virus",
  "Milne" = "Sheep tick",
  "Trigal and Ruete" = "Invasive flagellate",
  "Cohen I" = "Hokkaido vole - Group I", 
  "Cohen IV" = "Hokkaido vole - Group IV",
  "Cohen V" = "Hokkaido vole - Group V",
  "Cohen VII" = "Hokkaido vole - Group VII"
)

gamma2_vals <- c(
  "Beall" = 0.1558,
  "Kleczkowski" = 0.1116,
  "Milne" = 0.2481,
  "Trigal and Ruete" = 0.1208,
  "Cohen I" = 0.1953, 
  "Cohen IV" = 0.1778,
  "Cohen V" = 0.1357,
  "Cohen VII" = 0.2554
)

col_palette <- colorRampPalette(
  c("#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00")
)

plot_b <- function(data, curve_data) {
  
  ggplot(data, aes(x = q, y = b)) +
    geom_point(color = "#F8766D", size = 3) +
    geom_line(
      data = curve_data,
      aes(y = b_fit),
      color = "#00BFC4",
      linewidth = 1.5
    ) +
    facet_wrap(
      ~ sheet,
      ncol = 4,
      scales = "free",
      labeller = labeller(sheet = as_labeller(titles))
    ) +
    labs(x = expression(q), y = expression(b(q))) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16, face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      strip.text = element_text(size = 16)
    )
}