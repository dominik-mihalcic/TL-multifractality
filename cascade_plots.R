library(tidyverse)
library(plot3D)
source("helper.R")

set.seed(123)

# ==================================================
# 1D log-normal multiplicative cascade
# ==================================================

# ---- parameters ----
gamma2_values <- c(0.05, 0.1, 0.5, 1)
N <- 10000

# ---- grid ----
t <- seq(0, 1, length.out = N)

# ---- data generation ----
cascade_values_1d <- map_dfr(
  gamma2_values,
  function(gamma2){
    xi <- lognorm_casc_1d(gamma2 = gamma2, n_out = N)
    tibble(
      t = t, 
      xi = xi,
      gamma2 = gamma2,
      label = sprintf("gamma^2 == %.4f", gamma2)
    )
  }
)

# ---- plotting ----
ggplot(cascade_values_1d, aes(x = t, y = xi)) +
  geom_line(linewidth = 1, color = "#2C3E50") +
  facet_wrap(
    ~ label, 
    ncol = 2, 
    scales = "free_y", 
    labeller = label_parsed
  ) +
  labs(
    x = expression(t), 
    y = expression(xi(dt))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    strip.text = element_text(size = 16)
  )

# ==================================================
# 2D log-normal multiplicative cascade
# ==================================================

# ---- parameters ----
N_c <- 100
N_r <- 100

# ---- grid ----
x <- seq(0, 1, length.out = N_c)
y <- seq(0, 1, length.out = N_r)

# ---- data generation ----
cascade_values_2d <- lapply(
  gamma2_values, 
  function(gamma2){
    lognorm_casc_2d(gamma2 = gamma2, ncol_out = N_c, nrow_out = N_r)
  }  
)

# ---- plotting ----
plot_labels <- as.expression(
  lapply(
    gamma2_values,
    function(g) bquote(gamma^2 == .(g))
  )
)

par(mfrow = c(2, 2), mar = c(2, 2, 3, 2))

lapply(
  seq_along(cascade_values_2d),
  function(i){
    xi <- cascade_values_2d[[i]]
    persp3D(
      x = x, y = y, z = xi,
      col = col_palette(100),
      colkey = FALSE,
      theta = 30,
      phi = 20,
      expand = 0.6,
      xlab = "x", ylab = "y", zlab = "",
      ticktype = "detailed",
      border = NA,
      main = plot_labels[i],
      clab = NULL
    )
  }
) |> invisible()
