# =====================================================
# Analysis of b(q) for all datasets
# =====================================================
library(tidyverse)
library(readxl)

source("TL_helpers.R")

file   <- "TL_data.xlsx"
orders <- seq(0.25, 6, 0.25)

set.seed(1)
B <- 1000     # number of bootstrap resamples   

# ---- analyse all datasets at once ----
out <- list()
all_results <- NULL
for(sheet in names(titles)){
  message("analysing sheet: ", sheet)
  out[[sheet]] <- analyse_sheet(sheet)
  all_results <- bind_rows(all_results, out[[sheet]]$results)
}

# ---- Tables 1 and 2: parameter estimates and sup-t measures ----
tables <- list()
params <- list(
  `log-normal`  = c("gamma2", "Delta"),
  `log-Poisson` = c("alpha", "beta")
)

for(nm in names(models)){
  tables[[nm]] <- all_results |>
    filter(model == nm, quantity %in% params[[nm]]) |>
    mutate(
      value  = sprintf("%.3g [%.3g, %.3g]", estimate, lo, hi),
      S      = round(S, 2),
      c_star = round(c_star, 2),
      CR     = round(CR, 2)
    ) |>
    select(sheet, quantity, value, S, c_star, CR) |>
    pivot_wider(names_from = quantity, values_from = value) |>
    relocate(S, c_star, CR, .after = last_col())
}

tables


# ---- Table 3: exponent of the empirical TL and b(2) ----
b_table <- list()

for(sheet in names(titles)){
  
  # exponent of the empirical TL from sample means and variances
  dat <- read_excel(file, sheet = sheet) |>
    mutate(groups = as.character(groups)) |>
    filter(!is.na(counts))
  mv <- dat |>
    group_by(groups) |>
    summarise(
      n = n(), M = mean(counts), V = var(counts),
      .groups = "drop"
    ) |>
    filter(n >= 2, M > 0, V > 0)
  b_mv <- coef(lm(log(V) ~ log(M), data = mv))[[2]]
  
  # estimated b(2) and b(2) of the fitted models
  b2_hat <- out[[sheet]]$b_tbl |> filter(q == 2) |> pull(b)
  b2_til <- out[[sheet]]$results |> filter(quantity == "b2")
  
  b_table[[sheet]] <- tibble(
    sheet     = titles[[sheet]],
    b_mv      = b_mv,
    b2_hat    = b2_hat,
    b2_til_ln = b2_til$estimate[b2_til$model == "log-normal"],
    b2_til_lp = b2_til$estimate[b2_til$model == "log-Poisson"]
  )
}

bind_rows(b_table) |>
  mutate(across(where(is.numeric), function(x) round(x, 2)))


# ---- Figure 4: empirical moments of orders 2 to 6 against the mean ----
moment_values <- NULL
for(sheet in names(titles)){
  m <- out[[sheet]]$moments |>
    filter(M1 > 0) |>
    select(groups, M1, M2, M3, M4, M5, M6) |>
    mutate(sheet = sheet)
  moment_values <- bind_rows(moment_values, m)
}

# ---- order of data sets in the plot ----
moment_values <- moment_values |>
  mutate(sheet = factor(sheet, levels = names(titles)))

# ---- figure ----
moment_values |>
  pivot_longer(
    cols = -c(groups, M1, sheet),
    values_to = "value",
    names_to = "moment"
  ) |>
  mutate(moment = paste0("q = ", sub("M", "", moment))) |>
  ggplot(aes(x = log10(M1), y = log10(value), color = moment)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1.5) +
  facet_wrap(~ sheet, ncol = 4, scales = "free",
             labeller = labeller(sheet = as_labeller(titles))) +
  labs(
    title = NULL,
    x = expression(log[10](M^{(j)})),
    y = expression(log[10](M[q]^{(j)}))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    strip.text = element_text(size = 16)
  )


# ---- Figure 5: b(q) - q with the fitted models and the sup-t band ----
q_grid    <- seq(0.25, 6, length.out = 200)
points_df <- NULL
curves_df <- NULL
bands_df  <- NULL

for(sheet in names(titles)){
  
  # estimates of b(q)
  b_tbl <- out[[sheet]]$b_tbl
  points_df <- bind_rows(points_df, b_tbl |> mutate(sheet = sheet))
  
  # simultaneous bootstrap band
  band <- sup_t_band(b_tbl, out[[sheet]]$b_boot)
  bands_df <- bind_rows(bands_df, band |> mutate(sheet = sheet))
  
  # fitted curves
  for(nm in names(out[[sheet]]$fits)){
    b_fit <- predict(out[[sheet]]$fits[[nm]], newdata = tibble(q = q_grid))
    curve <- tibble(q = q_grid, sheet = sheet, model = nm, b_fit = b_fit)
    curves_df <- bind_rows(curves_df, curve)
  }
}

# ---- order of datasets and models in the plots ----
points_df <- points_df |>
  mutate(sheet = factor(sheet, levels = names(titles)))

curves_df <- curves_df |>
  mutate(
    sheet = factor(sheet, levels = names(titles)),
    model = factor(model, levels = names(models))
  )

bands_df <- bands_df |>
  mutate(sheet = factor(sheet, levels = names(titles)))

# ---- figure ----
ggplot(points_df, aes(q, b - q)) +
  geom_ribbon(
    data = bands_df,
    aes(q, ymin = lo - q, ymax = hi - q),
    inherit.aes = FALSE, fill = "grey85"
  ) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(data = curves_df,
            aes(q, b_fit - q, colour = model),
            linewidth = 1.1) +
  geom_point(colour = "grey25", size = 2.2) +
  facet_wrap(~ sheet, ncol = 4, scales = "free",
             labeller = labeller(sheet = as_labeller(titles))) +
  labs(x = expression(q), y = expression(b(q) - q),
       colour = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90"),
    panel.border = element_rect(color = "black", fill = NA,
                                linewidth = 0.8),
    strip.text = element_text(size = 16),
    legend.position = "bottom",
    legend.text = element_text(size = 14)
  )
