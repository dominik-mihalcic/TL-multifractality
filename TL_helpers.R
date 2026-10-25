# ---- sheet titles ----
titles <- c(
  "Beall" = "Colorado potato beetle",
  "Kleczkowski" = "Tobacco mosaic virus",
  "Milne" = "Sheep tick",
  "Trigal and Ruete" = "Invasive flagellate",
  "Cohen I" = "Hokkaido vole - I",
  "Cohen IV" = "Hokkaido vole - IV",
  "Cohen V" = "Hokkaido vole - V",
  "Cohen VII" = "Hokkaido vole - VII"
)

# ---- empirical moments of order q for every group ----
compute_moments <- function(data, orders){
  data <- data[!is.na(data$counts), ]
  x <- split(abs(data$counts), data$groups)
  out <- list(groups = names(x))
  for(q in orders){
    out[[paste0("M", q)]] <- unname(sapply(x, function(v) mean(v^q)))
  }
  as_tibble(out)
}

# ---- compute b(q) ----
compute_b <- function(moments){
  m <- moments[moments$M1 > 0, ]
  cols <- setdiff(names(m), c("groups", "M1"))
  slope <- lm.fit(
    cbind(1, log10(m$M1)),
    log10(as.matrix(m[cols]))
  )$coefficients[2, ]
  
  tibble(
    moment = cols, b = unname(slope),
    q = as.numeric(substring(cols, 2))
  ) |>
    bind_rows(tibble(moment = "M1", b = 1, q = 1)) |>
    arrange(q)
}

# ---- specifications of the models ----
models <- list(
  `log-normal` = list(
    formula = b ~ q - gamma2/2*(q^2-q) / 2,
    starts  = list(list(gamma2 = 0.1)),
    lower   = c(1e-8), upper = c(10)
  ),
  `log-Poisson` = list(
    formula = b ~ q - alpha*((beta^q-1)-q*(beta-1)) / 2,
    starts  = list(list(alpha = 0.4, beta = 0.5)),
    lower = c(1e-6, 1e-3), upper = c(1e4, 0.999)
  )
)

# ---- fit ----
fit_model <- function(spec, b_tbl, starts){
  for(s in starts){
    f <- try(nls(
      spec$formula, data = b_tbl, start = s,
      algorithm = "port",
      lower = spec$lower, upper = spec$upper,
      control = nls.control(maxiter = 1000)),
      silent = TRUE)
    if(class(f) != "try-error") return(f)
  }
  NULL
}

# ---- parameters, fitted b(2) and curvature difference ----
fit_summary <- function(fit, b_tbl){
  resid <- b_tbl$b - predict(fit)
  Delta <- coef(lm(resid ~ q + I(q^2), data = b_tbl))[[3]]
  as_tibble(as.list(coef(fit))) |>
    mutate(
      b2 = predict(fit, newdata = data.frame(q = 2)),
      Delta = Delta
    )
}

# ---- bootstrap replicates of b(q) ----
boot_b <- function(dat, B, orders) {
  b_all <- vector("list", B)
  for(i in 1:B){
    dat_i <- dat |>
      group_by(groups) |>
      slice_sample(prop = 1, replace = TRUE) |>
      ungroup()
    b_all[[i]] <- compute_b(compute_moments(dat_i, orders)) |> mutate(rep = i)
  }
  bind_rows(b_all)
}

# ---- sup-t band for b(q) ----
sup_t_band <- function(b_tbl, b_boot){
  se_tbl <- b_boot |>
    filter(q != 1) |>
    group_by(q) |>
    summarise(se = sd(b))
  
  t_boot <- b_boot |>
    filter(q != 1) |>
    inner_join(b_tbl |> select(q, b_hat = b), by = "q") |>
    inner_join(se_tbl, by = "q") |>
    group_by(rep) |>
    summarise(t = max(abs(b - b_hat) / se))
  c_star <- unname(quantile(t_boot$t, 0.95))
  
  b_tbl |>
    inner_join(se_tbl, by = "q") |>
    mutate(lo = b - c_star*se, hi = b + c_star*se, c_star = c_star)
}

# ---- measures of a fitted model ----
boot_measures <- function(fit, b_tbl, b_boot){
  band  <- sup_t_band(b_tbl, b_boot)
  model <- predict(fit, newdata = data.frame(q = band$q))
  tibble(
    S = max(abs(band$b - model) / band$se),
    c_star = band$c_star[1],
    CR = mean(abs(band$b - model) / band$se <= band$c_star)
  )
}

# ---- analysis of a single dataset ----
analyse_sheet <- function(sheet){
  
  # data
  dat <- read_excel(file, sheet = sheet) |>
    mutate(groups = as.character(groups)) |>
    filter(!is.na(counts))
  moments <- compute_moments(dat, orders)
  b_tbl   <- compute_b(moments)
  
  # fits and estimates on the original data
  fits <- list()
  est  <- list()
  for(nm in names(models)){
    fits[[nm]] <- fit_model(models[[nm]], b_tbl, models[[nm]]$starts)
    if(is.null(fits[[nm]])) next
    est[[nm]] <- fit_summary(fits[[nm]], b_tbl)
  }
  
  # bootstrap
  b_all   <- boot_b(dat, B, orders)
  est_all <- list()
  for(i in 1:B){
    b_i <- b_all |> filter(rep == i)
    for(nm in names(fits)){
      fit_i <- fit_model(models[[nm]], b_i, models[[nm]]$starts)
      if(is.null(fit_i)) next
      est_all[[length(est_all) + 1]] <- fit_summary(fit_i, b_i) |>
        mutate(model = nm)
    }
  }
  est_all <- bind_rows(est_all)
  
  # 95% percentile bootstrap intervals
  intervals <- est_all |>
    pivot_longer(-model, names_to = "quantity", values_drop_na = TRUE) |>
    group_by(model, quantity) |>
    summarise(
      lo = quantile(value, 0.025),
      hi = quantile(value, 0.975),
      .groups = "drop"
    )
  
  # measures of each model
  measures <- list()
  for(nm in names(fits)){
    measures[[nm]] <- boot_measures(fits[[nm]], b_tbl, b_all)
  }
  
  results <- bind_rows(est, .id = "model") |>
    pivot_longer(
      -model, names_to = "quantity",
      values_to = "estimate",
      values_drop_na = TRUE
    ) |>
    left_join(intervals, by = c("model", "quantity")) |>
    left_join(bind_rows(measures, .id = "model"), by = "model") |>
    mutate(sheet = sheet, .before = 1)
  
  list(
    results = results,
    moments = moments,
    b_tbl   = b_tbl,
    b_boot  = b_all,
    fits    = fits
  )
}