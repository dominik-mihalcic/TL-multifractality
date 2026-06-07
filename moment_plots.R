library(tidyverse)
source("helper.R")

# =====================================================
# log-moment scaling
# =====================================================

orders <- seq(1, 6, 1)

moment_values <- map(
  names(titles), 
  function(sheet){
    readxl::read_excel("TL_data.xlsx", sheet = sheet) |>
      mutate(groups = as.character(groups)) |>
      compute_moments(orders = orders) |>
      mutate(sheet = factor(sheet, levels = names(titles)))
    }
) |> bind_rows()

moment_values |> pivot_longer(
  cols = -c(groups, M1, sheet),
  values_to = "value",
  names_to = "moment"
) |> 
  mutate(moment = paste0("q = ", sub("M", "", moment))) |>
  ggplot(aes(x = log10(M1), y = log10(value), color = moment)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.5) +
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
