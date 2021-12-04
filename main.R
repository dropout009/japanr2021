source("setting.R")

fs::dir_delete("figure")
fs::dir_create("figure")

set.seed(42)

get_nls_formula = function(max_fq) {
  map_chr(0:max_fq, ~ glue("(1 - theta) * (1 - pi)^{.} * s_{.}")) %>% 
    c("awareness_prop ~ 1", .) %>% 
    glue_collapse(sep = " - ") %>% 
    as.formula()
}


simulate <- function(
  n_samples = 1000,
  n_times = 20,
  max_fq = 4,
  reach_prob_fn,
  awareness_prob_fn
) {
  
  df <- tibble(
    i = 1:n_samples,
    reach_prob = reach_prob_fn(n_samples)
  ) %>%
    crossing(t = 1:n_times) %>%
    mutate(is_reach = as.double(rbernoulli(n_samples * n_times, p = reach_prob))) %>%
    group_by(i) %>%
    mutate(
      fq = cumsum(is_reach),
      fq = if_else(fq >= max_fq, max_fq, fq)
    ) %>%
    ungroup() %>%
    mutate(
      awareness_prob = awareness_prob_fn(fq),
      is_aware = rbernoulli(n_samples * n_times, p = awareness_prob)
    )

  df_awareness <- df %>%
    group_by(t) %>%
    summarise(awareness_prop = mean(is_aware))

  df_share <- df %>%
    group_by(t, fq) %>%
    summarise(fq_prop = n() / n_samples) %>%
    pivot_wider(
      id_cols = t,
      names_from = fq,
      values_from = fq_prop,
      values_fill = 0,
      names_prefix = "s_"
    )
  
  df_aggregated <- df_awareness %>%
    left_join(df_share)
  
  # Single Source
  awareness_single_source <- df %>%
    group_by(fq) %>%
    summarise(mean(is_aware)) %>%
    pull()
  
  
  
  # Nonparametric estimation
  S <- df_aggregated %>%
    select(starts_with("s")) %>%
    as.matrix()
  
  y <- df_aggregated %>% pull(awareness_prop)
  
  bounds <- function(x) {
    d <- numeric(length(x) - 1)
    for (i in 2:length(x)) {
      d[i - 1] <- x[i] - x[i - 1]
    }
    return(d)
  }
  
  solution <- solnp(
    pars = seq(0.1, 0.9, length.out = max_fq + 1),
    fun = function(x) sum((y - S %*% x)^2),
    ineqfun = bounds,
    ineqLB = rep(0, max_fq),
    ineqUB = rep(1, max_fq),
    LB = rep(0, max_fq + 1),
    UB = rep(1, max_fq + 1),
    control = list(tol = 1e-10)
  )
  
  awareness_nonpara <- solution$pars
  
  # NLS 
  model_nls <- nls(
    formula = get_nls_formula(max_fq),
    start = c(theta = 0.1, pi = 0.1),
    lower = c(0, 0),
    upper = c(1, 1),
    algorithm = "port",
    control = list(maxiter = 100, tol = 1e-10),
    data = df_aggregated
  )
  
  c <- coef(model_nls)
  awareness_nls <- 1 - (1 - c["theta"]) * (1 - c["pi"])^(0:max_fq)
  
  # つっくけて完成
  levels <- c("Theoretical", "Single Source", "Nonpara", "NLS")
  
  df_result = tibble(
    fq = 0:max_fq,
    Theoretical = awareness_prob_fn(fq),
    `Single Source` = awareness_single_source,
    Nonpara = awareness_nonpara,
    NLS = awareness_nls
  ) %>%
    pivot_longer(!fq) %>%
    mutate(name = factor(name, levels = levels))
  
  
  return(list(
    y = y,
    S = S,
    result = df_result,
    raw = df,
    aggregated = df_aggregated
  ))
}


draw_lines  = function(df) {
  df %>% 
    ggplot(aes(fq, value, color = name, fill = name)) +
    geom_line() +
    geom_point(shape = 21, color = "white", size = 3) +
    scale_y_continuous(
      breaks = breaks_width(0.2),
      labels = label_percent(1),
      limits = c(0, 1)
    ) +
    scale_color_manual(
      name = "", 
      values = cols[c(8, 2, 4, 6)]
    ) +
    scale_fill_manual(
      name = "", 
      values = cols[c(8, 2, 4, 6)]
    ) +
    labs(
      title = "フリークエンシーとCM認知度の関係",
      x = "フリークエンシー", 
      y = "CM認知率"
    ) +
    theme_line()
}


result_list = simulate(
  n_samples = 10000,
  n_times = 10,
  max_fq = 5,
  reach_prob_fn = partial(runif, min = 0.4, max = 0.6),
  awareness_prob_fn = function(f) 1 - (1 - 0.2) * (1 - 0.3)^f
)

draw_lines(result_list[["result"]])
save_plot(plot = last_plot(), fname = "figure/result01.png")

result_list = simulate(
  n_samples = 10000,
  n_times = 20,
  max_fq = 5,
  reach_prob_fn = partial(runif, min = 0.4, max = 0.6),
  awareness_prob_fn = function(f) c(0.1, 0.8, 0.85, 0.9, 0.95, 1)[f + 1]
)

draw_lines(result_list[["result"]])
save_plot(plot = last_plot(), fname = "figure/result02.png")



# ベイズ ---------------------------------------------------------------------

fit_bayes = function(data, max_fq) {
  brm(
    bf(
      get_nls_formula(max_fq),
      theta ~ 1,
      pi ~ 1,
      nl = TRUE
    ),
    prior = c(
      prior(beta(1, 1), lb = 0, ub = 1, nlpar = "theta"),
      prior(beta(1, 1), lb = 0, ub = 1, nlpar = "pi")
    ),
    data = data,
    iter = 8000,
    cores = 4, 
    seed = 42
  )
}

draw_parameter_density = function(model) {
  model %>%
    gather_draws(b_theta_Intercept, b_pi_Intercept) %>% 
    mutate(
      .variable = .variable %>% str_remove("b_") %>% str_remove("_Intercept")
    ) %>% 
    ggplot(aes(.value)) +
    stat_halfeye() + 
    facet_wrap(~.variable, scales = "free_x", labeller = label_parsed) +
    labs(x = NULL, y = NULL, title = "パラメータの事後分布") + 
    theme_line()  
}

draw_frecuency_awareness_curve = function(model) {
  model %>%
    spread_draws(b_theta_Intercept, b_pi_Intercept) %>% 
    rename(theta = b_theta_Intercept, pi = b_pi_Intercept) %>%
    crossing(fq = 0:5L) %>%
    mutate(awareness_prop = 1 - (1 - theta) * (1 - pi)^fq) %>% 
    ggplot(aes(fq, awareness_prop)) +
    # geom_line(alpha = 0.1, color = cols[2]) +
    stat_lineribbon(
      .width = c(.95), 
      color = brewer.pal(5, "Blues")[[5]]
    ) + 
    scale_y_continuous(
      breaks = breaks_width(0.2),
      labels = label_percent(1),
      limits = c(0, 1)
    ) +
    scale_fill_brewer() + 
    labs(
      title = "フリークエンシーとCM認知度の関係",
      x = NULL, 
      y = NULL,
    ) +
    theme_line() +
    theme(legend.position = "none")
}


result_list1 = simulate(
  n_samples = 10000,
  n_times = 5,
  max_fq = 5,
  reach_prob_fn = partial(runif, min = 0.4, max = 0.6),
  awareness_prob_fn = function(f) 1 - (1 - 0.2) * (1 - 0.3)^f
)

model1 = fit_bayes(
  data = result_list1[["aggregated"]], 
  max_fq = 5
)

result_list2 = simulate(
  n_samples = 10000,
  n_times = 20,
  max_fq = 5,
  reach_prob_fn = partial(runif, min = 0.4, max = 0.6),
  awareness_prob_fn = function(f) 1 - (1 - 0.2) * (1 - 0.3)^f
)

model2 = fit_bayes(
  data = result_list2[["aggregated"]], 
  max_fq = 5
)



draw_parameter_density(model1)
save_plot(plot = last_plot(), fname = "figure/param01.png", width = 5)
draw_frecuency_awareness_curve(model1)
save_plot(plot = last_plot(), fname = "figure/curve01.png", width = 5)



draw_parameter_density(model2)
save_plot(plot = last_plot(), fname = "figure/param02.png", width = 5)
draw_frecuency_awareness_curve(model2)
save_plot(plot = last_plot(), fname = "figure/curve02.png", width = 5)


tibble(fq = 0:5, value = c(0.1, 0.5, 0.8, 0.9, 0.925, 0.95)) %>% 
  ggplot(aes(fq, value)) +
  geom_line(color = cols[2]) +
  geom_point(shape = 21, color = "white", fill = cols[2], size = 3) +
  scale_y_continuous(
    breaks = breaks_width(0.2),
    labels = label_percent(1),
    limits = c(0, 1)
  ) +
  labs(
    x = "フリークエンシー", 
    y = "CM認知率"
  ) +
  theme_line()

save_plot(plot = last_plot(), fname = "figure/image.png")
