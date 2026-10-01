# Dynamic Hedge-Ratio Pairs Trading: NVDA vs AVGO
# Rolling 60-day OLS beta | Rolling Z-score | Entry / Target / Hard Stop

required_pkgs <- c("quantmod", "PerformanceAnalytics", "xts", "zoo",
                   "tseries", "knitr", "tidyverse")

install_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(missing_pkgs) > 0) install.packages(missing_pkgs)
}
install_if_missing(required_pkgs)

suppressPackageStartupMessages({
  library(quantmod)
  library(PerformanceAnalytics)
  library(tidyverse)
})

# 1. Configuration
cfg <- list(
  tickers              = c("NVDA", "AVGO"),
  years_total          = 5,
  years_train          = 3,
  lookback             = 60,
  entry_z              = 2.0,
  exit_z               = 0.0,
  stop_z               = 3.5,
  cost_bps             = 5,      # per side, on gross notional
  rf_annual            = 0.0,
  trading_days         = 252,
  eg_crit_5pct         = -3.34,  # approx. Engle-Granger 5% critical value (2 vars, constant)
  require_cointegration = FALSE  # TRUE = abort if in-sample test fails
)

# 2. Data
get_prices <- function(cfg) {
  end_date   <- Sys.Date()
  start_date <- end_date - round(365.25 * cfg$years_total)
  
  px_list <- lapply(cfg$tickers, function(tk) {
    raw <- tryCatch(
      quantmod::getSymbols(tk, src = "yahoo", from = start_date, to = end_date,
                           auto.assign = FALSE, warnings = FALSE),
      error = function(e) stop("Download failed for ", tk, ": ", conditionMessage(e))
    )
    out <- quantmod::Ad(raw)
    colnames(out) <- tk
    out
  })
  
  px <- na.omit(do.call(merge, px_list))
  
  tibble(
    date = zoo::index(px),
    NVDA = as.numeric(px[, cfg$tickers[1]]),
    AVGO = as.numeric(px[, cfg$tickers[2]])
  )
}

split_sample <- function(px_tbl, years_train) {
  split_date <- min(px_tbl$date) + round(365.25 * years_train)
  list(
    train      = dplyr::filter(px_tbl, date <  split_date),
    test       = dplyr::filter(px_tbl, date >= split_date),
    split_date = split_date
  )
}

# 3. In-sample cointegration (Engle-Granger two-step with ADF)
test_cointegration <- function(train, cfg) {
  fit <- lm(NVDA ~ AVGO, data = train)
  adf <- tseries::adf.test(residuals(fit), alternative = "stationary")
  
  tibble(
    static_beta        = unname(coef(fit)["AVGO"]),
    static_alpha       = unname(coef(fit)["(Intercept)"]),
    adf_statistic      = unname(adf$statistic),
    adf_p_value        = adf$p.value,
    eg_crit_5pct       = cfg$eg_crit_5pct,
    reject_at_eg_5pct  = unname(adf$statistic) < cfg$eg_crit_5pct
  )
}

# 4. Rolling OLS hedge ratio, spread, Z-score
rolling_ols <- function(y, x, width) {
  n     <- length(y)
  alpha <- rep(NA_real_, n)
  beta  <- rep(NA_real_, n)
  
  if (n >= width) {
    for (i in width:n) {
      idx    <- (i - width + 1):i
      fit    <- stats::lm.fit(cbind(1, x[idx]), y[idx])
      alpha[i] <- fit$coefficients[1]
      beta[i]  <- fit$coefficients[2]
    }
  }
  tibble(alpha = alpha, beta = beta)
}

build_signal_frame <- function(px_tbl, split_date, cfg) {
  w <- cfg$lookback
  
  px_tbl %>%
    dplyr::bind_cols(rolling_ols(px_tbl$NVDA, px_tbl$AVGO, w)) %>%
    dplyr::mutate(
      spread    = NVDA - beta * AVGO,
      spread_mu = as.numeric(zoo::rollapplyr(spread, width = w, FUN = mean, fill = NA)),
      spread_sd = as.numeric(zoo::rollapplyr(spread, width = w, FUN = sd,   fill = NA)),
      z         = (spread - spread_mu) / spread_sd,
      is_oos    = date >= split_date
    )
}

# 5. Signal state machine
# +1 = long spread (long NVDA / short AVGO), -1 = short spread, 0 = flat.
# After a stop-loss, re-entry is blocked until |z| has fallen back inside the entry band, which prevents same-regime whipsaw re-entries.
generate_positions <- function(z, beta, is_oos, entry, exit, stop) {
  n      <- length(z)
  pos    <- numeric(n)
  reason <- character(n)
  cur    <- 0
  armed  <- TRUE
  
  for (i in seq_len(n)) {
    zi       <- z[i]
    tradable <- isTRUE(is_oos[i]) && !is.na(zi) && !is.na(beta[i]) && beta[i] > 0
    
    if (!tradable) {
      if (cur != 0) {
        reason[i] <- "data_invalid"
        cur <- 0
      }
      pos[i] <- cur
      next
    }
    
    if (cur == 0) {
      if (!armed && abs(zi) < entry) armed <- TRUE
      if (armed) {
        if (zi < -entry && zi > -stop) {
          cur <- 1
        } else if (zi > entry && zi < stop) {
          cur <- -1
        }
      }
    } else if (cur == 1) {
      if (zi <= -stop) {
        reason[i] <- "stop_loss";     cur <- 0; armed <- FALSE
      } else if (zi >= exit) {
        reason[i] <- "profit_target"; cur <- 0
      }
    } else {
      if (zi >= stop) {
        reason[i] <- "stop_loss";     cur <- 0; armed <- FALSE
      } else if (zi <= exit) {
        reason[i] <- "profit_target"; cur <- 0
      }
    }
    pos[i] <- cur
  }
  
  if (n > 0 && pos[n] != 0) {
    pos[n]    <- 0
    reason[n] <- "end_of_sample"
  }
  list(pos = pos, reason = reason)
}

# 6. Backtest engine
run_backtest <- function(sig, cfg) {
  sp <- generate_positions(sig$z, sig$beta, sig$is_oos,
                           cfg$entry_z, cfg$exit_z, cfg$stop_z)
  
  full <- sig %>%
    dplyr::mutate(
      signal_pos  = sp$pos,
      exit_reason = sp$reason,
      pos_lag     = dplyr::lag(signal_pos, default = 0),
      beta_lag    = dplyr::lag(beta),
      nvda_lag    = dplyr::lag(NVDA),
      avgo_lag    = dplyr::lag(AVGO),
      r_nvda      = NVDA / nvda_lag - 1,
      r_avgo      = AVGO / avgo_lag - 1,
      w_nvda      = nvda_lag / (nvda_lag + beta_lag * avgo_lag),
      w_avgo      = 1 - w_nvda,
      gross_ret   = dplyr::if_else(pos_lag == 0, 0,
                                   pos_lag * (w_nvda * r_nvda - w_avgo * r_avgo)),
      cost        = cfg$cost_bps / 1e4 * abs(signal_pos - dplyr::lag(signal_pos, default = 0)),
      strat_ret   = dplyr::coalesce(gross_ret - cost, 0)
    )
  
  bt <- full %>%
    dplyr::filter(is_oos) %>%
    dplyr::mutate(
      strat_equity = cumprod(1 + strat_ret),
      bench_equity = 0.5 * NVDA / dplyr::first(NVDA) + 0.5 * AVGO / dplyr::first(AVGO),
      bench_ret    = bench_equity / dplyr::lag(bench_equity, default = 1) - 1
    )
  
  trades <- bt %>%
    dplyr::mutate(
      entry_event = signal_pos != 0 & dplyr::lag(signal_pos, default = 0) == 0,
      trade_id    = cumsum(entry_event)
    ) %>%
    dplyr::filter(signal_pos != 0 | pos_lag != 0) %>%
    dplyr::group_by(trade_id) %>%
    dplyr::summarise(
      direction    = dplyr::if_else(dplyr::first(signal_pos[signal_pos != 0]) == 1,
                                    "Long spread", "Short spread"),
      entry_date   = min(date),
      exit_date    = max(date),
      holding_days = sum(pos_lag != 0),
      trade_return = prod(1 + strat_ret) - 1,
      exit_reason  = dplyr::last(exit_reason[exit_reason != ""]),
      .groups      = "drop"
    )
  
  list(bt = bt, trades = trades)
}

# 7. Performance analytics
compute_metrics <- function(bt, trades, cfg) {
  strat <- xts::xts(bt$strat_ret, order.by = bt$date)
  bench <- xts::xts(bt$bench_ret, order.by = bt$date)
  rf_d  <- cfg$rf_annual / cfg$trading_days
  
  tot_ret <- function(x) as.numeric(PerformanceAnalytics::Return.cumulative(x))
  sharpe  <- function(x) as.numeric(PerformanceAnalytics::SharpeRatio.annualized(
    x, Rf = rf_d, scale = cfg$trading_days))
  mdd     <- function(x) as.numeric(PerformanceAnalytics::maxDrawdown(x))
  pct     <- function(x) sprintf("%.2f%%", 100 * x)
  
  n_trades  <- nrow(trades)
  win_rate  <- if (n_trades > 0) mean(trades$trade_return > 0) else NA_real_
  n_stop    <- sum(trades$exit_reason == "stop_loss", na.rm = TRUE)
  n_stop_lo <- sum(trades$exit_reason == "stop_loss" & trades$direction == "Long spread",  na.rm = TRUE)
  n_stop_hi <- sum(trades$exit_reason == "stop_loss" & trades$direction == "Short spread", na.rm = TRUE)
  n_target  <- sum(trades$exit_reason == "profit_target", na.rm = TRUE)
  n_open    <- sum(trades$exit_reason == "end_of_sample", na.rm = TRUE)
  n_forced  <- sum(trades$exit_reason == "data_invalid", na.rm = TRUE)
  
  tibble(
    Metric = c("Total Return", "Annualised Sharpe Ratio", "Maximum Drawdown",
               "Win Rate (%)", "Total Trades",
               "Stop-Loss Triggers (+/-3.5)",
               "- Long-spread stops (z <= -3.5)",
               "- Short-spread stops (z >= +3.5)",
               "Profit-Target Exits (z -> 0)",
               "Closed at End of Sample",
               "Forced Exits (beta <= 0 / invalid data)"),
    Strategy = c(pct(tot_ret(strat)),
                 sprintf("%.2f", sharpe(strat)),
                 pct(-mdd(strat)),
                 if (is.na(win_rate)) "n/a" else sprintf("%.1f", 100 * win_rate),
                 as.character(n_trades),
                 as.character(n_stop), as.character(n_stop_lo), as.character(n_stop_hi),
                 as.character(n_target), as.character(n_open), as.character(n_forced)),
    `50/50 Buy & Hold` = c(pct(tot_ret(bench)),
                           sprintf("%.2f", sharpe(bench)),
                           pct(-mdd(bench)),
                           rep("n/a", 8))
  )
}

# 8. Charts
plot_dynamic_spread <- function(sig, split_date, cfg) {
  p1 <- "1. Adjusted close (USD)"
  p2 <- "2. Rolling hedge ratio (beta)"
  p3 <- "3. Rolling 60-day Z-score"
  
  long <- dplyr::bind_rows(
    sig %>% dplyr::select(date, NVDA, AVGO) %>%
      tidyr::pivot_longer(c(NVDA, AVGO), names_to = "series", values_to = "value") %>%
      dplyr::mutate(panel = p1),
    sig %>% dplyr::transmute(date, series = "Beta", value = beta, panel = p2),
    sig %>% dplyr::transmute(date, series = "Z-score", value = z, panel = p3)
  ) %>%
    dplyr::mutate(panel = factor(panel, levels = c(p1, p2, p3)))
  
  levels_df <- tibble(
    panel = factor(p3, levels = c(p1, p2, p3)),
    y     = c(cfg$entry_z, -cfg$entry_z, 0, cfg$stop_z, -cfg$stop_z),
    level = c("Entry (+/-2.0)", "Entry (+/-2.0)", "Target (0)",
              "Stop-loss (+/-3.5)", "Stop-loss (+/-3.5)")
  )
  
  ggplot(long, aes(date, value, colour = series)) +
    geom_line(linewidth = 0.45, na.rm = TRUE) +
    geom_hline(data = levels_df, aes(yintercept = y, linetype = level),
               colour = "grey30", linewidth = 0.4) +
    geom_vline(xintercept = split_date, linetype = "dotted", colour = "black") +
    facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = c(NVDA = "#76b900", AVGO = "#cc092f",
                                   Beta = "#1f4e79", `Z-score` = "#6a1b9a")) +
    scale_linetype_manual(values = c(`Entry (+/-2.0)` = "dashed",
                                     `Target (0)` = "solid",
                                     `Stop-loss (+/-3.5)` = "longdash")) +
    labs(title = "NVDA / AVGO: Dynamic Spread Diagnostics",
         subtitle = "Dotted vertical line = start of out-of-sample period",
         x = NULL, y = NULL, colour = NULL, linetype = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold", hjust = 0))
}

plot_equity_curves <- function(bt) {
  bt %>%
    dplyr::select(date,
                  `Pairs strategy`   = strat_equity,
                  `50/50 buy & hold` = bench_equity) %>%
    tidyr::pivot_longer(-date, names_to = "portfolio", values_to = "equity") %>%
    ggplot(aes(date, equity, colour = portfolio)) +
    geom_hline(yintercept = 1, colour = "grey60") +
    geom_line(linewidth = 0.7) +
    scale_colour_manual(values = c(`Pairs strategy` = "#1f4e79",
                                   `50/50 buy & hold` = "#e67e22")) +
    labs(title = "Out-of-Sample Performance Attribution",
         subtitle = "Growth of 1 unit of capital (net of transaction costs)",
         x = NULL, y = "Equity (start = 1)", colour = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

# 9. Pipeline
run_pipeline <- function(cfg, px_tbl = NULL) {
  if (is.null(px_tbl)) px_tbl <- get_prices(cfg)
  parts <- split_sample(px_tbl, cfg$years_train)
  
  coint <- test_cointegration(parts$train, cfg)
  cat("\n=== In-sample Engle-Granger / ADF test (residuals of NVDA ~ AVGO) ===\n")
  print(as.data.frame(coint), row.names = FALSE)
  if (!coint$reject_at_eg_5pct) {
    msg <- "In-sample residuals do NOT reject a unit root at the ~5% Engle-Granger level."
    if (cfg$require_cointegration) stop(msg) else warning(msg, call. = FALSE)
  }
  
  sig     <- build_signal_frame(px_tbl, parts$split_date, cfg)
  res     <- run_backtest(sig, cfg)
  metrics <- compute_metrics(res$bt, res$trades, cfg)
  
  cat("\n=== Out-of-sample analytics ===\n")
  print(knitr::kable(metrics, format = "markdown", align = c("l", "r", "r")))
  
  chart1 <- plot_dynamic_spread(sig, parts$split_date, cfg)
  chart2 <- plot_equity_curves(res$bt)
  print(chart1)
  print(chart2)
  ggsave("chart1_dynamic_spread.png",    chart1, width = 10, height = 9, dpi = 150)
  ggsave("chart2_equity_curves.png",     chart2, width = 10, height = 5, dpi = 150)
  
  invisible(list(prices = px_tbl, signals = sig, backtest = res$bt,
                 trades = res$trades, metrics = metrics, coint = coint,
                 chart1 = chart1, chart2 = chart2))
}

# 10. Run
results <- run_pipeline(cfg)