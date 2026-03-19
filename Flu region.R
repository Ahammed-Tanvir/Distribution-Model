library(readr)
library(data.table)
library(dplyr)
library(lubridate)
library(tidyr)
library(pracma)
library(Metrics)
library(ggplot2)

# ══════════════════════════════════════════════════════════════════════════════
# 1. IMPORT DATA
# ══════════════════════════════════════════════════════════════════════════════

MUSC_Weekly_Influenza_Region <- read_csv(
   '/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/MUSC/Infectious Disease EHR/Weekly Data/Latest Weekly Data/MUSC_Weekly_Influenza_Region_dx_cond_lab_Incident.csv'
)

Prisma_Weekly_Influenza_Region <- read_csv(
   '/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/Prisma Health/Infectious Disease EHR/Weekly Data/Latest Weekly Data/Prisma_Health_Weekly_Influenza_Region_dx_cond_lab_Incident.csv'
)

CDC <- read_csv("/Users/tanvirahammed/Downloads/target-hospital-admissions-8.csv") %>%
   filter(location == "45") %>%
   dplyr::select(Week = date, Weekly_Total = value)

state <- read_csv(
   '/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Tanvir/FluSight/Tanvir_QR/State/CDC_submission/2026-03-14-DMAPRIME-QR.csv'
) %>%
   filter(output_type_id == 0.5) %>%
   select(target_end_date, value) %>%
   rename(Week = target_end_date, Weekly_Total = value)

RFA_region <- read_csv(
   '/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/SC Health Records/RFA/Years 2017-2025/Weekly data/RFA_weekly_influenza_region_incident.csv'
) %>%
   select(Week = week_end, Region = region, Actual_Count = weekly_IP) %>%
   filter(format(Week, "%Y") %in% c("2021", "2022", "2023", "2024", "2025")) %>%
   group_by(Week) %>%
   mutate(Weekly_Total = sum(Actual_Count, na.rm = TRUE)) %>%
   ungroup()


# ══════════════════════════════════════════════════════════════════════════════
# 2. BUILD HEALTHSYSTEM DATA
# ══════════════════════════════════════════════════════════════════════════════

Healthsystem_Weekly_Influenza_Region <- full_join(
   MUSC_Weekly_Influenza_Region %>%
      dplyr::select(Region, Week, Weekly_Tests, Weekly_Positive_Tests, Weekly_Inpatient_Hospitalizations) %>%
      dplyr::rename(
         MUSC_Tests            = Weekly_Tests,
         MUSC_Positive_Tests   = Weekly_Positive_Tests,
         MUSC_Hosp             = Weekly_Inpatient_Hospitalizations
      ),
   Prisma_Weekly_Influenza_Region %>%
      dplyr::select(Region, Week, Weekly_Tests, Weekly_Positive_Tests, Weekly_Inpatient_Hospitalizations) %>%
      dplyr::rename(
         Prisma_Tests          = Weekly_Tests,
         Prisma_Positive_Tests = Weekly_Positive_Tests,
         Prisma_Hosp           = Weekly_Inpatient_Hospitalizations
      ),
   by = c("Region", "Week")
) %>%
   mutate(
      Weekly_Tests                      = rowSums(across(c(MUSC_Tests, Prisma_Tests)),                   na.rm = TRUE),
      Weekly_Positive_Tests             = rowSums(across(c(MUSC_Positive_Tests, Prisma_Positive_Tests)), na.rm = TRUE),
      Weekly_Inpatient_Hospitalizations = rowSums(across(c(MUSC_Hosp, Prisma_Hosp)),                    na.rm = TRUE)
   ) %>%
   dplyr::select(Region, Week, Weekly_Tests, Weekly_Positive_Tests, Weekly_Inpatient_Hospitalizations) %>%
   filter(year(Week) %in% 2021:2026)


# ══════════════════════════════════════════════════════════════════════════════
# 3. SMOOTHING FUNCTIONS
# ══════════════════════════════════════════════════════════════════════════════

smooth_fun <- function(data, window_size = 2) {
   n      <- length(data)
   result <- numeric(n)
   for (i in 1:n) {
      start     <- max(1, i - window_size)
      end       <- min(n, i + window_size)
      result[i] <- mean(data[start:end])
   }
   result
}

smooth_between_peaks <- function(y, peaks, window_size = 2) {
   n          <- length(y)
   smoothed_y <- y
   peaks      <- sort(unique(c(1, peaks, n)))
   for (i in seq_along(peaks[-1])) {
      start <- peaks[i]; end <- peaks[i + 1]
      if ((end - start + 1) > 3) {
         segment   <- y[start:end]
         segment_s <- smooth_fun(segment, window_size)
         segment[2:(length(segment) - 1)] <- segment_s[2:(length(segment_s) - 1)]
         smoothed_y[start:end] <- segment
      }
   }
   smoothed_y[peaks] <- y[peaks]
   smoothed_y
}

process_smooth <- function(y, distance = 8, prominence_percentile = 0.25, window_size = 2) {
   y <- as.numeric(y)
   if (all(is.na(y)) || length(na.omit(y)) < 5) return(rep(NA, length(y)))
   threshold    <- quantile(y, prominence_percentile, na.rm = TRUE)
   peaks_matrix <- tryCatch(
      pracma::findpeaks(y, minpeakdistance = distance, minpeakheight = threshold),
      error = function(e) NULL
   )
   if (is.null(peaks_matrix) || nrow(peaks_matrix) == 0) {
      smooth_fun(y, window_size)
   } else {
      peaks      <- peaks_matrix[, 2]
      y_smoothed <- smooth_between_peaks(y, peaks, window_size)
      y_smoothed[1:length(y)]
   }
}


# ══════════════════════════════════════════════════════════════════════════════
# 4. BUILD REGION_WEEKLY_COUNTS
# ══════════════════════════════════════════════════════════════════════════════

Region_weekly_counts <- RFA_region %>%
   group_by(Region) %>%
   arrange(Week) %>%
   mutate(Smoothed_Actual_Count = process_smooth(Actual_Count)) %>%
   ungroup() %>%
   mutate(Region = as.character(Region))

Region_weekly_counts <- Region_weekly_counts %>%
   full_join(
      Healthsystem_Weekly_Influenza_Region %>% mutate(Region = as.character(Region)),
      by = c("Region", "Week")
   )

Region_weekly_counts <- Region_weekly_counts %>%
   left_join(CDC, by = "Week", suffix = c("", "_cdc")) %>%
   mutate(Weekly_Total = coalesce(Weekly_Total, Weekly_Total_cdc)) %>%
   select(-Weekly_Total_cdc)

# ── Horizon definitions ───────────────────────────────────────────────────────
# horizon -1 = 0-wk ahead: target 2026-03-07, state total from CDC (observed)
# horizon  0 = nowcast:    target 2026-03-14, state total from state model
# horizon  1 = 1-wk ahead: target 2026-03-21, state total from state model
# horizon  2 = 2-wk ahead: target 2026-03-28, state total from state model
# horizon  3 = 3-wk ahead: target 2026-04-04, state total from state model

horizon_map <- tibble(
   horizon       = c(-1,             0,             1,             2,             3),
   forecast_week = as.Date(c(
      "2026-03-07", "2026-03-14", "2026-03-21", "2026-03-28", "2026-04-04"
   )),
   lag_n         = c(0, 1, 2, 3, 4),
   lag_col       = c("Weekly_Tests",      "Weekly_Tests_lag1", "Weekly_Tests_lag2",
                     "Weekly_Tests_lag3", "Weekly_Tests_lag4"),
   label         = c("0-wk (lag0)", "Nowcast (lag1)", "1-wk (lag2)",
                     "2-wk (lag3)", "3-wk (lag4)")
)

# State totals: 0-wk from CDC, rest from state model
state_totals <- bind_rows(
   CDC %>% filter(Week == as.Date("2026-03-07")) %>%
      rename(forecast_week = Week, state_total = Weekly_Total),
   state %>% rename(forecast_week = Week, state_total = Weekly_Total)
)

horizon_map <- horizon_map %>%
   left_join(state_totals, by = "forecast_week")

cat("\n── Horizon map ──\n")
print(horizon_map)

# Append all forecast weeks that aren't already in data
existing_weeks <- unique(Region_weekly_counts$Week)
new_weeks      <- horizon_map$forecast_week[!horizon_map$forecast_week %in% existing_weeks]

if (length(new_weeks) > 0) {
   new_rows <- expand_grid(
      Week   = new_weeks,
      Region = unique(Region_weekly_counts$Region)
   )
   Region_weekly_counts <- bind_rows(Region_weekly_counts, new_rows)
}

# Time features and all lags
Region_weekly_counts <- Region_weekly_counts %>%
   mutate(
      Year    = year(Week),
      Week_No = week(Week),
      Region  = as.factor(Region),
      Week_No = as.factor(Week_No),
      Year    = as.factor(Year)
   ) %>%
   mutate(Smoothed_Actual_Count = round(Smoothed_Actual_Count)) %>%
   arrange(Region, Week) %>%
   group_by(Region) %>%
   mutate(
      Weekly_Tests_lag1 = lag(Weekly_Tests, 1),
      Weekly_Tests_lag2 = lag(Weekly_Tests, 2),
      Weekly_Tests_lag3 = lag(Weekly_Tests, 3),
      Weekly_Tests_lag4 = lag(Weekly_Tests, 4)
   ) %>%
   ungroup() %>%
   mutate(Region = recode(Region, "Pee Dee" = "Pee_Dee")) %>%
   data.table()


# ══════════════════════════════════════════════════════════════════════════════
# 5. TRAIN / TEST SPLIT
# ══════════════════════════════════════════════════════════════════════════════

Region_weekly_counts_train <- Region_weekly_counts[
   Week >= as.Date("2022-09-01") & Week < as.Date("2023-04-01")
]

Region_weekly_counts_test <- Region_weekly_counts[
   Week >= as.Date("2024-11-01") & Week < as.Date("2025-02-15")
]


# ══════════════════════════════════════════════════════════════════════════════
# 6. COMPUTE RATE SUMMARIES FOR ALL LAGS (TRAIN ONLY)
# ══════════════════════════════════════════════════════════════════════════════

compute_rate_summary <- function(train_data, lag_col) {
   
   predictor_col <- lag_col   # e.g. "Weekly_Tests" or "Weekly_Tests_lag1"
   
   train_data %>%
      as_tibble() %>%
      filter(!is.na(Actual_Count), !is.na(.data[[predictor_col]]),
             .data[[predictor_col]] > 0) %>%
      mutate(Rate = Actual_Count / .data[[predictor_col]]) %>%
      group_by(Region) %>%
      summarise(
         Mean_Rate   = mean(Rate,   na.rm = TRUE),
         Median_Rate = median(Rate, na.rm = TRUE),
         .groups     = "drop"
      )
}

# Compute one rate table per lag
rate_summaries <- horizon_map %>%
   select(lag_n, lag_col) %>%
   distinct() %>%
   mutate(
      rate_table = purrr::map(
         lag_col,
         ~ compute_rate_summary(Region_weekly_counts_train, .x)
      )
   )

cat("\n── Rate summaries per lag ──\n")
for (i in seq_len(nrow(rate_summaries))) {
   cat("lag", rate_summaries$lag_n[i], "—", rate_summaries$lag_col[i], "\n")
   print(rate_summaries$rate_table[[i]])
}


# ══════════════════════════════════════════════════════════════════════════════
# 7. HELPER: DISTRIBUTE STATE TOTAL
# ══════════════════════════════════════════════════════════════════════════════

distribute_rate_test <- function(week_data, rate_table, lag_col, weekly_total) {
   
   week_data %>%
      mutate(Region = as.character(Region)) %>%
      left_join(rate_table %>% mutate(Region = as.character(Region)),
                by = "Region") %>%
      mutate(
         Tests_used        = .data[[lag_col]],
         Expected          = Tests_used * Mean_Rate,
         Share             = Expected / sum(Expected, na.rm = TRUE),
         Distributed_Count = round(Share * weekly_total)
      ) %>%
      select(Region, Tests_used, Expected, Share, Distributed_Count)
}


# ══════════════════════════════════════════════════════════════════════════════
# 8. RETROSPECTIVE VALIDATION — ALL HORIZONS ON TEST SET
# ══════════════════════════════════════════════════════════════════════════════

test_weeks <- sort(unique(Region_weekly_counts_test$Week))
retro_list <- list()

for (h in seq_len(nrow(horizon_map))) {
   
   lag_col   <- horizon_map$lag_col[h]
   lag_n     <- horizon_map$lag_n[h]
   label     <- horizon_map$label[h]
   rate_table <- rate_summaries$rate_table[[
      which(rate_summaries$lag_col == lag_col)
   ]]
   
   for (hw in as.list(test_weeks)) {
      hw <- as.Date(hw)
      
      week_data <- Region_weekly_counts_test %>%
         as_tibble() %>%
         filter(Week == hw)
      
      if (nrow(week_data) != 4 | any(is.na(week_data[[lag_col]]))) next
      
      wt     <- first(week_data$Weekly_Total)
      result <- distribute_rate_test(week_data, rate_table, lag_col, wt)
      
      retro_list[[paste(hw, lag_n)]] <- tibble(
         Week              = hw,
         horizon           = lag_n,
         label             = label,
         Region            = result$Region,
         Weekly_Total      = wt,
         Share             = result$Share,
         Distributed_Count = result$Distributed_Count
      ) %>%
         left_join(
            week_data %>% as_tibble() %>%
               mutate(Region = as.character(Region)) %>%
               select(Region, Actual_Count),
            by = "Region"
         )
   }
}

retro <- bind_rows(retro_list) %>%
   filter(!is.na(Actual_Count)) %>%
   mutate(
      Ai = ifelse(
         Actual_Count == 0 & Distributed_Count == 0,
         1,
         pmin(Actual_Count, Distributed_Count) /
            pmax(Actual_Count, Distributed_Count)
      )
   )


# ══════════════════════════════════════════════════════════════════════════════
# 9. VALIDATION METRICS — ALL HORIZONS
# ══════════════════════════════════════════════════════════════════════════════

# ── PA by horizon ─────────────────────────────────────────────────────────────
cat("\n── Overall Mean PA by Horizon ──\n")
retro %>%
   group_by(horizon, label) %>%
   summarise(Mean_Ai = mean(Ai, na.rm = TRUE), .groups = "drop") %>%
   arrange(horizon) %>%
   print()

# ── PA by region × horizon ────────────────────────────────────────────────────
cat("\n── Mean PA by Region × Horizon ──\n")
retro %>%
   group_by(Region, horizon, label) %>%
   summarise(Mean_Ai = mean(Ai, na.rm = TRUE), .groups = "drop") %>%
   arrange(Region, horizon) %>%
   pivot_wider(names_from = label, values_from = Mean_Ai) %>%
   mutate(across(where(is.numeric), ~ round(., 3))) %>%
   print()

# ── Full metrics by region × horizon ─────────────────────────────────────────
cat("\n── Full Metrics by Region × Horizon ──\n")
retro_metrics <- retro %>%
   group_by(Region, horizon, label) %>%
   summarise(
      Mean_Ai = mean(Ai,                                                na.rm = TRUE),
      MAE     = mae(Actual_Count, Distributed_Count),
      RMSE    = rmse(Actual_Count, Distributed_Count),
      MAPE    = mean(abs((Actual_Count - Distributed_Count) /
                            pmax(Actual_Count, 1))) * 100,
      Corr    = cor(Actual_Count, Distributed_Count, use = "complete.obs"),
      n       = n(),
      .groups = "drop"
   ) %>%
   arrange(Region, horizon)

print(retro_metrics)


# ══════════════════════════════════════════════════════════════════════════════
# 10. PLOT: ALL HORIZONS VS ACTUAL — FACETED BY REGION
# ══════════════════════════════════════════════════════════════════════════════

# Colour palette for 5 horizons
horizon_colors <- c(
   "0-wk (lag0)"    = "#2166ac",
   "Nowcast (lag1)" = "#4dac26",
   "1-wk (lag2)"    = "#f4a582",
   "2-wk (lag3)"    = "#d6604d",
   "3-wk (lag4)"    = "#92c5de"
)

horizon_linetypes <- c(
   "0-wk (lag0)"    = "solid",
   "Nowcast (lag1)" = "solid",
   "1-wk (lag2)"    = "solid",
   "2-wk (lag3)"    = "solid",
   "3-wk (lag4)"    = "solid"
)

# Actual counts — one row per week per region (no horizon dimension)
actuals <- retro %>%
   select(Week, Region, Actual_Count) %>%
   distinct()

ggplot() +
   # Actual line — black, thick
   geom_line(
      data = actuals,
      aes(x = Week, y = Actual_Count),
      color    = "black",
      linewidth = 1.2
   ) +
   # Predicted lines — one per horizon
   geom_line(
      data = retro,
      aes(x = Week, y = Distributed_Count,
          color    = label,
          linetype = label),
      linewidth = 0.8
   ) +
   facet_wrap(~ Region, scales = "free_y", ncol = 2) +
   scale_color_manual(
      name   = "Horizon",
      values = horizon_colors
   ) +
   scale_linetype_manual(
      name   = "Horizon",
      values = horizon_linetypes
   ) +
   labs(
      title    = "Regional Influenza Nowcast & Forecast: All Horizons vs Actual",
      subtitle = "Retrospective Validation Nov 2024 – Feb 2025 | Rate-by-Test Method | Train: Jan 2021 – Oct 2024",
      x        = "Week",
      y        = "Influenza Hospitalizations"
   ) +
   theme_bw(base_size = 12) +
   theme(
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold"),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank()
   ) +
   guides(
      color    = guide_legend(nrow = 1),
      linetype = guide_legend(nrow = 1)
   )


# ══════════════════════════════════════════════════════════════════════════════
# 11. PRODUCTION FORECASTS — ALL HORIZONS
# ══════════════════════════════════════════════════════════════════════════════

forecast_list <- list()

for (h in seq_len(nrow(horizon_map))) {
   
   fw         <- horizon_map$forecast_week[h]
   lag_col    <- horizon_map$lag_col[h]
   wt         <- horizon_map$state_total[h]
   label      <- horizon_map$label[h]
   rate_table <- rate_summaries$rate_table[[
      which(rate_summaries$lag_col == lag_col)
   ]]
   
   forecast_data <- Region_weekly_counts %>%
      as_tibble() %>%
      filter(Week == fw)
   
   cat("\n── Predictor check:", label, "(", as.character(fw), ") ──\n")
   forecast_data %>% select(Region, all_of(lag_col)) %>% print()
   
   result <- distribute_rate_test(forecast_data, rate_table, lag_col, wt)
   
   forecast_list[[h]] <- tibble(
      Week              = fw,
      horizon           = horizon_map$horizon[h],
      label             = label,
      Region            = result$Region,
      Weekly_Total      = wt,
      Share             = result$Share,
      Distributed_Count = result$Distributed_Count
   )
}

all_forecasts <- bind_rows(forecast_list)

cat("\n── Production Forecasts: All Horizons ──\n")
all_forecasts %>%
   select(Week, label, Region, Weekly_Total, Share, Distributed_Count) %>%
   mutate(Share = round(Share, 3)) %>%
   arrange(Week, Region) %>%
   print()

# ── Summary table: distributed counts by region and horizon ──────────────────
cat("\n── Distributed Count Summary ──\n")
all_forecasts %>%
   select(Week, label, Region, Distributed_Count) %>%
   pivot_wider(names_from = label, values_from = Distributed_Count) %>%
   arrange(Week, Region) %>%
   print()


# ══════════════════════════════════════════════════════════════════════════════
# 12. ATTACH ALL FORECASTS BACK TO MAIN DATA
# ══════════════════════════════════════════════════════════════════════════════
# Use horizon 0 (nowcast, lag1) as primary fill for Actual_Count
# Other horizons stored as separate columns for reference

nowcast_fill <- all_forecasts %>%
   filter(horizon == 0) %>%
   select(Week, Region, Distributed_Count) %>%
   rename(Distributed_Count_nowcast = Distributed_Count)

Region_weekly_counts <- Region_weekly_counts %>%
   as_tibble() %>%
   left_join(
      all_forecasts %>%
         select(Week, Region, horizon, Distributed_Count) %>%
         pivot_wider(names_from = horizon, values_from = Distributed_Count,
                     names_prefix = "Distributed_h"),
      by = c("Week", "Region")
   ) %>%
   left_join(nowcast_fill, by = c("Week", "Region")) %>%
   mutate(
      Actual_Count = coalesce(Actual_Count, as.numeric(Distributed_Count_nowcast))
   ) %>%
   data.table()

cat("\n── Final check: forecast weeks in Region_weekly_counts ──\n")
Region_weekly_counts[
   Week %in% horizon_map$forecast_week,
   .(Week, Region, Actual_Count,
     `Distributed_h-1`, Distributed_h0, Distributed_h1,
     Distributed_h2,    Distributed_h3)
] %>% print()
