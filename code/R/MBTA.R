################################################################################
# Loading Packages
################################################################################
library(AER)
library(quantmod)
library(tidyr)
library(ggplot2)
library(tseries)
library(dplyr)
library(forecast)
library(vars)
library(lubridate)
library(fastDummies)
library(purrr)
library(broom)

################################################################################
# Loading Data
################################################################################
# The line variable can be changed to any line of interest
# Additionally the start and end times here set the parameters of the training set
################################################################################

line <- "orange"  
start_time <- as.POSIXct("2023-01-01 00:00:00")
end_time <- as.POSIXct("2023-12-31 23:59:59")


line_url <- sprintf("https://raw.githubusercontent.com/zwinship/MBTA_Time_Series/refs/heads/main/data/processed/lines/%s_line.csv", line)
day_url <- sprintf("https://raw.githubusercontent.com/zwinship/MBTA_Time_Series/refs/heads/main/data/processed/line_day/%s_line_day.csv", line)


# This data frame is the main data used for all the following code
df <- read.csv(line_url)
# This data frame is used for day aggregated visualizations only, at the end of the code
dfday <- read.csv(day_url)


################################################################################
# Data Preprocessing and Cleaning
################################################################################
# Drop the index column from python (0 index)
df <- df[, !names(df) %in% "X"]

# Make the time column into R date time
df$time <- as.POSIXct(df$time, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

# Temp dropping glitched time rows
df <- df[!is.na(df$time), ]

# Creating complete time series with consistent intervals
seq_start_time <- min(df$time)
seq_end_time <- max(df$time)
complete_time_index <- seq(from = seq_start_time, to = seq_end_time, by = "30 min", tz = "UTC")
complete_time_df <- data.frame(time = complete_time_index)
df <- merge(complete_time_df, df, by = "time", all.x = TRUE)

# Forward filling these new rows with the last known value
df <- df %>% fill(everything(), .direction = "down")

# Adding date features
# First are the hour, weekday, and month
# The hours dummies never actually got used in the paper or further code
df <- df %>%
  mutate(
    date = as.Date(date),
    hour = as.numeric(format(time, "%H")),  # Extracts hour from datetime
    weekday = as.numeric(factor(weekdays(date), levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"))),
    month = as.numeric(factor(months(date), levels = month.name))
  ) %>%
  dummy_cols(select_columns = c("weekday", "month", "hour"), remove_selected_columns = TRUE)

# Hour dummies
df <- df %>%
  mutate(
    hour = as.numeric(format(time, "%H"))  # Extracts 0–23
  ) %>%
  dummy_cols(select_columns = "hour", remove_selected_columns = TRUE)

# Closes dummy
df <- df %>%
  mutate(
    closed = ifelse(format(time, "%H:%M") >= "00:00" & format(time, "%H:%M") < "05:00", 1, 0)
  )

# Rush hour dummy
df <- df %>%
  mutate(
    rush_hour = case_when(
      (
        (format(time, "%H:%M") >= "06:00" & format(time, "%H:%M") < "09:30") |
          (format(time, "%H:%M") >= "15:00" & format(time, "%H:%M") < "18:30")
      ) & weekday_6 == 0 & weekday_7 == 0 ~ 1,
      TRUE ~ 0
    )
  )


# GE 7 day lag
df <- df %>%
  mutate(gated_entries_lag7d = dplyr::lag(gated_entries, 336))

# GE 1 day lag
df <- df %>%
  mutate(gated_entries_lag1d = dplyr::lag(gated_entries, 48))

################################################################################
# Training Data Preparation
################################################################################
# Creating constants for the bounds of the model training
start_index <- which.min(abs(df$time - start_time))
end_index <- which.min(abs(df$time - end_time))

# Subsetting the df for the training time period
df_train <- df[start_index:(end_index - 1), ]
df_valid <- df[df$time > end_time, ]

################################################################################
# Time Series Creation
################################################################################
# Fixing log problems by filling in epsilon
epsilon <- 1e-6
df_train$gated_entries[df_train$gated_entries <= 0] <- epsilon

# Making the basic time series used in the code
ts <- ts(df_train$gated_entries, frequency = 48, start = c(2023, 1, 1))
xts <- xts(df_train$gated_entries, df_train$time)

# Creating differenced time series
xts_diff <- xts(diff(df_train$gated_entries), df_train$time[-1])
# frequency(xts_diff) <- 48

# Creating log time series with handling for zero values
xts_log <- xts(log(df_train$gated_entries), df_train$time)


################################################################################
# Visualization
################################################################################
# Basic time series plot
# These visualizations just end up looking very bad
# The frequency of the data makes it preety hard to actually tell whats going on
ggplot(df_train, aes(x = time, y = gated_entries)) +
  geom_line(color = "blue") +
  labs(title = "Daily Gated Entries Over Time", 
       x = "Date", 
       y = "Daily Gated Entries")

# Time Series Decomposition
decomposed_ts <- decompose(ts)
plot(decomposed_ts)

################################################################################
# Stationarity Testing
################################################################################
# ADF Tests on the 3 time series of interest
adf.test(xts, k = 48)

adf.test(xts_diff, k = 48)

adf.test(xts_log, k = 48)
# All 3 are concluded to be stationary


# Autocorrelation Function parameter
# The parameter can be changed to any of the transformed series
# Doing this instead of copy pasting the code 3 seperate time
################################################################################
# Set transformation parameter for autocorrelation analysis
# Choose one of: "level", "diff", "log"
autocorr_transformation <- "level"

################################################################################
# ACF and PACF Analysis
################################################################################
# Function to create ACF/PACF plots based on selected transformation
create_correlation_plot <- function(transformation, func_type) {
  # Select time series based on transformation parameter
  if (transformation == "level") {
    ts_data <- xts
    title_prefix <- "Gated Entries"
  } else if (transformation == "diff") {
    ts_data <- xts_diff
    title_prefix <- "Differenced Gated Entries"
  } else if (transformation == "log") {
    ts_data <- xts_log
    title_prefix <- "Log Gated Entries"
  } else {
    stop("Invalid transformation parameter. Choose 'level', 'diff', or 'log'.")
  }
  
  # Create ACF or PACF plot with capitalized function name
  if (func_type == "acf") {
    func_name <- "ACF"
    acf(as.numeric(ts_data), xaxt = "n", lag.max = 48 * 7, 
        main = paste(func_name, "for", title_prefix, "Lagged to One Week (Level)"))
  } else if (func_type == "pacf") {
    func_name <- "PACF"
    pacf(as.numeric(ts_data), xaxt = "n", lag.max = 48 * 7, 
         main = paste(func_name, "for", title_prefix, "Lagged to One Week (Level)", toupper(line), 'Line'))
  }
  
  # Adjusting x-axis and lines for each day
  # Currently the way that these functions work the x-axis doesn't look very clean
  # Therefore I replace the x-axis with the days of the week
  # Preety easy to do if you know the length of the x-axis and the time the plot strats at
  day_boundaries <- seq(0, 48 * 7, by = 48)
  day_midpoints <- day_boundaries[-length(day_boundaries)] + 24
  abline(v = day_boundaries, col = "gray", lty = 2)
  axis(1, at = day_midpoints, labels = paste0("Day ", 1:7))
}

# Generate ACF plot based on selected transformation
create_correlation_plot(autocorr_transformation, "acf")

# Generate PACF plot based on selected transformation
create_correlation_plot(autocorr_transformation, "pacf")



################################################################################
# Model Creation
################################################################################
# Creating a paramter of the exogenous variables used in the models
xreg_vars <- df_train %>%
  dplyr::select(closed, rush_hour, gated_entries_lag7d, starts_with("weekday_"), starts_with("month_")) %>%
  dplyr::select(-weekday_1, -month_1) %>% 
  as.matrix()


# Training the actual models
# trace = TRUE allows me to actually see the AIC as each model is tried
# Hence I can see models that end in INF AIC, meaning these models diverged
arimax_model <- auto.arima(xts, xreg = xreg_vars, trace = TRUE)
arimax_model_diff <- auto.arima(xts_diff, xreg = xreg_vars[-1, ], trace = TRUE)
arimax_model_log <- auto.arima(xts_log, xreg = xreg_vars, trace = TRUE)


# Checking and plotting residuals
level_res <- residuals(sarimax_model)
log_res <- residuals(sarimax_model_log)

checkresiduals(sarimax_model)
checkresiduals(sarimax_model_log)


################################################################################
# Creating the function to actually forecast these models
# Used AI for most of this
# Since the handout forecasting in class didnt really put these forecasts into a df
# I used this instead so the forecasts can be analyzed easier and plotted easier further into the code
################################################################################
forecast_day <- function(model, date, df_valid, model_type) {
  # Create sequence of times from 5:00 AM to midnight for the given day
  date_str <- format(date, "%Y-%m-%d")
  
  forecast_times <- seq(
    from = as.POSIXct(paste(date_str, "05:00:00"), tz = "UTC"),
    to = as.POSIXct(paste(date_str, "23:59:59"), tz = "UTC"),
    by = "30 min"
  )
  
  
  xreg_subset <- df_valid %>%
    filter(time %in% forecast_times) %>%
    dplyr::select(closed, rush_hour, gated_entries_lag7d, starts_with("weekday_"), starts_with("month_")) %>%
    dplyr::select(-weekday_1, -month_1) %>% 
    as.matrix()
  
  # Generate forecast for the number of periods needed
  n_periods <- length(forecast_times)
  forecast_result <- forecast(model, h = n_periods, xreg = xreg_subset)
  
  # Create a data frame with the forecast results
  if (grepl("log", model_type)) {
    forecast_df <- data.frame(
      timestamp = forecast_times,
      log_forecast = exp(as.numeric(forecast_result$mean)),
      log_lower_95 = exp(as.numeric(forecast_result$lower)),
      log_upper_95 = exp(as.numeric(forecast_result$upper))
    )
  } else {
    forecast_df <- data.frame(
      timestamp = forecast_times,
      level_forecast = as.numeric(forecast_result$mean),
      level_lower_95 = as.numeric(forecast_result$lower),
      level_upper_95 = as.numeric(forecast_result$upper)
    )
  }
  return(forecast_df)
}

################################################################################
# Forcasting the log model
################################################################################
unique_days <- unique(as.Date(df_valid$time))

df_forecast <- data.frame()

processed_days <- c()

for (day in unique_days) {
  day <- as.Date(day)
  
  # Check if we've already processed this day
  if (day %in% processed_days) {
    next
  }
  
  day_forecast <- forecast_day(sarimax_model_log, day, df_valid, "log")
  
  # Add to our processed days list
  processed_days <- c(processed_days, day)
  
  # Append to the all_forecasts dataframe
  df_forecast <- rbind(df_forecast, day_forecast)
  if (as.character(day) == "2025-02-28") break
}


df_forecast <- df_forecast %>%
  distinct(timestamp, .keep_all = TRUE)

df_forecast <- df_forecast %>%
  rename(time = timestamp)

df_valid <- df_valid %>%
  left_join(df_forecast %>% dplyr::select(time, log_forecast, log_lower_95, log_upper_95), by = "time")

df_valid$log_forecast_residual <- df_valid$gated_entries - df_valid$log_forecast


################################################################################
# Forecasting the standard level model
################################################################################

unique_days <- as.Date(unique(as.Date(df_valid$time)))

df_forecast <- data.frame()

for (day in unique_days) {
  # Generate forecast for this day using the pre-fitted model
  day <- as.Date(day)
  day_forecast <- forecast_day(sarimax_model, day, df_valid, "level")
  
  # Append to the all_forecasts dataframe
  df_forecast <- rbind(df_forecast, day_forecast)
  if (as.character(day) == "2025-02-28") break
}


df_forecast <- df_forecast %>%
  distinct(timestamp, .keep_all = TRUE)

df_forecast <- df_forecast %>%
  rename(time = timestamp)

df_valid <- df_valid %>%
  left_join(df_forecast %>% dplyr::select(time, level_forecast, level_lower_95, level_upper_95), by = "time")

df_valid$level_forecast_residual <- df_valid$gated_entries - df_valid$level_forecast



################################################################################
# Comparing accuracy
################################################################################
accuracy(df_valid$log_forecast, df_valid$gated_entries)
accuracy(df_valid$level_forecast, df_valid$gated_entries)


################################################################################
# Plotting Residuals (Normallity)
################################################################################
hist(df_valid$log_forecast_residual, breaks = 250, 
     main = 'Distribution of Forecast Residual (Log) Orange Line',
     xlab = 'Forecast Residual (Log)')
hist(df_valid$level_forecast_residual, breaks = 250,
     main = 'Distribution of Forecast Residual (Level) Orange Line',
     xlab = 'Forecast Residual (Level)')



################################################################################
# Model Tests for Normality
################################################################################
jarque.bera.test(na.omit(log_res))
jarque.bera.test(na.omit(level_res))


Box.test(na.omit(log_res, lag = 48))
Box.test(na.omit(level_res, lag = 48))


jarque.bera.test(na.omit(df_valid$log_forecast_residual))
jarque.bera.test(na.omit(df_valid$level_forecast_residual))



################################################################################
# DISCLAIMER
# Some of the following code for the supply optimization was done using AI
# I understand most of it was beyond the handouts in class
# But I wanted to expand the paper to actual use the forecasts I know how to do from class
# Even though most of the code isn't from class, it's actually nothing super complex
# Most of it just base R code and a lot of dplyr and some ggplot
# Anything that isn't dplyr below I did not personally write, although there isn't much
# I know some ggplot but I had to google a lot and use ai for a lot of the options in ggplot
################################################################################

################################################################################
# Lag of MAE from the log model
################################################################################
df_valid$mae <- abs(df_valid$gated_entries - df_valid$log_forecast)

df_valid <- df_valid %>%
  mutate(date = as.Date(time)) %>%
  group_by(date) %>%
  mutate(mae_day = mean(mae, na.rm = TRUE)) %>%
  ungroup()

df_valid <- df_valid %>%
  arrange(time) %>%
  mutate(mae_shifted = lag(mae_day, n = 7 * 48))


################################################################################
# Optimization for Log Model
################################################################################
# Don't want to include times where the system is closed
df_opt <- df_valid[df_valid$closed == 0, ]

# Current MBTA metrics observed
mbta_prop_full <- mean(df_opt$supply_gap_complex >= 0)
mbta_avg_supply_gap <- mean(df_opt$supply_gap_complex)
mbta_avg_interval_min <- mean((df_opt$line_length_sec / df_opt$traincnt) / 60)

# Selection lambda value, can be changed to whatever you want
lambda = 2

# Adjusting the demand and suppy gaps
df_opt$adjusted_demand <- lambda * (df_opt$log_forecast + df_opt$mae_shifted)
df_opt$forecast_traincnt <- ceiling(df_opt$adjusted_demand / df_opt$single_traincap)
df_opt$forecast_tot_linecap <- df_opt$forecast_traincnt * df_opt$single_traincap
df_opt$forecast_supply_gap_simple <- df_opt$gated_entries - df_opt$forecast_tot_linecap

df_opt <- df_opt %>%
  arrange(time) %>%
  mutate(
    lag_gap = lag(forecast_supply_gap_simple),
    forecast_prev_supply_gap_simple = ifelse(lag_gap > 0, lag_gap, 0)
  )
# Further adjusting the demand and suppy gaps
df_opt$forecast_gated_entires_complex <- df_opt$gated_entries + df_opt$forecast_prev_supply_gap_simple
df_opt$forecast_supply_gap_complex <- df_opt$forecast_gated_entires_complex - df_opt$forecast_tot_linecap
df_opt$train_interval_sec <- df_opt$line_length_sec / df_opt$traincnt
df_opt$forecast_train_interval_sec <- df_opt$line_length_sec / df_opt$forecast_traincnt

# Simple hisotgrams of the supply gaps
hist(df_opt$supply_gap_complex, breaks = 250,
     main = 'Supply Gap Distribution of MBTA Orange Line (Lambda = 2)',
     xlab = 'Supply Gap')
hist(df_opt$forecast_supply_gap_complex, breaks = 250,
     main = 'Supply Gap Distribution of Log Forecast Orange Line (Lambda = 2)',
     xlab = 'Supply_Gap')

# Choosing lambda values for the for loop below, and the steps between each value
# I found 0.01 seemed to be the best, any lambda past 4 also just doesn't make sense
lambda_values <- seq(0.01, 4, by = 0.01)

results <- data.frame(lambda = numeric(),
                      avg_supply_gap = numeric(),
                      pct_supply_gap_gt_0 = numeric())


# For loop to get the data for the lambda visual
for (lambda in lambda_values) {
  df_temp <- df_opt %>%
    mutate(
      adjusted_demand = lambda * (log_forecast + mae_shifted),
      forecast_traincnt = ceiling(adjusted_demand / single_traincap),
      forecast_tot_linecap = forecast_traincnt * single_traincap,
      forecast_supply_gap_simple = gated_entries - forecast_tot_linecap,
      lag_gap = lag(forecast_supply_gap_simple),
      forecast_prev_supply_gap_simple = ifelse(lag_gap > 0, lag_gap, 0),
      forecast_gated_entires_complex = gated_entries + forecast_prev_supply_gap_simple,
      forecast_supply_gap_complex = forecast_gated_entires_complex - forecast_tot_linecap,
      forecast_train_interval_sec = line_length_sec / forecast_traincnt
    )
  
  avg_gap <- mean(df_temp$forecast_supply_gap_complex, na.rm = TRUE)
  pct_gap_pos <- mean(df_temp$forecast_supply_gap_complex > 0, na.rm = TRUE) * 100
  avg_interval_min <- mean(df_temp$forecast_train_interval_sec[is.finite(df_temp$forecast_train_interval_sec)], na.rm = TRUE) / 60
  
  results <- rbind(results, data.frame(lambda = lambda,
                                       avg_supply_gap = avg_gap,
                                       pct_supply_gap_gt_0 = pct_gap_pos,
                                       avg_train_interval_min = avg_interval_min))
}

# Gap ranges, min and max used for the visual
gap_range <- range(results$avg_supply_gap, na.rm = TRUE)
gap_min <- gap_range[1]
gap_max <- gap_range[2]

# Scale factor again for the visual
scale_factor_pct <- (gap_max - gap_min) / 100

interval_range <- range(results$avg_train_interval_min, na.rm = TRUE)
interval_min <- interval_range[1]
interval_max <- interval_range[2]
scale_factor_interval <- (gap_max - gap_min) / (interval_max - interval_min)

results <- results %>%
  mutate(
    pct_gap_scaled = pct_supply_gap_gt_0 * scale_factor_pct + gap_min,
    interval_scaled = (avg_train_interval_min - interval_min) * scale_factor_interval + gap_min
  )

# Actual visualization of the data obtained above
ggplot(results, aes(x = lambda)) +
  geom_line(aes(y = avg_supply_gap, color = "Average Supply Gap"), size = 1.2) +
  geom_line(aes(y = pct_gap_scaled, color = "Percent Times Train at Capacity"), size = 1.2) +
  geom_line(aes(y = interval_scaled, color = "Train Frequency (Minutes)"), size = 1.2) +
  geom_hline(yintercept = mbta_avg_supply_gap, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = mbta_prop_full * scale_factor_pct + gap_min, linetype = "dashed", color = "red") +
  geom_hline(yintercept = mbta_avg_interval_min * scale_factor_interval + gap_min, linetype = "dashed", color = "green") +
  scale_y_continuous(
    name = "Avg Supply Gap",
    sec.axis = sec_axis(
      ~ (. - gap_min) / scale_factor_pct, 
      name = "Percent Gap > 0 (%) / Train Frequency (Minutes)"
    )
  ) +
  scale_color_manual(values = c(
    "Average Supply Gap" = "blue", 
    "Percent Times Train at Capacity" = "red", 
    "Train Frequency (Minutes)" = "green"
  )) +
  labs(x = "Lambda", color = "Metric", title = "Supply Gap Metrics vs Lambda (Log Model) Orange Line") +
  theme_minimal()



################################################################################
################################################################################
################################################################################


################################################################################
# Lag of MAE from the level model
################################################################################
df_valid$mae <- abs(df_valid$gated_entries - df_valid$level_forecast)

df_valid <- df_valid %>%
  mutate(date = as.Date(time)) %>%
  group_by(date) %>%
  mutate(mae_day = mean(mae, na.rm = TRUE)) %>%
  ungroup()

df_valid <- df_valid %>%
  arrange(time) %>%
  mutate(mae_shifted = lag(mae_day, n = 7 * 48))

################################################################################
# Optimization for Level Model
################################################################################
# Don't want to include times where the system is closed
df_opt <- df_valid[df_valid$closed == 0, ]

# Current MBTA metrics observed
mbta_prop_full <- mean(df_opt$supply_gap_complex >= 0)
mbta_avg_supply_gap <- mean(df_opt$supply_gap_complex)
mbta_avg_interval_min <- mean((df_opt$line_length_sec / df_opt$traincnt) / 60)

# Selection lambda value, can be changed to whatever you want
lambda = 2

# Adjusting the demand and suppy gaps
df_opt$adjusted_demand <- lambda * (df_opt$level_forecast + df_opt$mae_shifted)
df_opt$forecast_traincnt <- ceiling(df_opt$adjusted_demand / df_opt$single_traincap)
df_opt$forecast_tot_linecap <- df_opt$forecast_traincnt * df_opt$single_traincap
df_opt$forecast_supply_gap_simple <- df_opt$gated_entries - df_opt$forecast_tot_linecap

df_opt <- df_opt %>%
  arrange(time) %>%
  mutate(
    lag_gap = lag(forecast_supply_gap_simple),
    forecast_prev_supply_gap_simple = ifelse(lag_gap > 0, lag_gap, 0)
  )

# Further adjusting the demand and suppy gaps
df_opt$forecast_gated_entires_complex <- df_opt$gated_entries + df_opt$forecast_prev_supply_gap_simple
df_opt$forecast_supply_gap_complex <- df_opt$forecast_gated_entires_complex - df_opt$forecast_tot_linecap
df_opt$train_interval_sec <- df_opt$line_length_sec / df_opt$traincnt
df_opt$forecast_train_interval_sec <- df_opt$line_length_sec / df_opt$forecast_traincnt


# Simple hisotgrams of the supply gaps
hist(df_opt$supply_gap_complex, breaks = 250,
     main = 'Supply Gap Distribution of MBTA Orange Line (Lambda = 2)',
     xlab = 'Supply Gap')
hist(df_opt$forecast_supply_gap_complex, breaks = 250,
     main = 'Supply Gap Distribution of Level Forecast Orange Line (Lambda = 2)',
     xlab = 'Supply_Gap')

# Choosing lambda values for the for loop below, and the steps between each value
# I found 0.01 seemed to be the best, any lambda past 4 also just doesn't make sense
lambda_values <- seq(0.01, 4, by = 0.01)

results <- data.frame(lambda = numeric(),
                      avg_supply_gap = numeric(),
                      pct_supply_gap_gt_0 = numeric())

# For loop to get the data for the lambda visual
for (lambda in lambda_values) {
  df_temp <- df_opt %>%
    mutate(
      adjusted_demand = lambda * (level_forecast + mae_shifted),
      forecast_traincnt = ceiling(adjusted_demand / single_traincap),
      forecast_tot_linecap = forecast_traincnt * single_traincap,
      forecast_supply_gap_simple = gated_entries - forecast_tot_linecap,
      lag_gap = lag(forecast_supply_gap_simple),
      forecast_prev_supply_gap_simple = ifelse(lag_gap > 0, lag_gap, 0),
      forecast_gated_entires_complex = gated_entries + forecast_prev_supply_gap_simple,
      forecast_supply_gap_complex = forecast_gated_entires_complex - forecast_tot_linecap,
      forecast_train_interval_sec = line_length_sec / forecast_traincnt
    )
  
  avg_gap <- mean(df_temp$forecast_supply_gap_complex, na.rm = TRUE)
  pct_gap_pos <- mean(df_temp$forecast_supply_gap_complex > 0, na.rm = TRUE) * 100
  avg_interval_min <- mean(df_temp$forecast_train_interval_sec[is.finite(df_temp$forecast_train_interval_sec)], na.rm = TRUE) / 60
  
  results <- rbind(results, data.frame(lambda = lambda,
                                       avg_supply_gap = avg_gap,
                                       pct_supply_gap_gt_0 = pct_gap_pos,
                                       avg_train_interval_min = avg_interval_min))
}


# Gap ranges, min and max used for the visual
gap_range <- range(results$avg_supply_gap, na.rm = TRUE)
gap_min <- gap_range[1]
gap_max <- gap_range[2]

# Scale factor again used for the visual
scale_factor_pct <- (gap_max - gap_min) / 100  # since percent is 0 to 100


interval_range <- range(results$avg_train_interval_min, na.rm = TRUE)
interval_min <- interval_range[1]
interval_max <- interval_range[2]
scale_factor_interval <- (gap_max - gap_min) / (interval_max - interval_min)


results <- results %>%
  mutate(
    pct_gap_scaled = pct_supply_gap_gt_0 * scale_factor_pct + gap_min,
    interval_scaled = (avg_train_interval_min - interval_min) * scale_factor_interval + gap_min
  )


# Actual visualization of the data obtained above
ggplot(results, aes(x = lambda)) +
  geom_line(aes(y = avg_supply_gap, color = "Average Supply Gap"), size = 1.2) +
  geom_line(aes(y = pct_gap_scaled, color = "Percent Times Train at Capacity"), size = 1.2) +
  geom_line(aes(y = interval_scaled, color = "Train Frequency (Minutes)"), size = 1.2) +
  geom_hline(yintercept = mbta_avg_supply_gap, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = mbta_prop_full * scale_factor_pct + gap_min, linetype = "dashed", color = "red") +
  geom_hline(yintercept = mbta_avg_interval_min * scale_factor_interval + gap_min, linetype = "dashed", color = "green") +
  scale_y_continuous(
    name = "Avg Supply Gap",
    sec.axis = sec_axis(
      ~ (. - gap_min) / scale_factor_pct, 
      name = "Percent Gap > 0 (%) / Train Frequency (Minutes)"
    )
  ) +
  scale_color_manual(values = c(
    "Average Supply Gap" = "blue", 
    "Percent Times Train at Capacity" = "red", 
    "Train Frequency (Minutes)" = "green"
  )) +
  labs(x = "Lambda", color = "Metric", title = "Supply Gap Metrics vs Lambda (Level Model) Orange Line") +
  theme_minimal()

################################################################################
################################################################################
################################################################################





################################################################################
# The following code I think could have been its own script
# Because it used the daily aggregated data instead of the 30-minute frequency
# But to keep it simple I included it all in 1 script
# The following code is only used for visualizations for the paper
# None of it is for the models, forecasting, summary statistics, etc.
################################################################################
# Daily Time Series Analysis
################################################################################
# Data Processing
################################################################################
# Creating R date time
dfday$date <- as.POSIXct(dfday$date, format = "%Y-%m-%d")

# Training Data Preparation
################################################################################
# Index for start and end time
start_index <- which.min(abs(dfday$date - start_time))
end_index <- which.min(abs(dfday$date - end_time))

# Subsetting the df for the training time period
dfday_train <- dfday[start_index:end_index, ]

################################################################################
# Time Series Creation
################################################################################
# Creating the daily time series with weekly seasonality (frequency=7)
tsday <- ts(dfday$gated_entries, frequency = 365, start = c(2022, 1), end = c(2025, 10))
tsday <- window(tsday, start = c(2023, 1), end = c(2025, 1))
xtsday <- xts(dfday_train$gated_entries, dfday_train$date)

################################################################################
# Visualization
################################################################################
# Creating basic daily TS plot
ggplot(dfday_train, aes(x = date, y = gated_entries)) +
  geom_line(color = "darkorange", size = 0.7) +  # Line plot
  labs(title = "Daily Gated Entries Over Time (Orange Line)", 
       x = "Date", 
       y = "Daily Gated Entries")

################################################################################
# Time Series Decomposition
################################################################################
# Plotting the decomposition of the daily TS with weekly seasonality
decomposed_ts_day <- decompose(tsday)
plot(decomposed_ts_day)

################################################################################
################################################################################
################################################################################


################################################################################
# Average day Visualization
################################################################################
# Loading in the data from the GitHub repo
# This is created from the python code also on the repo
line <- "blue"  
line_url <- sprintf("https://raw.githubusercontent.com/zwinship/MBTA_Time_Series/refs/heads/main/data/processed/lines/%s_line.csv", line)
dfblue <- read.csv(line_url)
line <- "red"  
line_url <- sprintf("https://raw.githubusercontent.com/zwinship/MBTA_Time_Series/refs/heads/main/data/processed/lines/%s_line.csv", line)
dfred <- read.csv(line_url)


# Making sure time is in the right format
df$time <- as.POSIXct(df$time, format = "%Y-%m-%d %H:%M:%S")
dfblue$time <- as.POSIXct(dfblue$time, format = "%Y-%m-%d %H:%M:%S")
dfred$time <- as.POSIXct(dfred$time, format = "%Y-%m-%d %H:%M:%S")


# Creating the time of day
df$time_of_day <- format(df$time, format = "%H:%M:%S")
dfblue$time_of_day <- format(dfblue$time, format = "%H:%M:%S")
dfred$time_of_day <- format(dfred$time, format = "%H:%M:%S")

# Average gated entries by time for each line
avg_by_time_df <- df %>%
  group_by(time_of_day) %>%
  summarise(avg_value = mean(gated_entries, na.rm = TRUE)) %>%
  arrange(time_of_day)

avg_by_time_blue <- dfblue %>%
  group_by(time_of_day) %>%
  summarise(avg_value = mean(gated_entries, na.rm = TRUE)) %>%
  arrange(time_of_day)

avg_by_time_red <- dfred %>%
  group_by(time_of_day) %>%
  summarise(avg_value = mean(gated_entries, na.rm = TRUE)) %>%
  arrange(time_of_day)



# Graphing them all together
# My favorite plot for the whole paper
# Helps shows the variation that we actually need to capture in the study
# Unfortunately we did come up short, as seen in the paper

ggplot() +
  geom_line(data = avg_by_time_df, aes(x = time_of_day, y = avg_value, group = 1), color = "darkorange", size = 1) +
  geom_line(data = avg_by_time_blue, aes(x = time_of_day, y = avg_value, group = 1), color = "blue", size = 1) +
  geom_line(data = avg_by_time_red, aes(x = time_of_day, y = avg_value, group = 1), color = "red", size = 1) +
  labs(title = "Average Gated Entries by Time of Day for All Heavy Rail Line (Colored By Line)",
       x = "Time of Day",
       y = "Average Gated Entries") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5)) +  # Rotate labels and adjust size
  scale_x_discrete(breaks = function(x) x[seq(1, length(x), by = 5)]) +  # Skip every 10th value on x-axis
  theme_minimal()



