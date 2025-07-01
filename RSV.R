
library(readr)
library(data.table)
library(dplyr)
library(lubridate)
library(tidyr)
library(caret)
library(doParallel)
library(randomForest)
library(rBayesianOptimization)
library(Metrics)  # For evaluation metrics



## --- Import MUSC Weekly Influenza Data ---
MUSC_Weekly_Influenza_County <- read_csv(
   '/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/MUSC/Infectious Disease EHR/Weekly Data/RSV/06_26_2025/MUSC_Weekly_RSV_County_06_26_2025.csv'
   )

## --- Import Prisma Health Weekly Influenza Data ---
Prisma_Weekly_Influenza_County <- read_csv(
   '/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/Prisma Health/Infectious Disease EHR/Weekly Data/RSV/06_26_2025/Prisma_Health_RSV_County.csv'
   )

Healthsystem_Weekly_Influenza_County <- full_join(
   MUSC_Weekly_Influenza_County %>% select(County, Week, Weekly_Encounters) %>% 
      rename(MUSC_Encounters = Weekly_Encounters),
   Prisma_Weekly_Influenza_County %>% select(County, Week, Weekly_Encounters) %>% 
      rename(Prisma_Encounters = Weekly_Encounters),
   by = c("County", "Week")
) %>%
   mutate(
      Weekly_Encounters = rowSums(across(c(MUSC_Encounters, Prisma_Encounters)), na.rm = TRUE)
   ) %>%
   select(County, Week, Weekly_Encounters)




#RFA -------------
hosp_data_rfa <- read_csv('/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/SC Health Records/RFA/Years 2019-2023/Clean Data/data_rsv.csv',
                          col_select = c(ZIP, COUNTY,
                                         DISYEAR, ADMYEAR, DISMTH, DISDAY, ADMMTH,
                                         ADMDAY, RFA_ID, ADMD_new, DISD_new, COUNTY))

##Remove duplicate (the difference is less than 14 days)
hosp_data_rfa <- hosp_data_rfa %>%
   distinct(RFA_ID, ADMD_new, .keep_all = TRUE)





# Find rows where RFA_ID is similar and difference in ADMD_new is less than 2 weeks
hosp_data_rfa <- as.data.table(hosp_data_rfa) # Convert to data.table
hosp_data_rfa$ADMD_new <- as.Date(hosp_data_rfa$ADMD_new, format = "%m/%d/%Y")# Ensure ADMD_new is in Date format
hosp_data_rfa <- hosp_data_rfa[order(RFA_ID, ADMD_new)] # Sort data by RFA_ID and ADMD_new
# Calculate the difference between consecutive admission dates for each RFA_ID
hosp_data_rfa[, admission_diff := ADMD_new - shift(ADMD_new), by = RFA_ID]
# Filter rows where the difference is less than 14 days and the difference is not NA
close_admissions <- hosp_data_rfa[admission_diff < 14 & !is.na(admission_diff)]

# Identify rows to remove based on RFA_ID and ADMD_new
rows_to_remove <- hosp_data_rfa[close_admissions, on = .(RFA_ID, ADMD_new), which = TRUE]
hosp_data_rfa <- hosp_data_rfa[-rows_to_remove]# Remove these rows from hosp_data_rfa_dt


RFA <- hosp_data_rfa %>%
   dplyr::select(RFA_ID, COUNTY, ADMD_new)

# Create RFA$Week based on RFA$ADMD_new 
RFA$Week <- floor_date(RFA$ADMD_new, unit = "week", week_start = 1)

# Filter rows where EncounterDate is in the years 2022, or 2023
RFA <- RFA[format(RFA$Week, "%Y") %in% c("2021", "2022", "2023"), ]

# Define a named vector mapping county codes to county names
county_mapping <- c(
   "1" = "Abbeville", "2" = "Aiken", "3" = "Allendale", "4" = "Anderson", "5" = "Bamberg",
   "6" = "Barnwell", "7" = "Beaufort", "8" = "Berkeley", "9" = "Calhoun", "10" = "Charleston",
   "11" = "Cherokee", "12" = "Chester", "13" = "Chesterfield", "14" = "Clarendon", "15" = "Colleton",
   "16" = "Darlington", "17" = "Dillon", "18" = "Dorchester", "19" = "Edgefield", "20" = "Fairfield",
   "21" = "Florence", "22" = "Georgetown", "23" = "Greenville", "24" = "Greenwood", "25" = "Hampton",
   "26" = "Horry", "27" = "Jasper", "28" = "Kershaw", "29" = "Lancaster", "30" = "Laurens",
   "31" = "Lee", "32" = "Lexington", "33" = "McCormick", "34" = "Marion", "35" = "Marlboro",
   "36" = "Newberry", "37" = "Oconee", "38" = "Orangeburg", "39" = "Pickens", "40" = "Richland",
   "41" = "Saluda", "42" = "Spartanburg", "43" = "Sumter", "44" = "Union", "45" = "Williamsburg", "46" = "York"
)

# Convert County_Code to character and replace County_Code with County_Name using the mapping
RFA$COUNTY <- county_mapping[as.character(RFA$COUNTY)]
RFA <- RFA %>% rename(County = COUNTY)
RFA <- RFA[!is.na(County)]




#keep only the counties where the total rows for 2023 exceed 100
# Filter rows for 2023 and count the number of rows for each county
county_counts_2023 <- RFA[format(ADMD_new, "%Y") == "2023", .N, by = County]

# Find counties where the count is greater than 500
counties_to_keep <- county_counts_2023[N > 100, County]

# Filter the RFA dataset to keep only rows for these counties
RFA <- RFA[County %in% counties_to_keep]




#compare the distributed counts with actual counts
# Calculate Weekly Row Count for Each County
County_weekly_counts <- RFA[, .N, by = .(Week, County)]
setnames(County_weekly_counts, "N", "Actual_Count")


# Define all possible combinations
all_combinations <- expand.grid(
   County = unique(RFA$County),
   Week = unique(RFA$Week)
)

# Left join to fill missing combinations
County_weekly_counts <- left_join(all_combinations, County_weekly_counts, by = c("County", "Week"))

# Replace NA values in occurrences with 0
County_weekly_counts <- County_weekly_counts %>%
   mutate(Actual_Count = replace_na(Actual_Count, 0))



# Calculate Weekly Total Row Count Across Counties
County_weekly_counts <- as.data.table(County_weekly_counts)
weekly_totals <- County_weekly_counts[, .(Weekly_Total = sum(Actual_Count)), by = Week]

# Merge weekly totals back to the County-week data
County_weekly_counts <- merge(County_weekly_counts, weekly_totals, by = "Week")






# Moving average smoother
smooth_fun <- function(data, window_size = 2) {
   n <- length(data)
   result <- numeric(n)
   
   for (i in 1:n) {
      start <- max(1, i - window_size)
      end <- min(n, i + window_size)
      result[i] <- mean(data[start:end])
   }
   
   return(result)
}

# Smoothing between peaks (preserve peaks)
smooth_between_peaks <- function(y, peaks, window_size = 2) {
   n <- length(y)
   smoothed_y <- y
   peaks <- sort(unique(c(1, peaks, n)))
   
   for (i in seq_along(peaks[-1])) {
      start <- peaks[i]
      end <- peaks[i + 1]
      
      if ((end - start + 1) > 3) {
         segment <- y[start:end]
         segment_s <- smooth_fun(segment, window_size)
         segment[2:(length(segment) - 1)] <- segment_s[2:(length(segment_s) - 1)]
         smoothed_y[start:end] <- segment
      }
   }
   
   smoothed_y[peaks] <- y[peaks]
   return(smoothed_y)
}

# Wrapper function
process_smooth <- function(y, distance = 8, prominence_percentile = 0.25, window_size = 2) {
   y <- as.numeric(y)
   
   if (all(is.na(y)) || length(na.omit(y)) < 5) {
      return(rep(NA, length(y)))
   }
   
   threshold <- quantile(y, prominence_percentile, na.rm = TRUE)
   peaks_matrix <- tryCatch({
      pracma::findpeaks(y, minpeakdistance = distance, minpeakheight = threshold)
   }, error = function(e) NULL)
   
   if (is.null(peaks_matrix) || nrow(peaks_matrix) == 0) {
      return(smooth_fun(y, window_size))
   } else {
      peaks <- peaks_matrix[, 2]
      y_smoothed <- smooth_between_peaks(y, peaks, window_size)
      return(y_smoothed[1:length(y)])
   }
}




County_weekly_counts <- County_weekly_counts %>%
   group_by(County) %>%
   arrange(Week) %>%
   mutate(Smoothed_Actual_Count = process_smooth(Actual_Count)) %>%
   ungroup()


ggplot(filter(County_weekly_counts, County == "Greenville"), aes(x = Week)) +
   geom_line(aes(y = Actual_Count), color = "gray") +
   geom_line(aes(y = Smoothed_Actual_Count), color = "blue") +
   labs(title = "Actual vs Smoothed Weekly Count - Greenville",
        y = "Encounters", x = "Week")










#Incorporate new predicted data--------------------------
State_prediction <- read_csv('/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Tanvir/Flu/Distribution/Data from Dr. Minjae/04:08:25/covid-19_state_data_04_14_25.csv')
State_prediction <- State_prediction %>%
   mutate(Date = mdy(Date))


Counties <- counties_to_keep


State_prediction_24_25 <- State_prediction %>%
   filter(lubridate::year(Date) %in% c(2024, 2025)) %>%
   select(Week=Date, Weekly_Total=`Projected Cases(post training)`)


# Expand the dataset and add Actual_Count column with NA values
State_prediction_24_25 <- State_prediction_24_25 %>%
   mutate(County = list(Counties)) %>%  # Create a list-column of Counties
   unnest(County) #%>%  # Expand each week into four rows, one per County
#mutate(Actual_Count = 4)  # Add Actual_Count with NA values

#State_prediction_24_25$Actual_Count<-round(State_prediction_24_25$Actual_Count)

# Combine the datasets row by row
County_weekly_counts <- bind_rows(County_weekly_counts, State_prediction_24_25)




# creating variables ----------
# Create Year and Week_No variables
County_weekly_counts <- County_weekly_counts %>%
   mutate(
      Year = year(Week),
      Week_No = week(Week)
   )



# Ensure County columns are the same type (character recommended)
County_weekly_counts <- County_weekly_counts %>%
   mutate(County = as.character(County))

Healthsystem_Weekly_Influenza_County <- Healthsystem_Weekly_Influenza_County %>%
   mutate(County = as.character(County))

# Left join on County and Week
County_weekly_counts <- County_weekly_counts %>%
   left_join(Healthsystem_Weekly_Influenza_County, by = c("County", "Week"))



# Convert County, Week_No, and Year to factors
County_weekly_counts <- County_weekly_counts %>%
   mutate(
      County = as.factor(County),
      Week_No = as.factor(Week_No),
      Year = as.factor(Year)
      
   )







County_weekly_counts$Smoothed_Actual_Count <- round(County_weekly_counts$Smoothed_Actual_Count)

# Ensure County_weekly_counts is a data.table
County_weekly_counts <- data.table(County_weekly_counts)

# Split datasets
County_weekly_counts_train <- County_weekly_counts[Week >= as.Date("2020-01-01") & Week < as.Date("2023-01-01")]
County_weekly_counts_valid <- County_weekly_counts[Week >= as.Date("2023-01-01") & Week < as.Date("2023-07-01")]
County_weekly_counts_test  <- County_weekly_counts[Week >= as.Date("2023-07-01") & Week <= as.Date("2023-12-25")]
County_weekly_counts_predict <- County_weekly_counts[Week >= as.Date("2024-01-01") & Week <= as.Date("2025-04-14")]

# Convert categorical variables to numeric (County must be one-hot encoded or label encoded)
County_levels <- unique(County_weekly_counts$County)
County_weekly_counts_train[, County := as.integer(factor(County, levels = County_levels))]
County_weekly_counts_valid[, County := as.integer(factor(County, levels = County_levels))]
County_weekly_counts_test[, County := as.integer(factor(County, levels = County_levels))]
County_weekly_counts_predict[, County := as.integer(factor(County, levels = County_levels))]


# One-hot encode all predictors including County
train_matrix <- model.matrix(Smoothed_Actual_Count ~ Weekly_Total + County + Week_No + Weekly_Encounters, 
                             data = County_weekly_counts_train)[, -1]

valid_matrix <- model.matrix(Smoothed_Actual_Count ~ Weekly_Total + County + Week_No + Weekly_Encounters, 
                             data = County_weekly_counts_valid)[, -1]

# One-hot encode test dataset
test_matrix <- model.matrix(~ Weekly_Total + County + Week_No + Weekly_Encounters, 
                            data = County_weekly_counts_test)[, -1]

predict_matrix <- model.matrix(~ Weekly_Total + County + Week_No + Weekly_Encounters, 
                               data = County_weekly_counts_predict)[, -1]


# Create DMatrix objects
dtrain <- xgb.DMatrix(data = train_matrix,
                      label = County_weekly_counts_train$Smoothed_Actual_Count)

dvalid <- xgb.DMatrix(data = valid_matrix,
                      label = County_weekly_counts_valid$Smoothed_Actual_Count)

dtest <- xgb.DMatrix(data = test_matrix,
                     label = County_weekly_counts_test$Smoothed_Actual_Count)

dpredict <- xgb.DMatrix(data = predict_matrix)  # ✅ no label




# Parallel setup
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Objective function
objective_function <- function(eta, max_depth, subsample, colsample_bytree, nrounds) {
   param <- list(
      objective = "reg:squarederror",
      eta = eta,
      max_depth = as.integer(round(max_depth)),
      subsample = subsample,
      colsample_bytree = colsample_bytree,
      verbosity = 0
   )
   
   model <- xgb.train(
      params = param,
      data = dtrain,
      nrounds = as.integer(round(nrounds)),
      watchlist = list(valid = dvalid),
      early_stopping_rounds = 10,   # stop if no improvement for 10 rounds
      eval_metric = "rmse",
      verbose = 0
   )
   
   preds <- predict(model, dvalid)
   rmse_score <- RMSE(preds, County_weekly_counts_valid$Smoothed_Actual_Count)
   
   list(Score = -rmse_score)
}

# Bayesian Optimization
set.seed(1234)
opt_results <- bayesOpt(
   FUN = objective_function,
   bounds = list(
      eta = c(0.01, 0.3),
      max_depth = c(2L, 15L),
      subsample = c(0.5, 1.0),
      colsample_bytree = c(0.5, 1.0),
      nrounds = c(50L, 500L)
   ),
   initPoints = 10,
   iters.n = 40,
   acq = "ucb",
   kappa = 2.576,
   verbose = 0
)

# Best parameters
best_params <- getBestPars(opt_results)

# Final model
final_model <- xgb.train(
   params = list(
      objective = "reg:squarederror",
      eta = best_params$eta,
      max_depth = as.integer(round(best_params$max_depth)),
      subsample = best_params$subsample,
      colsample_bytree = best_params$colsample_bytree,
      verbosity = 0
   ),
   data = dtrain,
   nrounds = as.integer(round(best_params$nrounds)),
   eval_metric = "rmse",
   watchlist = list(valid = dvalid),
   early_stopping_rounds = 10,
   verbose = 0
)

# Stop parallel cluster
stopCluster(cl)


# # Predict using XGBoost model
# dtest <- xgb.DMatrix(data = test_matrix,
#                      label = County_weekly_counts_test$Smoothed_Actual_Count)

County_weekly_counts_test$Predicted_Count <- predict(final_model, newdata = dtest)
County_weekly_counts_test$Predicted_Count <- pmax(0, County_weekly_counts_test$Predicted_Count)

County_weekly_counts_test <- as.data.table(County_weekly_counts_test)

# Now you can apply the weighting and distribution
County_weekly_counts_test[, Weight := Predicted_Count / sum(Predicted_Count), by = Week]
County_weekly_counts_test[, Distributed_Count := Weekly_Total * Weight]


County_weekly_counts_test <- County_weekly_counts_test %>%
   mutate(County = as.character(County)) %>%  # if needed
   group_by(County) %>%
   arrange(Week) %>%
   mutate(Smoothed_Distributed_Count = process_smooth(Distributed_Count)) %>%
   ungroup()

County_weekly_counts_test$Smoothed_Distributed_Count <- as.numeric(County_weekly_counts_test$Smoothed_Distributed_Count)
County_weekly_counts_test$Smoothed_Distributed_Count <- round(County_weekly_counts_test$Smoothed_Distributed_Count)



County_weekly_counts_test$Ai <- ifelse(
   County_weekly_counts_test$Smoothed_Actual_Count == 0 & County_weekly_counts_test$Smoothed_Distributed_Count == 0,
   1,
   pmin(County_weekly_counts_test$Smoothed_Actual_Count, County_weekly_counts_test$Smoothed_Distributed_Count) /
      pmax(County_weekly_counts_test$Smoothed_Actual_Count, County_weekly_counts_test$Smoothed_Distributed_Count)
)





mean_ai_by_County <- County_weekly_counts_test %>%
   group_by(County) %>%
   summarise(Mean_Ai = mean(Ai, na.rm = TRUE))

print(summary(mean_ai_by_County$Mean_Ai))
print(summary(County_weekly_counts_test$Ai))







# Predict on validation set
County_weekly_counts_valid$Predicted_Count <- predict(final_model, newdata = dvalid)

County_weekly_counts_valid[, Weight := Predicted_Count / sum(Predicted_Count), by = Week]
County_weekly_counts_valid[, Distributed_Count := Weekly_Total * Weight]
County_weekly_counts_valid$Distributed_Count <- round(County_weekly_counts_valid$Distributed_Count)

County_weekly_counts_valid$Ai <- ifelse(
   County_weekly_counts_valid$Smoothed_Actual_Count == 0 & County_weekly_counts_valid$Distributed_Count == 0,
   1,
   pmin(County_weekly_counts_valid$Smoothed_Actual_Count, County_weekly_counts_valid$Distributed_Count) /
      pmax(County_weekly_counts_valid$Smoothed_Actual_Count, County_weekly_counts_valid$Distributed_Count)
)

mean_ai_by_County <- County_weekly_counts_valid %>%
   group_by(County) %>%
   summarise(Mean_Ai = mean(Ai, na.rm = TRUE))

print(summary(mean_ai_by_County$Mean_Ai))
print(summary(County_weekly_counts_valid$Ai))




# Final predictions
County_weekly_counts_predict$Predicted_Count <- predict(final_model, newdata = dpredict)

County_weekly_counts_predict[, Weight := Predicted_Count / sum(Predicted_Count), by = Week]
County_weekly_counts_predict[, Distributed_Count := Weekly_Total * Weight]
County_weekly_counts_predict$Distributed_Count <- round(County_weekly_counts_predict$Distributed_Count)




# Convert County to integer
County_weekly_counts_train <- County_weekly_counts_train %>%
   mutate(County = as.integer(County))

County_weekly_counts_test <- County_weekly_counts_test %>%
   mutate(County = as.integer(County))

County_weekly_counts_valid <- County_weekly_counts_valid %>%
   mutate(County = as.integer(County))

# Combine datasets
df_combined <- bind_rows(
   County_weekly_counts_train,
   County_weekly_counts_test,
   County_weekly_counts_valid
)

# Retrieve original county names from training data
County_levels <- unique(County_weekly_counts$County)  # original names in correct order

# Re-assign factor with correct labels
df_combined$County <- factor(
   df_combined$County,
   levels = 1:length(County_levels),
   labels = County_levels
)

# Alphabetize County levels
df_combined$County <- factor(df_combined$County, levels = sort(unique(as.character(df_combined$County))))

# Define test start date
test_start_week <- as.Date("2023-07-01")

# Plot
ggplot(df_combined, aes(x = Week)) +
   geom_line(aes(y = Distributed_Count, color = "Distributed Encounters"), size = 0.5) +
   geom_line(aes(y = Smoothed_Actual_Count, color = "Smoothed Observed Encounters"), size = 0.5) +
   geom_line(aes(y = Actual_Count, color = "Observed Encounters"), size = 0.5) +
   geom_vline(xintercept = as.numeric(test_start_week), linetype = "dotted", color = "black", size = 0.7) +
   facet_wrap(~ County, scales = "free") +
   labs(
      title = "RSV: Observed, Smoothed Observed, and Distributed Encounters Over Weeks",
      x = "Week",
      y = "Count",
      color = ""
   ) +
   theme_minimal() +
   scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
   scale_color_manual(values = c(
      "Distributed Encounters" = "blue",
      "Observed Encounters" = "red",
      "Smoothed Observed Encounters" = "black"
   )) +
   theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      strip.text = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 12),
      legend.position = c(1, 0),
      legend.justification = c(1, 0)
   ) 



selected_counties <- c("Aiken", "Charleston", "Horry", "Greenville", "Richland", "Spartanburg")

ggplot(df_combined %>% filter(County %in% selected_counties), aes(x = Week)) +
   geom_line(aes(y = Distributed_Count, color = "Distributed Encounters"), size = .5) +
   geom_line(aes(y = Smoothed_Actual_Count, color = "Smoothed Observed Encounters"), size = .5) +
   geom_line(aes(y = Actual_Count, color = "Observed Encounters"), size = .5) +
   geom_vline(xintercept = as.numeric(test_start_week), linetype = "dotted", color = "black", size = 0.7) +
   facet_wrap(~ County, scales = "free") +
   labs(
      title = "Flu: Observed vs. Distributed Encounters (Selected Counties)",
      x = "Week",
      y = "Count",
      color = ""
   ) +
   theme_minimal() +
   scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
   scale_color_manual(values = c("Distributed Encounters" = "blue",
                                 "Observed Encounters" = "red",
                                 "Smoothed Observed Encounters" = "black")) +
   theme(
      axis.text.x = element_text(angle = 60, hjust = 1), 
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom"
   )
