library(readr)
library(data.table)
library(dplyr)
library(lubridate)
library(tidyr)
library(caret)
library(doParallel)
library(randomForest)
library(rBayesianOptimization)
library(Metrics)  



## --- Import MUSC Weekly Influenza Data ---
MUSC_Weekly_Influenza_Region <- read_csv(
   '/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/MUSC/Infectious Disease EHR/Weekly Data/Influenza/06_26_2025/MUSC_Weekly_Influenza_Region_06_26_2025.csv'
   )

## --- Import Prisma Health Weekly Influenza Data ---
Prisma_Weekly_Influenza_Region <- read_csv(
   '/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/Prisma Health/Infectious Disease EHR/Weekly Data/Influenza/06_26_2025/Prisma_Health_Influenza_Region.csv'
   )

Healthsystem_Weekly_Influenza_Region <- full_join(
   MUSC_Weekly_Influenza_Region %>% select(Region, Week, Weekly_Encounters) %>% 
      rename(MUSC_Encounters = Weekly_Encounters),
   Prisma_Weekly_Influenza_Region %>% select(Region, Week, Weekly_Encounters) %>% 
      rename(Prisma_Encounters = Weekly_Encounters),
   by = c("Region", "Week")
) %>%
   mutate(
      Weekly_Encounters = rowSums(across(c(MUSC_Encounters, Prisma_Encounters)), na.rm = TRUE)
   ) %>%
   select(Region, Week, Weekly_Encounters)




#RFA ------------
hosp_data_rfa <- read_csv('/Users/tanvirahammed/Library/CloudStorage/Box-Box/BoxPHI-PHMR Projects/Data/SC Health Records/RFA/Years 2019-2023/Clean Data/data_influenza.csv',
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





# Create a new column 'Region' based on 'County'
RFA <- RFA %>%
   mutate(Region = case_when(
      County %in% c("Abbeville", "Anderson", "Cherokee", "Greenville", "Greenwood",
                    "Laurens", "McCormick", "Oconee", "Pickens", "Spartanburg", "Union") ~ "Upstate",
      County %in% c("Chesterfield", "Clarendon", "Darlington", "Dillon", "Florence",
                    "Georgetown", "Horry", "Lee", "Marion", "Marlboro", "Sumter", "Williamsburg") ~ "Pee Dee",
      County %in% c("Aiken", "Barnwell", "Chester", "Edgefield", "Fairfield", "Kershaw",
                    "Lancaster", "Lexington", "Newberry", "Richland", "Saluda", "York") ~ "Midlands",
      County %in% c("Allendale", "Bamberg", "Beaufort", "Berkeley", "Calhoun", "Charleston",
                    "Colleton", "Dorchester", "Hampton", "Jasper", "Orangeburg") ~ "Lowcountry",
      TRUE ~ "Other"  # For counties not in the specified lists (if any)
   ))








# Filter rows for 2023 and count the number of rows for each county
county_counts_2023 <- RFA[format(ADMD_new, "%Y") == "2023", .N, by = Region]

# Find counties where the count is greater than 500
counties_to_keep <- county_counts_2023[N > 10, Region]

# Filter the RFA dataset to keep only rows for these counties
RFA <- RFA[Region %in% counties_to_keep]




#compare the distributed counts with actual counts
# Calculate Weekly Row Count for Each Region
Region_weekly_counts <- RFA[, .N, by = .(Week, Region)]
setnames(Region_weekly_counts, "N", "Actual_Count")


# Define all possible combinations
all_combinations <- expand.grid(
   Region = unique(RFA$Region),
   Week = unique(RFA$Week)
)

# Left join to fill missing combinations
Region_weekly_counts <- left_join(all_combinations, Region_weekly_counts, by = c("Region", "Week"))

# Replace NA values in occurrences with 0
Region_weekly_counts <- Region_weekly_counts %>%
   mutate(Actual_Count = replace_na(Actual_Count, 0))



# Calculate Weekly Total Row Count Across Counties
Region_weekly_counts <- as.data.table(Region_weekly_counts)
weekly_totals <- Region_weekly_counts[, .(Weekly_Total = sum(Actual_Count)), by = Week]

# Merge weekly totals back to the Region-week data
Region_weekly_counts <- merge(Region_weekly_counts, weekly_totals, by = "Week")






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




Region_weekly_counts <- Region_weekly_counts %>%
   group_by(Region) %>%
   arrange(Week) %>%
   mutate(Smoothed_Actual_Count = process_smooth(Actual_Count)) %>%
   ungroup()


ggplot(filter(Region_weekly_counts, Region == "Greenville"), aes(x = Week)) +
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
   mutate(Region = list(Counties)) %>%  # Create a list-column of Counties
   unnest(Region) #%>%  # Expand each week into four rows, one per Region
#mutate(Actual_Count = 4)  # Add Actual_Count with NA values

#State_prediction_24_25$Actual_Count<-round(State_prediction_24_25$Actual_Count)

# Combine the datasets row by row
Region_weekly_counts <- bind_rows(Region_weekly_counts, State_prediction_24_25)




# creating variables ----------
# Create Year and Week_No variables
Region_weekly_counts <- Region_weekly_counts %>%
   mutate(
      Year = year(Week),
      Week_No = week(Week)
   )



# Ensure Region columns are the same type (character recommended)
Region_weekly_counts <- Region_weekly_counts %>%
   mutate(Region = as.character(Region))

Healthsystem_Weekly_Influenza_Region <- Healthsystem_Weekly_Influenza_Region %>%
   mutate(Region = as.character(Region))

# Left join on Region and Week
Region_weekly_counts <- Region_weekly_counts %>%
   left_join(Healthsystem_Weekly_Influenza_Region, by = c("Region", "Week"))



# Convert Region, Week_No, and Year to factors
Region_weekly_counts <- Region_weekly_counts %>%
   mutate(
      Region = as.factor(Region),
      Week_No = as.factor(Week_No),
      Year = as.factor(Year)
      
   )







Region_weekly_counts$Smoothed_Actual_Count <- round(Region_weekly_counts$Smoothed_Actual_Count)

# Ensure Region_weekly_counts is a data.table
Region_weekly_counts <- data.table(Region_weekly_counts)

# Split datasets
Region_weekly_counts_train <- Region_weekly_counts[Week >= as.Date("2020-01-01") & Week < as.Date("2023-01-01")]
Region_weekly_counts_valid <- Region_weekly_counts[Week >= as.Date("2023-01-01") & Week < as.Date("2023-07-01")]
Region_weekly_counts_test  <- Region_weekly_counts[Week >= as.Date("2023-07-01") & Week <= as.Date("2023-12-25")]
Region_weekly_counts_predict <- Region_weekly_counts[Week >= as.Date("2024-01-01") & Week <= as.Date("2025-04-14")]

# Convert categorical variables to numeric (Region must be one-hot encoded or label encoded)
Region_levels <- unique(Region_weekly_counts$Region)
Region_weekly_counts_train[, Region := as.integer(factor(Region, levels = Region_levels))]
Region_weekly_counts_valid[, Region := as.integer(factor(Region, levels = Region_levels))]
Region_weekly_counts_test[, Region := as.integer(factor(Region, levels = Region_levels))]
Region_weekly_counts_predict[, Region := as.integer(factor(Region, levels = Region_levels))]


# One-hot encode all predictors including Region
train_matrix <- model.matrix(Smoothed_Actual_Count ~ Weekly_Total + Region + Week_No + Weekly_Encounters, 
                             data = Region_weekly_counts_train)[, -1]

valid_matrix <- model.matrix(Smoothed_Actual_Count ~ Weekly_Total + Region + Week_No + Weekly_Encounters, 
                             data = Region_weekly_counts_valid)[, -1]

# One-hot encode test dataset
test_matrix <- model.matrix(~ Weekly_Total + Region + Week_No + Weekly_Encounters, 
                            data = Region_weekly_counts_test)[, -1]

predict_matrix <- model.matrix(~ Weekly_Total + Region + Week_No + Weekly_Encounters, 
                               data = Region_weekly_counts_predict)[, -1]


# Create DMatrix objects
dtrain <- xgb.DMatrix(data = train_matrix,
                      label = Region_weekly_counts_train$Smoothed_Actual_Count)

dvalid <- xgb.DMatrix(data = valid_matrix,
                      label = Region_weekly_counts_valid$Smoothed_Actual_Count)

dtest <- xgb.DMatrix(data = test_matrix,
                     label = Region_weekly_counts_test$Smoothed_Actual_Count)

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
   rmse_score <- RMSE(preds, Region_weekly_counts_valid$Smoothed_Actual_Count)
   
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
#                      label = Region_weekly_counts_test$Smoothed_Actual_Count)

Region_weekly_counts_test$Predicted_Count <- predict(final_model, newdata = dtest)
Region_weekly_counts_test$Predicted_Count <- pmax(0, Region_weekly_counts_test$Predicted_Count)

Region_weekly_counts_test <- as.data.table(Region_weekly_counts_test)

# Now you can apply the weighting and distribution
Region_weekly_counts_test[, Weight := Predicted_Count / sum(Predicted_Count), by = Week]
Region_weekly_counts_test[, Distributed_Count := Weekly_Total * Weight]


Region_weekly_counts_test <- Region_weekly_counts_test %>%
   mutate(Region = as.character(Region)) %>%  # if needed
   group_by(Region) %>%
   arrange(Week) %>%
   mutate(Smoothed_Distributed_Count = process_smooth(Distributed_Count)) %>%
   ungroup()

Region_weekly_counts_test$Smoothed_Distributed_Count <- as.numeric(Region_weekly_counts_test$Smoothed_Distributed_Count)
Region_weekly_counts_test$Smoothed_Distributed_Count <- round(Region_weekly_counts_test$Smoothed_Distributed_Count)



Region_weekly_counts_test$Ai <- ifelse(
   Region_weekly_counts_test$Smoothed_Actual_Count == 0 & Region_weekly_counts_test$Smoothed_Distributed_Count == 0,
   1,
   pmin(Region_weekly_counts_test$Smoothed_Actual_Count, Region_weekly_counts_test$Smoothed_Distributed_Count) /
      pmax(Region_weekly_counts_test$Smoothed_Actual_Count, Region_weekly_counts_test$Smoothed_Distributed_Count)
)





mean_ai_by_Region <- Region_weekly_counts_test %>%
   group_by(Region) %>%
   summarise(Mean_Ai = mean(Ai, na.rm = TRUE))

print(summary(mean_ai_by_Region$Mean_Ai))
print(summary(Region_weekly_counts_test$Ai))







# Predict on validation set
Region_weekly_counts_valid$Predicted_Count <- predict(final_model, newdata = dvalid)

Region_weekly_counts_valid[, Weight := Predicted_Count / sum(Predicted_Count), by = Week]
Region_weekly_counts_valid[, Distributed_Count := Weekly_Total * Weight]
Region_weekly_counts_valid$Distributed_Count <- round(Region_weekly_counts_valid$Distributed_Count)

Region_weekly_counts_valid$Ai <- ifelse(
   Region_weekly_counts_valid$Smoothed_Actual_Count == 0 & Region_weekly_counts_valid$Distributed_Count == 0,
   1,
   pmin(Region_weekly_counts_valid$Smoothed_Actual_Count, Region_weekly_counts_valid$Distributed_Count) /
      pmax(Region_weekly_counts_valid$Smoothed_Actual_Count, Region_weekly_counts_valid$Distributed_Count)
)

mean_ai_by_Region <- Region_weekly_counts_valid %>%
   group_by(Region) %>%
   summarise(Mean_Ai = mean(Ai, na.rm = TRUE))

print(summary(mean_ai_by_Region$Mean_Ai))
print(summary(Region_weekly_counts_valid$Ai))




# Final predictions
Region_weekly_counts_predict$Predicted_Count <- predict(final_model, newdata = dpredict)

Region_weekly_counts_predict[, Weight := Predicted_Count / sum(Predicted_Count), by = Week]
Region_weekly_counts_predict[, Distributed_Count := Weekly_Total * Weight]
Region_weekly_counts_predict$Distributed_Count <- round(Region_weekly_counts_predict$Distributed_Count)




Region_weekly_counts_train <- Region_weekly_counts_train %>%
   mutate(Region = as.integer(Region))

Region_weekly_counts_test <- Region_weekly_counts_test %>%
   mutate(Region = as.integer(Region))

Region_weekly_counts_valid <- Region_weekly_counts_valid %>%
   mutate(Region = as.integer(Region))

# Combine datasets
df_combined <- bind_rows(
   Region_weekly_counts_train,
   Region_weekly_counts_test,
   Region_weekly_counts_valid
)


# Retrieve original county names from training data
Region_levels <- unique(Region_weekly_counts$Region)  # original names in correct order

# Re-assign factor with correct labels
df_combined$Region <- factor(
   df_combined$Region,
   levels = 1:length(Region_levels),
   labels = Region_levels
)

# Alphabetize County levels
df_combined$Region <- factor(df_combined$Region, levels = sort(unique(as.character(df_combined$Region))))

# Define test start date
test_start_week <- as.Date("2023-07-01")

# Plot
ggplot(df_combined, aes(x = Week)) +
   geom_line(aes(y = Distributed_Count, color = "Distributed Encounters"), size = 0.5) +
   geom_line(aes(y = Smoothed_Actual_Count, color = "Smoothed Observed Encounters"), size = 0.5) +
   geom_line(aes(y = Actual_Count, color = "Observed Encounters"), size = 0.5) +
   geom_vline(xintercept = as.numeric(test_start_week), linetype = "dotted", color = "black", size = 0.7) +
   facet_wrap(~ Region, scales = "free") +
   labs(
      title = "COVID-19: Observed, Smoothed Observed, and Distributed Encounters Over Weeks",
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

df_combined <- df_combined %>%
   mutate(Region = factor(Region))

test_start_week <- as.Date("2023-07-01")

# Plot
ggplot(df_combined, aes(x = Week)) +
   geom_line(aes(y = Distributed_Count, color = "Distributed Encounters"), size = 0.5) +
   geom_line(aes(y = Smoothed_Actual_Count, color = "Smoothed Observed Encounters"), size = 0.5) +
   geom_line(aes(y = Actual_Count, color = "Observed Encounters"), size = 0.5) +
   geom_vline(xintercept = as.numeric(test_start_week), linetype = "dotted", color = "black", size = 0.7) +
   facet_wrap(~ Region, scales = "free") +
   labs(
      title = "Influenza: Observed vs. Distributed Encounters Over Weeks",
      x = "Week",
      y = "Count",
      color = ""
   ) +
   theme_minimal() +
   scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
   scale_color_manual(values = c("Distributed Encounters" = "blue", "Observed Encounters" = "red", "Smoothed Observed Encounters" = "black")) +
   theme(
      axis.text.x = element_text(angle = 60, hjust = 1), 
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "bottom"
   )



