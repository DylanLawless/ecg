
library(dplyr)
library(lubridate)
library(ggplot2)

df <- read.table("../data/TGA_003.csv", sep = ",", header =TRUE)

# remove empty colunms
df <- df |> select(Subject, RR_interval, RR_value_ms, ECG_date,  RR_time, Quality)

df_events <- read.table("../data/TGA_003_Events.csv", sep = ",", header =TRUE)

head(df)
head(df_events)

names(df)
names(df_events)

# Type `einzelne SVES` = bad
# and `length RR: * ms` = bad

# rename headers to match between datasets
# colnames(df)[colnames(df) == 'oldName'] <- 'newName'
colnames(df)[colnames(df) == 'Subject'] <- 'Patient_Number'
colnames(df)[colnames(df) == 'RR_time'] <- 'Start_time'
colnames(df_events)[colnames(df_events) == 'Date_'] <- 'ECG_date'
df$RR_time <- df$Start_time # but keep 

df_merge <- merge(df, df_events)

# fix the dates
df$ECG_date <- dmy(df$ECG_date)
df_events$ECG_date <- dmy(df_events$ECG_date)

# Combine the date and time into a single datetime column for df_events
df_events$Datetime <- as.POSIXct(paste(df_events$ECG_date, df_events$Start_time), format = "%Y-%m-%d %H:%M:%S")
df$Datetime <- as.POSIXct(paste(df$ECG_date, df$Start_time), format = "%Y-%m-%d %H:%M:%S")

head(df)
head(df_events)

df_merge <- merge(df, df_events, by = c("Patient_Number", "Datetime"))
head(df_merge)

# THESE CANNOT MEGE AS SUBJECT IDs to not overlap: subject 0-1, and 2-4, respectively. 
set.seed(123) # for reproducibility
df <- df |>
  filter(RR_value_ms < 100000) # remove the obvious error outlier
  
# distribution
# Histogram of RR_value_ms
p <- ggplot(df, aes(x = RR_value_ms)) +
  geom_histogram(binwidth = 10, # Adjust binwidth to suitable scale
                 fill="blue", color="black") +
  labs(x = "RR Interval in ms", y = "Frequency", title = "Histogram of RR_value_ms") +
  theme_minimal()

p
ggsave("../images/plot_hist.pdf", plot = p) 


# Points and Violin plot of RR_value_ms
ggplot(df,  aes(x = "", y = RR_value_ms)) +
  geom_violin(trim=FALSE, fill="skyblue") +
  geom_point(stat = "identity", position = position_jitter(width = 0.1), alpha = 0.5) +
  labs(y = "RR Interval in ms", title = "Violin Plot with Points of RR_value_ms") +
  theme_minimal()

p
ggsave("../images/plot_point.pdf", plot = p) 


# Filter the data based on the initial visual cutoff
# Subset to 10% of the data
set.seed(123) # for reproducibility
df_subset <- df |>
  filter(RR_value_ms < 1000) |> # remove the obvious error outlier
  sample_frac(.2) # sample a smaller %


# Histogram of RR_value_ms
p <- df_subset |> ggplot( aes(x = RR_value_ms)) +
  geom_histogram(binwidth = 10, # Adjust binwidth to suitable scale
                 fill="blue", color="black") +
  labs(x = "RR Interval in ms", y = "Frequency", title = "Histogram of RR_value_ms") +
  theme_minimal()

p
ggsave("../images/plot_iqr_hist.pdf", plot = p) 


# Points and Violin plot of RR_value_ms
p <- df_subset |> ggplot(  aes(x = "", y = RR_value_ms)) +
  geom_violin(trim=FALSE, fill="skyblue") +
  geom_point(stat = "identity", position = position_jitter(width = 0.1), alpha = 0.5) +
  labs(y = "RR Interval in ms", title = "Violin Plot with Points of RR_value_ms") +
  theme_minimal()

p
ggsave("../images/plot_iqr_dot.pdf", plot = p) 


# Figure out cut-off ---
# IQR method ----
# 
# # Recalculate IQR and determine refined outlier thresholds
iqr_val_refined <- IQR(df_subset$RR_value_ms)
q1_refined <- quantile(df_subset$RR_value_ms, 0.25)
q3_refined <- quantile(df_subset$RR_value_ms, 0.75)
lower_bound_refined <- q1_refined - 1.5 * iqr_val_refined
upper_bound_refined <- q3_refined + 1.5 * iqr_val_refined

# Output the refined thresholds
print(paste("Refined Lower Bound:", lower_bound_refined))
print(paste("Refined Upper Bound:", upper_bound_refined))

# Filter out values outside the IQR boundaries
df_final <- df_subset %>%
  filter(RR_value_ms >= lower_bound_refined & RR_value_ms <= upper_bound_refined)

# It's useful to compare how many data points remain
cat("Original Data Points: ", nrow(df), "\n")
cat("Data Points after Initial Filtering: ", nrow(df_subset), "\n")
cat("Data Points after Outlier Removal: ", nrow(df_final), "\n")

# subjects
p <- df_subset |> ggplot(  aes(x = Patient_Number, y = RR_value_ms)) +
  geom_violin(trim=FALSE, fill="skyblue") +
  geom_point(stat = "identity", position = position_jitter(width = 0.1), alpha = 0.5) +
  theme_minimal() + 
  geom_hline(aes(yintercept = upper_bound_refined), linetype = "dashed", color = "red") + 
  geom_hline(aes(yintercept = lower_bound_refined), linetype = "dashed", color = "blue")

p
ggsave("../images/plot_iqr_dot_subject.pdf", plot = p) 


# plot time course
p <- df_subset |> ggplot(  aes(x = Datetime, y = RR_value_ms)) +
  geom_point(stat = "identity", position = position_jitter(width = 0.1), alpha = 0.5) +
  theme_minimal() + 
  geom_hline(aes(yintercept = upper_bound_refined), linetype = "dashed", color = "red") + 
  geom_hline(aes(yintercept = lower_bound_refined), linetype = "dashed", color = "blue")

p
ggsave("../images/plot_iqr_time.pdf", plot = p) 



p <- df_subset |> 
  filter(Datetime >= as.POSIXct("2015-04-01 00:00:00") & 
           Datetime <= as.POSIXct("2015-04-01 00:15:00")) |>
  ggplot(  aes(x = Datetime, y = RR_value_ms)) +
  geom_point(stat = "identity", position = position_jitter(width = 0.1), alpha = 0.5) +
  geom_line() +
  labs(title = "15 minutes") +
  theme_minimal() + 
  geom_hline(aes(yintercept = upper_bound_refined), linetype = "dashed", color = "red") + 
  geom_hline(aes(yintercept = lower_bound_refined), linetype = "dashed", color = "blue")

p
ggsave("../images/plot_iqr_time_15min.pdf", plot = p) 


  
# Modified Z-score method ----

# Constants
median_rr <- median(df_subset$RR_value_ms, na.rm = TRUE)
mad_rr <- mad(df_subset$RR_value_ms, constant = 1)  # standard deviation estimator
threshold <- 3.5  # Adjust this as needed

# Calculate thresholds
upper_threshold_z <- median_rr + threshold * mad_rr / 0.6745
lower_threshold_z <- median_rr - threshold * mad_rr / 0.6745

# Plotting
p <- ggplot(df_subset, aes(x = Datetime, y = RR_value_ms)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  geom_hline(yintercept = upper_threshold_z, linetype = "dashed", color = "red") +
  geom_hline(yintercept = lower_threshold_z, linetype = "dashed", color = "blue") +
  labs(title = "RR_value_ms with Modified Z-score Thresholds") +
  theme_minimal()

p
ggsave("../images/plot_zscore.pdf", plot = p) 


# Percentile Thresholds method ----

# Calculate percentiles
p1 <- quantile(df_subset$RR_value_ms, 0.01)
p99 <- quantile(df_subset$RR_value_ms, 0.99)

# Plotting
p <- ggplot(df_subset, aes(x = Datetime, y = RR_value_ms)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  geom_hline(yintercept = p99, linetype = "dashed", color = "red") +
  geom_hline(yintercept = p1, linetype = "dashed", color = "blue") +
  labs(title = "RR_value_ms with Percentile Thresholds") +
  theme_minimal()
p
ggsave("../images/plot_percentile.pdf", plot = p) 


# KDE-Based Outlier Detection method ----

# Calculate the density
dens <- density(df_subset$RR_value_ms)
cutoff_density <- 0.01  # Low-density cutoff, adjust based on plot(dens)

# Find cutoffs
low_cutoff <- min(dens$x[dens$y > cutoff_density])
high_cutoff <- max(dens$x[dens$y > cutoff_density])

# Plotting
p <- ggplot(df_subset, aes(x = Datetime, y = RR_value_ms)) +
  geom_point(position = position_jitter(width = 0.1), alpha = 0.5) +
  geom_hline(yintercept = high_cutoff, linetype = "dashed", color = "red") +
  geom_hline(yintercept = low_cutoff, linetype = "dashed", color = "blue") +
  labs(title = "RR_value_ms with KDE-Based Thresholds") +
  theme_minimal()
p
ggsave("../images/plot_kde.pdf", plot = p) 


# K means method ----

library(dplyr)
library(ggplot2)
library(cluster)  # For clustering

# Assuming df_subset has the relevant 'RR_value_ms' data
set.seed(123)  # For reproducibility
kmeans_result <- kmeans(df_subset$RR_value_ms, centers = 3)

# Add cluster results to the dataframe
df_subset$cluster <- as.factor(kmeans_result$cluster)

# Calculate the mean of each cluster
cluster_means <- aggregate(RR_value_ms ~ cluster, data = df_subset, mean)
sorted_clusters <- cluster_means[order(cluster_means$RR_value_ms),]

# Assuming clusters are sorted by their mean RR_value_ms values
threshold_lower <- (sorted_clusters$RR_value_ms[1] + sorted_clusters$RR_value_ms[2]) / 2
threshold_upper <- (sorted_clusters$RR_value_ms[2] + sorted_clusters$RR_value_ms[3]) / 2

p <- ggplot(df_subset, aes(x = Datetime, y = RR_value_ms, color = cluster)) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1)) +
  geom_hline(yintercept = threshold_lower, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = threshold_upper, linetype = "dashed", color = "red") +
  labs(title = "RR_value_ms with k-Means Clustering") +
  theme_minimal() +
  scale_color_brewer(palette = "Dark2")  # Color coding clusters

p
ggsave("../images/plot_kmean.pdf", plot = p) 


#  dbscan ----
# When clusters vary significantly in size, standard k-means clustering, which minimizes variance within clusters, might not perform optimally. This could potentially skew the cluster centers towards the larger cluster. To handle this, you might consider using a clustering method that can cope better with such discrepancies.
# 
# Alternatives and Adjustments
# Weighted K-Means: One approach is to use a weighted version of k-means that accounts for cluster size by assigning weights inversely proportional to the number of points in each cluster. This can help in balancing the influence of each cluster.
# 
# Density-Based Clustering (e.g., DBSCAN): This method does not require specifying the number of clusters beforehand and can handle clusters of varying sizes and shapes. DBSCAN groups together closely packed points and marks points in low-density regions as outliers.
# 
# Trimmed K-Means: This method involves running k-means while ignoring a certain percentage of the most extreme observations during the calculation of centroids, potentially leading to more robust cluster centers.

# Applying DBSCAN; you may need to tune eps and minPts
# eps is the maximum distance between two samples for them to be considered as in the same neighborhood
# minPts is the number of samples in a neighborhood for a point to be considered as a core point

library(dbscan)

df_subset2 <- df_subset |> sample_frac(0.2) # subset if slow, e.g.: 0.2

# Reshape RR_value_ms into a matrix
rr_matrix <- matrix(df_subset2$RR_value_ms, ncol = 1)

# Create a range for eps and minPts
eps_range <- seq(50, 500, by = 50)  # Adjust this range based on your data scale
minPts_range <- seq(5, 20, by = 5)

# Placeholder to store the best result
best_result <- NULL
found_clusters <- FALSE

for (eps in eps_range) {
  for (minPts in minPts_range) {
    # Apply DBSCAN
    dbscan_result <- dbscan(rr_matrix, eps = eps, minPts = minPts)
    
    # Check the number of clusters (excluding noise)
    num_clusters <- length(unique(dbscan_result$cluster)) - 1
    
    # Check if more than one cluster (excluding noise)
    if (num_clusters > 1) {
      # Store the result and break out of loops if clusters are found
      best_result <- list(eps = eps, minPts = minPts, clusters = dbscan_result$cluster)
      found_clusters <- TRUE
      print(best_result)
      break
    }
  }
  if (found_clusters) break
}

# Check if clusters were found and add results to the dataframe
if (!is.null(best_result)) {
  df_subset2$cluster <- as.factor(best_result$clusters)
  
  # Plotting
  p <- ggplot(df_subset2, aes(x = Datetime, y = RR_value_ms, color = cluster)) +
    geom_point(alpha = 0.5, position = position_jitter(width = 0.1)) +
    labs(title = paste("RR_value_ms with DBSCAN Clustering - eps:", best_result$eps, "minPts:", best_result$minPts)) +
    theme_minimal() +
    scale_color_brewer(palette = "Set1")  # Color coding clusters
  
  ggsave("../images/plot_dbscan.pdf", plot = p) 
  p
} else {
  print("No suitable clustering found with the specified range of parameters.")
}

# Filter to keep only cluster 2
df_cluster2 <- df_subset2 %>% 
  filter(cluster == 1)

library(ggplot2)

# Histogram of RR_value_ms for cluster 2
p <- ggplot(df_cluster2, aes(x = RR_value_ms)) +
  geom_histogram(binwidth = 10,  # Adjust binwidth to suitable scale
                 fill="blue", color="black") +
  labs(x = "RR Interval in ms", y = "Frequency", title = "Histogram of RR_value_ms for Cluster 2") +
  theme_minimal()

p
ggsave("../images/plot_dbscan_hist.pdf", plot = p)



