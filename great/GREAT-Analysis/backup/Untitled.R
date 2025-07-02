# --- Categorize distances ---
all_data <- all_data %>%
  mutate(TSS_distance_bin = case_when(
    abs(distance2nearestTSS) <= 1000 ~ "promoter (+/-1kb)",
    abs(distance2nearestTSS) <= 5000 ~ "proximal (1-5kb)",
    abs(distance2nearestTSS) <= 20000 ~ "distal (5-20kb)",
    TRUE ~ "intergenic (>20kb)"
  ))

# Convert to factor to control plotting order
all_data$TSS_distance_bin <- factor(all_data$TSS_distance_bin, levels = c("promoter (+/-1kb)", "proximal (1-5kb)", "distal (5-20kb)", "intergenic (>20kb)"))

# --- Boxplot: JSD at 22C by distance ---
ggplot(all_data, aes(x = TSS_distance_bin, y = JSD_bit_22C, fill = context)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "JSD at 22°C by Distance to Nearest TSS",
       x = "Distance to TSS", y = "JSD at 22°C")

# --- JSD Difference (22C - 10C) by TSS distance ---
all_data <- all_data %>%
  mutate(jsd_diff_22_10 = JSD_bit_22C - JSD_bit_10C)

ggplot(all_data, aes(x = TSS_distance_bin, y = jsd_diff_22_10, fill = context)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(title = "JSD Difference (22°C - 10°C) by TSS Distance",
       x = "Distance to TSS", y = "Δ JSD (22C - 10C)")

# --- Methylation direction (hyper/hypo) analysis ---
all_data <- all_data %>%
  mutate(meth_change = meth_22C - meth_10C,
         meth_direction = case_when(
           meth_change > 0.1 ~ "hypermethylated",
           meth_change < -0.1 ~ "hypomethylated",
           TRUE ~ "stable"
         ))

# Bar plot: Direction of methylation change by distance
ggplot(all_data, aes(x = TSS_distance_bin, fill = meth_direction)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Methylation Change Direction by Distance to TSS",
       x = "Distance to TSS", y = "Proportion of Sites")
