# --- Subset TEs ---
te_data <- all_data %>% filter(grepl("TE", class, ignore.case = TRUE))

# Optional: extract TE family from gene_id if applicable
te_data <- te_data %>%
  mutate(te_family = ifelse(grepl("AT[1-5]TE\\d+", nearestTSS.gene_id),
                            substr(nearestTSS.gene_id, 7, 9), "unknown"))

# --- JSD changes ---
te_data <- te_data %>%
  mutate(jsd_diff_16_10 = JSD_bit_16C - JSD_bit_10C,
         jsd_diff_22_16 = JSD_bit_22C - JSD_bit_16C)

# --- Methylation changes ---
te_data <- te_data %>%
  mutate(meth_diff_16_10 = meth_16C - meth_10C,
         meth_diff_22_16 = meth_22C - meth_16C)

# --- Boxplot: JSD in TE classes ---
ggplot(te_data, aes(x = te_family, y = jsd_diff_22_16, fill = context)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(title = "JSD Difference (22C - 16C) across TE families",
       x = "TE Family", y = "Δ JSD (22C - 16C)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- Compare JSD in TEs vs. flanking (non-TE) regions ---
all_data <- all_data %>%
  mutate(is_TE = grepl("TE", class, ignore.case = TRUE))

ggplot(all_data, aes(x = is_TE, y = JSD_bit_22C, fill = context)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(title = "JSD at 22°C: TEs vs Non-TEs",
       x = "Is Transposable Element", y = "JSD at 22°C")

ggplot(all_data, aes(x = is_TE, y = JSD_bit_22C - JSD_bit_10C, fill = context)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(title = "JSD Change (22C - 10C): TEs vs Non-TEs",
       x = "Is Transposable Element", y = "Δ JSD (22C - 10C)")
