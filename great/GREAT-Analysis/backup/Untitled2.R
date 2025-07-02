# --- Add methylation differences ---
all_data <- all_data %>%
  mutate(
    meth_diff_16_10 = meth_16C - meth_10C,
    meth_diff_22_16 = meth_22C - meth_16C,
    abs_meth_diff_16_10 = abs(meth_diff_16_10),
    abs_meth_diff_22_16 = abs(meth_diff_22_16)
  )

# --- Correlation between absolute meth change and JSD ---
cor_16_10 <- cor.test(all_data$abs_meth_diff_16_10, all_data$JSD_bit_16C, use = "complete.obs")
cor_22_16 <- cor.test(all_data$abs_meth_diff_22_16, all_data$JSD_bit_22C, use = "complete.obs")

cat("Correlation between abs(methylation diff 16-10) and JSD at 16C:\n")
print(cor_16_10)

cat("\nCorrelation between abs(methylation diff 22-16) and JSD at 22C:\n")
print(cor_22_16)



ggplot(all_data, aes(x = abs_meth_diff_16_10, y = JSD_bit_16C, color = context)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_minimal() +
  labs(title = "JSD vs. Abs Methylation Change (16C - 10C)",
       x = "Absolute Methylation Difference", y = "JSD at 16°C")

ggplot(all_data, aes(x = abs_meth_diff_22_16, y = JSD_bit_22C, color = context)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  theme_minimal() +
  labs(title = "JSD vs. Abs Methylation Change (22C - 16C)",
       x = "Absolute Methylation Difference", y = "JSD at 22°C")
