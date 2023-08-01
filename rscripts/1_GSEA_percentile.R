# Split the data by clustering codification <- ADD IN LINE 147
# Absolute zscore value for percentile calculation
zscore_k1$zscore_abs <- abs(zscore_k1$zscore)
zscore_k2$zscore_abs <- abs(zscore_k2$zscore)
zscore_k3$zscore_abs <- abs(zscore_k3$zscore)
zscore_k4$zscore_abs <- abs(zscore_k4$zscore)
zscore_k5$zscore_abs <- abs(zscore_k5$zscore)
zscore_k6$zscore_abs <- abs(zscore_k6$zscore)

# Percentile based on abolute zscore value
zscore_k1$zscore_per <- percent_rank(zscore_k1$zscore_abs)
zscore_k2$zscore_per <- percent_rank(zscore_k2$zscore_abs)
zscore_k3$zscore_per <- percent_rank(zscore_k3$zscore_abs)
zscore_k4$zscore_per <- percent_rank(zscore_k4$zscore_abs)
zscore_k5$zscore_per <- percent_rank(zscore_k5$zscore_abs)
zscore_k6$zscore_per <- percent_rank(zscore_k6$zscore_abs)

# Filter by percentile
zscore_k1 <- zscore_k1 %>% 
  filter(zscore_per > 0.5)
zscore_k2 <- zscore_k2 %>% 
  filter(zscore_per > 0.5)
zscore_k3 <- zscore_k3 %>% 
  filter(zscore_per > 0.5)
zscore_k4 <- zscore_k4 %>% 
  filter(zscore_per > 0.5)
zscore_k5 <- zscore_k5 %>% 
  filter(zscore_per > 0.5)
zscore_k6 <- zscore_k6 %>% 
  filter(zscore_per > 0.5)
