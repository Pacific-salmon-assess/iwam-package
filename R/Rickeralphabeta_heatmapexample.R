# HEATMAP of alpha and beta parameter values on SMSY and SREP ####
library(ggplot2)
library(reshape2)
library(pracma)  # For Lambert W function

# Define the range for a.par and b.par
a_par_values <- seq(1, 10, length.out = 20)
b_par_values_new <- seq(0.001, 0.002, length.out = 20)

# Create an empty matrix to store SMSY values
SMSY_matrix_new <- matrix(0, nrow = length(a_par_values), ncol = length(b_par_values_new))
SREP_matrix_new <- matrix(0, nrow = length(a_par_values), ncol = length(b_par_values_new))

# Compute SMSY values
for (i in seq_along(a_par_values)) {
  for (j in seq_along(b_par_values_new)) {
    a_par <- a_par_values[i]
    b_par <- b_par_values_new[j]
    
    W0_term <- gsl::lambert_W0(exp(1 - log(a_par)))
    SMSY_matrix_new[i, j] <- (1 - W0_term) / b_par
    SREP_matrix_new[i, j] <- log(a_par) / b_par
  }
}

# Convert matrix to data frame for ggplot
df <- expand.grid(a_par = a_par_values, b_par = b_par_values_new)
df$SMSY <- as.vector(SMSY_matrix_new)

df2 <- expand.grid(a_par = a_par_values, b_par = b_par_values_new)
df2$SREP <- as.vector(SREP_matrix_new)

# HEATMAP: SMSY
ggplot(df, aes(x = b_par, y = a_par, fill = SMSY)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  labs(title = "Heatmap of SMSY values with Adjusted b.par Range",
       x = "b.par values (0.00001 to 0.002)",
       y = "a.par values",
       fill = "SMSY") +
  theme_minimal()

# HEATMAP: SREP
ggplot(df2, aes(x = b_par, y = a_par, fill = SREP)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  labs(title = "Heatmap of SREP values with Adjusted b.par Range",
       x = "b.par values (0.00001 to 0.002)",
       y = "a.par values",
       fill = "SREP") +
  theme_minimal()

# LINEAR COMPARISON
ggplot(df, aes(x = b_par, y = a_par)) +
  geom_line() +
  theme_minimal()
