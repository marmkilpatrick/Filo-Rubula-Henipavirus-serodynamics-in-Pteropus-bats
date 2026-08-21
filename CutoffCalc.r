# Define parameters for the two normal distributions
# m1 <- 0;    s1 <- 1   # Distribution 1: Mean = 0, SD = 1
# m2 <- 1.5;  s2 <- 1.5 # Distribution 2: Mean = 1.5, SD = 1.5

find_normal_intersections <- function(m1, s1, m2, s2) {
  # Formulate the quadratic equation coefficients: Ax^2 + Bx + C = 0
  # Derived from setting dnorm(x, m1, s1) == dnorm(x, m2, s2)
  A <- (1 / (2 * s1^2)) - (1 / (2 * s2^2))
  B <- (m2 / s2^2) - (m1 / s1^2)
  C <- (m1^2 / (2 * s1^2)) - (m2^2 / (2 * s2^2)) - log(s2 / s1)
  
  # Find roots using polyroot (expects coefficients from lowest degree to highest)
  roots <- polyroot(c(C, B, A))
  
  # Filter for real roots (discard negligible imaginary components)
  real_roots <- Re(roots[abs(Im(roots)) < 1e-6])
  return(sort(real_roots))
}

# # 1. Compute the switching point(s)
# intersections <- find_normal_intersections(m1, s1, m2, s2)
# print("The intersection/switching point(s) are located at:")
# print(intersections)
# 
# # 2. Visualize the crossover using ggplot2
# library(ggplot2)
# 
# # Create a sequence of x values spanning across both distributions
# x_vals <- seq(min(m1 - 4*s1, m2 - 4*s2), max(m1 + 4*s1, m2 + 4*s2), length.out = 1000)
# df <- data.frame(
#   x = rep(x_vals, 2),
#   Density = c(dnorm(x_vals, m1, s1), dnorm(x_vals, m2, s2)),
#   Distribution = rep(c("Dist 1", "Dist 2"), each = length(x_vals))
# )
# 
# # Plot the densities and highlight the intersection lines
# ggplot(df, aes(x = x, y = Density, color = Distribution)) +
#   geom_line(size = 1) +
#   geom_vline(xintercept = intersections, linetype = "dashed", color = "darkgray", size = 0.8) +
#   geom_point(data = data.frame(x = intersections, y = dnorm(intersections, m1, s1)), 
#              aes(x = x, y = y), color = "red", size = 3) +
#   labs(title = "Normal Distribution Crossover Points",
#        subtitle = paste("Intersections at x =", paste(round(intersections, 4), collapse = ", ")),
#        x = "Value", y = "Probability Density") +
#   theme_minimal()
