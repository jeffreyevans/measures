# Load the necessary libraries
library(rgl)

# Create 10 time series
set.seed(123)
n <- 1000
time <- 1:n
time_series_data <- lapply(1:100, function(i) {
  cumsum(rnorm(n))
})

# Create a grid of time and series indices
time_grid <- rep(time, times = 100)
series_indices <- rep(1:100, each = n)

# Create a matrix to store the values for the surface plot
surface_matrix <- matrix(unlist(time_series_data), ncol = 100, byrow = TRUE)

# Create the 3D surface plot
plot3d(time_grid, series_indices, surface_matrix, type = "n", col = "blue", xlab = "Time", ylab = "Series", zlab = "Value")

# Add the 10 time series as lines on the surface
for (i in 1:100) {
  lines3d(time, rep(i, n), time_series_data[[i]], col = i, size = 2)
}
