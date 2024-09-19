# Function to calculate a1
calc_a1 <- function(n) {
  a1 <- sum(1 / (1:(n-1)))
  return(a1)
}

# Function to calculate a2
calc_a2 <- function(n) {
  a2 <- sum(1 / (1:(n-1))^2)
  return(a2)
}

# Function to calculate b1
calc_b1 <- function(n) {
  b1 <- (n + 1) / (3 * (n - 1))
  return(b1)
}

# Function to calculate b2
calc_b2 <- function(n) {
  b2 <- (2 * (n^2 + n + 3)) / (9 * n * (n - 1))
  return(b2)
}

# Function to calculate c1
calc_c1 <- function(n, a1) {
  b1 <- calc_b1(n)
  c1 <- b1 - (1 / a1)
  return(c1)
}

# Function to calculate c2
calc_c2 <- function(n, a1, a2) {
  b2 <- calc_b2(n)
  c2 <- b2 - (n + 2) / (a1 * n) + (a2 / a1^2)
  return(c2)
}

# Function to calculate e1
calc_e1 <- function(c1, a1) {
  e1 <- c1 / a1
  return(e1)
}

# Function to calculate e2
calc_e2 <- function(c2, a1, a2) {
  e2 <- c2 / (a1^2 + a2)
  return(e2)
}