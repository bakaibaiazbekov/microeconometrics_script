###############################################################################
# Introduction to Microeconometrics
# SoSe 2019
###############################################################################
# Class 02
###############################################################################

required_packages <- c("haven", "estimatr", "stats4", "maxLik")	
if (!all(required_packages %in% installed.packages()[,1])) {
	install.packages(required_packages[!required_packages %in%
									 installed.packages()[,1]])
}

# Exercise 2a: Monte Carlo illustration

P <- 0.3		# probability of pulling out a white ball
N <- 100		# number of trials in each simulation run
R <- 1e+4 	# number of replications

one_sim <- function(p, n) {
	y <- sample(c(1, 0), n, prob = c(p, 1-p), replace = TRUE)
	mean(y)
}

# Run R simulations
p_ <- sapply(1:R, function(x) one_sim(P, N))

# Look at the theoretical and simulated std.dev.'s
cat("Theoretical std.dev. = ", sqrt(P*(1-P)/N), "\n")
cat("Simulated std.dev. = ", sd(p_), "\n")

# Plot the densities
plot(density(p_), main = "Density of the estimated parameter")
curve(dnorm(x, P, sqrt(P*(1-P)/N)), add = TRUE, col = "red")

# Exercise 2c: R's stats4::mle function

# The log likelihood function
mybinom <- function(y, p) {
	sum(y * log(p) + (1 - y) * log(1 - p))
}

# Generate the sample
y <- sample(c(1, 0), N, prob = c(P, 1-P), replace = TRUE)

# See ?stats4::mle for details
stats4::mle(function(p) -mybinom(y, p), start = list(p = 0.5))

# Another option is to use the maxLik package
maxLik::maxLik(function(p) mybinom(y, p), start = 0.5)

# Exercise 3b: Newton-Raphson algorirthm
f <- function(theta) log(theta) - exp(theta)
theta <- 2
crit <- 10
iter <- 1

while(crit > 1e-8) {
	d1 <- 1/theta - exp(theta)
	d2 <- -1/theta^2 - exp(theta)
	crit <- -d1^2/d2
	theta <- theta - d1/d2
	sprintf("theta = %8.6g, f = %8.6g, d1 = %8.6g, d2 = %8.6g, Crit =  %8.6g",
					theta, f(theta), d1, d2, crit)
}

# Numerical derivatives
h <- 1e-6
theta <- 2
crit <- 10
iter <- 1
while(crit > 1e-8) {
	d1 <- (f(theta + h) - f(theta))/h
	d2 <- (f(theta + h) - 2*f(theta) + f(theta - h)) / h^2
	crit <- -d1^2/d2
	theta <- theta - d1/d2
	sprintf("theta = %8.6g, f = %8.6g, d1 = %8.6g, d2 = %8.6g, Crit =  %8.6g",
					theta, f(theta), d1, d2, crit)
}
