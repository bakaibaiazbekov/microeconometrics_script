###############################################################################
# Introduction to Microeconometrics
# SoSe 2019
###############################################################################
# Class 04
###############################################################################

# Install requred packages unless they are already installed
required_packages <- c("haven", "maxLik", "numDeriv", "mfx", "tidyverse",
											 "magrittr")
if (!all(required_packages %in% installed.packages()[,1])) {
	install.packages(required_packages[!required_packages %in%
									 installed.packages()[,1]])
}
library(magrittr)

# Exercise 1
# (b)
d <- haven::read_dta("binlfp2.dta")

# Logistic CDF
logistic <- function(x) {
	exp(x)/(1 + exp(x))
}

# Conditional log likelihood
logit_lf <- function(coefs, formula, data) {
	X <- model.matrix(formula, data = data)
	Y <- model.frame(formula, data = data)[[1]]
	attributes(Y) <- NULL
	as.numeric(Y*log(logistic(X %*% coefs)) +
						 (1 - Y)*log(logistic(-X %*% coefs)))
}

# Model formula
frm <- lfp ~ k5 + k618 + age + wc + inc
# Log likelihood for our sample
loglik <- function(coefs) logit_lf(coefs, frm, d)

# Maximization
ml <- maxLik::maxLik(loglik, start = rep(0, 6))
# The Hessian calculated by maxLik is rubbish, so use the numDeriv package
H <- numDeriv::hessian(function(b) sum(loglik(b)), ml$estimate)

# A standard matrix reporting estimation results
coefs <- cbind(`Estimate` = ml$estimate,
							 `Std.Err.` = sqrt(diag(-solve(H))))
coefs <- cbind(coefs, `z value` = coefs[,"Estimate"]/coefs[,"Std.Err."])
coefs <- cbind(coefs, `p-value` = 2*pnorm(abs(coefs[,"z value"]),
																					lower.tail = FALSE))
rownames(coefs) <- c("(Intercept)", all.vars(frm)[-1])
coefs

# Stata output to compare
#-------------------------------------------------------------------------
#    lfp |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
#--------+----------------------------------------------------------------
#  _cons |   3.799878   .6228765     6.10   0.000     2.579062    5.020693
#     k5 |  -1.463093   .1946494    -7.52   0.000    -1.844599   -1.081587
#   k618 |  -.0875886   .0672691    -1.30   0.193    -.2194336    .0442565
#    age |  -.0636601     .01258    -5.06   0.000    -.0883165   -.0390037
#     wc |   1.062299   .1987034     5.35   0.000      .672847     1.45175
#    inc |  -.0307647    .007685    -4.00   0.000     -.045827   -.0157024
#-------------------------------------------------------------------------

# Compare to the R built-in command
# The std.errors are a bit differemt since a diferrent method is used
m <- glm(frm, d, family = binomial(link = "logit"))
summary(m)

# (c)
# For this question, we will use the results from the GLM function. They are
# close enough to those obtained by implementing the MLE from scratch.

# Marginal effects at the mean (MEM)
#------------------------------------
# We will use the mfx package which is specifically designed to calculate
# marginal effects in glm models
# Note that unlike Stata, R figures out what variables are dummies without
# hints from the user
frm2 <- lfp ~ k5 + k618 + age + wc + hc + lwg + inc
me_mfx <- mfx::logitmfx(frm2, d, atmean = TRUE)
me_mfx

# Now let us try to implement this manually (including the delta method)
d_m <- dplyr::summarize_all(d, mean)
m2 <- glm(frm2, data = d, family = binomial(link = "logit"))
d_m$xb <- predict(m2, newdata = d_m, type = "link")
# The marginal effects are
me_man <- exp(d_m$xb)/(1 + exp(d_m$xb))^2 * m2$coefficients[-1]
me_man
# For continuous variables, the manual calculation yields the same
# results


# Of course, the wc and hc variables have to be handled differently
d_m <- rbind(d_m, d_m, d_m, d_m)
d_m[1, "wc"] <- 1
d_m[3, "hc"] <- 1
d_m[2, "wc"] <- 0
d_m[4, "hc"] <- 0
d_m$prob <- predict(m2, newdata = d_m, type = "response")
# wc
sprintf("The ME for wc is %.8f", d_m[1, "prob"] - d_m[2, "prob"])
# hc
sprintf("The ME for hc is %.8f", d_m[3, "prob"] - d_m[4, "prob"])
# Again, the same figures as using the logitmfx function

# The delta method
d$prob <- predict(m2, newdata = d, type = "response")
R <- function(x, m) { # here x is a vector of values of variables for 1 obs
											# (including the predicted probability)
											# and m is the model
	k <- length(all.vars(m$formula))
	p <- x["prob"]
	c1 <- p * (1 - p) * (1 - 2 * p)
	c2 <- p * (1 - p)
	x <- matrix(c(1, x[all.vars(m$formula)[-1]]), ncol = 1)
	b <- matrix(m$coefficients[-1], ncol = 1)
	c1 * (b %*% t(x)) + c2 * cbind(rep(0, k-1), diag(rep(1, k-1)))
}

d_m <- dplyr::summarize_all(d, mean)
d_m$prob <- predict(m2, newdata = d_m, type = "response")

Rhat <- R(d_m %>% unlist(), m2)
Vme <- Rhat %*% summary(m2)$cov.unscaled %*% t(Rhat)
std_errs <- sqrt(diag(Vme))
names(std_errs) <- all.vars(frm2)[-1]
std_errs
# compare to mfx's output
me_mfx
# For continuous variables the results are identical


# Note that for the dummy variables the output is incorrect
# Let us handle them separately
d_m_d <- rbind(d_m, d_m, d_m, d_m)
d_m_d[1, "wc"] <- 1
d_m_d[2, "wc"] <- 0
d_m_d[3, "hc"] <- 1
d_m_d[4, "hc"] <- 0

d_m_d$prob <- predict(m2, newdata = d_m_d, type = "response")
# wc
p1 <- d_m_d[1, "prob"] %>% as.numeric()
p0 <- d_m_d[2, "prob"] %>% as.numeric()
x <- model.matrix(frm2, d_m_d[1:2,])
R <- t(p1 * (1 - p1) * x[1,] - p0 * (1 - p0) * x[2,])
V <- summary(m2)$cov.unscaled
me_wc <- p1 - p0
se_wc <- sqrt(R %*% V %*% t(R)) %>% as.numeric()
# hc
p1 <- d_m_d[3, "prob"] %>% as.numeric()
p0 <- d_m_d[4, "prob"] %>% as.numeric()
x <- model.matrix(frm2, d_m_d[3:4,])
R <- t(p1 * (1 - p1) * x[1,] - p0 * (1 - p0) * x[2,])
V <- summary(m2)$cov.unscaled
me_hc <- p1 - p0
se_hc <- sqrt(R %*% V %*% t(R)) %>% as.numeric()

sprintf("wc: %.8f (%.8f)", me_wc, se_wc)
sprintf("hc: %.8f (%.8f)", me_hc, se_hc)
me_mfx
# We see that we managed to reproduce the results of the logitmfx function

# Average marginal effects
#--------------------------
# It's as simple as
mfx::logitmfx(frm2, d, atmean = FALSE)
# The results are close to MEM but not the same


# Odds ratios
#-------------
# Here we simply calculate the exponent of the coefficients
m2$coefficients %>% exp()


# Exercise 2: Small-sample analytics
# (a)
set.seed(10101)
rm(list = ls())
simprobit <- function(n, trueb2) {
	x <- rnorm(n)
	y_  <- trueb2*x + rnorm(n)
	y <- ifelse(y_ > 0, 1, 0)
	m <- glm(y ~ x, family = binomial(link = "probit"))
	b2hat <- m$coefficients[2]
	seb2hat <- summary(m)$coefficients[2, 2]
	ztestforb2 <- (b2hat - 1)/seb2hat
	c(b2hat = b2hat, seb2hat = seb2hat, ztestforb2 = ztestforb2)
}
# number of MC replications
R <- 1e+4
# Run the simulation
d <- sapply(1:R, function(i) simprobit(40, 1))
d <- as.data.frame(t(d))
# 1. Compare the mean of b2hat with the true value (1)
# 2. Compare the s.e. of b2hat with the mean estimate of this s.e. (seb2hat)
summary(d)
sd(d$b2hat.x)

d$reject0.01 <- ifelse(abs(d$ztestforb2.x) > qnorm(1 - 0.01/2), 1, 0)
d$reject0.05 <- ifelse(abs(d$ztestforb2.x) > qnorm(1 - 0.05/2), 1, 0)
d$reject0.10 <- ifelse(abs(d$ztestforb2.x) > qnorm(1 - 0.10/2), 1, 0)
d$reject0.20 <- ifelse(abs(d$ztestforb2.x) > qnorm(1 - 0.20/2), 1, 0)
d %>% dplyr::summarize(size0.01 = mean(reject0.01),
								 size0.05 = mean(reject0.05),
								 size0.10 = mean(reject0.10),
								 size0.20 = mean(reject0.20))
hist(d$ztestforb2.x, freq = FALSE, xlim = c(-4, 4),
	main = expression(paste("Histogram of ", z(hat(beta)[2]))),
	xlab = "z")
curve(dnorm(x), add = TRUE, lwd = 2, col = "blue")
# Summary: The distribution of the Wald test statistic is not staggeringly close
# to normal.

# (b)
R <- 1e+4
d <- sapply(1:R, function(i) simprobit(40, 2))
d <- as.data.frame(t(d))
# 1. Compare the mean of b2hat with the true value (2)
# 2. Compare the s.e. of b2hat with the mean estimate of this s.e. (seb2hat)
summary(d)
sd(d$b2hat.x)
# Somathing weird is going on. Due to outliers, the mean of our estimates
# is noticeably different from the median

d$reject0.01 <- ifelse(abs(d$ztestforb2.x) > qnorm(1 - 0.01/2), 1, 0)
d$reject0.05 <- ifelse(abs(d$ztestforb2.x) > qnorm(1 - 0.05/2), 1, 0)
d$reject0.10 <- ifelse(abs(d$ztestforb2.x) > qnorm(1 - 0.10/2), 1, 0)
d$reject0.20 <- ifelse(abs(d$ztestforb2.x) > qnorm(1 - 0.20/2), 1, 0)
d %>% dplyr::summarize(size0.01 = mean(reject0.01),
								 size0.05 = mean(reject0.05),
								 size0.10 = mean(reject0.10),
								 size0.20 = mean(reject0.20))
# Reject rate grows as we increase the significance level
