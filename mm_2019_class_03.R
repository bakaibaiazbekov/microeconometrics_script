###############################################################################
# Introduction to Microeconometrics
# SoSe 2019
###############################################################################
# Class 03
###############################################################################

required_packages <- c("haven", "car", "lmtest", "nlWaldTest",
												"sampleSelection")
if (!all(required_packages %in% installed.packages()[,1])) {
	install.packages(required_packages[!required_packages %in%
									 installed.packages()[,1]])
}

d <- haven::read_dta("mus10cs.dta")
mus10cs$visit <- ifelse(mus10cs$docvis > 0, 1, 0)

# The glm function uses iteratively reweighted least squares to estimate
# the probit model
m <- glm(visit ~ private + chronic + female + income, data = mus10cs,
				 family = binomial(link = "probit"))
# The results are slightly different from Stata since a different method is used
summary(m)

# The probit function from sampleSelection package uses ML and Newton-Raphson
# algorithm as Stata does
m <- sampleSelection::probit(visit ~ private + chronic + female + income,
												data = mus10cs)
# Here we get the same result as in Stata since the same algorithm is used
# We will use these results
summary(m)

# Exercise 2b: Wald test
# The following three are equivalent
lmtest::waldtest(m, "female")
car::linearHypothesis(m, "female")
nlWaldTest::nlWaldtest(m, "b[4]=0")

# As are the following two
car::linearHypothesis(m, "private + chronic = 5")
nlWaldTest::nlWaldtest(m, "b[2] + b[3] = 5")

# Testing multiple hypotheses
nlWaldTest::nlWaldtest(m, "b[2] + b[3] = 5; b[4] = 0")

# Nonlinear restriction
nlWaldTest::nlWaldtest(m, "b[4]/b[2] = 1")
nlWaldTest::nlWaldtest(m, "b[4] = b[2]")

# Exercise 2c: LR test
lmtest::lrtest(m, . ~ . -private - chronic)

# Manually
mres <- update(m, . ~ .-private -chronic)
LR <- -2*(logLik(mres) - logLik(m))
cat("LR = ", as.numeric(LR), "\n")
cat("p-value = ", 1 - pchisq(LR, 2), "\n")
