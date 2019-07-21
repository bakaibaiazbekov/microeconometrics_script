###############################################################################
# Introduction to Microeconometrics
# SoSe 2019
###############################################################################
# Class 06
###############################################################################

# Install requred packages unless they are already installed
required_packages <- c("haven", "mlogit", "magrittr", "dplyr", "stargazer")
if (!all(required_packages %in% installed.packages()[,1])) {
  install.packages(required_packages[!required_packages %in%
                   installed.packages()[,1]])
}
library(magrittr)

# Question 1
#-----------
d <- haven::read_dta("mus15data.dta")
# Look at variable labels
sapply(d, function(x) attr(x, "label")) %>% as.matrix()
# Summary of the data
d %>% dplyr::summarize_all(list(min, max, mean, sd)) %>%
  matrix(nrow = ncol(d),
         dimnames = list(colnames(d), c("min", "max", "mean", "sd")))
# Labels for the mode variable
haven::print_labels(d$mode)

# Case-specific (or individual specific) variables
d %>% dplyr::group_by(mode) %>%
  dplyr::summarize(N = length(income), mean(income), sd(income))

# Alternative-specific variables
d %>% dplyr::group_by(mode) %>%
  dplyr::summarize(mean(pbeach), mean(ppier), mean(pprivate), mean(pcharter))
d %>% dplyr::group_by(mode) %>%
  dplyr::summarize(mean(qbeach), mean(qpier), mean(qprivate), mean(qcharter))


# Question 2
#-----------
# Read the "formula.data" vignette for the mlogit package
help(package = "mlogit")

# Multinomial logit: case-specific variable only
# First, prepare the data set
# Geberate id's
d$id <- 1:nrow(d)
# Rename the alternative-specific variables
# E.g., pbeach -> p_beach
for (i in 4:15) {
  name <- colnames(d)[i]
  varname <- substr(name, 1, 1)
  alt <- substr(name, 2, nchar(name))
  colnames(d)[i] <- paste(varname, alt, sep = "_")
}
rm(i, name, varname, alt)

# Transform the mode variable into a factor
d$mode <- names(attr(d$mode, "labels"))[d$mode] %>% as.factor()
# Transform the data
# the command does pretty much the same thing as the reshape Stata command
d <- mlogit::mlogit.data(d, shape = "wide", choice = "mode", varying = 4:15,
                         sep = "_")
# Define the formula of the model
# No alternative-specific variables
frm_ml <- mlogit::mFormula(mode ~ 0 | income)

# Estimate the ML model
m_ml <- mlogit::mlogit(frm_ml, data = d)
summary(m_ml)

# Question 3
#-----------
# To get rrr, just exponentiate the coefficients
sum_ml_rrr <- summary(m_ml)
sum_ml_rrr$CoefTable[,1] <- exp(sum_ml_rrr$coefficients)
colnames(sum_ml_rrr$CoefTable)[1] <- "RRR"
# Use delta method to get the standard errors
R <- diag(exp(sum_ml_rrr$coefficients))
V <- -solve(sum_ml_rrr$hessian)
sum_ml_rrr$CoefTable[,2] <- sqrt(diag(R %*% V %*% R))
# Look at the results
sum_ml_rrr
rm(R, V)

# Question 4
#-----------
# First, calculate the means of variables
md <- dplyr::group_by(d, alt) %>%
  dplyr::summarize(d = mean(d), p = mean(p), q = mean(q), income = mean(income))
# MEM
effects(m_ml, covariate = "income", data = md)
# Unfortunately, the mlogit package does not produce s.e. for marginal effects
# AME
# Implement AME manually (allowing for alternative-specific variables to use
# later)
ame <- function(alt, dataset, model, variable) {
         d1 <- dataset
         delta <- sd(d1[[variable]])/1e+5
         if (!is.null(alt)) {
           delta <- ifelse(d1$alt == alt, delta, 0)
         }
         d1[[variable]] <- d1[[variable]] + delta
         ((predict(model, newdata = d1) -
           predict(model, newdata = dataset)) / max(delta)) %>%
         colMeans()
}
ame(alt = NULL, dataset = d, model = m_ml, variable = "income")

# Question 5
#-----------
m_asc <- mlogit::mlogit(mode ~ p + q | income, data = d)
summary(m_asc)
# MEM
effects(m_asc, covariate = "p", data = md)
# AME
sapply(unique(d$alt), FUN = ame, dataset = d, model = m_asc, variable = "p")
sapply(unique(d$alt), FUN = ame, dataset = d, model = m_asc, variable = "q")
# The matrices are slightly asymmetric because we are using numeric derivatives


# Question 7
#-----------
# A subset model
m_asc2 <- mlogit::mlogit(mode ~ p + q | income,
                         alt.subset = c("beach", "private", "charter"),
                         data = d)
# The hmftest function from the mlogit package is wrong since it includes
# the intercepts which clearly have different interpretations depending on the
# alternatives we consider: p(x = 0) will be mechanically different with
# different subsets of alternatives. And the test is built in such a way that
# under null the two vectors of coefficients should be asymptotically equal.
# If we run hmftest nevertheless, the H statistic is negative!
# (Recall that under null it is asymptotically chi-squared)
# So shame on you, Yves Croissant (the author of mlogit)!

# Calculate the test statistic "manually"
b1 <- m_asc$coefficients
b2 <- m_asc2$coefficients
b1 <- b1[!grepl("intercept", names(b1))]
b2 <- b2[!grepl("intercept", names(b2))]
b1 <- b1[names(b2)]
V1 <- -solve(summary(m_asc)$hessian)
V2 <- -solve(summary(m_asc2)$hessian)
V1 <- V1[names(b1), names(b1)]
V2 <- V2[names(b1), names(b1)]
H <- t(b2 - b1) %*% solve(V2 - V1) %*% (b2 - b1) %>% as.numeric()
# Display the result
cat("\nHausman test\n-------------\nH =", H, "(df =", paste0(length(b1), ")\n"),
  "\rp-value =", pchisq(H, length(b1), lower = FALSE), "\n\n")

rm(b1, b2, V1, V2)  

# Question 9
#-----------
m_nest <- mlogit::mlogit(mode ~ p + q | income,
                         nests = list(shore = c("beach", "pier"),
                                      boat = c("private", "charter")),
                         data = d, tol = 1e-16, ftol = 1e-16, steptol = 1e-16)
summary(m_nest)
# Note that here R is actually doing a better job than Stata does.
# Comparing the log likelihood,
# Stata ll = -1192.41156191258
# R ll = -1192.41155943438
m_nest$logLik %>% as.numeric() %>% format(digits = 16)
# It also takes Stata forever for some reason.
# (To be fair, I tweaked the tolerance options quite a bit here but doing the
# same with Stata would make it run even longer).
# The main conclusion here is that the likelihood of a nested logit is not
# very well-behaved.

# Test IIA using the LR test
mlogit:::lrtest.mlogit(m_asc, m_nest)
# The result is fairly close to Stata

# Question 10
#-----------
# AME of prices
sapply(unique(d$alt), FUN = ame, dataset = d, model = m_nest, variable = "p") 
# AME of catch rates
sapply(unique(d$alt), FUN = ame, dataset = d, model = m_nest, variable = "q") 
# Note that now the signs of the effects of prices/catch rates of correlated
# modes are the same

# Question 11
#-----------
# Define a function to calculate the BIC
BIC <- function(model, dataset) {
  n <- nrow(dataset)
  k <- logLik(model) %>% attr("df")
  log(n) * k - 2 * as.numeric(logLik(model))
}

# Compare the models
stargazer::stargazer(m_asc, m_nest, type = "text")
matrix(c(AIC(m_asc), BIC(m_asc, d), AIC(m_nest), BIC(m_nest, d)),
       nrow = 2,
       dimnames = list(c("AIC", "BIC"),
                       c("Conditional logit", "Nested logit")))

