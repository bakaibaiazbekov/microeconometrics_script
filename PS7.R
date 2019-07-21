# PS 7 

library(haven)
library(AER)

sah <- read_dta("sah.dta")

# ex 1 
head(sah$sah)
hist(sah$sah)

# ex 2
Y <- cbind(sah$sah)
X <- cbind(sah$age, sah$educ, sah$income, sah$female, sah$black, sah$prdrugs)
Xvar <- c("age", "educ", "income", "female", "black", "prdrugs")

# Descriptive statistics 
summary(Y)
summary(X)
table(Y)

ols <- lm(Y ~ X, data = sah)
summary(ols)

# create dummies 
sah$good1 <- as.numeric(sah$sah == 1 | sah$sah == 2)
sah$good2 <- as.numeric(sah$sah == 1 | sah$sah == 2 | sah$sah == 3)

# estimate the simple probit model with good1
Y <- cbind(sah$good1)
probit1 <- glm(Y ~ X, 
                  family = binomial(link = "probit"), 
                  data = sah)

coeftest(probit1, vcov. = vcovHC, type = "HC1")

# Now estimate the simple probit model with good2
Y <- cbind(sah$good2)
probit2 <- glm(Y ~ X, 
              family = binomial(link = "probit"), 
              data = sah)
coeftest(probit2, vcov. = vcovHC, type = "HC1")
stargazer::stargazer(probit2, type = 'text')

library(MASS)
# Ordered probit
Y <- cbind(sah$sah)
oprob <- polr(as.factor(Y) ~ X, data = sah, Hess = TRUE)
summary(oprob)
stargazer::stargazer(oprob, type = 'text')

# Calculate marginal effects for probit
library(mfx)
probitmfx(probit2, atmean = TRUE, data = sah)
probitmfx(probit2, atmean = FALSE, data = sah)

# Calculate marginal effects for oprobit
library(erer)
attach(sah)
oprob <- polr(factor(sah) ~ age + educ + income + female + black + prdrugs , data = sah, Hess = TRUE)

ocME(oprob)

