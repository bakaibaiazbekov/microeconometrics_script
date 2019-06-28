rm(list=ls()) #Removes all items in Environment!
library(nlWaldTest) # for the `nlWaldtest()` function
library(lmtest) #for `coeftest()` and `bptest()`.
library(broom) #for `glance(`) and `tidy()`
library(PoEdata) #for PoE4 datasets
library(car) #for `hccm()` robust standard errors
library(sandwich)
library(knitr) #for `kable()`
library(forecast) 
library(AER)
library(xtable)
library(devtools)
library(dplyr)
library(mlogit)
################################################################################################

data("TravelMode", package = "AER") #load the data
head(TravelMode) #see the head of the data
nrow(TravelMode) #number of rows
length(unique(TravelMode$individual)) #how many individuals?
table(is.na(TravelMode)) #is there any missing data
table(filter(TravelMode, choice =="yes")$mode)/(nrow(TravelMode)/4) #average frequency of choices

# Set Reference Group: car,constant normalize to 0
# Parameters: , length of 6
# constants: 3
# weight on independent variables: 3
# Data: TravelMode(mode,choice,gcost,wait,travel)
# 1. Utility for each alternative given X,Θ
# 2. Choice probability for one single individual given utility of four alternatives
# 3. Choice probabilites for all the individuals
# 4. Sum of LogLikelihood

# 1 Utility for each alternative given X,Θ

theta = c(1,2,3,-0.001,-0.003,-0.005) #generate one set of parameters
sample <- filter(TravelMode, individual == 1)  #filter one individual
sample$constant <- 0
sample$constant[sample$mode=="air"] = theta[1]
sample$constant[sample$mode=="train"] = theta[2]
sample$constant[sample$mode=="bus"] = theta[3]

sample$utility = theta[4]*sample$gcost + theta[5]*sample$wait + theta[6]*sample$travel + sample$constant 
sample

# 2. Choice probability for each individual given utility of four alternatives

real_choice  <- filter(sample,choice=="yes")
real_choice


choice_prob = exp(real_choice$utility) / sum(exp(sample$utility))
choice_prob


# 3. Write functions to repeat step 1 & step 2 for all individuals

# Choice Probability Function
# Input: sample, a data.frame for one individual,contains utility for each alternative and real choice
# Output: choice probability for the individual

# Choice Probability Function for each individual
# Suppose we have the utility for each alternative
choice.prob <- function(sample){
  x = filter(sample,choice=="yes") #filter the choice 
  prob = exp(x$utility)/sum(exp(sample$utility)) #caculate probability
  return(prob)
}


# 4. Grouping, Operations and Aggregations(Within-across)
# Group the dataset by individual
# Compute choice probability within group
# Summarize choice probabilities across group
# group_by(),summarise(),do() function in package dplyr

group = group_by(TravelMode,individual) #group the data by individual
summarise(group,avg_travel = mean(travel)) # calculate the avg travel time for each individual

TravelMode$utility <- runif(nrow(TravelMode),0,1) #generate some random utilities
Probability <- TravelMode  %>% group_by(individual)  %>% do(data.frame( prob= choice.prob(.))) #apply choice.prob() function on each individuals

Probability # a data.frame, each row contains the choice probability for each individual

sum(log(Probability$prob)) #sum up the log likelihoods

# Likelihood function
# Input: parameter
# Intermediate: utility,choice probability
# Output: Loglikelihood
# Likelihood function for the whole data
Likelihood<-function(theta){
  # caclulating utility for each alternative for all the individuals
  TravelMode$constant=0
  TravelMode$constant[TravelMode$mode=="air"] = theta[1]
  TravelMode$constant[TravelMode$mode=="train"] = theta[2]
  TravelMode$constant[TravelMode$mode=="bus"] = theta[3]
  TravelMode$utility = theta[4]*TravelMode$gcost+theta[5]*TravelMode$wait
  +theta[6]*TravelMode$travel+TravelMode$constant 
  
  # caclulating choice probability for each individual
  Probability=TravelMode %>% 
    group_by(individual) %>% # group the data by individual
    do(data.frame(prob=choice.prob(.)))  # use Choice() function on each individual
  return(-sum(log(Probability$prob))) # return the sum log likelihood 
}


# MLE 

# optim()
# Starting value
# Function to optimize
# optimizing method

est <- optim(c(4.0,4.0,4.0,-0.01,-0.01,-0.01),Likelihood,method = "BFGS") 
est$par

# Using mlogit Package

# formatting data
TM <- mlogit.data(TravelMode, choice = "choice", shape = "long", 
                  chid.var = "individual", alt.var = "mode", drop.index = TRUE)
head(TM)

# estimate with mlogit
ml.TM <- mlogit(choice ~ gcost +wait +travel, TM, reflevel = "car")
#show results
summary(ml.TM)

#how fitted choice probability match with data
apply(fitted(ml.TM, outcome=FALSE), 2, mean) # fitted mean choice probability

ml.TM$freq/sum( ml.TM$freq) # mean choice probability in data


################################################################################################
install_git("https://github.com/ccolonescu/PoEdata")

auto.probit <- glm(auto~dtime, family=binomial(link="probit"), 
                   data=transport)

kable(tidy(auto.probit), digits=4, align='c', caption="Transport example, estimated by probit")

# calculate this effect at  dtime=2(time difference of 20 minutes)
xdtime <- data.frame(dtime=2)
predLinear <- predict(auto.probit, xdtime, 
                      data=transport, type="link")
DpDdtime <- coef(auto.probit)[[2]]*dnorm(predLinear)
DpDdtime

# calculate the predicted probability of choosing  auto  when the time difference is 30 minutes ( dtime=3 ):
xdtime <- data.frame(dtime=3)
predProbit <- predict(auto.probit, xdtime, 
                      data=transport, type="response")

# The marginal effect at the average predicted value
avgPredLinear <- predict(auto.probit, type="link")
avgPred <- mean(dnorm(avgPredLinear))
AME <- avgPred*coef(auto.probit)
AME


# COKE
coke.logit <- glm(coke~pratio+disp_coke+disp_pepsi, 
                  data=coke, family=binomial(link="logit"))
kable(tidy(coke.logit), digits=5, align="c",
      caption="Logit estimates for the 'coke' dataset")

coke.LPM <- lm(coke ~ pratio + disp_coke + disp_pepsi, 
               data=coke)
coke.probit <- glm(coke~pratio+disp_coke+disp_pepsi, 
                   data=coke, family=binomial(link="probit"))
stargazer(coke.LPM, coke.probit, coke.logit,
          header=FALSE, 
          title="Three binary choice models for the 'coke' dataset",
          type=.stargazertype,
          keep.stat="n",digits=4, single.row=FALSE,
          intercept.bottom=FALSE,
          model.names=FALSE,
          column.labels=c("LPM","probit","logit"))

# Prediction and marginal effects 
tble <- data.frame(table(true=coke$coke, 
                         predicted=round(fitted(coke.logit))))
kable(tble, align='c', caption="Logit prediction results")

# hypothesis test

Hnull <- "disp_coke+disp_pepsi=0"
linearHypothesis(coke.logit, Hnull)

# Multinomial Logit
library(nnet)

nels.multinom <- multinom(psechoice~grades, data=nels_small)

summary(nels.multinom)

medGrades <- median(nels_small$grades)

fifthPercentileGrades <- quantile(nels_small$grades, .05)

newdat <- data.frame(grades=c(medGrades, fifthPercentileGrades))

pred <- predict(nels.multinom, newdat, "probs")

pred

# conditional multinomial logit 
library(MCMCpack)

N <- nrow(cola)

N3 <- N/3

price1 <- cola$price[seq(1,N,by=3)]
price2 <- cola$price[seq(2,N,by=3)]
price3 <- cola$price[seq(3,N,by=3)]

bchoice <- rep("1", N3)

for (j in 1:N3){
  if(cola$choice[3*j-1]==1) bchoice[j] <- "2"
  if(cola$choice[3*j]==1) bchoice[j] <- "3"
}

cola.clogit <- MCMCmnl(bchoice ~
                         choicevar(price1, "b2", "1")+
                         choicevar(price2, "b2", "2")+
                         choicevar(price3, "b2", "3"),
                       baseline="3", mcmc.method="IndMH")
sclogit <- summary(cola.clogit)
tabMCMC <- as.data.frame(sclogit$statistics)[,1:2]

row.names(tabMCMC)<- c("b2","b11","b12")
kable(tabMCMC, digits=4, align="c",
      caption="Conditional logit estimates for the 'cola' problem")

pPepsi <- 1
pSevenup <- 1.25
pCoke <- 1.10
b13 <- 0
b2  <- tabMCMC$Mean[1]
b11 <- tabMCMC$Mean[2]
b12 <- tabMCMC$Mean[3]

# The probability that individual i chooses Pepsi:
PiPepsi <- exp(b11+b2*pPepsi)/
  (exp(b11+b2*pPepsi)+exp(b12+b2*pSevenup)+
     exp(b13+b2*pCoke))
# The probability that individual i chooses Sevenup:
PiSevenup <- exp(b12+b2*pSevenup)/
  (exp(b11+b2*pPepsi)+exp(b12+b2*pSevenup)+
     exp(b13+b2*pCoke))
# The probability that individual i chooses Coke:
PiCoke <- 1-PiPepsi-PiSevenup


# Ordered Choice Models 
library(MCMCpack)
nels.oprobit <- MCMCoprobit(psechoice ~ grades, data = nels_small, mcmc = 10000)
sOprobit <- summary(nels.oprobit)

tabOprobit <- sOprobit$statistics[, 1:2]

kable(tabOprobit, digits=4, align="c", caption="Ordered probit estimates for the 'nels' problem")

# The probabilities for each choice
mu1 <- -tabOprobit[1] # mu1 = -(Intercept)
b <- tabOprobit[2] # beta = grades
mu2 <- tabOprobit[3]-tabOprobit[1] # gamma2 - (Intercept)
xGrade <- c(mean(nels_small$grades), quantile(nels_small$grades, 0.05))

# Probabilities:
prob1 <- pnorm(mu1-b*xGrade)
prob2 <- pnorm(mu2-b*xGrade)-pnorm(mu1-b*xGrade)
prob3 <- 1-pnorm(mu2-b*xGrade)

# Marginal effects:
Dp1DGrades <- -pnorm(mu1-b*xGrade)*b
Dp2DGrades <- (pnorm(mu1-b*xGrade)-pnorm(mu2-b*xGrade))*b
Dp3DGrades <- pnorm(mu2-b*xGrade)*b

# the marginal effect of grades on the prob of attending a 4 year college for a student with average grade 
# and for student in the top 5 %. 
average(grade)     top 5% 
-0.14322299       -0.03074883

# Models for Count Data
olympics.count <- glm( medaltot ~ log(pop) + log(gdp), 
                      family = "poisson", 
                      na.action = na.omit,
                      data = olympics)
kable(tidy(olympics.count), digits = 4, align = 'c',
      caption = "Poisson model for the 'olympics' problem")


library(AER)
dispersiontest(olympics.count)

## the Tobit or Cencored Data model
head(mroz)
hist(mroz$hours, breaks=20, col="grey")

mroz.tobit <- tobit(hours ~ educ + exper + age + kidsl6, 
                    data = mroz)
sMrozTobit <- summary(mroz.tobit)
sMrozTobit


# the marginal effect of education on hours for some given values of the regressors.
xEduc <- 12.29
xExper <- 10.63
xAge <- 42.54
xKids <- 1
bInt <- coef(mroz.tobit)[[1]]
bEduc <- coef(mroz.tobit)[[2]]
bExper <- coef(mroz.tobit)[[3]]
bAge <- coef(mroz.tobit)[[4]]
bKids <- coef(mroz.tobit)[[5]]
bSigma <- mroz.tobit$scale
Phactor <- pnorm((bInt + bEduc*xEduc + bExper*xExper + bAge *xAge + bKids*xKids) / bSigma)
DhoursDeduc <- bEduc*Phactor

DhoursDeduc

library(sampleSelection)
wage.heckit <- selection(lfp ~ age + educ + I(kids618+kidsl6) + mtr, log(wage) ~ educ + exper, 
                         data=mroz, method="ml")
summary(wage.heckit) 
