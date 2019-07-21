install.packages("faraway")

library(faraway)

data(gavote)

# Data Exloration
head(gavote)

# create undercount of voters
gavote$undercount <- (gavote$ballots-gavote$votes)/gavote$ballots

summary(gavote$undercount)

# see the hist
hist(gavote$undercount, main = "Undercount", xlab = "Percent Undercount")

# see the pie chart & barplot
pie(table(gavote$equip), col=gray(0:4/4))
barplot(sort(table(gavote$equip),decreasing=TRUE),las=2)

# proportion voting for Gore relates to the proportion of African Americans
gavote$pergore <- gavote$gore/gavote$votes

plot(pergore ~ perAA, gavote, xlab="Proportion African American", ylab="Proportion for Gore")


