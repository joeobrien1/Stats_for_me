data("iris")
?iris
head(iris)
fitPetalSpecies<- lm(Petal.Length~Species,iris)
summary(fitPetalSpecies)
# tells us average length longer than setosa


# does species matter

anova(fitPetalSpecies)

# must make lm as input to anova

#SS shows small variance w in groups, large btw groups

# F value: high F value -> low p value -> highly significant