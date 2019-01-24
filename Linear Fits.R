# linear models demo

par(mfrow=c(1,1))

data(cars)

?cars
head(cars)
cars[1:5,]
plot(cars)
y<- cars[,"dist"]
x<- cars[,"speed"]
n<- nrow(cars)
n
beta1.hat <- cor(x,y) * sd(y)/sd(x)
beta0.hat<- mean(y) - beta1.hat*mean(x)
y.hat<- beta0.hat+beta1.hat*x
resids<-y-y.hat
rrs<- sum(resids^2)
abline(a=beta0.hat, b=beta1.hat, col=2)


#testing H0: beta1 = 0
t.stat<- cor(x,y)/ sqrt(1-cor(x,y)^2)* sqrt(n-2)
t.stat
lm(dist ~ speed, data = cars)
fit<-lm(dist ~ speed, data = cars)
summary(fit)
#check normality
hist(fit)
hist(resids)
qqnorm(resids)
plot(fit)


# residuals vs fitted

plot(anscombe)


data("mtcars")
?mtcars


# convert transmisson to say auto or stick instead of 0 or 1
mtcars$trans<- factor(c("auto", "stick")[mtcars$am+1])


#"factor" tells R to treat something as a categorical variable... 
#assigns levels or categories to numerical variables

plot(mtcars$hp,mtcars$mpg)

summary(lm(mpg~hp,data=mtcars))

boxplot(mtcars$mpg~mtcars$trans)

#lm() is linear model


#interaction? two ways to write it
summary(lm(mpg~trans+hp+trans:hp,data=mtcars))
summary(lm(mpg~trans*hp,data=mtcars))

summary(lm(mpg~trans+hp,data=mtcars))

fitMPG<-(lm(mpg~trans+hp+mtcars$trans:hp,data=mtcars))

fitMPG2<-(lm(mpg~trans*hp,data=mtcars))
fitMPG2
#interaction is not significantly significant on mpg
#--> but if interaction is taken out, both effectors are highly significant effectors of mpg


summary(fit)

summary(fitMPG)


# ask R to do reduction for you

# R will perform AIC--> pick lowest AIC model

AIC(fitMPG)

AIC(fitMPG2)

#stepwise progression
step(fitMPG)

#gives back the best linear fit model

fitMGPreduced<-lm(formula = mpg ~ trans + hp, data = mtcars)

summary(fitMGPreduced)


# model tumor size. throw in all variables... Then perform step() or "stepwise regression"

#^^^^linear Regression^^^^




# in general in logistic regression, beta1 represents increase in log odds of outcome per unit increase in predictor
# and that means e^beta1 represents odds ratio for unit increase




admissions <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
head(admissions)
summary(admissions)

# generalized linear model (glm)

fitGPA<-glm(admit ~ gpa, family = "binomial", data = admissions)
fitGPA

lm(fitGPA)
fitGPA
summary(fitGPA)
pairs(admissions)

# estimate standard gives gives the difference in log odds... here gpa estimate std. = 1.05 ... so for every unit of GPA(1.0)
# the fold increase = exp(1.05) ... this is a 2.86 fold increase.






