data(iris)
#1. Explore the iris data:
#(a) What variables are present? Are they continuous or categorical?

iris
# The variables are: Sepal.Length, Sepal.Width, Petal.Length, Petal.Width, and Species

# All are continuous except for "Species", which is categorical.




#(b) For the categorical variables, provide a table of how many samples there are in each category. Show the R code you used.

cat.var<- table(iris$"Species")

summary(cat.var)
cat.var

#(c) Provide the range (minimum, maximum) of the continuous variables.

summary(iris)
#Sepal.Length: min: 4.3, max: 7.9
#Sepal.Width: min: 2.0, max: 4.4
#Petal.Length: min: 1.0, max: 6.9
#Petal.Width: min: .1, max: 2.5

#(d) Produce a scatterplot matrix (a pairs plot) of the continuous variables.

cont.data<-  iris[ -c(5) ]
cont.data
plot(cont.data)

#2. Do iris flower sizes depend on species?
#(a) Produce boxplots of sepal length broken down by iris species.

boxplot(iris$Sepal.Length~iris$Species)



#(b) For each iris species, compute the mean sepal length. (Show your work.)

seto<-(iris[ iris$"Species" == "setosa", ])
vers<-(iris[ iris$Species == "versicolor", ])
virg<-(iris[ iris$Species == "virginica",])


seto<-seto[complete.cases(seto),]
summary(seto)
# setosa mean sepal length = 5.00
summary(vers)
# versicolor mean sepal length = 5.94
summary(virg)
# virginica mean sepal length = 6.59


#(c) For each iris species, compute the sample SD of the sepal length. (Show your work.)

setosl<- sd(seto$Sepal.Length)
setosl
#setosa sd= .35

verssl<- sd(vers$Sepal.Length)
verssl
#versicolor sd=.52

virgsl<- sd(virg$Sepal.Length)
virgsl
#virginica sd=.63

nrow(seto)



#(d) Using your answers above, fill out a complete ANOVA table for sepal length vs. species. (Show your work.)

Y <- iris[,"Sepal.Length"]
Y

ybar<- mean(iris$Sepal.Length)
ybar
setoi <- seto[,"Sepal.Length"]
versi<- vers[,"Sepal.Length"]
virgi<-virg[,"Sepal.Length"]
n <- nrow(seto)


#  SSb = 63.9

setossb<-((5-ybar)^2)*50
versssb<- ((5.94-ybar)^2)*50
virgssb<-((6.59-ybar)^2)*50

SSb.1<-setossb+versssb+virgssb
SSb.1
#  SSw = 6.86
seto<-na.omit(seto)
setossw<- (sum(setoi-5)^2)*(nrow(seto)-1)
setossw

versssw<-(sum(versi-5.94)^2)*(nrow(vers)-1)
versssw
virgssw<-(sum(virgi-6.59)^2)*(nrow(virg)-1)
virgssw
SSw.1<- setossw+versssw+virgssw
SSw.1
SStot<-SSw.1+SSb.1
#  SStotal  = 6.86+63.9 = 70.7
SStot

var(iris$Sepal.Length)


# DFb = K-1 = 2
K<-3
DFb.1<- K-1



# DFw = N - K = 147
N<- 150
DFw.1<- N-K

# DFtotal  = N - 1 = 149
DFtot<- 149


# MSb 
setomsb.1<- (sum(setoi-5)^2)*(1/(3-1))
setomsb.1
versmsb.1<-(sum(versi-5.94)^2)*(1/(3-1))
versmsb.1
virgmsb.1<- (sum(virgi-6.59)^2)*(1/(3-1))
virgmsb.1


MSb.1<- SSb.1/2


MSb.1
#MSb = 31.95



#MSw = MStot-MSb = 0.04666667
MSw.1<- SSw.1/147
MSw.1




MSb.1+MSw.1
# MStotal = 31.99833


MStot.1<-var(iris$Sepal.Length)


# F value = 684.68

f.val <-MSb.1/MSw.1


f.val

#(e) Based on your table, and using the “p*” family of functions (such as pnorm and pt), compute a p value for the ANOVA. Also provide, in words, an interpretation of the result.

?pnorm
pnorm(f.val)
pt(f.val, DFtot)
?pt
# Based on the table that I put together, the p value = 1. With this value, there is no way that the species will effect the sepal length.
# I believe that I did this wrong.



#(f) In class, we learned about using “anova(lm(Y~X,data=myData))” to obtain an ANOVA analysis of Y versus the categorical variable X using the dataframe “myData”.
#Use this syntax to perform the ANOVA analysis you did by hand above. Do the results agree with yours?
(anova(lm(iris$Sepal.Length~iris$Species,data=iris)))
pf(119.26, 2, 147)

# The anova table agreed with my DF calculations and my sum of squares calculations. 
#It also agreed with my mean sum of squares between groups.
# The table did not agree with my Mean sum of squares within groups, or my F value. 
#It also did not agree with my p value. For the p value, it is close to 0, while mine was 1.




#(g) Repeat part (f) to perform three more ANOVAs for the sepal width, petal length, and petal width.
anova(lm(iris$Sepal.Width~iris$Species,data=iris))
anova(lm(iris$Petal.Length~iris$Species,data=iris))
anova(lm(iris$Petal.Width~iris$Species,data = iris))

#(h) At this point, we have conducted multiple tests of whether the flower measurements depend on species. 
#What would our p-value threshold need to be to ensure a FWER (family-wise error rate) of 0.05?
#Which tests meet that threshold?

p.adjust(2.2e-16,n = length(iris$Sepal.Width~iris$Species),"bonferroni")

?p.adjust()

# 6.6*10^-16


#3. ANOVA assumes that the groups have equal variances. Let’s check that assumption.

#(a) Consider the boxplots of sepal length broken down by iris species that you produced in part (a) of the previous problem. By eye, do you believe the variances are equal?

# by eye, I do not believe that the variences are equal. I looks to me the viginica has the largest varience and setosa has the smallest.



#(b) Using your answers to part (c) of the previous problem, compute a statistic to test whether the variance of the setosa and virginica sepal lengths are equal. 
#To what distribution would we compare this test statistic? Obtain the associated p-value and provide an interpretation (in words) of the result.


plot(var(seto$Sepal.Length)~var(virg$Sepal.Length))
setosl1<- complete.cases(setosl)
virgsl1<- complete.cases(virgsl)

var.test(seto$Sepal.Length, virg$Sepal.Length, alternative = "two.sided")
#F value: .30
#p value: 6.366e-05

# There is a 95% confidence interval between .174 and .541 on in an F distribution. This is where we would compare the F statistic
# the p value is low, the allows me to reject the null hypothesis that the variances are equal in the two groups.



#(c) In light of your result from part (b), does your interpretation of the ANOVA analysis in the previous problem change?
#If your scientific question was whether sepal length depends on species and you found that the ANOVA equal-variance assumption was violated, 
#how else might you analyze the data to answer your question? Carry out that analysis and provide a written interpretation of your results.

# Yes this does change my interpretation of the ANOVA analysis. The ANOVA test assumes equal variance, meaning that the answer and p value attained from that table is likely incorrect.
# I can analyze my data using the kruskal wallis test, which does not assume equal variance.

kruskal.test(seto$Sepal.Length ~ virg$Sepal.Length, data = iris)

# The kruskal.test returned a p value of .7177 and chi-squared result of 16. The p value of this test tells me that I cannot reject the null hypothesis that the medians of these two groups are not equal.

#(d) Is the test used in part (b) parametric or non-parametric? What assumptions are being made for that test?

# The variance test is a parametric test. There is an assumption that the populations are normally distributed.

#(e) Suppose the assumptions required by the test in (b) are violated.
#Describe in words how you could test for the equality of variances using resampling. 
#You may (optionally) provide R code to carry it out (hint: use the replicate() and sample() functions).


seto.samp<- rep(seto$Sepal.Length, sample, 10000)
hist(seto.samp)
var(seto.samp)

virg.samp<- rep(virg$Sepal.Length, sample, 10000)
hist(virg.samp)
var(virg.samp)

var.test(seto.samp, virg.samp)

# in order to mimic a normal distribution with an un normal dataset, repeat the sampling many times.
#With repeted sampling the data will tend to go toward a gaussian curve and appear normally distributed,


#4. Next we’ll consider the relationship between sepal length and sepal width.
#(a) Refer to the scatterplot matrix you produced in problem 1(d). 
#Do you expect sepal length and sepal width to be positively correlated, negatively correlated, or uncorrelated?

# Refering to the scatter plots from question 1(d), I do not expect the sepal.length and sepal.width to be correlated.

plot(iris)

#(b) Compute the correlation between sepal length and sepal width. Are your expectations confirmed?

cor.test(iris$Sepal.Length,iris$Sepal.Width)

# using person's product-moment correlation, the correlation returned -.118. 
# this correlation was not significant though, because the p value returned .152

#(c) Now compute the correlation between sepal length and sepal width separately for each of the iris species. 
#How do these compare to the overall correlation you observed in part (b)?

cor.test(seto$Sepal.Length, seto$Sepal.Width)
# cor: .743
# p value: 6.71e-10
cor.test(vers$Sepal.Length, vers$Sepal.Width)
# cor: .526
# p value: 8.772e-05
cor.test(virg$Sepal.Length, virg$Sepal.Width)
# cor: .457
# p value.000844

# Each of these taken as individuals are positivly correlated, and very significant based on the p value. 
# This is compared to all of them taken at once, which showed a negative correlation that was insignificant.

#(d) Write down a linear model (ie, an equation of the form Yi = β0 + β1Xi + ··· + εi) 
# that predicts sepal length as a function only of sepal width.

sepal.lenght.bar<-mean(iris$Sepal.Length)

sepal.width.bar<-mean(iris$Sepal.Width)
y.int <- sepal.width.bar*-0.118 -  sepal.lenght.bar
y.int


seto.cor<- cor(seto$Sepal.Length,seto$Sepal.Width)
vers.cor<- cor(vers$Sepal.Length,vers$Sepal.Width)
virg.cor<-cor(virg$Sepal.Length,virg$Sepal.Width)

plot(iris)
sep.wid<- iris$Sepal.Width
sep.len<- iris$Sepal.Length
n<- nrow(iris)
n
beta1.hat <- cor(sep.wid,sep.len) * sd(sep.len)/sd(sep.wid)
beta0.hat<- mean(sep.len) - beta1.hat*mean(sep.wid)
y.hat<- beta0.hat+beta1.hat*sep.wid

beta1.hat
beta0.hat

abline(a=beta0.hat, b=beta1.hat, col=2)


#y.int = 6.53

# Yi = 6.53 + beta1.hat(x) + εi



#(e) Using R’s lm() function, fit the model you wrote in part (d) and assign it to fit1.


fit1<- (lm(sep.len~sep.wid))
plot(fit1)
summary(fit1)
?lm()

# the y int returned is 3.42, the same as I got from my model in part (d)
# from this fit, we have a p value of .152 and an F statistic o 2.07

#(f) Repeat part (d), but now include an additional term to model the difference in average
#sepal length between species. Then repeat part (e), assigning it to fit2.

# # Yi = 3.42 + beta1.hat(x1) + εi

fit2<- (lm(sep.len~sep.wid:iris$Species))
fit2

#(g) As above, but now include the interaction between sepal width and species. Assign this to
#fit3.


fit3<- lm(Sepal.Length ~ Sepal.Width:Species+Sepal.Length:Species, data = iris)
fit3
# # Yi = 3.42 + beta1.hat(x1)+ beta2.hat(x2)+ εi



#(h) Print out the summary of fit3 and interpret the results. What are the results telling you? 
#Relate the numbers there to the β’s in your mathematical model. 
#Finally, write down equations for the best-fit line for sepal width vs length in each of the three species 
#(ie, three different line equations) with the slopes and intercepts from your model.

summary(fit3)
#based on the results in the table below, the p values for the sepal width in all species were not significant, 
#and we cannot reject the null hypothesis that the sepal width has no effect on the sepal length for any species. 
# the sepal length/species interaction did have an extremely low p value, meaning that it was very significant in relation
# to overall sepal length

#(i) predict(fit1) will compute the estimated sepal lengths (yˆ) for your data based on the i
#model you fit in fit1. Mathematically, how are these related to the residuals you’d get from residuals(fit1)? Check this using R.

summary(fit1)
# The residuals equal the obserbed - fit. The residuals are the amount of error that is between the observed and the predicted fit data.


predict(fit1)
residuals(fit1)
residuals(fit1)==(iris$Sepal.Length~iris$Sepal.Width - predict(fit1))

# the data returned showed that this was not the case in my equation. This may be because the predicted fit did not take into
# account the different subset of species, making it unreliabe

#(j) Ideally, the estimated sepal lengths (yˆ) and the observed sepal lengths y should be corre- ii
#lated. Check this by plotting the true sepal lengths vs. the predicted lengths from fits 1-3.
#We can use the squared correlation between them as a measure of how well our model fits.
#Using R, compute the squared correlation between y and yˆ for each of the three models. ii
#Which gives the best R2?

plot(iris$Sepal.Length,predict(fit1))
plot(iris$Sepal.Length,predict(fit2))
plot(iris$Sepal.Length,predict(fit3))

summary(lm(iris$Sepal.Length~predict(fit1)))# R2: .014
summary(lm(iris$Sepal.Length~predict(fit2)))# R2: 0.723
summary(lm(iris$Sepal.Length~predict(fit3)))# R2: 1.00
# from the plots shown, clearly fit3 gives the best R2. The R2 value is 1.

#(k) Look at the summaries for fit1, fit2, and fit3. The “Multiple R-squared” value gives the
#square of the correlation between the observed dependent variable yi (in this case, sepal
#length) and the estimated yˆ from your model. Are the ones you obtained in part (j) the i
#same as what R gives?

summary(fit1)# R2:.014
summary(fit2)# R2:.723
summary(fit3)# R2: 1.00

# the summary R2 values are all the exact same as teh ones from the linear model.

cor.test(iris$Sepal.Length, predict(fit1))
cor.test(iris$Sepal.Length, predict(fit2))
cor.test(iris$Sepal.Length, predict(fit3))


#(l) Generally, adding more covariates will tend to produce a higher R2. 
#That doesn’t necessarily mean it’s a better model, however – you may be overfitting the data! 
#A better way to select a model is to test if the variance explained by the larger model compensates for the degrees of freedom you introduce by adding covariates. 
#Try anova(fit1,fit2). Does fit2 account for significantly more variance than fit1? What about fit3? 
#Based on the ANOVA comparisons, which model would you choose?

anova(fit1,fit2)

# fit two does account for significantly higher variance compared to fit1. The F value returned is much higher, and the p value is extremely low.


anova(fit1,fit3)
# Fit 3 adds very little variance to fit 1. The F value is returned at 2.07 and hte p value at .15. 
# I would choose fit3 to compare to because it would add almost no variance to my sample.


#
#5. Here we’ll consider a model for predicting wide sepals.
#(a) In the iris dataset, what are the odds that a given specimen has a sepal width greater than 3cm?

quantile(iris$Sepal.Width)
up.3<-sep.wid>3
sum(up.3)
prob.g.3<-sum(up.3)/150
# there is a .44 chance that a given specimen has a sepal width greater than 3




#(b) Suppose the fraction of irises with sepal widths > 3cm in the wild is the same as it is in our sample. 
#If I were to pick a dozen irises at random, what is the probability that at least 8 would have sepals wider than 3cm? 
#(Show your work.)

(prob.g.3)^8

dbinom(8, size=12, prob=0.44)
?dbinom
# there is a .068 probability that you would pick 8/12 sepals wider than 3 cm.


#(c) Test the hypothesis that having a sepal width greater than 3cm is independent of species. 
#Show any relevant tables and justify your choice of test.

t.test(seto$Sepal.Width, iris$Sepal.Width, alternative = "greater", .22)# p value = .011
t.test(vers$Sepal.Width, iris$Sepal.Width, alternative = "greater", .22)# p value = 1.00
t.test(virg$Sepal.Width, iris$Sepal.Width, alternative = "greater", .22) # p value = .00071
mean(seto$Sepal.Width)
mean(vers$Sepal.Width)
mean(virg$Sepal.Width)

# through performing t.tests on each of the species, setosa is the only species with a significant p value. This means
# setosa is the only speices that we can reject the null hypothesis that the mean that that the sepal width is dependent on the species

#(d) Now we’ll fit a logistic regression model. In R, this is done with the glm() (“generalized linear model”) function. 
#Because there are many types of generalized linear models, 
#we need to tell it the model family to use; in the case of logistic regression, that’s a binomial family.
#Consider the following regression:
# WideSepalFit <- glm( (Sepal.Width>3)~Petal.Width+Species, family=binomial, data=iris) 
#In this model, what is the outcome (dependent) variable? What are the independent covariates/predictors?

WideSepalFit<- glm( (Sepal.Width>3)~Petal.Width+Species, family=binomial, data=iris) 
iris$Species
# The dependent variable is the Sepal width. The independent variables are the species and the petal width

#(e) Look at summary(WideSepalFit). Describe how you interpret the table of coefficients that it spits out.

summary(WideSepalFit)

# when plotted in a linear model, the Estimate value represents the y int for each subset listed. the standard error is a measure 
# of how far away the subset point is typically away from the "fit" or predicted line of the slope. The p values are the level of significance each 
# group has on the dependent variable, and the z value is the point on a gaussian model where the sample population has shifted from the mean.

#(f) Based on your output, what is the odds ratio for wide sepals in versicolor vs. setosa?
exp(coef(WideSepalFit))

#(versicolor) 1.68*10^-5: (setosa) 1.20

#(g) Based on your output, what is the odds ratio for wide sepals in versicolor vs. virginica?

# (versicolor)1.681064e-05 : (virginica) 4.959749e-07 
vers.o<- 1.681064e-05

virg.o<-4.959749e-07 

#(h) Based on your output, by how much would one need to increase the petal width to double the odds 
#that the sepal width would be greater than 3cm?





#(i) What is the probability that a virginica iris with a petal width of 2cm has sepals wider than 3cm?

sep.w.pet.w.fit<- glm( (virg$Sepal.Width>3)~virg$Petal.Width==2, family=binomial) 

summary(sep.w.pet.w.fit)

# there is a .04 probability that the sepal width would be greater than 3 cm.

#6) Finally, we’ll try to predict species based on sepal and petal widths.
#(a) Create a scatter plot of sepal width vs. petal width, colored by species.


plot(iris$Sepal.Width,iris$Petal.Width, col = (iris$Species))
                                        

#(b) Using logistic regression, fit a model that will allow you to estimate the probability that an iris with a 
#3cm wide sepal and a 1.5cm wide petal is species virginica (vs any other species).




#(c) What do you obtain for the estimated probability? Refer to the graph you created in (a); 
#does your probability estimate make sense? Why or why not?


