

#1 Download from Canvas the file “prestige.txt” and save it somewhere on your computer. 
#Examine it using a text editor or spreadsheet program, being very careful not to modify it. 
#Describe, in words, how the data is organized. How many observations (samples) are there? 

#The data presented in the text file is organized in a txt file. There are 102 different samples, and they are looking at 5 different categories with respect to a particular profession.


#What are the variables (measurements) named in the data file? Are they continuous or categorical?

# The variables measured are job, prestige, education, income, type, and maleDominated.
#the continuous measurments are: prestige, education, and income.
#the categorical measurments are: job, type and maleDominated.

jobSurvey <- read.csv("/Users/josephobrien/Desktop/BioStats/prestige.txt",row.names=1)


nrow(jobSurvey)
colnames(jobSurvey)
(jobSurvey$"job")



#2. Suppose your data file were tab-delimitted instead of comma-delimitted. How would you change the read.csv line above to read it in?

#In order to read tab-deliminatted data you would need to use the read.table() function rather than the read.csv


#3. Using which() and rownames(), find out which are the 9 “professional” occupations that were not male–dominated in Canada in 1971.

class(jobSurvey[,"maleDominated"])
class(jobSurvey[,"type"])

levels(jobSurvey[,"type"])
table(jobSurvey[,"type"],jobSurvey[,"maleDominated"])
prof1<-which(jobSurvey$"type"=="prof")
male1<-which(jobSurvey$"maleDominated"==FALSE)
jobSurvey$"type"=="prof" & jobSurvey$"maleDominated"==FALSE
h<-which(jobSurvey$"type"=="prof" & jobSurvey$"maleDominated"==FALSE)
h
g<- (jobSurvey[c(h),] )
g


rownames(g,)



#4 Recall what we learned about the χ2 test.

#(a) Using a hand-held calculator or computer (not R), compute the expected counts for each of the cells in the type vs. maleDominated table under the null hypothesis that job type is independent of gender balance.


mean(jobSurvey$"income")

jobSurvey$income
nrow(jobSurvey)
ncol(jobSurvey)

prof1<-which(jobSurvey$"type"=="prof")

table(jobSurvey[,"type"],jobSurvey[,"maleDominated"])

# The expected x^2 score for each test would be 
#bc: 4.6
#prof: 11.8
#wc: 8.7

#(b) Using a hand-held calculator or computer (not R), compute the χ2 statistic to test that hypothesis.

#bc:11.9
#prof: .664
#wc: 2.6
11.9+.664+2.6
#x^2 statistic = 15.164

#(c) Using your answer to (b) and the pchisq() function, test the hypothesis that the occupation type is independent of gender balance. Give your code, the p-value, and a written interpretation of your results (be sure to discuss any assumptions you have made).
?pchisq

pchisq(15.164,2)

1-pchisq(15.164,2) # pvalue = .00051, meaning we can reject the null hypothesis using a .95 confidence interval


#(d) Read the help page for the function chisq.test() and apply it to test the hypothesis above. How does it compare to your answer in (c)?

?chisq.test()

chisq.test(table(jobSurvey[,"type"],jobSurvey[,"maleDominated"]))

# this test gave an x-squared value of 10.515 and a p-value of .0052. This is different from the value that I calculated by hand, the p value is greater and x-squared is less

#(e) We learned in class that Fisher’s Exact Test is more accurate (albeit more computation- ally expensive) than the χ2-test for count data. Read the help page for the function fisher.test() and apply it here. How does it compare to your answers from (b) and (c)?

fisher.test(table(jobSurvey[,"type"],jobSurvey[,"maleDominated"]))

# this is a two sided t-test, which has a slightly greater pvalue than the chisq.test run through R, and is much greater than the .00051 value which I calculated by hand.

#5. In this problem, we’ll compare the job type designation to average education in each profession by categorizing the education data into low, medium, and high groups.



#(a) Using the quantile() command, find the 33% and 67% quantiles of education in the job- Survey data.

?quantile
quantile(jobSurvey$education, probs = seq(0,  1, .33))
# 33% = 8.9
# 66% = 11.5

#(b) Using R, create a new variable called eduCats that categorizes the average educational attainment for each job as “low”, “med”, or “hi”. Tabulate eduCats by itself. Are the counts what you expected? Why/why not?




eduCats <- factor(c(jobSurvey$education,"low","med","hi"))

?factor
"low"<-jobSurvey$education<8.9
"med"<-jobSurvey$education>=8.9 & jobSurvey$education<=11.5
"hi"<-jobSurvey$education>11.5

?rep
#eduAboveMean <- rep(NA,nrow(jobSurvey))
#eduAboveMean[jobSurvey$education<mean(jobSurvey$education)] <- "below"                    
#eduAboveMean[jobSurvey$education>=mean(jobSurvey$education)] <- "above"      

eduCats<- rep(NA,nrow(jobSurvey))

eduCats[jobSurvey$education<8.9]<- "low"

eduCats[jobSurvey$education>=8.9 & jobSurvey$education<=11.5] <- "med"

eduCats[jobSurvey$education>11.5]<- "hi"

table(eduCats)

table(eduCats, jobSurvey$type)

# The counts are each 34. This is what I expected because they are separated into thirds based on the 33% quartile that I set. There for they have to be equal numbers because the total sample population is divisible by 3.


#(c) Tabulate eduCats vs job type. Give both the R output and a version of the table that would be suitable for publication.



table(eduCats, jobSurvey$type)




#(d) Would a χ2 test be appropriate to test the hypothesis that the level of education (low, med, hi) is independent of job classification (blue collar, white collar, professional)? Why or why not?

# the x^2 test would be appropriate to use in this situation because it test the relationship between two catagorical variables. This test could easily be used on the table presented from the code, and would give insite between education and profession.


# (e) How would you describe the relationship between job type and education level? Justify your answer with the data.

# The relationship between highly educated individuals and job placement is apparent from the table. There is a higher proportion of professors with high education(29/31), and bc workers with low education (33/44) than would be found in a random sampling.

#6. Create the following boxplots and comment on any patterns you notice: 

# boxplot of educational attainment vs. job type using "~"

boxplot(jobSurvey$education~jobSurvey$type)

# let's add axis labels:

boxplot(jobSurvey$education~jobSurvey$type,
        ylab="average education (yrs)",xlab="job type")

#(a) Education vs job-type
boxplot(jobSurvey$education~jobSurvey$type,
        ylab="average education (yrs)",xlab="job type")

# there is a strong correlation between the type of profession and average years of education. For example professors have on about average 5 more years of education than bc workers.

#(b) Income vs job-type

boxplot(jobSurvey$income~jobSurvey$type,
        ylab="Income", xlab="job type")
#Professors look to have the highest income, and have the a higher income than the other two professions. Between the bc and wc workers, however, education does not seem to make any difference.
# I say this because the wc workers are typically more educated, and make almost the same income as bc workers.

#(c) Education vs maleDominated

boxplot(jobSurvey$education~jobSurvey$maleDominated,
        ylab="average education (yrs)", xlab= "Male Dominated" )
# Male dominated professions have a mean less education than non Male dominated professions. They also have a much wider range of the education that they have.


#(d) Income vs maleDominated

boxplot(jobSurvey$income~jobSurvey$maleDominated,
        ylab="Income", xlab= "Male Dominated" )
# Male dominated professions have a mean income higher than non male dominated professions. They also have a maximum income far greater than non male dominated professions.


#7

# plot of income vs education
plot(jobSurvey$income~jobSurvey$education)
# color points by jobSurvey$type:
plot(jobSurvey$income~jobSurvey$education,col=jobSurvey$type) 
# change the shape ("pch") to reflect maleDominated:
# we need to add 1 to the boolean to make the pch an integer 
plot(jobSurvey$income~jobSurvey$education,col=jobSurvey$type,
  pch=jobSurvey$maleDominated+1)
# First, load the ggplot2 package you installed:
library(ggplot2) 
# I am using an old R version, so I get a warning.

# Now, we'll tell qplot to plot income vs education, using the jobSurvey data
qplot(x=education, y=income, data=jobSurvey)

# Same thing, but using colour (note the u!) to indicate the job "type"
qplot(x=education, y=income, data=jobSurvey, colour=type)

# Let's also use shape to indicate male-dominated jobs:
qplot(x=education, y=income, data=jobSurvey, colour=type, shape=maleDominated)

#7. That last plot is a bit hard to see; let’s improve it.
#(a) Using what you learned in the first lab, how would you select the rows of jobSurvey that
#had only blue-collar jobs?

bc <- jobSurvey[which( jobSurvey$type == "bc"),]


#  (b) Use your answer to (a) to plot income vs education for blue-collar jobs only. Color the
#points according to whether or not the job was male-dominated.

qplot(x=income, y=education, data=bc, colour=maleDominated)

#(c) Do the same for white-collar and professional jobs. Are the plots directly comparable? Why or why not?

wc<- jobSurvey[which( jobSurvey$type == "wc"),]
prof<- jobSurvey[which( jobSurvey$type == "prof"),]

qplot(x=income, y=education, data=wc, colour=maleDominated)
qplot(x=income, y=education, data=prof, colour=maleDominated)
# From plotting the data, there are some similar trends: male dominated fields have greater incomes. WC workers have a similar education level in male dominated than otherwise. 
#In BC work, education looks to have a greater impact on income than WC however. There also does not seem to be a great difference between the average mean income of bc, and wc workers.
# The male professors appear to make far more income regardless of the amount of education they have received. They also have a mean income of nearly double that of wc or bc workers. 


#  (d) What is the minimum income of any job? What is the maximum?

summary(jobSurvey)
summary(bc)
summary(wc)
summary(prof)

# the Minimum income of any job type is 611, and the maximum is 25879

#  (e) In both the plot() and qplot() functions, you can change the range of the axes using xlim=c(xmin,xmax) and ylim=c(ymin,ymax) inside the plot command, where xmin, xmax, ymin, and ymax are numbers that you specify. For example, plot(1,1,ylim=c(0,20)) will plot the point (1,1) with a y axis that goes from 0 to 20. Using this information and your answer to part (d), replot the plots from (b) and (c) on the same scale. Compare the plots and describe anything you notice.
?plot
plot(1,1,ylim=c(0,20))

qplot(x=income, y=education, data=bc, colour=maleDominated, xlim=c(611,25879))

qplot(x=income, y=education, data=wc, colour=maleDominated, xlim=c(611,25879))

qplot(x=income, y=education, data=prof, colour=maleDominated, xlim=c(611,25879))

# it appears from these plots that in all cases MaleDominated professions make a greater income than non-maleDominated professions. 
#Also there appears to be a strong difference in the amount of education and income between professors and both wc and bc workers.
# In general it appears that wc workers have more education than bc workers, but do not make a greater income.

#8 Another useful exploratory tool is the pairs() function, which plots the columns of a data frame against each other. Try pairs(jobSurvey).

#(a) Which panel shows the relationship between job prestige and income? How would you describe that relationship in words? Does moving into a higher-income profession also increase prestige?
pairs(jobSurvey)


# The panel in row 1, column 3 shows the relationship between prestige and income. It does appear that they are directly correlated, as in, people with higher prestige also have a higher income.

# it also appears that prestige is correlated with job type as well. This is shown from having one catagorical group much higher on the plot than the other two. 


#(b) Based on the pairs plot, which do you think have the largest correlation: prestige and income, prestige and education, or education and income? Which do you expect will be least correlated? Why?

# prestige v income, prestige v education, education v income

#it appears as if they all share some correlation based on the plots. It looks to me as if the prestige and education share the strongest correlation bc the slope looks as if the slope is closest to 1. Prestige and income look to share the weakest correlation. 

#(c) Read the help page for the cor() function and compute the correlations r and rank corre- lations ρ for the three pairs of variables in (b). Do they confirm your intuitions?

?cor()
cor(jobSurvey$prestige,jobSurvey$income)
cor(jobSurvey$prestige,jobSurvey$education)
cor(jobSurvey$education,jobSurvey$income)
# prestige, education: .85
# prestige, income: .71
# education, income: .58

# I was correct in saying that prestige and education are the most closely related because of the highest correlation factor.
# I was incorrect in saying that prestige and income have the least correlation however, the least correlation is between education and income.

#9 Modify the example above to test the following hypotheses. For each test, explain in words your interpretation of the output.

# test differences in mean income as a function of maleDominated
t.test(jobSurvey$income~jobSurvey$maleDominated)

#(a) The null hypothesis that average income is the same in maleDominated professions vs the alternative hypothesis that maleDominated professions have higher average income.

?t.test
md<- jobSurvey[which( jobSurvey$maleDominated == TRUE),]

md$income

t.test( md$income, jobSurvey$income, alternative="greater")

# the one sided t.test shows that we can reject null hypothesis. The true mean is greater than 0.

#(b) The null hypothesis that average education doesn’t differ between maleDominated and non- maleDominated blue collar professions, vs the alternative hypothesis that they are unequal in blue collar professions.

md<- jobSurvey[which( jobSurvey$maleDominated == TRUE),]

md.bc <- md[which(md$type == "bc"),]

md.bc

nrow(md.bc)
nrow(md)
md
bc
nrow(bc)
nmd.bc<-(bc[ bc$maleDominated==FALSE, ])

t.test( md.bc$education,bc$education)
t.test( nmd.bc$education, bc$education)

# The t test tells us that we can reject the null hypothesis and that the education does differ in Male Dominated blue collar professions vs the all male dominated professions.


#(c) The null hypothesis that average education is the same for blue collar vs other workers (white collar and professional combined).



wc.prof<-(jobSurvey[ jobSurvey$type != "bc", ])
summary(wc.prof)

wc.prof<-wc.prof[complete.cases(wc.prof),]

wc.prof


t.test( wc.prof$education,bc$education)

# based on this t test we can reject the null hypothesis that the average education is the same for blue collar workers as it is for all other workers.

#(d) The null hypothesis that there is an average education difference of four years between blue collar and other workers (white collar and professional combined), vs a null hypothesis that the average education difference is greater than four years.

?t.test

t.test(wc.prof$education, bc$education,"greater" ,4)

# I do not understand the question. I do not see how there can be 2 null hypotheses. I am assumming that there was an error and that the alternative hypothesis is that there is a greater than four year difference. 

# This t test shows that we cannot reject the null hypothesis the the mean difference of education is 4 years or less because the p-values is .095. 


#(e) We just conducted four hypothesis tests. Are they independent? Why or why not?



# they are not all independent of one another because they take from the same data sets and only switch variables. 


#10 Consider the relationship between prestige and education.
#(a) Plot prestige vs education. Would a linear model be appropriate here? Why or why not?

qplot(x=prestige, y=education, data=jobSurvey)
qqplot(jobSurvey$prestige,jobSurvey$education)

# I believe a linear model would work in this instance because the two samples are directly correlated. The plot given shows that there could be a slope drawn that would represent the correlation between the two variables. 
# this is further supported by a normal looking qqplot that follows a linear trend.

#(b) By eye, sketch a best-fit line on the graph and estimate its slope and intercept from your sketch. Show your work.

# From my personal sketch, I estimate that the slope is approximatly .1, and that the y-int is 4.2. The line is attached to an additional PDF.

# y= (.1) x + 4.2 


#(c) Using the formulae we learned in class, compute the estimates βˆ and βˆ for the linear
# model prestigei = β0 + β1educationi + εi. Show your work.



#prestigei = β0 + β1educationi + εi

# we want to predict prestige as a function of education

#εi->0

ep<-cor(jobSurvey$prestige,jobSurvey$education)

#β1ˆ = cor(x,y)(sdx/sdy)

pres<-sd(jobSurvey$prestige)
pres
edu<-sd(jobSurvey$education)
edu

edu.pres<- edu/pres
edu.pres

m <- ep*edu.pres
m

mean.pres<-mean(jobSurvey$prestige)

mean.edu<-mean(jobSurvey$education)

b<- mean.edu - mean.pres*m

b

# Using the formulas given in class, I estimate that the βˆ1 is 0.13. I also estimate that the βˆ0 is 4.4.  

#(d) How close were the estimates in (c) to your by-eye guesses in (b)?

# part B: slope = .10 , y-int = 4.2
# part C: slope = 1.3, y-int = 4.4
# I would trust the answers taken in part C as more accurate. They are close, but in part B I approximated from what I saw, and it was not exact.


#(e) Using the formulae we learned in class, compute a test statistic to test H0 : β1 = 0 versus HA :β1 ̸=0

r<- edu/pres
r

n<- nrow(jobSurvey)
n

t.b1 <- (r/(sqrt(1-r^2)))*(sqrt(n-2))
t.b1

#test statistic = 1.6

#(f) To what distribution should you compare that test statistic? Do it and report the p value you obtain. Show your work!

# y = mx+b ------- dist^ = β1ˆx + β0ˆ

#t-stat = 1.6, df = 102-2 = 100

# two-sided pvalue

?pt

ht <- pt(1.6, 100,,lower.tail = FALSE)
ht
lt<- pt(-1.6,100,,lower.tail = TRUE)
lt
twotails<- ht+lt
twotails

# p value = .11


y<- jobSurvey[,"education"]
x<- jobSurvey[,"prestige"]
n<- nrow(jobSurvey)
beta1.hat <- cor(x,y) * sd(y)/sd(x)
beta0.hat<- mean(y) - beta1.hat*mean(x)
y.hat<- beta0.hat+beta1.hat*x
resids<-y-y.hat
rrs<- sum(resids^2)
rrs

se<- sqrt(rrs/100)
se

sx<- sd(jobSurvey$prestige)
sx


high.95<- m + qt (1 - (1.96/2), df = 100)*(se/(sx*sqrt(101)))
high.95
low.95 <- m - qt (1 - (1.96/2), df = 100)*(se/(sx*sqrt(101)))
low.95

high.95
low.95

# The t statistic should be compared in the 95% confidence interval. This would allow the distribution to be between .12 and .15 range of the t distribution. The p value returned is .11.

#11 Here we’ll examine the fit obtained in the above example.
# Construct a linear fit of prestige vs education.
# Using explicit naming of the data frame variables:
fit <- lm(jobSurvey$education ~ jobSurvey$prestige)
# Using just the column names & specifying the dataframe with data= 
fit <- lm(prestige ~ education, data=jobSurvey)

#(a) Use the summary() command to obtain βˆ0 and βˆ1 . Does R’s estimate agree with your 10
#estimate in the previous problem? What about the p-value for H0 : β1 = 0

summary(lm(jobSurvey$education ~ jobSurvey$prestige))

# y-int: 4.4
# slope: .13
#p.value: < 2.2e-16
# The y-int and the slope both agreeed with my estimates in the previous problem. The p value is much lower in the lm model.


#(b) Use the residuals() command to obtain the model residuals. How do we expect the residuals to be distributed? Plot a histogram of the residuals to check.

lm.resids<-residuals(lm(jobSurvey$education ~ jobSurvey$prestige))
# we would expect the residuals to be evenly distributed around the mean (0).
hist(lm.resids)
# The histagram shows that my assumption was correct.

#(c) “Statistician” was not listed amongst the professions in the survey. Suppose that in 1971 in Canada, the average statistician had 15 years of education, on par with other scientists. Based on your model, what do you think the Pineo-Porter prestige score for statisticians would have been?

# y = mx + b

# educationi = (.13)*prestigei + (4.4)

#prestigei = (educationi - 4.4) / .13

(15-4.4)/.13
# The prestige score of a statistician would have been approximatly prestige 81.5

summary(lm(jobSurvey$education ~ jobSurvey$prestige))



#(d) Suppose we wanted to examine the relationship between prestige and education only amongst blue- and white-collar occupations (excluding professional occupations). 
#How would we modify the lm() command above? 
#Produce another fit with this model. How does it compare to the fit with all data?

# In order to exclude the professional occupations, we would modify the data set to exclude those entries.

wc.bc1<-(jobSurvey[ jobSurvey$type != "prof", ])

wc.bc1<-wc.bc1[complete.cases(wc.bc1),]

summary(lm(wc.bc1$education ~ wc.bc1$prestige))
summary(lm(jobSurvey$education ~ jobSurvey$prestige))

cor(jobSurvey$prestige,jobSurvey$education)

cor(wc.bc1$prestige,wc.bc1$education)

# The y-int increased in the new sample, and the slope decreased. 
#For these two subgroups of the jobSurvey data, there is less of a correlation between education and prestige (.59) than in the jobSurvey data as a whole(.85).






