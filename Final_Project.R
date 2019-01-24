data(package = .packages(all.available = TRUE))


data(package = .packages(all.available = TRUE), lung)

data(package = .packages(all.available = TRUE), myeloid)

View(myeloid)

qqline(lo,log.time)
?qqline

hist(lung$ph.ecog)

male<- lung[which(lung$sex==2),]

qplot(lung$time, lung$ph.ecog)

hist(lung$pat.karno)

nrow(lung)

hist(lung$meal.cal)
hist(lung$wt.loss)

time.c<- complete.cases(lung$time) 
ph.ecog.c<-complete.cases(lung$ph.ecog)
wtl<-complete.cases(lung$wt.loss)
meal.cal.c<- complete.cases(lung$meal.cal)

cor(meal.cal.c,lung$time)

cor(wt.loss.c, lung$time)

cor(lung$time,lung$ph.ecog)

cor(lung$status, lung$time)

cor(lung$sex,lung$time)
pairs(lung)

cor(time.c, ph.ecog.c)
pair(lung$time)
cor(lung$time,lung$time^2)
cor(ph.ecog.c, lung$time)
cor(pat.k.c, lung$time)
cor(lung$pat.karno,lung$time)


wtl<-complete.cases(lung$wt.loss)
cor( pat.k.c,ph.k.c)
ph.k.c<-complete.cases(lung$ph.karno)
pat.k.c<- complete.cases(lung$pat.karno)
cor(wtl,time.c)
plot(lung$time, lung$wt.loss)


lung1<- lung[complete.cases(lung),]

lung1

boxplot(lung$time~lung$status==2)

cor.test(lung$time,lung$wt.loss)

cor.test(lung$time, lung$sex)

cor.test(lung$time, lung$status)

qqplot(lung$time, lung$status)

cor.test(lung$time,lung$age)

cor.test(lung$time, lung$meal.cal)

cor.test(lung$time, lung$ph.ecog)

cor.test(lung$time, lung$ph.karno)

cor.test(lung$time, lung$pat.karno)

wilcox.test(lung$time, lung$pat.karno)

wilcox.test(lung$time, lung$status, paired = TRUE)
?wilcox.test
shapiro.test(lung$time)
shapiro.test(lung$age)
hist(lung$wt.loss)

hist()


hist(lung$time)
qqplot(lung$time,lung$time)
qqplot(lung$wt.loss,lung$time)
qplot(lung$time,lung$wt.loss)


cor(lung$time,lung$meal.cal, method = "pearson")

time.o<-na.omit(lung$time)
meal.cal.o<-na.omit(lung$meal.cal)
cor(time.o, meal.cal.o)

boxplot(lung$time~lung$sex)

boxplot(lung$time~lung$status)

hist(lung$age)

chisq.test(lung$meal.cal,lung$wt.loss)
?chisq.test()
nrow(lung)

qqplot(lung$time, lung$age)

wilcox.test(lung1$wt.loss, lung1$time)






plot(lung$time, lung$meal.cal)

plot(lung$time, lung$wt.loss)

nrow(lung1)

qqplot(lung$time, lung$status)
shapiro.test(lung$time)
nrow(lung)
lung$age
wtl
wilcox.test(lung$time~lung$age)
wilcox.test( lung$time ~ lung$sex)
wilcox.test(lung$time~lung$status)
wilcox.test()
mean(lung1$time)

wilcox.test(lung1$age,309.9341)
wilcox.test(lung1$wt.loss,309.9341)
wilcox.test(lung1$meal.cal,309.9341)
wilcox.test(lung1$time,309.9341 )
wilcox.test(lung1$ph.ecog,309.9341)
wilcox.test(lung1$time~lung1$sex)
boxplot(lung$time~ lung$status)

wilcox.test(lung$time,lung$wt.loss, conf.int = .95)

wilcox.test(lung$wt.loss, mean(lung$time), conf.int = .95)

cor.test(lung$time,lung$wt.loss)

cor.test(lung$time,-lung$status)



quantile(lung1$meal.cal, probs = seq(0,  1, .33))
#low - <825
#med - >=825, <= 1068.4
#high - > 1068.4

meal.cal.b <- factor(c(jobSurvey$education,"low","med","hi"))
"low"<-lung1$meal.cal<825
"med"<-lung1$meal.cal>=825 & jobSurvey$education<=1068.4
"hi"<-jobSurvey$education>1068.4

meal.cal.b<- rep(NA,nrow(lung1))

meal.cal.b[lung1$meal.cal<825]<- "low"

meal.cal.b[lung1$meal.cal>=825 & lung1$meal.cal<=1068.4] <- "med"

meal.cal.b[lung1$meal.cal>1068.4]<- "hi"

table(meal.cal.b, lung1$status,header=T)





boxplot(meal.cal.b~lung1$time)


wilcox.test(lung1$status, mean(lung1$time), paired = FALSE, conf.int = .95)


wilcox.test(lung$sex, mean(lung$time), paired = FALSE)
wilcox.test(lung$age,mean(lung$time),paired = FALSE)
wilcox.test(lung$wt.loss,mean(lung$time),paired = FALSE)
wilcox.test(lung$meal.cal,mean(lung$time),paired = FALSE)
boxplot(lung$time~lung$status)
shapiro.test(lung$time)

qqplot(lung$time,lung$wt.loss)
cor.test(lung$time,lung$status, method = "spearman")
?t.test
t.test(lung1$wt.loss,lung1$time,alternative = "two.sided",0,paired = TRUE)

?barplot()
?plot

qqplot(lung$time, lung$meal.cal)
t.test(lung["time"], lung["wt.loss"])

summary(lm(lung$time~lung$meal.cal+ lung$wt.loss))
fit.lung.time<-lm(formula = lung$time ~ lung$meal.cal + lung$wt.loss + lung$status, data = lung)

?boxplot

wilcox.test(lung$status,mean(lung$time),alternative = "less")

cor.test(lung$time,lung$status)

log.time.meal<-summary(lm(log.time~log.meal))

log.time.meal


plot(height, width = 1, space = NULL,
        names.arg = NULL, legend.text = NULL, beside = FALSE,
        horiz = FALSE, density = NULL, angle = 45,
        col = NULL, border = par("fg"),
        main = NULL, sub = NULL, xlab = NULL, ylab = NULL,
        xlim = NULL, ylim = NULL, xpd = TRUE, log = "",
        axes = TRUE, axisnames = TRUE,
        cex.axis = par("cex.axis"), cex.names = par("cex.axis"),
        inside = TRUE, plot = TRUE, axis.lty = 0, offset = 0,
        add = FALSE, args.legend = NULL, ...)



summary(log.time.meal)

kruskal.test(lung$time, lung$status)

cor.test(lung$time,lung$status)

?lm()

cor.test(log.time, log.meal, alternative = "greater")


boxplot(lung$time~lung$status)

summary(boxplot(lung$time~lung$status))

boxplot(lung$time~lung$status,
        ylab="Time(Days)",xlab="Anaplasia Status")
qqplot(lung$time, lung$wt.loss,
       xlab="Time(Days)",ylab="Weight loss (lbs.)", col())
qqplot(lung$time, lung$meal.cal,
       xlab="Time(Days)",ylab="Calories per Day")


quantile(lung$time, probs = seq(0,  1, .25))

cor.test(lung$time, lung$age, method = "kendall")

cor.test(lung$time, lung$sex, method = "kendall")

cor.test(lung$time, lung$status, method = "kendall")

ana.1<- lung[which(lung$status==1),]
ana.2<- lung[which(lung$status==2),]
quantile(ana.1$time, probs = seq(0,  1, .25))
quantile(ana.2$time, probs = seq(0,  1, .25))

nrow(lung)

log.meal<-log(lung$meal.cal)

qqplot(lung$time, log.meal)

log.time<- log(lung$time)
qqplot(log.time, log.meal)
qqline(log.time, log.meal)
qqplot(log.time, log.meal,
       xlab="log(Time(Days))",ylab="log(Calories per Day)")
plot(log.time, log.meal,
     xlab="log(Time(Days))",ylab="log(Calories per Day)")

qqplot(lung$time)
qqnorm(log.time)
qqnorm(log.meal)
log.meal.rep<- rep(log.meal, sample, 1000000)
hist(log.meal.rep)

qqnorm(lung$time)
lung$time
?sample()

cor.test(log.meal,log.time)

qqnorm(lung$time)
qqnorm((lung$meal.cal))
qqnorm(lung$wt.loss)
?barplot()
barplot(c(4,3,2), width = 1, space = NULL,
        names.arg = NULL, legend.text = NULL, beside = FALSE,
        horiz = FALSE, density = NULL, angle = 45,
        col = NULL, border = par("fg"),
        main = NULL, sub = NULL, xlab = "hi", ylab = NULL,
        xlim = NULL, ylim = NULL, xpd = TRUE, log = "",
        axes = TRUE, axisnames = TRUE,
        cex.axis = par("cex.axis"), cex.names = par("cex.axis"),
        inside = TRUE, plot = TRUE, axis.lty = 0, offset = 0,
        add = FALSE, args.legend = NULL)
