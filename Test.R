train <- read.table(file = "D:/Кости/R/R Script/Task_BigDataSchool/train.txt",header = T)
test <- read.table(file = "D:/Кости/R/R Script/Task_BigDataSchool/test.txt",header = T,sep = "\t")
base_1 <- read.table(file = "D:/Кости/R/R Script/Task_BigDataSchool/Base1.txt",header = T,sep = "\t")
base_2 <- read.table(file = "D:/Кости/R/R Script/Task_BigDataSchool/Base2.txt",header = T,sep = "\t")

library(dplyr)
library(psych)

magb1 <- merge(base_1,train,by="ID")
magb2 <- merge(base_2,train,by="ID")
magb3 <- merge(agb1,magb2,by="ID")
magb5 <- merge(magb1,magb2,by='ID')
magb6 <- na.omit(magb3)
magb7 <- magb6[,c(3:44,46:50)]
#-------------------------------------------------------------------------------
agb1  <- aggregate(magb1[,-2],by=list(magb1$ID),FUN=mean,na.rm=T)
agb2 <- aggregate(agb1[,-c(1:2)],by=list(agb1$TARGET),FUN=mean,na.rm=T,digits=1)
agb3 <- aggregate(magb1[,-2],by=list(magb1$TARGET),FUN = median,na.rm=T)
#-------------------------------------------------------------------------------
levels(magb2$T1) <-c(1:25) 
levels(magb2$T2) <- c(1:7)
levels(magb2$T3) <- c(1:152)
levels(magb2$T4) <- c(1:1394)
levels(magb2)
str(magb2)
dim(magb2)

#-------------------------------------------------------------------------------
t1 <- table(magb2$TARGET,magb2$T2)
t1
prop.table(t1)
barplot(t1,legend.text = T,args.legend = list(x="topright"),beside = T)
#chisq.test(t1)
#fisher.test(t1)
boxplot(V40~TARGET,magb1)
#-------------------------------------------------------------------------------
str(magb1)
summary(magb1)
magb1$TARGET <- as.factor(magb1$TARGET)
df1 <- subset(magb1, TARGET!= "0"& TARGET!= "2")
hist(df1$V40)
library(ggplot2)
agb1$TARGET <- as.factor(agb1$TARGET)
ggplot(agb1,aes(TARGET,V40))+
        geom_boxplot()
ggplot(df1,aes(TARGET,V40))+
        stat_summary(fun.data = mean_cl_normal,
                     geom = "errorbar",with=0.1
        )+
        stat_summary(fun.y = mean, geom="point", size=2)
ggplot(agb1,aes(TARGET,V40))+
        stat_summary(fun.data = mean_cl_boot,
                     geom = "errorbar"
        )+
        stat_summary(fun.y = mean, geom="point", size=2)
mean_cl_boot(agb1$V40)
summary(df1$V40[df1$TARGET=="1"],na.action=na.pass)
boxplot(df1$V40[df1$TARGET=="1"])
hist(agb1$V40[agb1$V40 != 0],breaks = 500)
hist(agb1$V40[df1$TARGET=="1"])
shapiro.test(sample(agb1$V40,4000))
shapiro.test(sample(agb1$V40[magb1$TARGET=="3"],4000)) #normalnost raspredelenia
bartlett.test(V40~TARGET,magb1) #Gomogennost despersij
t.test(V40~TARGET,df1)
t.test(df1$V40,df1$V35,paired = T,na.action=na.pass)
wilcox.test(V40~TARGET,df1) #Kriterij Vilkoxona,MonaVitni
ggplot(df1,aes(TARGET,V40))+
        geom_boxplot()
ggplot(df1,aes(TARGET,V40))+geom_boxplot()
fit <- aov(V40~TARGET,agb1)
summary(fit)

fit1 <- aov(V40~TARGET.x+T2,agb4)
summary(fit1)
model.tables(fit1,"means")

library(ggplot2)
pd = position_dodge(0.1)
ggplot(agb4,aes(x=T2, y=V40,color=TARGET.x,group=TARGET.x))+
        stat_summary(fun.data = mean_cl_boot,geom = 'errorbar',position = pd)+
        stat_summary(fun.data = mean_cl_boot,geom = 'line',position = pd)+
        stat_summary(fun.data = mean_cl_boot,geom = 'point',size=2,position = pd,pch=15)+
        theme_bw()
fit3 <- aov(V40~TARGET.x*T2,agb4)
summary(fit3)
fit4 <- aov(V40~TARGET.x,agb4)
summary(fit4)
TukeyHSD(fit4)

ggplot(agb4,aes(TARGET.x,V40))+geom_boxplot()
ggplot(agb4,aes(y=V40,x=TARGET.x))+
        stat_summary(fun.data = mean_cl_boot,geom = 'errorbar')


magb5$ID <- as.factor(magb5$ID)
magb5$TARGET.y <- as.numeric(magb5$TARGET.y)
magb5$MONTH <- as.factor(magb5$MONTH)
fit5 <- aov(TARGET.y~T2,data=magb5)
summary(fit5)
fit5 <- aov(V40~TARGET.y+Error(ID/V40),data=magb5)
summary(fit5)
ggplot(magb5,aes(y=TARGET.y,x=T2))+
        stat_summary(fun.data = mean_cl_boot,geom = 'errorbar')+
        stat_summary(fun.y = mean, geom="point", size=1)
ggplot(magb5,aes(x=T2,y=TARGET.y))+geom_boxplot()+
        facet_grid(~MONTH)

fit6 <- aov(TARGET.y~T2*MONTH,data=magb5)
summary(fit6)
fit6 <- aov(TARGET.y~T2*MONTH+Error(ID/(T3*MONTH)),data=magb5)
summary(fit6)



#-------------------------------------------------------------------------------
my_na_rm <- function(x,y){
        x[is.na(x[,y])] <- a <- mean(x[,y],na.rm = T)
        print(a)
}
my_na_rm(magb5,3)
magb5[,3]
sum(is.na(magb5[,3]))

my_na_rm2 <- function(data,col=1:ncol(data)){
                for(i in col){
                data[is.na(i)] <- mean(i,na.rm = T)
        }
}
magb5[,3] <- my_na_rm2(magb5,magb5[,3])
sum(is.na(magb5[,3]))

fit <- cor.test(agb1$V40,agb1$TARGET)
summary(fit)
cor.test(~V40+as.numeric(TARGET),agb1)
plot(y=agb1$V40,x=agb1$V35)
ggplot(agb1, aes(V40,V35,col=factor(TARGET)))+geom_point(size=0.5)
df <- magb1[,c("V40","V35","TARGET")]
pairs(df) #
df <- magb1[,c(3:45)]
cor(df)
library(psych)
fit <- corr.test(df,adjust = "holm",method = "spearman")
fit
fit <- lm(as.numeric(TARGET.x)~V40,magb6)
summary(fit)

fit <- lm(V40~TARGET.x,magb6)
summary(fit)
confint(fit) #doveritelnie intervali (esli peresekajut 0 - ne predskazivajut)

ggplot(magb6,aes(V40,as.numeric(TARGET.x)))+
        geom_point()+
        geom_smooth(method = "lm")+
        facet_grid(.~T2)
ggplot(magb6,aes(TARGET.x,V40))+
        geom_point()+
        geom_smooth(method = "lm")

fitted.val <- data.frame(V40=magb6$V40, fitted=fit$fitted.values)  

newdf <- as.data.frame(magb6$V40)
newdf$newtarg <- predict(fit,newdf)
newdf$TARGET <- magb6$TARGET

hist(as.numeric(magb6$TARGET.x))
fit <- lm(as.numeric(magb6$TARGET.x)~as.factor(magb6$T2)+as.factor(magb6$T3),magb6)
summary(fit)

magb7 <- magb6[,c(3:44,46:50)]
labels(magb6$TARGET.y,"TARGET")
fit <- lm(TARGET.y~.,data = magb6)
summary(fit) 
fit2 <- lm(TARGET.y~V1+V25+V28+V35+V39+V41+V42+T1+T3,data = magb6)
summary(fit2)
anova(fit,fit2)
fit3 <- lm(TARGET.y~V1+V25+V28+V35+V39+V41+V42,data = magb6)
summary(fit3)

#Poisk optimalnoj modeli predskazanija
optimalfit <- step(fit,direction = 'backward')
optimalfit <-lm(formula = TARGET.y ~ V1 + V2 + V4 + V5 + V8 + V10 + V11 + 
           V13 + V15 + V16 + V19 + V21 + V22 + V25 + V26 + V27 + V28 + 
           V31 + V32 + V34 + V35 + V38 + V39 + V40 + V41 + V42 + T1, 
   data = magb6)
summary(optimalfit)

#Start analising Regression model
pairs(magb6[,c(5:10,47)]) #proverka koreljacii
#Postroenie linii trenda
ggplot(magb6,aes(x=V40,y=TARGET.y))+
        geom_point()+
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10,face="bold"))+
        geom_smooth(method='lm')
ggplot(magb6,aes(x=V40,y=TARGET.y))+
        geom_point()+
        theme(axis.text=element_text(size=10),
              axis.title=element_text(size=10,face="bold"))+
        geom_smooth()
#Proverka normalnosti raspredelenija
ggplot(magb6,aes(x=V40))+
        geom_histogram()

ggplot(magb6,aes(x=TARGET.y))+
        geom_histogram()

ggplot(magb6,aes(x=log(V40)))+
        geom_histogram()
#Sozdanie modeli

lm1 <- lm(TARGET.y~V40,magb6)
summary(lm1)
magb6$V40_sq <- (magb6$V40)^2
lm2 <- lm(TARGET.y~V40+V40_sq,magb6)
summary(lm2)
anova(lm1,lm2) #Sravnenie dvuh modelej, luchshe li vtoraja pervoj

magb6$lm1 <- lm1$fitted.values
magb6$lm2 <- lm2$fitted.values
magb6$res1 <- lm1$residuals
magb6$res2 <- lm2$residuals
magb6$number <- 1:nrow(magb6)
#Sravnenie 2-h modelej (linejnaja i kvadratichnaja)
ggplot(magb6,aes(V40,TARGET.y))+
        geom_point()+
        geom_line(aes(x=V40,y=lm1),col='red',lwd=1)+
        geom_line(aes(x=V40,y=lm2),col='blue',lwd=1)
#Proverka normalnosti predskazanija modeli
#Postroenie grafika ostatkov i znachenij predskazannoj modeli
ggplot(magb6,aes(x=lm1,y=res1))+
        geom_point(size=3)+
        geom_hline(yintercept = 0,col='red',lwd=1)
ggplot(magb6,aes(x=lm2,y=res2))+
        geom_point(size=3)+
        geom_hline(yintercept = 0,col='red',lwd=1)

#Logisticheskaja regressija
magb8 <-subset(magb7,TARGET.y==2|TARGET.y==3)
magb8$TARGET.y <- as.factor(magb8$TARGET.y)  
ggplot(magb8, aes(V40,V35,col=T2))+
        geom_point()+
        facet_grid(.~TARGET.y)+
        theme(axis.title = element_text(),
              axis.text = element_text())
fit <- glm(TARGET.y~V40+V35+T2,magb8,family = "binomial")
summary(fit)
exp(fit$coefficients) #coeficienti bez logarifma
head(predict(object = fit)) #natur logarifm ot odds
head(predict(object = fit,type = "response")) #Veojatnost togo, chto budet odin iz factorov (opredelennij)
magb8$prob <- predict(object = fit,type ="response") 

library(ROCR)
pred_fit <- prediction(magb8$prob,magb8$TARGET.y)
perf_fit <- performance(pred_fit,"tpr","fpr")
plot(perf_fit)
auc <- performance(pred_fit,measure = "auc")
str(auc)
perf2 <- performance(pred_fit,x.measure = "cutoff",measure = "spec")
perf3 <- performance(pred_fit,x.measure = "cutoff",measure = "sens")
perf4 <- performance(pred_fit,x.measure = "cutoff",measure = "acc")
plot(perf2,col="red")
plot(add=T,perf3,col="blue")
plot(add=T,perf4,col="green")
legend(x=0.6,y=0.5,c("spec","sens","acc"),lty=1,col = c('red','blue','green'),bty='n',cex=1,lwd=1)
abline(v=0.31,lwd=2)
magb8$pred_resp <- factor(ifelse(magb8$prob>0.31,1,0),labels = c("2","3"))
magb8$correct <- ifelse(magb8$pred_resp==magb8$TARGET.y,1,0)
ggplot(magb8,aes(prob,fill=factor(correct)))+
        geom_dotplot()+
        theme(axis.title = element_text(),
              axis.text = element_text())
mean(magb8$correct)

testdf <- subset(magb7,TARGET.y==2|TARGET.y==3)
testdf$TARGET.y <- NA
testdf$TARGET.y <- predict(fit,newdata = testdf,type="response")

#Export rezultatov analiza
library(xtable)
library(stargazer)
magb8$TARGET.y <- as.numeric(magb8$TARGET.y)
fit1 <- lm(TARGET.y~V40+T2,magb8)
summary(fit1)
fit2 <- aov(TARGET.y~V40+V1,data=magb8)
summary(fit2)

fit_table1 <- xtable(fit1)
fit_table2 <- xtable(fit2)

print(fit_table1,type="html",file="fit_table1.html")
print(fit_table2,type="html",file="fit_table2.html")

stargazer(fit1,type="html",dep.var.labels = "Target",covariate.labels = c("V40","T2"),out = "models1.html")
#-------------------------------------------------------------------------------



#-------------------------------------------------------------------------------

