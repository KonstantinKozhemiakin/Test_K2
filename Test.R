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
magb8 <-subset(magb7,TARGET.y==2|TARGET.y==3)
magb9 <- magb7[,1:42]
magb10 <- magb7[sapply(magb7,is.numeric)]
#-------------------------------------------------------------------------------
agb1  <- aggregate(magb1[,-2],by=list(magb1$ID),FUN=mean,na.rm=T)
agb2 <- aggregate(agb1[,-c(1:2)],by=list(agb1$TARGET),FUN=mean,na.rm=T,digits=1)
agb3 <- aggregate(magb1[,-2],by=list(magb1$TARGET),FUN = median,na.rm=T)
agb4 <- na.omit(agb1[,c(3:ncol(agb1))])
agb5 <- log(agb4[,c(1:42)])
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
T1 <- table(magb7[,c(44,47)])
chisq.test(T1)

patients <- rbind(c(18, 7), c(6, 13))
#подпишем строки и столбцы
colnames(patients) <- c("Yes", "No")
rownames(patients) <- c("Placebo", "Aspirin")
chisq.test(patients)
#Proverka ostatkov shitesta
mosaicplot(patients, color=T, shade=T, ylab="Thrombosis", xlab="Group")

fisher.test(cbind(c(1,3),c(3,1))) 
#-------------------------------------------------------------------------------
#Klasterizacija K-means
library(ggplot2)
d <- iris[, c("Sepal.Length", "Petal.Width")]

fit <- kmeans(d, 3)
d$clusters <- factor(fit$cluster)

ggplot(d, aes(Sepal.Length, Petal.Width, col = clusters))+
        geom_point(size = 2)+
        theme_bw() 

#Klasterizacija ierarhicheskaja
library(ggplot2) 
library(ggrepel) # для симпатичной подписи точек на графике

x <- rnorm(10)
y <- rnorm(10)
test_data <- data.frame(x, y)
test_data$labels <- 1:10

ggplot(test_data, aes(x, y, label = labels))+
        geom_point()+
        geom_text()

d = dist(test_data)
fit <- hclust(d, method = "single")
plot(fit, labels = test_data$labels)
rect.hclust(fit, 3) # укажите желаемое число кластеров, сейчас стоит 2

test_data <- read.csv("https://stepic.org/media/attachments/course/524/test_data_hclust.csv")
str(test_data)
smart_hclust(test_data, 3) # выделено три кластера
#В этой и следующей задаче на кластерный анализ предполагается, что мы используем функцию hclust() для 
#кластеризации данных с параметрами �по умолчанию:
hclust(d, method = "complete", members = NULL)
#Для расчета матрицы расстояний предполагается, что используется функция dist() также с параметрами по умолчанию:
dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
#Для выделения желаемого числа кластеров по результатам иерархической кластеризации воспользуйтесь функцией cutree().
#Иными словами, для кластеризации данных swiss на три кластера мы бы использовали команды:
dist_matrix <- dist(swiss) # расчет матрицы расстояний
fit <- hclust(dist_matrix) # иерархическая кластеризация 
cluster <- cutree(fit, 3) # номер кластера для каждого наблюдения

d = dist(agb1[,c(3:42)],method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
fit <- hclust(d, method = "complete", members = NULL)
plot(fit, labels = agb1$TARGET)
rect.hclust(fit, 4) # укажите желаемое число кластеров

dist_matrix <- dist(agb1[,c(3:44)]) 
fit <- hclust(dist_matrix)
cluster <- cutree(fit, 4)
#-------------------------------------------------------------------------------
#Metod analiza glavnh komponent
test_data <- read.csv("https://stepic.org/media/attachments/course/524/pca_test.csv")
get_pc <- prcomp(test_data)
#-------------------------------------------------------------------------------
#Problema s plavajushchej zapjatoj
0.1+0.05==0.15
#[1] FALSE
isTRUE(all.equal(0.1+0.05, 0.15))
#-------------------------------------------------------------------------------
#Proverka geteroskedastichnosti
qplot(x=carat,y=price,data=diamonds)+
        geom_smooth(method = lm)
qplot(x=log(carat),y=log(price),data=diamonds)+
        geom_smooth(method = lm)
library(lmtest)
bptest(lm(carat~price,diamonds))
bptest(lm(log(carat)~log(price),diamonds))
fit6 <- lm(carat~price,diamonds)
fit6 <- lm(log(carat)~log(price),diamonds)
shapiro.test(fit6$residuals)
plot(fit6)
#-------------------------------------------------------------------------------               
#Test
hist(magb7$V40)
hist(log(magb7$V40))
qqplot(log(magb7$V40),magb1$TARGET)

fit <- glm(TARGET.y~.,data = magb7)
optimalfit <- step(fit,direction = 'backward')
summary(optimalfit)
summary(fit)
plot(optimalfit)
anova(fit,optimalfit)
fit2 <- lm(as.numeric(TARGET)~.,data=agb4)
summary(fit2)
fit3 <- aov(as.numeric(TARGET)~.,data=agb4)
summary(fit3)
anova(fit2,fit3) 

library("car") 
vif(optimalfit) #Proverka na nalichie multikolinearnosti

fit.res <- lm(I(optimalfit$residuals^2)) 
summary(fit.res)
#-------------------------------------------------------------------------------
#Smeshannie regressionnie modeli
library(lme4)
library(mlmRev)
library(lmerTest)
#-------------------------------------------------------------------------------
magb9 <- magb7[,1:42]


NA_position <- function(x,y){
        all(is.na(x)==is.na(y))
}

sum(is.na(magb9))
magb9[is.na(magb9)]
sum(apply(magb9,2,is.na))

outliers_count <- function(x){
        otliers <- x[abs(x - mean(x)) > 2 * sd(x)]
        if (length(otliers) > 0){
                return(otliers)
        } else {
                return("There are no otliers")
        }
}

magb9_outliers <- apply(magb9, 2, outliers_count)
str(magb9_outliers)

sum_isnan <- function(x){
        sum(is.nan(x))
}
sum_isnan_m <- apply(magb9,2,sum_isnan)

sum_0 <- function(x){
        sum(x==0)
}
sum_0_m <- apply(magb9,2,sum_0)

log(magb9$V1[magb9$V1!="0"])

apply(magb9,2,log_m)

log_m <- function(x){
        log_mm <- log(x[x!=0])
}
log_m(magb9)


aov_res <- apply(magb7[,1:42],2,function(x) aov(x~magb7$TARGET.y))
summary(aov_res$V40)
norm_test <- apply(magb7[,1:42],2,function(x) shapiro.test(x[sample(x,4000)])) 
#-------------------------------------------------------------------------------
#Vivod negativnih znachenij tablici bez NA
test_data <- as.data.frame(list(V1 = c(-9.7, -10, -10.5, -7.8, -8.9), V2 = c(NA, -10.2, -10.1, -9.3, -12.2), V3 = c(NA, NA, -9.3, -10.9, -9.8)))
test_data <- as.data.frame(list(V1 = c(NA, -0.5, -0.7, -8), V2 = c(-0.3, NA, -2, -1.2), V3 = c(1, 2, 3, NA)))

get_negative_values <- function(X){
        res <- apply(X[apply(X, 2, function(X) length(X[X<0 & !is.na(X)])>0)], 2, function(X) X[X<0 & !is.na(X)])
        return(res)
}

get_negative_values <- function(test_data){    
        negative_col <- apply(test_data, 2, function(x) any(x[!is.na(x)] < 0))    
        return(apply(test_data[negative_col], 2, function(x) x[!is.na(x) & x <0]))}

get_negative_values <- function(df) {
        d = apply(df, 2, function (x) x[ x < 0  &  !is.na(x) ])
        d[ lengths(d) > 0 ]
}

get_negative_values(test_data)
#-------------------------------------------------------------------------------
#Delete NA in dataframe
test_data <- as.data.frame(list(V1 = c(NA, NA, NA, NA, 13, 12, 9, 10, 8, 9, 11, 11, 10, 12, 9), V2 = c(NA, 12, 8, NA, 11, 11, 9, 8, 8, 10, 10, 11, 10, 10, 10), V3 = c(NA, 5, NA, 13, 12, 11, 11, 14, 8, 12, 8, 8, 10, 10, 8), V4 = c(10, 10, 10, 10, 13, 10, 11, 7, 12, 10, 7, 10, 13, 10, 9)))
na_rm <- function(x){
        na <- function(x){x[is.na(x)] <- mean(x,na.rm=T);x}
        test1 <- apply(x,2,na)
        return(as.data.frame(test1))
}

na_rm <- function(x){
        return(as.data.frame(apply(x,2,function(x){ x[is.na(x)] <- mean(x,na.rm=T);x})))
}

na_rm  <- function(x){
        as.data.frame(
                mapply(
                        function(vec, val) {vec[is.na(vec)] <- val; vec},
                        x,
                        colMeans(x, na.rm = T)
                )
        )
}

na_rm  <- function(x) {
        if (anyNA(x)) {
                f <- function(xx) {
                        nas <- is.na(xx)
                        xx[nas] <- sum(xx, na.rm = TRUE) / sum(!nas)
                        xx
                }
                x <- as.data.frame(apply(x, 2, f))
        }
        x
}

na_rm  <- function(x){
        as.data.frame(apply(x, 2, function(y) replace (y, is.na(y), mean(y, na.rm = T))))
        
}

na_rm  <- function(x){
        setMean <- function(y){
                mean_val <- mean(y[!is.na(y)])
                y[is.na(y)] <- mean_val
                return(y)
        } 
        
        as.data.frame(apply(x,2,setMean))
}


na_rm(test_data)
#-------------------------------------------------------------------------------
# Logarifmiruem tablici
log_m <- function(x){
        pr <- function(x){x[x<0.01] <- 0.01;x}
        pr1 <- apply(x,2,function(x)pr(x))
        return(as.data.frame(apply(pr1,2,function(x) log(x))))
}

x <- magb9
logmagb <- log_m(x)

log_m <- function(x){
        pr <- apply(x,2,function(x){x[x<0.01] <- 0.01;x})
        return(as.data.frame(apply(pr1,2,function(x) log(x))))
}

log_m <- function(x){
        return(as.data.frame(apply(apply(x,2,function(x){x[x<0.01] <- 0.01;x}),2,function(x) log(x))))
}
#-------------------------------------------------------------------------------
#Summa polozhitelnih znachenij list 
d <- data.frame(X1 = c(-1, -2, 0), X2 = c(10, 4, NA), X3 = c(-4, NA, NA))
positive_sum <-  function(test_data){
        lapply(test_data,function(x)sum(x[!is.na(x)& x>0]))
        
}

positive_sum <-  function(d){
        lapply(d, FUN = function(x) sum(subset(x, x > 0)))
}

positive_sum <- function(df) {
        lapply(lapply(df, function(x) subset(x, x>=0)), sum)
}

positive_sum(d)
#-------------------------------------------------------------------------------
m_names <- mapply(paste, list("row", "col"), list(1:100, 1:200), sep = "_")
str(m_names)

# Poluchit srednee otklonenie
get_sd <- function(x){
        num_var <- sapply(x, is.numeric)
        sapply(x[num_var], sd)
}

# C��������� ��������� ������� ��������� � �������� dataframe:
# my_df[1] - ������� dataframe
# my_df[[1]] - ������� ������
# my_df[, 1] - ������� ������
#-------------------------------------------------------------------------------
# Otbor dannih po chasti nazvanija peremennoj
test_data <- as.data.frame(list(name = c("p4@HPS1", "p7@HPS2", "p4@HPS3", "p7@HPS4", "p7@HPS5", "p9@HPS6", "p11@HPS7", "p10@HPS8", "p15@HPS9"), expression = c(118.84, 90.04, 106.6, 104.99, 93.2, 66.84, 90.02, 108.03, 111.83)))
names = c("HPS5", "HPS6", "HPS9", "HPS2", "HPS3", "HPS7", "HPS4", "HPS8")
names = c("HPS5")
dataset <- as.data.frame(list(name = c("p4@HPS1", "p7@HPS2", "p4@HPS3", "p7@HPS4", "p7@HPS5", "p9@HPS6", "p11@HPS7", "p10@HPS8", "p15@HPS9"), expression = c(118.84, 90.04, 106.6, 104.99, 93.2, 66.84, 90.02, 108.03, 111.83)))


my_names <- function (dataset, names){ 
       a <- as.data.frame(sapply(names,function(x) grepl(x,dataset[[1]])))
       b <- as.logical(rowSums(a))
       c <- test_data$name
       d <- as.character(c[b])
       e <- as.data.frame(sapply(d,function(x)x==dataset[[1]]))
       f <- as.logical(rowSums(e))
       g <- dataset[f,seq(dataset)]
       return(g)
       
}

my_names <- function (dataset, names){    
        dataset[as.numeric(lapply(names, function(x) which(grepl(x, dataset$name)))),]}
#����� ������ names � ��� ������� �������� ���� ��������� � ������� x$mane � ������� grepl
#which - ������� ������ ������ (�.�. ��� ��������� ���� �������)


my_names <- function (dataset, names) {
        dataset[sapply(names, grep, dataset$name), ]
}

# s uchetom povtorenij
my_names <- function (dataset, names) dataset[grepl(paste(names,collapse = "|"), dataset$name),]

# Bez grep
my_names <- function (dataset, names) {
        dataset[gsub(".*@", "", dataset$name) %in% unique(names), ]
}

my_names <- function (dataset, names){
        f <- function (y) any(sapply(names, function (x) grepl(x, y)))
        dataset[c(apply(dataset, 1, f)),]
}
        
my_names <- function (dataset, names){
        t<-sapply(names,function(x) grepl(x,dataset[,1]))
        dataset[apply(t,1,any),]
}

my_names <- function(test_data, names){
        e <- sapply(names, function(x) grepl(x, test_data[,'name']))
        q <- which(e == T, arr.ind = T)
        return(test_data[q[,1],])
}

my_names <- function (dataset, names){
        dataset[as.logical(apply(sapply(names, function(x) grepl(x,dataset$name)), 1, sum)),]
}

my_names <- function (dataset, names){
        m1 <- sapply(names, function(x) grepl(x, dataset[,1]))
        m2 <- apply(m1,1,any)
        dataset[m2,]
}

my_names <- function (dataset, names){ 
        ans<- dataset[grepl(names[1], dataset$name), ]
        for(i in 1:length(names)){
                ans[i, ] <- dataset[grepl(names[i], dataset$name), ]
        }
        return(ans)
}

my_names <- function(df, names) {
        is_in_names <- function(cell) any(sapply(names, function(n) (grepl(n, cell))))
        df[sapply(df$name, is_in_names),]
}

my_names(test_data,names)
#-------------------------------------------------------------------------------
# Postroenie regressionnoj modeli na peremennih s norm raspredeleniem
smart_lm(swiss)

test_data <- read.csv("https://stepik.org/media/attachments/course/724/test.csv")
smart_lm(test_data)

test_data <- data.frame(x = 1:100, y = 1:100, z = 1:100)
smart_lm(test_data)

smart_lm <- function(x){
        col <- x[2:ncol(x)]
        y <- x[(sapply(col[sapply(col, is.numeric)],function(x)shapiro.test(x)$p.value))>0.05]
        return(if(ncol(y)<1) {paste("There are no normal variables in the data")}
                  else{as.vector((lm(paste(paste((labels(y)[[2]])[1]), paste((labels(y)[[2]])[2:ncol(y)], collapse=" + ") , sep=" ~ "),data=y))$coefficients)
                  })
}

smart_lm <- function(x){
        y <- x[(sapply(x[-1][sapply(x[-1], is.numeric)],function(x)shapiro.test(x)$p.value))>0.05]
        return(if(ncol(y)<1) {paste("There are no normal variables in the data")}
               else{lm(y)$coefficients
               })
}

smart_lm <- function(df){
        is.normal <- sapply(df[-1], function(x) shapiro.test(x)$p > 0.05)
        if (sum(is.normal) > 0) {
                lm(df[c(T, is.normal)])$coeff
        } else {
                "There are no normal variables in the data"
        }
}

smart_lm <- function(test_data){
        filtered_data <- cbind(test_data[1],test_data[-1][sapply(test_data[-1],function(x) shapiro.test(x)$p.value) > 0.05])
        if (length(filtered_data) == 1) return("There are no normal variables in the data")
        fit <- lm(filtered_data)
        fit$coefficients
}

smart_lm <- function(x){
        idx <- names(which(sapply(lapply(x[-1], shapiro.test), '[[', "p.value") > 0.05))
        if(length(idx)>0) lm(x[,1] ~ ., data = x[idx])$coef else print("There are no normal variables in the data") 
}
#-------------------------------------------------------------------------------
#Vitaskivaem p-znachenie iz otcheta shapiro-test
x <- rnorm(30, 10, 1)
t.test(x, mu=10)

test_1 <- shapiro.test(iris$Sepal.Length)
test_2 <- shapiro.test(iris$Sepal.Width)
test_3 <- shapiro.test(iris$Petal.Length)
test_4 <- shapiro.test(iris$Petal.Width)

normality_tests <- list(test_1, test_2, test_3, test_4)
normality_tests <- lapply(iris[, 1:4], shapiro.test)

get_p_value <- function(test_list){
        return(sapply(test_list,function(x) x[2]))
        
}

get_p_value <- function(test_list){
        return(lapply(test_list, function(x) x$p.value))
}

get_p_value = function(test_data){    
        sapply(test_data, '[', 2)}

get_p_value <- function(test_list){
        lapply(test_list, .subset2, "p.value")
}

x <- normality_tests
x[[1]]$p.value
sapply(x,function(x) x[2])

#-------------------------------------------------------------------------------
#Sravnenie viborochnih srednih so srednim generalnoj sovokupnosti 
x <- rnorm(30, 10, 1)
t.test(x, mu=10)
one_sample_t(iris[, 1:4], 4)

one_sample_t <- function(test_data, general_mean <- 4){
        test_data <- sapply(test_data,function(x)x[is.numeric(x)])
        test_data[,!is.numeric(test_data)]
        sapply(test_data,function(x)t.test(x,mu=4))
        
}
#-------------------------------------------------------------------------------
aggregate(mpg ~ cyl + am, mtcars, mean) 
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
