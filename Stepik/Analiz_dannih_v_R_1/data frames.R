
# Reading data

?read.table
?read.csv


mydata <- read.csv('evals.csv')


# Summaries

head(mydata, 3)
tail(mydata)

View(mydata)

str(mydata)

a <- names(mydata)

summary(mydata)




# Variables

b <- mydata$score

mean(mydata$score)

summary(mydata$score)

mydata$score * 2

mydata$ten_point_scale <- mydata$score * 2



summary(mydata$ten_point_scale)

mydata$new_varible <- 0
mydata$number <- 1:nrow(mydata)
summary(mydata$number)

nrow(mydata)
ncol(mydata)





# Subsetting

mydata$score[1:10]

mydata[1,1]
mydata[c(2,193,225),1]
mydata[101:200,1]

mydata[5,]
mydata[,1] == mydata$score

mydata[,2:5]
head(mydata[,2:5])

##


# Subsetting with condition

mydata$gender
mydata$gender == 'female'
head(mydata[mydata$gender == 'female',1:3])

head(subset(mydata, gender == 'female'))
head(subset(mydata, score > 3.5))



# rbind, cbind

mydata2 <- subset(mydata, gender == 'female')
mydata3 <- subset(mydata, gender == 'male')
mydata4 <- rbind(mydata2, mydata3)

mydata5 <- mydata[,1:10]
mydata6 <- mydata[,11:24]
mydata7 <- cbind(mydata6, mydata5)


#ѕам€тка
#1. „тобы изучить структуру данных воспользуйтесь командой str()

str(mtcars)
#2. „тобы отобрать нужные колонки (переменные) в данных вы можете:
        
#- использовать номера колонок (не забудьте обернуть индексы в вектор):
mtcars[, c(1, 3, 4)] 

#- использовать имена колонок:
mtcars[, c("mpg", "hp")]

#3. „тобы отобрать нужные строки в данных:
mtcars[c(1, 5, 7), ]

#Ёти приемы можно комбинировать:
mtcars[c(1, 4, 5), c(1, 4)] 

#«апомните, сначала идут индексы строк, потом индексы колонок! “акже обратите внимание, что мы можем использовать отрицательную индексацию, чтобы отобрать все колонки/строки кроме указанных:
mtcars[, -c(3, 4)] # отберем все строчки и все колонки кроме 3 и 4. 

#4. ƒл€ более сложных запросов используйте функцию subset():
subset(mtcars, hp > 100 & am == 1)

#5. ƒобавить новую переменную можно при помощи оператора $
mtcars$new_var <- 1:32

#6. „тобы удалить переменную из данных, используйте такую конструкцию:
mtcars$new_var <- NULL