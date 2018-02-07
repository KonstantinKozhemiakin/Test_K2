#Список установленных пакетов
search()

# Список доступных пакетов для установки
a <- available.packages()
#проверка установлен ли пакет 
find.package("devtools")
#проверка работает ли ртулз
find_rtools()

#узнать рабочую директорию
getwd()
#список файлов рабочей директории
dir()
#установить рабочую директорию
setwd()
#задать ресурс данных
source("mydata.R")
#перечень таблиц и переменных в рабочей директории 
ls()
# Увидеть аргументы
args()
#Create a file in your working directory called "mytest.R" 
file.create("mytest.R")
#Use dir.create() to create a directory in the current working directory called "testdir"
dir.create(path = old.dir)
#Check to see if "mytest.R" exists in the working directory using the 
file.exists("mytest.R")
#Access information about the file
file.info()
#You can use the $ operator --- e.g., file.info("mytest.R")$mode --- to grab specific items.
#
#Change the name of the file "mytest.R" to "mytest2.R" by using file.rename().
file.rename("mytest.R","mytest2.R")
#Deleted files
file.remove()
#copy of "mytest2.R" called "mytest3.R" using file.copy()
file.copy("mytest2.R","mytest3.R")
#Provide the relative path to the file
file.path()
#Создать папку testdir3 в папке testdir2
dir.create(file.path('testdir2', 'testdir3'), recursive = TRUE)

##SQWIRL
install.packages("swirl") #интерактивное обучение
packageVersion("swirl")
library(swirl)
install_from_swirl("R Programming")
swirl()
#You can exit swirl and return to the R prompt (>) at any time by pressing the Esc key. If
#| you are already at the prompt, type bye() to exit and save your progress. When you exit
#| properly, you'll see a short message letting you know you've done so.

#| When you are at the R prompt (>):
#| -- Typing skip() allows you to skip the current question.
#| -- Typing play() lets you experiment with R on your own; swirl will ignore what you do...
#| -- UNTIL you type nxt() which will regain swirl's attention.
#| -- Typing bye() causes swirl to exit. Your progress will be saved.
#| -- Typing main() returns you to swirl's main menu.
#| -- Typing info() displays these options again.

#If at any point you'd like more information on a particular topic related to R, you can
#| type help.start() at the prompt, which will open a menu of resources (either within
#| RStudio or your default web browser, depending on your setup). Alternatively, a simple
#| web search often yields the answer you're looking for.

# Справка по символу ":"
?`:`
seq(1,10)
seq(0, 10, by=0.5)
seq(5, 10, length=30)
length(my_seq)
1:length(my_seq)
seq(along.with=my_seq)
seq_along(my_seq)
rep(0, times= 40)
rep(c(0, 1, 2), times = 10)
rep(c(0, 1, 2), each = 10)

##Vectors
my_char <- c("My", "name","is")
paste(my_char, collapse = " ")
c(my_char, "your_name_here")
my_name <- c(my_char, "Konstantin")
paste("Hello", "world!", sep = " ")
paste(1:3, c("X", "Y", "Z"), sep = "")
paste(LETTERS, 1:4, sep = "-")

y <- rnorm(1000)
z <- rep(NA, 1000)
my_data <- sample(c(y,z),100)
my_na <- is.na(my_data)

x[1:10]
#Index vectors come in four different flavors -- logical vectors, vectors of positive, integers, vectors of negative integers, and vectors of character strings
x[is.na(x)]
y <- x[!is.na(x)]
y[y > 0]
x[!is.na(x)&x>0]
x[c(3,5,7)]
x[c(2, 10)]
#вектор без 2-го и 10-го значения
x[c(-2, -10)]
x[-c(2, 10)]
vect <- c(foo = 11, bar = 2, norf = NA)
names(vect)
vect2 <- c(11,2,NA)
names(vect2) <-c("foo", "bar", "norf")
#Одинаковые ли векторы?
identical(vect,vect2)
vect["bar"]
vect[c("foo","bar")]

##Matrices and Data Frames
dim(my_vector)
length()
dim(my_vector) <- c(4,5)
dim(my_vector)
attributes(my_vector)
class(my_vector)
my_matrix2 <- matrix(1:20,4,5)
patients <- c("Bill","Gina","Kelly","Sean")
cbind(patients,my_matrix)
my_data<- data.frame(patients, my_matrix)
class(my_data)
cnames <- c("patient", "age", "weight", "bp", "rating","test")
colnames(my_data) <- cnames

##Logic
TRUE == TRUE
(FALSE == TRUE) == FALSE
6<7
10<=10
TRUE != FALSE
5!=7 #5 не равно 7
(TRUE != FALSE) == !(6 ==7)
!FALSE
FALSE & FALSE
TRUE & c(TRUE, FALSE, FALSE) #сравнивает все значения по очереди
c(TRUE,TRUE, TRUE) & c(TRUE, FALSE, FALSE)
TRUE && c(TRUE, FALSE, FALSE) #сравнивает только первые значения
TRUE | c(TRUE, FALSE, FALSE)#сравнивает все значения по очереди
TRUE||c(TRUE,FALSE,FALSE)#сравнивает только первые значения
6 != 10 && FALSE && 1 >= 2  
TRUE || 5 < 9.3 || FALSE
5>8||6!=8 && 4>3.9

#рандомный вектор со 100 значениями с нормальным распределением
rnorm(100)
#функция
myfunction <- function() {
        x <- rnorm(100)
        mean(x)
}
ls()
myfunction()
#функция
myfunction <- function(x) {
        y <- rnorm(length(y))
        mean(y)}

x <- vector("numeric",length = 10)
x <- c(1.7,"a") #character
x <- c(TRUE,2) #numeric
x <- c("a",TRUE) # character

x <- 0:6
class(x)
as.numeric(x)
as.logical(x)
as.character(x)
as.complex(x)

x <- c("A","B","C")
class(x)
as.numeric(x)
as.logical(x)
as.character(x)
as.complex(x)

X <- list(1,"A",TRUE,1+4i)
class(x)

m <- matrix(nrow = 2,ncol = 3)
dim(m)
attributes(m)

m <- matrix(1:6, nrow = 2, ncol = 3)

m <- 1:10
dim(m) <- c(2,5)
m

x <- 1:3
y <- 10:12
cbind(x,y)
rbind(x,y)

x <- as.factor(c("yes","yes","no","yes","no"))
table(x)
unclass(x)

x <- factor(c("yes","yes","no","yes","no"), levels=c("yes","no"))

x <- c(1,2,NA,NaN,4)
is.na(x)
is.nan(x)

x <- data.frame(foo=1:4,bar=c(T,T,F,F))
nrow(x)
ncol(x)

x <- 1:3
names(x)
names(x) <- c("foo","bar","norf")

x <- list(a=1,b=2,c=3)

m <- matrix(1:4,nrow = 2,ncol = 2)
dimnames(m) <- list(c("a","b"),c("c","d"))

read.csv()
read.table()
source()
dget()
load()
unserialize()

write.table()
write.csv()
dump()
save()
serialize()

read.table(file = ,sep = ,colClasses = ,nrows = , comment.char = ,skip = ,
           stringsAsFactors = )

getwd()
setwd("C:/Users/User/Documents/R")

file.show("Base1.txt")
count.fields("Base1.txt")
scan(file = "")
str("C:/Users/User/Documents/R/Base1.txt") 
summary("C:/Users/User/Documents/R/Base1.txt")
class("C:/Users/User/Documents/R/Base1.txt")

initial <- read.table("Base1.txt",header = T,nrows = 100,sep="\t")
classes <- sapply(initial,class)
tabAll <- read.table("Base1.txt",header = T,sep="\t")
str(tabAll) 
summary(tabAll)
class(tabAll)

getwd()
setwd("C:/Users/User/Documents/R")

ya <-  data.frame(a=1, b="a")
dput(ya)
dput(ya,file="ya.R")
new.y <- dget("ya.R")

x <- "foo"
y <- data.frame(a=1,b="a")
dump(c("x","y"),file="data.R")
rm(x,y)
source("data.R")

file()
gzfile()
bzfile()

str(file)
con <- file("foo.txt","r")
data <- read.csv(con)
close(con)

data <- read.csv("foo.txt")

con <- gzfile("words.gz")
x <- readline(con,10)

x <- c("a","b","c","c","d","a")
x[1]
x[2]
x[1:4]
x[x>"a"]
u <- x>"a"
x[u]

x <- list(foo=1:4,bar=0.6)
x[1] #list
x[[1]]
x$bar 
x["bar"] #list
x[["bar"]]

x <- list(foo=1:4,bar=0.6,baz="hello")
x[c(1,3)] #list
x[[c(1,3)]] 
x[[1]][[3]]

x <- list(foo=1:4,bar=0.6,baz="hello")
name <- "foo"
x[(name)]
x$name
x$foo

x <- list(a=list(10,12,14),b=c(3.14,2.81))
x[(c(1,3))] #list
x[[c(1,3)]]
x[[1]][[3]]
x[[c(2,1)]]

x <- matrix(1:6,2,3)
x[1,2] 
x[2,1]
x[1,2,drop=F]

x <- matrix(1:6,2,3)
x[1,]
x[1,,drop=F]

x <- list(adrwark=1:5)
x$a
x[["a"]]
x["a"] 
x[["a",exact=F]]

x <- c(1,2,NA,4,5)
bad <- is.na(x)
x[!bad]
x[bad]
y <- c("a","b",NA,"d","f")
good <- complete.cases(x,y)
x[good]

x <- airquality[1:6,]
good <- complete.cases(airquality)
airquality[good,][1:6,]

x <- 1:4; y <- 6:9
x>2
x>=2
y==8 #знак равно

x <- list(hi=1:9,3,3)
y <- list(by=2:10,3,3)
z <- x$hi+y$by

x <- matrix(1:4,2,2);y <- matrix(rep(10,4),2,2)
x+y
x/y
x*y
x%*%y # умножение матриц


file.remove('testdir')

if(x>3){
        y<- 10
}else{
        y <- 0
}

y <- if(x>3) {
        10
} else {
        0
}

for (i in 1:10){
        print(i)
}

x <- c("a","b","c","d")
for(i in 1:4){
        print(x[i])
}

for(i in seq_along(x)){
        print(x[i])
}

for(i in 1:4)print(x[i])

x <- matrix(1:6,2,3)
for(i in seq_len(nrow(x))){
        for(j in seq_len(ncol(x))){
                print(x[i,j])
        }
}

count <-0
while(count<10){
        print(count)
        count <- count+1
}

z <- 5
while(z >=3 && z <= 10){
        print(z)
        coin <- rbinom(1,1,0.5)
        if (coin==1){ #random walk
                z <- z+1
        } else{
                z <- z-1
        }
}

x <- 1
tol <- 1e-8

repeat{
        x1 <- computeEstimate()
        if(abs(x1-x0)<tol){
                break
        }else{
                x0 <- x1
        }
}

install.packages("sos")
require("sos")
findFn("computeEstimate")
findFn('multiply',maxPages = 1)

for(i in 1:100){
        if(i<=20){
                ##Skip first 20 iterations
                next
        }
        ##do something here
}

add2 <- function(x,y){
        x+y
}
add2(3,5)

above10 <- function(x){
        use <- x>10
        x[use]
}
x <- c(1:100)
above10(x)


above <- function(x,n){
        use <- x>n
        x[use]
}
x <- 1:20
above(x,10)

above <- function(x,n=10){ #по умолчанию 10, но можно менять
        use <- x>n
        x[use]
}
x <- 1:20
above(x)

above <- function(x,n=10){
        use <- x>n
        x[use]
}
x <- 1:20
above(x,15)

columnmeen <- function(y){
        nc <- ncol(y)
        means <- numeric(nc)
        for(i in 1:nc){
                means[i] <- mean(y[,i])
        }
}

columnmeen(airquality)

columnmeen <- function(y){
        nc <- ncol(y)
        means <- numeric(nc)
        for(i in 1:nc){
                means[i] <- mean(y[,i])
        }
        means #выводит данные
}
columnmeen(airquality)

#уберём NA
columnmeen <- function(y,removeNA=TRUE){ #добавили аргумент ремувНА
        nc <- ncol(y)
        means <- numeric(nc)
        for(i in 1:nc){
                means[i] <- mean(y[,i],na.rm = removeNA) #вставили его в формулу и прировняли на.рм к ремувНА
        }
        means
}
columnmeen(airquality)

#-------------------------------------------------------------------------------

setwd("D:/Кости/R/R Script/")
filelist <- list.files(path = "D:/Кости/R/R Script/specdata/",pattern = ".csv",full.names = T)

complete <- function(directory,id=1:332){
        filelist <- list.files(path = directory,pattern = ".csv",full.names = T)
        nobs <- numeric()
        for(i in id){
                data <- read.csv(filelist[i])
                nobs <- c(nobs,sum(complete.cases(data)))
        }
        data.frame(id,nobs)
}
nobs <- c(sum(complete.cases(filelist)))
nobs
complete.cases(filelist)
sum(complete.cases(filelist))
complete.cases(filelist)
nobs <- complete.cases(filelist)

corr <- function(directory,threshold=0){
        filelist <- list.files(path = directory,pattern = ".csv",full.names = T)
        df <-  complete(directory)
        ids <-  df[df["nobs"] > threshold, ]$id
        corrr <-  numeric()
        for(i in ids){
                data <- read.csv(filelist[i])
                nobs <- c(nobs,sum(complete.cases(data)))
                dff <- data[complete.cases(data), ]
                corrr = c(corrr, cor(dff$sulfate, dff$nitrate))
        }
        corrr
}
corr("specdata/",1:10)

#-------------------------------------------------------------------------------
reading_file <- function(directory,column,names_col,col_m_1,col_m_2) {
        catalog <- list.files(path = directory,
                              pattern = ".csv",
                              full.names = TRUE)
        read_catalog <- function(catalog) {
                dat<- read.csv(file = catalog, header = T)
        }
        names <- list(NULL,names_col)
        fin_dat <-matrix(ncol=column,byrow=T,dimnames = names)
        for (i in catalog) {
                temp_dat<- read_catalog(i)
                fin_dat <- rbind(temp_dat,fin_dat)
        }
        fin_dat
        #fin_mean <-c(mean(fin_dat[,col_m_1],na.rm = T),mean(fin_dat[,col_m_2],na.rm = T))
        #fin_mean
        #mean(fin_dat[["nitrate"]],na.rm = T)
}
reading_file(directory = "C:/Users/User/Documents/R/specdata",column = 4,names_col = c("Date","sulfate","nitrate","ID"),col_m_1 = 2,col_m_2 = 3)
#-------------------------------------------------------------------------------


x < list(a = 1:5, b = rnorm(10))
x < list(a = 1:5, b = rnorm(10))
lapply(x, mean)
x < list(a = 1:5, b = rnorm(10))
x <-  list(a = 1:5, b = rnorm(10))
lapply(x, mean)
x <- list(a = 1:4, b = rnorm(10), c = rnorm(20, 1), d = rnorm(100, 5))
lapply(x, mean)
x <- 1:4
lapply(x, runif)
x <- 1:4
lapply(x, runif, min = 0, max = 10)
lapply(x, function(elt) elt[,1])
x <- list(a = matrix(1:4, 2, 2), b = matrix(1:6, 3, 2))
x
lapply(x, function(elt) elt[,1])
x <- list(a = 1:4, b = rnorm(10), c = rnorm(20, 1), d = rnorm(100, 5))
lapply(x, mean)
sapply(x, mean)
x <- matrix(rnorm(200), 20, 10)
apply(x, 2, mean)
x
apply(x, 2, mean)
apply(x, 1, mean)
rowSums = apply (x, 1, sum)
x <- matrix(rnorm(200), 20, 10)
apply(x, 1, quantile, probs = c(0.25, 0.75))
a <- array(rnorm(2 * 2 * 10), c(2, 2, 10))
apply(a, c(1, 2), mean)
mapply(rep, 1:4, 4:1)
list(rep(1, 4), rep(2, 3), rep(3, 2), rep(4, 1))
noise <- function(n, mean, sd) {
        rnorm(n, mean, sd)
}
noise(5, 1, 2)
noise <- function(n, mean, sd) {
        rnorm(n, mean, sd)
}
noise(5, 1, 2)
noise(1:5, 1:5, 2)
mapply(noise, 1:5, 1:5, 2)
x <- c(rnorm(10), runif(10), rnorm(10, 1))
f <- gl(3, 10)
f
tapply(x, f, mean)
x
f <- gl(3, 10)
f
tapply(x, f, mean)
tapply(x, f, mean, simplify = FALSE)
tapply(x, f, range)
x <- c(rnorm(10), runif(10), rnorm(10, 1))
f <- gl(3, 10)
split(x, f)
x <- c(rnorm(10), runif(10), rnorm(10, 1))
f <- gl(3, 10)
f
tapply(x, f, mean)
lapply(split(x, f), mean)

library(datasets)
head(airquality)
s <- split(airquality, airquality$Month)
s
head(airquality)
s <- split(airquality, airquality$Month)
lapply(s, function(x) colMeans(x[, c("Ozone", "Solar.R", "Wind")]))
sapply(s, function(x) colMeans(x[, c("Ozone", "Solar.R", "Wind")]))
sapply(s, function(x) colMeans(x[, c("Ozone", "Solar.R", "Wind")], na.rm  = TRUE))
x <- rnorm(10)
x <- c(rnorm(10), runif(10), rnorm(10, 1))
x
x <- rnorm(10) #случайная переменная с нормальным распределением с 10 наблюдениями
x
f1 <- gl(2, 5) #фактор с двумя уровнями, каждый с которых повторяеться 5 раз
f1
f2 <- gl(5, 2) #фактор с пятью уровнями, каждый с которых повторяеться 2 раза
f2
interaction(f1, f2) #комбинирует все уровни первого фактора со всеми факторами второго второго
str(split(x, list(f1, f2)))
split(x, list(f1, f2))
x <- rnorm(10) #случайная переменная с нормальным распределением с 10 наблюдениями
f1 <- gl(2, 5) #фактор с двумя уровнями, каждый с которых повторяеться 5 раз
f2 <- gl(5, 2) #фактор с пятью уровнями, каждый с которых повторяеться 2 раза
f1
f2
interaction(f1, f2) #комбинирует все уровни первого фактора со всеми факторами второго второго
str(split(x, list(f1, f2))) #разбить вектор x в соответствии с двумя факторами (f1, f2)
x <- rnorm(10) #случайная переменная с нормальным распределением с 10 наблюдениями
str(split(x, list(f1, f2))) #разбить вектор x в соответствии с двумя факторами (f1, f2)
split(x, list(f1, f2))
str(split(x, list(f1, f2))) #разбить вектор x в соответствии с двумя факторами (f1, f2)

printmessage <- function(x) {
        if(x > 0)
                + 		print("x is greater than zero")
        else
                + 		print("x is less than or equal to zero")
        invisible(x)
}

printmessage <- function(x) {
        if(x > 0)
                print("x is greater than zero")
        else
                print("x is less than or equal to zero")
        invisible(x) #предотвращает автоматический вывод результатов по окончанию фунции.
        #Можно пользоваться функцией load, которая работает с объектами простраства
}
printmessage(1)
printmessage <- function(x) {
        if(x > 0)
                print("x is greater than zero")
        else
                print("x is less than or equal to zero")
        #invisible(x) #предотвращает автоматический вывод результатов по окончанию фунции.
        #Можно пользоваться функцией load, которая работает с объектами простраства
}
printmessage(1)
printmessage2 <- function(x) {
        if(is.na(x))
                print("x is a missing value!")
        else if(x > 0)
                print("x is greater than zero")
        else
                print("x is less than or equal to zero")
        invisible(x)
}
