#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## R-Programming
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

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

## Write a function that reads a directory full of files and reports the number of completely observed cases in each data file.
## The function should return a data frame where the first column is the name of the file and the second column is the number of complete cases.
## A prototype of this function follows
complete <- function(directory, id = 1:332) {
        
        ## Get a list of filenames
        filenames <- list.files(path=directory, pattern="*.csv")
        
        ## Initialize variables
        ids <-vector()
        counts = vector()
        
        ## Loop over the passed id's
        for(i in id) {
                
                ## Pad the i to create a filename
                filename <- sprintf("%03d.csv", i)
                filepath <- paste(directory, filename, sep="/")
                
                ## Load the data
                data <- read.csv(filepath)
                
                ## Store the id
                ids <- c(ids, i)
                
                ## Calculate and store the count of complete cases
                completeCases <- data[complete.cases(data),]
                counts <- c(counts, nrow(completeCases))
        }
        
        ## Return the data frame
        data.frame(id=ids, nobs=counts)
}

#-------------------------------------------------------------------------------

## Write a function that takes a directory of data files and a threshold for
## complete cases and calculates the correlation between sulfate and nitrate for
## monitor locations where the number of completely observed cases (on all
## variables) is greater than the threshold.
## The function should return a vector of correlations for the monitors that
## meet the threshold requirement. If no monitors meet the threshold
## requirement, then the function should return a numeric vector of length 0.
corr <- function(directory, threshold = 0) {
        
        completes <- complete(directory, 1:332)
        completes <- subset(completes, nobs > threshold )
        
        ## Initialize variables
        correlations <- vector()
        
        ## Loop over the passed id's
        for(i in completes$id ) {
                
                ## Pad the i to create a filename
                filename <- sprintf("%03d.csv", i)
                filepath <- paste(directory, filename, sep="/")
                
                ## Load the data
                data <- read.csv(filepath)
                
                ## Calculate and store the count of complete cases
                completeCases <- data[complete.cases(data),]
                count <- nrow(completeCases)
                
                ## Calculate and store the count of complete cases
                ## if threshhold is reached
                if( count >= threshold ) {
                        correlations <- c(correlations, cor(completeCases$nitrate, completeCases$sulfate) )
                }
        }
        
        ## Return a numeric vector of correlations
        correlations
}

#-------------------------------------------------------------------------------

## Calculates the mean of a pollutant (sulfate or nitrate) across a specified list of monitors.
## The function 'pollutantmean' takes three arguments: 'directory', 'pollutant', and 'id'.
## Given a vector monitor ID numbers, 'pollutantmean' reads that monitors' particulate matter data
## from the directory specified in the 'directory' argument and returns the mean of the pollutant
## across all of the monitors, ignoring any missing values coded as NA.
pollutantmean <- function(directory, pollutant, id = 1:332) {
        
        ## Get a list of filenames
        filenames <- list.files(path=directory, pattern="*.csv")
        
        ## Initialize a vector to hold values
        vals <- vector()
        
        ## Loop over the passed id's
        for(i in id) {
                
                ## Pad the i to create a filename
                filename <- sprintf("%03d.csv", i)
                filepath <- paste(directory, filename, sep="/")
                
                ## Load the data
                data <- read.csv(filepath)
                
                ## Select our column
                d <- data[,pollutant]
                
                ## Ignore NAs
                d <- d[!is.na(d)]
                
                ## append to our vector
                vals <- c(vals, d)
        }
        
        ## Return the value rounded to 3 dec places
        round(mean(vals), 3)
}

#-------------------------------------------------------------------------------

checkPkgs <- function() {
        pkg.inst <- installed.packages()
        pkgs <- c("RCurl", "digest")
        have.pkg <- pkgs %in% rownames(pkg.inst)
        
        if(any(!have.pkg)) {
                cat("Some packages need to be installed\n")
                r <- readline("Install necessary packages [y/n]? ")
                if(tolower(r) == "y") {
                        need <- pkgs[!have.pkg]
                        message("installing packages ",
                                paste(need, collapse = ", "))
                        install.packages(need)
                }
        }
}

checkPkgs()

CLASS <- "rprog-002"
challenge.url <- paste("http://class.coursera.org", CLASS,
                       "assignment/challenge", sep = "/")
submit.url <- paste("http://class.coursera.org", CLASS,
                    "assignment/submit", sep = "/")

loginPrompt <- function() {
        email <- readline("Submission login (email): ")
        passwd <- readline("Submission  password: ")
        r <- list(email = email, passwd = passwd)
        assign(".CourseraLogin", r, globalenv())
        invisible(r)
}

submit <- function(manual = FALSE, resetLogin = FALSE) {
        library(RCurl)
        library(digest)
        if(!manual) {
                if(exists(".CourseraLogin") && !resetLogin)
                        cred <- get(".CourseraLogin")
                else
                        cred <- loginPrompt()
                if(!is.list(cred) || !(names(cred) %in% c("email", "passwd")))
                        stop("problem with login/password")
                email <- cred$email
                password <- cred$passwd
        }
        ## Prompt Submission Part
        sid <- partPrompt()
        
        ## Get output
        output <- getOutput(sid)        
        
        if(!manual) {
                ## Get challenge
                ch <- getChallenge(email)
                
                ## Attempt submission with challenge
                ch.resp <- challengeResponse(password, ch$ch.key)
                results <- submitSolution(email, ch.resp, sid, output, ch$state)
                if(!length(results))
                        results <- "Incorrect!"
                cat("Result: ", results, "\n")
        }
        else {
                outfile <- paste(sid, "output.txt", sep = "-")
                writeLines(as.character(output), outfile)
                cat(sprintf("Please upload the file '%s' to Coursera\n",
                            outfile))
        }
        invisible()
}

getOutput <- function(sid) {
        ## JUST FOR TESTING
        ### sid <- sub("-dev", "", sid, fixed = TRUE)
        if(sid == "pollutantmean-1") {
                source("pollutantmean.R")
                pollutantmean("specdata", "sulfate", 1:10)
        }
        else if(sid == "pollutantmean-2") {
                source("pollutantmean.R")
                pollutantmean("specdata", "nitrate", 70:72)
        }
        else if(sid == "pollutantmean-3") {
                source("pollutantmean.R")
                pollutantmean("specdata", "sulfate", 34)
        }
        else if(sid == "pollutantmean-4") {
                source("pollutantmean.R")
                pollutantmean("specdata", "nitrate")
        }
        else if(sid == "complete-1") {
                source("complete.R")
                cc <- complete("specdata", c(6, 10, 20, 34, 100, 200, 310))
                paste(cc$nobs, collapse = "\n")
        }
        else if(sid == "complete-2") {
                source("complete.R")
                cc <- complete("specdata", 54)
                cc$nobs
        }
        else if(sid == "complete-3") {
                source("complete.R")
                set.seed(42)
                cc <- complete("specdata", 332:1)
                use <- sample(332, 10)
                paste(cc[use, "nobs"], collapse = "\n")                
        }
        else if(sid == "corr-1") {
                source("corr.R")
                cr <- corr("specdata")
                cr <- sort(cr)
                set.seed(868)
                out <- round(cr[sample(length(cr), 5)], 4)
                paste(out, collapse = "\n")
        }
        else if(sid == "corr-2") {
                source("corr.R")
                cr <- corr("specdata", 129)
                cr <- sort(cr)
                n <- length(cr)
                set.seed(197)
                out <- c(n, round(cr[sample(n, 5)], 4))
                paste(out, collapse = "\n")
                
        }
        else if(sid == "corr-3") {
                cr <- corr("specdata", 2000)
                n <- length(cr)
                cr <- corr("specdata", 1000)
                cr <- sort(cr)
                paste(c(n, round(cr, 4)), collapse = "\n")
        }
        else {
                stop("invalid part number")
        }
}

partPrompt <- function() {
        sid <- c("pollutantmean-1",
                 "pollutantmean-2",
                 "pollutantmean-3",
                 "pollutantmean-4",
                 "complete-1",
                 "complete-2",
                 "complete-3",
                 "corr-1",
                 "corr-2",
                 "corr-3"
        )
        ## JUST FOR TESTING
        ## sid <- paste(sid, "dev", sep = "-")
        
        sidname <- c("'pollutantmean' part 1",
                     "'pollutantmean' part 2",
                     "'pollutantmean' part 3",
                     "'pollutantmean' part 4",
                     "'complete' part 1",
                     "'complete' part 2",
                     "'complete' part 3",
                     "'corr' part 1",
                     "'corr' part 2",
                     "'corr' part 3"                     
        )
        numparts <- length(sid)
        cat(paste(paste("[", seq_len(numparts), "]", sep = ""), sidname),
            sep = "\n")
        partnum <- readline(sprintf("Which part are you submitting [1-%d]? ",
                                    numparts))
        partnum <- as.integer(partnum)
        if(is.na(partnum))
                stop("please specify part number")
        if(partnum > numparts)
                stop("invalid part number")
        sid[partnum]
}

getChallenge <- function(email) {
        params <- list(email_address = email, response_encoding = "delim")
        result <- getForm(challenge.url, .params = params)
        s <- strsplit(result, "|", fixed = TRUE)[[1]]
        list(ch.key = s[5], state = s[7])
}

challengeResponse <- function(password, ch.key) {
        x <- paste(ch.key, password, sep = "")
        digest(x, algo = "sha1", serialize = FALSE)
}

submitSolution <- function(email, ch.resp, sid, output, signature, src = "",
                           http.version = NULL) {
        output <- as.character(base64(output))
        src <- as.character(base64(src))
        params <- list(assignment_part_sid = sid,
                       email_address = email,
                       submission = output,
                       submission_aux = src,
                       challenge_response = ch.resp,
                       state = signature)
        params <- lapply(params, URLencode)
        result <- postForm(submit.url, .params = params)
        s <- strsplit(result, "\\r\\n")[[1]]
        tail(s, 1)
}

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
reading_file(directory = "D:/Кости/R/R Script/specdata",column = 4,names_col = c("Date","sulfate","nitrate","ID"),col_m_1 = 2,col_m_2 = 3)
findat <- reading_file(directory = "D:/Кости/R/R Script/specdata",column = 4,names_col = c("Date","sulfate","nitrate","ID"),col_m_1 = 2,col_m_2 = 3)
findat

#-------------------------------------------------------------------------------

## A pair of functions that cache the inverse of a matrix


## Creates a special matrix object that can cache its inverse
makeCacheMatrix <- function( m = matrix() ) {
        
        ## Initialize the inverse property
        i <- NULL
        
        ## Method to set the matrix
        set <- function( matrix ) {
                m <<- matrix
                i <<- NULL
        }
        
        ## Method the get the matrix
        get <- function() {
                ## Return the matrix
                m
        }
        
        ## Method to set the inverse of the matrix
        setInverse <- function(inverse) {
                i <<- inverse
        }
        
        ## Method to get the inverse of the matrix
        getInverse <- function() {
                ## Return the inverse property
                i
        }
        
        ## Return a list of the methods
        list(set = set, get = get,
             setInverse = setInverse,
             getInverse = getInverse)
}


## Compute the inverse of the special matrix returned by "makeCacheMatrix"
## above. If the inverse has already been calculated (and the matrix has not
## changed), then the "cachesolve" should retrieve the inverse from the cache.
cacheSolve <- function(x, ...) {
        
        ## Return a matrix that is the inverse of 'x'
        m <- x$getInverse()
        
        ## Just return the inverse if its already set
        if( !is.null(m) ) {
                message("getting cached data")
                return(m)
        }
        
        ## Get the matrix from our object
        data <- x$get()
        
        ## Calculate the inverse using matrix multiplication
        m <- solve(data) %*% data
        
        ## Set the inverse to the object
        x$setInverse(m)
        
        ## Return the matrix
        m
}

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
                print("x is greater than zero")
        else
                print("x is less than or equal to zero")
        invisible(x)
}

printmessage <- function(x) {
        if(x > 0)
                print("x is greater than zero")
        else
                print("x is less than or equal to zero")
        invisible(x) 
        #предотвращает автоматический вывод результатов по окончанию фунции.
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


#rnorm (): генерировать случайные нормальные переменные с заданным средним и стандартным отклонением
#dnorm (): оценить плотность нормальной вероятности (с заданным значением / SD) в точке (или векторе точек)
#pnorm (): оценить кумулятивную функцию распределения для нормального распределения
#rpois (): генерирует случайные вариации Пуассона с заданной скоростью
#Функции распределения вероятности обычно имеют связанные с ними функции foud. Функции имеют префикс с ...

#d для плотности
#r для генерации случайных чисел
#p для кумулятивного распределения
#q для квантильной функции
#Работа с нормальным распределением требует использования этих функций foud

dnorm(x, mean = 0, sd = 1, log = FALSE)
pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
qnorm(p, mean = 0, sd = 1, lower.tail = TRUE, loq.p = FALSE)
rnorm(n, mean = 0, sd = 1)
#Если phi - кумулятивная функция распределения для стандартного нормального распределения, то pnorm (q) = phi (q) и qnorm = (phi ^ -1) * p.

x <- rnorm(10)
x

x <- rnorm(10, 20, 2)
x
summary(x)
set.seed(5) #Начальное число для генератора случайных чисел 
rnorm(5)

rpois(10, 1)
rpois(10, 2)
rpois(10, 20)
ppois(2, 2)
ppois(4, 2)
ppois(6, 2)


set.seed(20)
x <- rnorm(100)
e <- rnorm(100, 0, 2)
y <- 0.5 + 2 * x + e
summary(y)
plot(x, y)

#x является двоичным
set.seed(10)
x <- rbinom(100, 1, 0.5)
e <- rnorm(100, 0, 2)
y <- 0.5 + 2 * x + e
summary(y)
plot(x, y)

#имитировать модель Пуассона
set.seed(1)
x <- rnorm(100)
log.mu <- 0.5 + 0.3 * x
y <- rpois(100, exp(log.mu))
summary(y)
plot(x, y)


set.seed(1)
sample(1:10, 4)
sample(1:10, 4)
sample(letters, 5)
sample(1:10)
sample(1:10, replace = TRUE)


system.time(readLines("http://www.jhsph.edu")) #для определения времени на работу программы

hilbert <- function(n) {
        i <- 1:n
        1 / outer(i - 1, i, "+")
}
x <- hilbert(1000)
system.time(svd(x))

system.time({
        n <- 1000
        r <- numeric(n)
        for (i in 1:n) {
                x <- rnorm(n)
                r[i] <- mean(x)
        }
})

#Rprof () запускает профилировщик в R
#summaryRprof () табулирует вывод R-профайлера и вычисляет, сколько времени тратится на какую функцию
#by.total» делит время, затрачиваемое на каждую функцию, на общее время выполнения
#«by.self» делает то же самое, но сначала вычитает время, затраченное на выполнение функций выше в стеке вызовов

library(swirl)
swirl()

object.size()
names()
tail() #to view the last 6 rows.
summary()
table(plants$Active_Growth_Period)

#simulate rolling four six-sided dice: 
sample(1:6, 4, replace = TRUE)
#The probability of rolling the exact same result is (1/6)^4 = 0.00077, which is pretty small!

LETTERS
sample(LETTERS)
# Подброс несправедливой монетки
flips <- sample(c(0,1), 100, replace = TRUE, prob = c(0.3, 0.7))
rbinom(1, size = 100, prob = 0.7)
flips2 <- rbinom(100, size = 1, prob = 0.7)

#Each probability distribution in R has an r*** function (for "random"), 
#a d*** function (for "density"), a p*** (for "probability"), and q*** (for "quantile").

rnorm(10,mean = 100,sd = 25)
cm <- colMeans(my_pois)

data(cars)
plot(x = cars$speed,y = cars$dist)
plot(x = cars$speed, y = cars$dist, ylab = "Stopping Distance",xlab = "Speed")
plot(cars,main= "My Plot")
plot(cars,sub ="My Plot Subtitle")
plot(cars,xlim = c(10,15))
plot(cars,pch = 2)

data(mtcars)
boxplot(formula = mpg ~ cyl, data = mtcars)
hist(mtcars$mpg)

#-------------------------------------------------------------------------------

setwd("D:/Кости/Coursera/Coursera Data Science/2. R Programming/4W/4. Programming Assignment/1/rprog%2Fdata%2FProgAssignment3-data/")
dir()

best <- function(state, outcome){
        ## Validating the outcome string
        ## Creating a vector the diseases whose outcome we may want
        outcoames = c("heart attack", "heart filure", "pneumonia")
        if(outcome %in% outcomes == FALSE) stop("invalid outcome")
        
        ## Read outcome data
        data <- read.csv("outcome-of-care-measures.csv", colClasses = "character")
        
        ## Filter and simplify the column names
        data <- data[c(2, 7, 11, 17, 23)]
        names(data)[1] <- "name"
        names(data)[2] <- "state"
        names(data)[3] <- "heart attack"
        names(data)[4] <- "heart failure"
        names(data)[5] <- "pneumonia"
        
        ## Validating the state string
        states <- data[, 2]
        states <- unique(states)
        if(state %in% states == FALSE) stop("invalid state")
        
        ## Take only those rows with have the required state value	
        data <- data[data$state==state & data[outcome] != 'Not Available', ]
        vals <- data[, outcome]
        rowNum <- which.min(vals)
        ## Return hospital name in that state with lowest 30-day death rate
        data[rowNum, ]$name
}


best(state = "KS",outcome = "heart attack")

#-------------------------------------------------------------------------------

rankhospital <- function(state, outcome, num = "best") {
        ## Reading the outcome data
        data <- read.csv("outcome-of-care-measures.csv", colClasses = "character")
        ## Selecting the columns that are reqiured and naming them
        data <- data[c(2, 7, 11, 17, 23)]
        names(data)[1] <- "name"
        names(data)[2] <- "state"
        names(data)[3] <- "heart attack"
        names(data)[4] <- "heart failure"
        names(data)[5] <- "pneumonia"
        
        ## Validating the outcome string
        outcomes = c("heart attack", "heart failure", "pneumonia")
        if(outcome %in% outcomes == FALSE) stop("invalid outcome")
        
        ## Validating the state string
        states <- data[, 2]
        states <- unique(states)
        if(state %in% states == FALSE) stop("invalid state")
        
        ## Validating the num value
        if(num != "best" && num != "worst" && num%%1 != 0) stop("invalid num")
        
        ## Grab only those rows which matches the reqiured state value and 
        ## whose data is available    
        data <- data[data$state==state & data[outcome] != 'Not Available', ]
        
        ## Ordering the data in ascending order, first according to the names 
        ## column and then according to their ranks for the specific outcome column
        data[outcome] <- as.data.frame(sapply(data[outcome], as.numeric))
        data <- data[order(data$name, decreasing = FALSE), ]
        data <- data[order(data[outcome], decreasing = FALSE), ]
        
        ## Processing the num argument for various conditions
        vals <- data[, outcome]
        if(num == "best"){
                rowNum <- which.min(vals)
        }else if(num == "worst"){
                rowNum <- which.max(vals)
        }else{
                rowNum <- num
        }
        
        ## Return hospital name in that state with lowest 30-day death rate
        data[rowNum, ]$name
}

rankhospital(state = "KS", outcome = "heart attack")

#-------------------------------------------------------------------------------

rankall <- function(outcome, num = "best"){
        ## Reading the outcome data
        data <- read.csv("outcome-of-care-measures.csv", colClasses = "character")
        ## Selecting the columns that are reqiured and naming them
        data <- data[c(2, 7, 11, 17, 23)]
        names(data)[1] <- "name"
        names(data)[2] <- "state"
        names(data)[3] <- "heart attack"
        names(data)[4] <- "heart failure"
        names(data)[5] <- "pneumonia"
        
        ## Validating the outcome string
        outcomes = c("heart attack", "heart failure", "pneumonia")
        if(outcome %in% outcomes == FALSE) stop("invalid outcome")
        
        ## Validating the num value
        if(num != "best" && num != "worst" && num%%1 != 0) stop("invalid num")
        
        ## Grab only those rows whose data is reqiured.
        data <- data[data[outcome] != 'Not Available', ]
        
        ## Ordering the data in ascending order, first according to the names 
        ## column and then according to their ranks for the specific outcome column
        data[outcome] <- as.data.frame(sapply(data[outcome], as.numeric))
        data <- data[order(data$name, decreasing = FALSE), ]
        data <- data[order(data[outcome], decreasing = FALSE), ]
        
        ## Helper functiont to process the num argument
        getHospByRank <- function(df, s, n) {
                df <- df[df$state==s, ]
                vals <- df[, outcome]
                if(n == "best"){
                        rowNum <- which.min(vals)
                }else if(n == "worst" ){
                        rowNum <- which.max(vals)
                }else{
                        rowNum <- n
                }
                df[rowNum, ]$name
        }
        
        ## For each state, find the hospital of the given rank
        states <- data[, 2]
        states <- unique(states)
        newdata <- data.frame("hospital"=character(), "state"=character())
        for(st in states) {
                hosp <- getHospByRank(data, st, num)
                newdata <- rbind(newdata, data.frame(hospital=hosp, state=st))
        }
        
        ## Return a data frame with the hospital names and the (abbreviated) 
        ## state name
        newdata <- newdata[order(newdata['state'], decreasing = FALSE), ]
        newdata
}

rankall(outcome = "heart attack", num = "worst")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
## Getting and Cleaning Data
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

setwd("d:/Кости/R/R Script/")

# Проверка на существование папки с конкретным именем и создание, если таковой не было
if(!file.exists("data")) {
        dir.create("data")}

fileUrl <- "https://data.baltimorecity.gov/api/views/dz54-2aru/rows.csv?accessType=DOWNLOAD"
download.file(fileUrl,destfile = "D:/Кости/R/R Script/data/cameras.csv")
list.files("D:/Кости/R/R Script/")
dateDownloadet <- date()

data <- read.table(file = "data/cameras.csv",header = T,sep = ",")
data <- read.csv("data/cameras.csv")
head(data)


library(xlsx)
fileUrl <- "https://cdn.online-convert.com/example-file/spreadsheet/xlsx/example.xlsx"
download.file(fileUrl,destfile = "D:/Кости/R/R Script/data/example.xlsx")
list.files("D:/Кости/R/R Script/")

colIndex <- 2:3
rowIndex <- 1:4
data <- read.xlsx("data/cameras_2.xlsx",sheetIndex = 1,header = T,colIndex = colIndex,rowIndex = rowIndex)

write.xlsx()
read.xlsx2() # much fuster then previous but reading subset may be slightly unstable 

#XML

# Start tags <section>
# End tags </section>
# Empry tags <line-break/>

# Example:
#<Greeting> Hello,world </greeting>
# Atributes are component of the lable
#<img src="jef.jpg" alt="instructor"/>
#<step number="3"> Connect A to B, </step>

## read xml file using library(xml)
library(XML)
library(RCurl)
curlVersion()$features
curlVersion()$protocol

fileUrl <- "https://www.w3schools.com/xml/simple.xml" # for some reason Do not work !!??
fileUrl <- getURL("https://www.w3schools.com/xml/simple.xml", ssl.verifyPeer=FALSE)
doc <- xmlTreeParse(fileUrl,useInternal=TRUE)
rootNode <- xmlRoot(doc) #Корневой элемент

xmlName(rootNode) 
names(rootNode)
rootNode[[1]]
rootNode[[1]][[1]]
rootNode[[1]][[2]]
xmlSApply(rootNode,xmlValue)

# XPath 
# /node Top level node (node is the name of tags)
# //node Node at any level
# node[@attr-name] Node with an attribute name
# node[@attr-name='bob'] Node with attribute name attr-name='bob'
# stat.berkeley.edu/https://www.stat.berkeley.edu/~statcur/Workshop2/Presentations/XML.pdf

## Using xpathSApply to parse data
xpathSApply(rootNode,"//name",xmlValue) 
xpathSApply(rootNode,"//price",xmlValue) 

fileUrl <- "http://www.espn.com/nfl/team/_/name/bal/baltimore-ravens"
doc <- htmlTreeParse(fileurl, useInternal=TRUE)
teams <- xpathSApply(doc, "//div[@class='game-info']", xmlValue)
teams
score <- xpathSApply(doc, "//div[@class='score']",xmlValue)
score

## Using jsonlite to parse JSON file

library(jsonlite)
jsonData <- fromJSON("https://api.github.com/users/jtleek/repos")
names(jsonData)
names(jsonData$owner)
jsonData$owner
jsonData$owner$login
myjson <- toJSON(iris,pretty=TRUE)
cat(myjson)
iris2 <- fromJSON(myjson)
iris2
head(iris2)
#r-bloggers.com/new-package-jsonlite-a-smarter-json-encoderdecoder/

## Create data tables just like data frames

library(data.table)
DF = data.frame(x=rnorm(9),y=rep(c("a","b","c"),each=3),z=rnorm(9))
head(DF)

DT = data.table(x=rnorm(9),y=rep(c("a","b","c"),each=3), z=rnorm(9))
head(DT,3)
tables() #View all tables with data which saves in memory in now time
DT[2,]

DT[DT$y=="a",]
#Выделяем таким образом вторую и третью строку (работает не как в DataFrame)
DT[c(2,3)] 
# Колонки выделяем так же иначе

{
        x = 1
        y = 2
}
k = {print(10);5}
print(k)

# Среднее по колонкам "x" и "z"
DT[,list(mean(x),sum(z))] 
DT[,table(y)]
DT[,w:=z^2]
DT2 <- DT
DT[,y:=2]
DT[,m:={tmp <- (x+z);log2(tmp+5)}]
DT[,a:=x>0]
DT[,b:=mean(x+w),by=a]

set.seed(123);
DT <- data.table(x=sample(LETTERS[1:3],1E5,TRUE))
DT[,.N,by=x] #.N подсчитывает количество раз втречаемых аргументов
DT <- data.table(x=rep(c("a","b","c"),each=100),y=rnorm(300))
setkey(DT,x) 
DT['a']

DT1 <- data.table(x=c('a','a','b','dt1'),y=1:4)
DT2 <- data.table(x=c('a','b','dt2'),z=5:7)
setkey(DT1,x);setkey(DT2,x)
merge(DT1,DT2)

big_df <- data.frame(x=rnorm(1E6),y=rnorm(1E6))
file <- tempfile()
write.table(big_df,file = file,row.names = F,col.names = T,sep = "\t",quote = F)
system.time(fread(file))
system.time(read.table(file,header = T,sep = "\t"))

ls()
rm(list=ls())
library(dplyr)
mydf <- read.csv(path2csv, stringsAsFactors = FALSE)
cran <- tbl_df(mydf)
rm("mydf")
5:20
select(cran, r_arch:country)
select(cran, country:r_arch)
select(cran, -time)
-5:20
-(5:20)
select(cran, -(X:size))
filter(cran, package == "swirl")
filter(cran,r_version == "3.1.1", country == "US")
filter(cran, r_version <= "3.0.2", country == "IN")
filter(cran, country == "US" | country == "IN")
filter(cran, size > 100500, r_os == "linux-gnu")
is.na(c(3, 5, NA, 10))
!is.na(c(3, 5, NA, 10))
filter(cran, !is.na(r_version))
cran2 <- select(cran,size:ip_id)
arrange(cran2, ip_id)
arrange(cran2, desc(ip_id))
arrange(cran2, package, ip_id) #will first arrange by package names (ascending alphabetically), then by ip_id. This means that if there are multiple rows with the same value for package, they will be sorted by ip_id (ascending numerically)
arrange(cran2,country,desc(r_version),ip_id)
cran3 <- select(cran,ip_id,package,size)
mutate(cran3, size_mb = size / 2^20)
mutate(cran3, size_mb = size / 2^20,size_gb = size_mb/2^10)
mutate(cran3, correct_size = size + 1000)
summarize(cran, avg_bytes = mean(size))

cran <- tbl_df(mydf)
rm("mydf")
by_package <- group_by(cran,package)
summarize(by_package, mean(size))

pack_sum <- summarize(by_package,
                      count = n(),
                      unique = n_distinct(ip_id),
                      countries = n_distinct(country),
                      avg_bytes = mean(size))

quantile(pack_sum$count, probs = 0.99)
top_counts <- filter(pack_sum,count>679)
View(top_counts)
top_counts_sorted <- arrange(top_counts,desc(count))
View(top_counts_sorted)
quantile(pack_sum$unique, probs = 0.99)
top_unique <- filter(pack_sum,unique>465)
top_unique_sorted <- arrange(top_unique,desc(unique))


top_countries <- filter(pack_sum, countries > 60)
result1 <- arrange(top_countries, desc(countries), avg_bytes)
print(result1)


result2 <-
        arrange(
                filter(
                        summarize(
                                group_by(cran,
                                         package
                                ),
                                count = n(),
                                unique = n_distinct(ip_id),
                                countries = n_distinct(country),
                                avg_bytes = mean(size)
                        ),
                        countries > 60
                ),
                desc(countries),
                avg_bytes
        )

print(result2)


result3 <-
        cran %>%
        group_by(package) %>%
        summarize(count = n(),
                  unique = n_distinct(ip_id),
                  countries = n_distinct(country),
                  avg_bytes = mean(size)
        ) %>%
        filter(countries > 60) %>%
        arrange(desc(countries), avg_bytes)

# Print result to console
print(result3)

cran %>%
        select(ip_id, country, package, size) %>%
        print

cran %>%
        select(ip_id, country, package, size) %>%
        mutate(size_mb = size / 2^20) %>%
        print

cran %>%
        select(ip_id, country, package, size) %>%
        mutate(size_mb = size / 2^20) %>%
        filter(size_mb <= 0.5)

cran %>%
        select(ip_id, country, package, size) %>%
        mutate(size_mb = size / 2^20) %>%
        filter(size_mb <= 0.5) %>%
        arrange(desc(size_mb))

library(tidyr)
library(dplyr)

# Tidy data
#1: Each variable forms a column
#2: Each observation forms a row
#3: Each type of observational unit forms a table

# messy data
#1: Column headers are values, not variable names
#2: Variables are stored in both rows and columns
#3: A single observational unit is stored in multiple tables
#4: Multiple types of observational units are stored in the same table
#5: Multiple variables are stored in one column

gather(students, sex, count, -grade)

res <- gather(students2, sex_class, count, -grade)
separate(data = res,col = sex_class, into = c("sex","class"))

students2 %>%
        gather(sex_class, count, -grade) %>%
        separate(sex_class, c("sex", "class")) %>%
        print

students3 %>%
        gather(class, grade, class1:class5, na.rm = TRUE) %>%
        print

library(readr)
parse_number("class5")
students3 %>%
        gather(class, grade, class1:class5, na.rm = TRUE) %>%
        spread(test, grade) %>%
        mutate(class = parse_number(class)) %>%
        print

student_info <- students4 %>%
        select(id, name, sex) %>%
        print

student_info <- students4 %>%
        select(id, name, sex) %>%
        unique %>%
        print
gradebook <- students4 %>%
        select(id, class, midterm, final) %>%
        print

passed <- passed %>% mutate(status = "passed")
failed <- failed %>% mutate(status = "failed")
bind_rows(passed, failed)

sat %>%
        select(-contains("total")) %>%
        gather(key = part_sex, value = count, -score_range) %>%
        separate(part_sex,c("part", "sex")) %>%
        print
sat %>%
        select(-contains("total")) %>%
        gather(part_sex, count, -score_range) %>%
        separate(part_sex, c("part", "sex")) %>%
        group_by(part, sex) %>%
        mutate(total = sum(count),
               prop = count / total
        ) %>% print


## Reading from MySQL
# https://class.coursera.org/getdata-017/lecture/21
library(RMySQL)
ucscDb <- dbConnect(MySQL(),user="genome",host="genome-mysql.cse.ucsc.edu")
result <- dbGetQuery(ucscDb,"show databases;"); dbDisconnect(ucscDb);

hg19 <- dbConnect(MySQL(),user="genome", db="hg19",host="genome-mysql.cse.ucsc.edu")
allTables <- dbListTables(hg19)
length(allTables)
allTables[1:5]
# Get dimensions of a specific table
dbListFields(hg19,"affyU133Plus2")  # Get columns
dbGetQuery(hg19, "select count(*) from affyU133Plus2")  # Get rows
# Read from the table
affyData <- dbReadTable(hg19, "affyU133Plus2")
head(affyData)
# Select a specific subset
query <- dbSendQuery(hg19, "select * from affyU133Plus2 where misMatches between 1 and 3")
affyMis <- fetch(query); quantile(affyMis$misMatches)

# limit n = 10, top 10 , and need clear query
affyMisSmall <- fetch(query,n=10); dbClearResult(query);
dim(affyMisSmall)
# Don't forget to close the connection!
dbDisconnect(hg19)

ucscDb <- dbConnect(MySQL(),user="ASUS",host="C:/Program Files/Microsoft SQL Server/MSSQL13.MSSQLSERVER/MSSQL/DATA/")
result2 <- dbGetQuery(ucscDb,"show databases;"); dbDisconnect(ucscDb)

## Reading from HDF5
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")

library(rhdf5)
created = h5createFile("example.h5")
created

created = h5createGroup("example.h5","foo")
created = h5createGroup("example.h5","baa")
created = h5createGroup("example.h5","foo/foobaa")
h5ls("example.h5")

A = matrix(1:10,nr=5,nc=2)
h5write(A,"example.h5","foo/A")
B = array(seq(0.1,2.0,by=0.1),dim = c(5,2,2))
attr(B,"scale") <- "liter"
h5write(B,"example.h5","foo/foobaa/B")
h5ls("example.h5")

df = data.frame(1L:5L,seq(0,1,length.out = 5),
                c("ab","cde","fghi","a","s"),stringsAsFactors = F)
h5write(df,"example.h5","df")
h5ls("example.h5")

readA = h5read("example.h5","foo/A")
readA = h5read("example.h5","foo/foobaa/B")
readdf = h5read("example.h5","df")
readA

h5write(c(12,13,14),"example.h5","foo/A",index = list(1:3,1))
h5read("example.h5","foo/A")

## Reading from Web

fileURL <- "https://github.com/KonstantinKozhemiakin/Test_K2"
xData <- getURL(fileURL)
con = xData
htmlCode = readLines(con)
close(con)
htmlCode

#Почему то не работает:
library(XML)
url <- "https://github.com/KonstantinKozhemiakin/Test_K2"
html <- htmlTreeParse(url, useInternalNodes=T)

#А это работает:
library(XML)
library(RCurl)
fileURL <- "https://github.com/KonstantinKozhemiakin/Test_K2"
xData <- getURL(fileURL)
html <- htmlTreeParse(xData, useInternalNodes=T)
xpathSApply(html, "//title", xmlValue)
xpathSApply(html, "//td[@class='content']", xmlValue)

library(httr);html2=GET(url)
content2=content(html2,as="text")
parsedHtml=htmlParse(content2,asText = T)
xpathSApply(parsedHtml,"//title",xmlValue)

pg1 = GET("http:/httpbin.org/basic-auth/user/passwd")
pg1

pg2 = GET("http://httpbin.org/basic-auth/user/passwd",
          authenticate("user","passwd"))
pg2
names(pg2)

google = handle("http://google.com")
pg1 = GET(handle=google,path="/")
pg2 = GET(handle=google,path="search")

## Reading From APIs

#Accessing Twitter from R
library(httr)
myapp = oauth_app("twitter", #name of application
                  key="yourConsumerKeyHere",secret="yourConsumerSecretHere") #from application
sig = sign_oauth1.0(myapp,
                    token = "yourTokenHere",
                    token_secret = "yourTokenSecretHere")
homeTL = GET("https://api.twitter.com/1.1/statuses/home_timeline.json", sig)
#Converting the json object
json1 = content(homeTL)
json2 = jsonlite::fromJSON(toJSON(json1))
json2[1,1:4]

#Subsetting - quick review
set.seed(13435)
X <- data.frame("var1"=sample(1:5),"var2"=sample(6:10),"var3"=sample(11:15))
X <- X[sample(1:5),]; X$var2[c(1,3)] = NA
X

#Subsetting - quick review
X[,1]
X[,"var1"]
X[1:2,"var2"]

#Logicals ands and ors
X[(X$var1 <= 3 & X$var3 > 11),]
X[(X$var1 <= 3 | X$var3 > 15),]

#Dealing with missing values
X[which(X$var2 > 8),]

#Sorting
sort(X$var1)
sort(X$var1,decreasing=TRUE)
sort(X$var2,na.last=TRUE)

#Ordering
X[order(X$var1),]
X[order(X$var1,X$var3),]

#Ordering with plyr
library(plyr)
arrange(X,var1)
arrange(X,desc(var1))

#Adding rows and columns
X$var4 <- rnorm(5)
X

Y <- cbind(X,rnorm(5))
Y

## summarizing Data

#Getting the data from the web

if(!file.exists("./data")){dir.create("./data")}
fileUrl <- "https://data.baltimorecity.gov/api/views/k5ry-ef3g/rows.csv?accessType=DOWNLOAD"
download.file(fileUrl,destfile="./data/restaurants.csv") #удалил- ,method="curl" -не работало
restData <- read.csv("./data/restaurants.csv")

#Look at a bit of the data
head(restData,n=3)
tail(restData,n=3)
summary(restData)
str(restData)
quantile(restData$councilDistrict,na.rm=TRUE)
quantile(restData$councilDistrict,probs=c(0.5,0.75,0.9))

#Make table
table(restData$zipCode,useNA="ifany")
table(restData$councilDistrict,restData$zipCode)

#Check for missing values
sum(is.na(restData$councilDistrict))
any(is.na(restData$councilDistrict))
all(restData$zipCode > 0)

#Row and column sums
colSums(is.na(restData))
all(colSums(is.na(restData))==0)

#Values with specific characteristics
table(restData$zipCode %in% c("21212"))
table(restData$zipCode %in% c("21212","21213"))
restData[restData$zipCode %in% c("21212","21213"),]

#Cross tabs
data(UCBAdmissions)
DF = as.data.frame(UCBAdmissions)
summary(DF)

xt <- xtabs(Freq ~ Gender + Admit,data=DF)
xt

#Flat tables
warpbreaks$replicate <- rep(1:9, len = 54)
xt = xtabs(breaks ~.,data=warpbreaks)
xt
ftable(xt)

#Size of a data set
fakeData = rnorm(1e5)
object.size(fakeData)
print(object.size(fakeData),units="Mb")

## creatingNewVariables

#Getting the data from the web
if(!file.exists("./data")){dir.create("./data")}
fileUrl <- "https://data.baltimorecity.gov/api/views/k5ry-ef3g/rows.csv?accessType=DOWNLOAD"
download.file(fileUrl,destfile="./data/restaurants.csv")
restData <- read.csv("./data/restaurants.csv")

#Creating sequences
#Sometimes you need an index for your data set
s1 <- seq(1,10,by=2) ; s1
s2 <- seq(1,10,length=3); s2
x <- c(1,3,8,25,100); seq(along = x)

#Subsetting variables
restData$nearMe = restData$neighborhood %in% c("Roland Park", "Homeland")
table(restData$nearMe)

#Creating binary variables
restData$zipWrong = ifelse(restData$zipCode < 0, TRUE, FALSE)
table(restData$zipWrong,restData$zipCode < 0)

#Creating categorical variables
restData$zipGroups = cut(restData$zipCode,breaks=quantile(restData$zipCode))
table(restData$zipGroups)
table(restData$zipGroups,restData$zipCode)

#Easier cutting
library(Hmisc)
restData$zipGroups = cut2(restData$zipCode,g=4)
table(restData$zipGroups)

#Creating factor variables
restData$zcf <- factor(restData$zipCode)
restData$zcf[1:10]
class(restData$zcf)

#Levels of factor variables
yesno <- sample(c("yes","no"),size=10,replace=TRUE)
yesnofac = factor(yesno,levels=c("yes","no"))
relevel(yesnofac,ref="yes")
as.numeric(yesnofac)

#Cutting produces factor variables
library(Hmisc)
restData$zipGroups = cut2(restData$zipCode,g=4)
table(restData$zipGroups)

#Using the mutate function
library(Hmisc); library(plyr)
restData2 = mutate(restData,zipGroups=cut2(zipCode,g=4))
table(restData2$zipGroups)

#Common transforms
abs(x) #absolute value
sqrt(x) #square root
ceiling(x) #ceiling(3.475) is 4
floor(x) #floor(3.475) is 3
round(x,digits=n) #roun(3.475,digits=2) is 3.48
signif(x,digits=n) #signif(3.475,digits=2) is 3.5
cos(x) sin(x) #etc.
log(x) #natural logarithm
log2(x) log10(x) #other common logs
exp(x) #exponentiating x

## reshapingData

#Each variable forms a column
#Each observation forms a row
#Each table/file stores data about one kind of observation (e.g. people/hospitals).

#Start with reshaping
library(reshape2)
head(mtcars)

#Melting data frames (Урезать датаcет)
mtcars$carname <- rownames(mtcars)
carMelt <- melt(mtcars,id=c("carname","gear","cyl"),measure.vars=c("mpg","hp"))
head(carMelt,n=3)
tail(carMelt,n=3)

#Casting data frames
cylData <- dcast(carMelt, cyl ~ variable)
cylData
cylData <- dcast(carMelt, cyl ~ variable,mean)
cylData

#Averaging values
head(InsectSprays)
tapply(InsectSprays$count,InsectSprays$spray,sum)

#Another way - split
spIns =  split(InsectSprays$count,InsectSprays$spray)
spIns

#Another way - apply
sprCount = lapply(spIns,sum)
sprCount

#Another way - combine
unlist(sprCount)
sapply(spIns,sum)

#Another way - plyr package
library(plyr)
ddply(InsectSprays,.(spray),summarize,sum=sum(count))

#Creating a new variable
spraySums <- ?ddply(InsectSprays,.(spray),summarize,sum=ave(count,FUN=sum))
dim(spraySums)
head(spraySums)

#also the functions
#acast - for casting as multi-dimensional arrays
#arrange - for faster reordering without using order() commands
#mutate - adding new variables

## dplyr

library(dplyr)
# select - return a subset of the column of data frame
# filter - extract a sabset of rows from a data frame based on logical condition
# arrange - reorder roes of DF
# rename - rename variables in DF
# mutate - add new variables/columns or transform existing variables
# summarise/summarize - generate summary statistics of different variables in the DF,
#possibly with strata

cars2 <- as.data.frame(mtcars)
names(cars2)
head(select(cars2,wt:carb))
head(select(cars2,-(wt:carb)))

i <- match("wt", names(cars2))
j <- match("carb",names(cars2))
head(cars2[,-(i:j)])

car.f <- filter(cars2,disp>200);car.f 
car.f <- filter(cars2,disp>200 & cyl<8);car.f 

cars2 <- arrange(cars2,mpg); head(cars2);tail(cars2)
cars2 <- arrange(cars2,desc(mpg)); head(cars2);tail(cars2)

cars2 <- rename(cars2,disps = disp, wts = wt);head(cars2)

cars2 <- mutate(cars2, disps2 = disps-mean(disps,na.rm = T))
head(select(cars2,disps,disps2))

cars2 <- as.data.frame(mtcars) ; head(cars2)
cars2 <- mutate(cars2, cil_6 = factor(1*(cyl>6),labels = c("down","up")));head(cars2)

carscil <- group_by(cars2,cil_6); head(carsdisp)
summarise(carscil,disp=mean(disp),mpgs=max(mpg),wts=median(wt))

cars2 %>% mutate(cil_6 = factor(1*(cyl>6),labels = c("down","up"))) %>% 
        group_by(cil_6) %>%
        summarise(disp=mean(disp),mpgs=max(mpg),wts=median(wt)) %>%
        print

## merging Data

#Peer review data
if(!file.exists("./data")){dir.create("./data")}
fileUrl1 = "https://github.com/jtleek/modules/tree/master/03_GettingData/03_05_mergingData/data/reviews.csv"
fileUrl2 = "https://github.com/jtleek/modules/tree/master/03_GettingData/03_05_mergingData/data/solutions.csv"

download.file(fileUrl1,destfile="./data/reviews.csv")
download.file(fileUrl2,destfile="./data/solutions.csv")
reviews = read.csv("./data/reviews.csv"); solutions <- read.csv("./data/solutions.csv") # Не рабочая ссылка
reviews = read.csv("./data/reviews2.csv",sep = ";"); solutions <- read.csv("./data/solutions2.csv",sep = ";")
head(reviews,2)
head(solutions,2)

#Merging data - merge()
#Merges data frames
#Important parameters: x,y,by,by.x,by.y,all
names(reviews)
names(solutions)

#Merging data - merge()
mergedData = merge(reviews,solutions,by.x="solution_id",by.y="id",all=TRUE)
head(mergedData)

#Default - merge all common column names
intersect(names(solutions),names(reviews))
mergedData2 = merge(reviews,solutions,all=TRUE)
head(mergedData2)

#Using join in the plyr package
#Faster, but less full featured - defaults to left join, see help file for more

library(plyr)
df1 = data.frame(id=sample(1:10),x=rnorm(10))
df2 = data.frame(id=sample(1:10),y=rnorm(10))
arrange(join(df1,df2),id)

#If you have multiple data frames
df1 = data.frame(id=sample(1:10),x=rnorm(10))
df2 = data.frame(id=sample(1:10),y=rnorm(10))
df3 = data.frame(id=sample(1:10),z=rnorm(10))
dfList = list(df1,df2,df3)
join_all(dfList)

## editing Text Variables

#Fixing character vectors - tolower(), toupper()
if(!file.exists("./data")){dir.create("./data")}
fileUrl <- "https://data.baltimorecity.gov/api/views/dz54-2aru/rows.csv?accessType=DOWNLOAD"
download.file(fileUrl,destfile="./data/cameras.csv")
cameraData <- read.csv("./data/cameras.csv")
names(cameraData)

toupper(names(cameraData))
tolower(names(cameraData)) 

#Fixing character vectors - strsplit()
#Good for automatically splitting variable names
#Important parameters: x, split
splitNames = strsplit(names(cameraData),"\\.")
splitNames[[5]]
splitNames[[6]]

#Quick aside - lists
mylist <- list(letters = c("A", "b", "c"), numbers = 1:3, matrix(1:25, ncol = 5))
head(mylist)

#Quick aside - lists
mylist[1]
mylist$letters
mylist[[1]]

#Fixing character vectors - sapply()
#Applies a function to each element in a vector or list
#Important parameters: X,FUN
splitNames[[6]][2]

firstElement <- function(x){x[1]}
sapply(splitNames,firstElement)

#Peer review data
reviews = read.csv("./data/reviews2.csv",sep = ";"); solutions <- read.csv("./data/solutions2.csv",sep = ";")
head(reviews,2);head(solutions,2)

#Fixing character vectors - sub()
#Important parameters: pattern, replacement, x
names(reviews)
sub("_","",names(reviews),)

#Fixing character vectors - gsub()
testName <- "this_is_a_test"
sub("_","",testName)
gsub("_","",testName)

#Finding values - grep(),grepl()
grep("Alameda",cameraData$intersection)
table(grepl("Alameda",cameraData$intersection))
cameraData2 <- cameraData[!grepl("Alameda",cameraData$intersection),]

#More on grep()
grep("Alameda",cameraData$intersection,value=TRUE)
grep("JeffStreet",cameraData$intersection)
length(grep("JeffStreet",cameraData$intersection))

#More useful string functions
library(stringr)
nchar("Jeffrey Leek")
substr("Jeffrey Leek",1,7)
paste("Jeffrey","Leek")

#More useful string functions
paste0("Jeffrey","Leek")
str_trim("Jeff      ")

## working With Dates

#Starting simple
d1 = date()
d1
class(d1)

#Date class
d2 = Sys.Date()
d2class(d2)

#Formatting dates
#%d = day as number (0-31), %a = abbreviated weekday,%A = unabbreviated weekday, %m = month (00-12), %b = abbreviated month, %B = unabbrevidated month, %y = 2 digit year, %Y = four digit year
format(d2,"%a %b %d")

#Creating dates
x = c("1січ1960", "2січ1960", "31бер1960", "30вер1960"); z = as.Date(x, "%d%b%Y")
z
z[1] - z[2]
as.numeric(z[1]-z[2])

#Converting to Julian
weekdays(d2)
months(d2)
julian(d2)

#Lubridate
library(lubridate); ymd("20140108")
mdy("08/04/2013")
dmy("03-04-2013")

#Dealing with times
ymd_hms("2011-08-03 10:15:03")
ymd_hms("2011-08-03 10:15:03",tz="Pacific/Auckland")
?Sys.timezone

#Some functions have slightly different syntax
x = dmy(c("1jan2013", "2jan2013", "31mar2013", "30jul2013"))
wday(x[1])
wday(x[1],label=TRUE)

## Coursera Test run_analysis.R #1
setwd("D:/Кости/Coursera/Coursera Data Science/3. Getting and Cleaning Data/4W/4. Course Project/1/UCI HAR Dataset/")
library(reshape2)

# Load the various datasets
test.subject <- read.table("./test/subject_test.txt")
test.x <- read.table("./test/X_test.txt")
test.y <- read.table("./test/y_test.txt")

train.subject <- read.table("./train/subject_train.txt")
train.x <- read.table("./train/X_train.txt")
train.y <- read.table("./train/y_train.txt")

features <- read.table("./features.txt")
activity.labels <- read.table("./activity_labels.txt")

# Merge the test and train subject datasets
subject <- rbind(test.subject, train.subject)
colnames(subject) <- "subject"

# Merge the test and train labels, applying the textual labels
label <- rbind(test.y, train.y)
label <- merge(label, activity.labels, by=1)[,2]

# Merge the test and train main dataset, applying the textual headings
data <- rbind(test.x, train.x)
colnames(data) <- features[, 2]

# Merge all three datasets
data <- cbind(subject, label, data)

# Create a smaller dataset containing only the mean and std variables
search <- grep("-mean|-std", colnames(data))
data.mean.std <- data[,c(1,2,search)]

# Compute the means, grouped by subject/label
melted = melt(data.mean.std, id.var = c("subject", "label"))
means = dcast(melted , subject + label ~ variable, mean)

# Save the resulting dataset
write.table(means, file="./data/tidy_data.txt")

# Output final dataset
means

## Coursera Test run_analysis.R #2

setwd("D:/Кости/Coursera/Coursera Data Science/3. Getting and Cleaning Data/4W/4. Course Project/1/UCI HAR Dataset/")

# Import the dplyr library
library(dplyr)

# Read in the X test dataset
x.test <- read.csv("test/X_test.txt", sep="",
                   header=FALSE)

# Read in the test labels
y.test <- read.csv("test/y_test.txt", sep="",
                   header=FALSE)

# Rest in the test subject dataset
subject.test <- read.csv("test/subject_test.txt",
                         sep="", header=FALSE)

# Merge the test datasets into a single dataframe
test <- data.frame(subject.test, y.test, x.test)

# Read in the X training dataset
x.train <- read.csv("train/X_train.txt", sep="",
                    header=FALSE)

# Read in the training labels
y.train <- read.csv("train/y_train.txt", sep="",
                    header=FALSE)

# Read in the training subject dataset
subject.train <- read.csv("train/subject_train.txt",
                          sep="", header=FALSE)

# Merge test training datasets into a single dataframe
train <- data.frame(subject.train, y.train, x.train)

# Combine the training and test running datasets
run.data <- rbind(train, test)

# Remove the files we don't need anymore from
# the environment.
remove(subject.test, x.test, y.test, subject.train,
       x.train, y.train, test, train)

# Read in the measurement labels dataset
features <- read.csv("features.txt", sep="", header=FALSE)
# Convert the 2nd column into a vector
column.names <- as.vector(features[, 2])
# Apply the measurement labels as column names to the combined
# running dataset
colnames(run.data) <- c("subject_id", "activity_labels", column.names)

# Select only the columns that contain mean or standard deviations.
# Make sure to bring along the subject and label columns.
# Exclude columns with freq and angle in the name.
run.data <- select(run.data, contains("subject"), contains("label"),
                   contains("mean"), contains("std"), -contains("freq"),
                   -contains("angle"))

# Read in the activity labels dataset
activity.labels <- read.csv("activity_labels.txt", 
                            sep="", header=FALSE)

# Replace the activity codes in the trimmed down running
# dataset with the labels from the activity labels dataset.
run.data$activity_labels <- as.character(activity.labels[
        match(run.data$activity_labels, activity.labels$V1), 'V2'])

# Clean up the column names. Remove parantheses and hyphens
# from column names, both of which are invalid characters in
# column names. Also fix a set of columns that repeat the
# word "Body".
setnames(run.data, colnames(run.data), gsub("\\(\\)", "", colnames(run.data)))
setnames(run.data, colnames(run.data), gsub("-", "_", colnames(run.data)))
setnames(run.data, colnames(run.data), gsub("BodyBody", "Body", colnames(run.data)))

# Group the running data by subject and activity, then
# calculate the mean of every measurement.
run.data.summary <- run.data %>%
        group_by(subject_id, activity_labels) %>%
        summarise_each(funs(mean))

# Write run.data to file
write.table(run.data.summary, file="run_data_summary.txt", row.name=FALSE)

setwd("D:/Кости/R/Git_R/Test_K2")

