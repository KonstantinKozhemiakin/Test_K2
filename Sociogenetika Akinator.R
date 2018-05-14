
# Entropy

library(entropy)


info <- function(CLASS.FREQ){
        freq.class <- CLASS.FREQ
        info <- 0
        for(i in 1:length(freq.class)){
                if(freq.class[[i]] != 0){ # zero check in class
                        entropy <- -sum(freq.class[[i]] * log2(freq.class[[i]]))  #I calculate the entropy for each class i here
                }else{ 
                        entropy <- 0
                } 
                info <- info + entropy # sum up entropy from all classes
        }
        return(info)
}


buys <- c("no", "no", "yes", "yes", "yes", "no", "yes", "no", "yes", "yes", "yes", "yes", "yes", "no")
credit <- c("fair", "excellent", "fair", "fair", "fair", "excellent", "excellent", "fair", "fair", "fair", "excellent", "excellent", "fair", "excellent")
student <- c("no", "no", "no","no", "yes", "yes", "yes", "no", "yes", "yes", "yes", "no", "yes", "no")
income <- c("high", "high", "high", "medium", "low", "low", "low", "medium", "low", "medium", "medium", "medium", "high", "medium")
age <- c(25, 27, 35, 41, 48, 42, 36, 29, 26, 45, 23, 33, 37, 44) # we change the age from categorical to numeric



buys <- c("no", "no", "yes", "yes", "yes", "no", "yes", "no", "yes", "yes", "yes", "yes", "yes", "no")
freqs <- table(buys)/length(buys)
info(freqs)


entropy.empirical(freqs, unit="log2")

#------------------------------------------------------------------------

# Bayesian method

library(rjags)


#-------------------------------------------------------------------------

# Permutation 1 Col

Soc <- c("Analitic","Upravlenec","Ispolnitel","Komunikator")
Intel <- c("Funkcional","Praktik","Emocii","Sensori")

permutations <- function(n){
        if(n==1){
                return(matrix(1))
        } else {
                sp <- permutations(n-1)
                p <- nrow(sp)
                A <- matrix(nrow=n*p,ncol=n)
                for(i in 1:n){
                        A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
                }
                return(A)
        }
}

Perm_1 <- matrix(Soc[permutations(4)],ncol=4)
Perm_2 <- matrix(Intel[permutations(4)],ncol=4)
Tipazh <- rbind(Perm_1,Perm_2)

#-------------------------------------------------------------------------

## Permutation 2 Col

Soc <- c("Analitic","Upravlenec","Ispolnitel","Komunikator")
Intel <- c("Funkcional","Praktik","Emocii","Sensori")

library(combinat)
M <- matrix(c(Soc,Intel), nrow=2, byrow=TRUE)
pcM <- permn(ncol(M))
expP <- expand.grid(1:length(pcM), 1:length(pcM))

Tipazh_2 <- Map(
        function(a,b) rbind( M[1, pcM[[a]]], M[2, pcM[[a]]] ),
        expP[,1],
        expP[,2]
)

Tipazh_2

#Vsego perestanovok   lapply(permn(M),matrix,nrow=2)

#freqs <- table(Tipazh)/length(Tipazh)
#entropy.empirical(freqs, unit="log2")

#-------------------------------------------------------------------------

# Classes

Class <- c("1_Analitic","1_Upravlenec","1_Ispolnitel","1_Komunikator","1_Funkcional","1_Praktik","1_Emocii","1_Sensori",
           "2_Analitic","2_Upravlenec","2_Ispolnitel","2_Komunikator","2_Funkcional","2_Praktik","2_Emocii","2_Sensori",
           "3_Analitic","3_Upravlenec","3_Ispolnitel","3_Komunikator","3_Funkcional","3_Praktik","3_Emocii","3_Sensori",
           "4_Analitic","4_Upravlenec","4_Ispolnitel","4_Komunikator","4_Funkcional","4_Praktik","4_Emocii","4_Sensori")

#-------------------------------------------------------------------------

#Table of classes and Permutation

z <- as.table(Tipazh)
z[c(1,0)]
str(z)

Class_tipazh <- function(c){
        for (i in c){
                
        }
        
}
