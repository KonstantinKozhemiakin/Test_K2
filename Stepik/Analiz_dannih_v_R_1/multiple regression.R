#
# multiple linear regression
#




# numeric predictors

fit <- lm(Fertility ~ Examination + Catholic, data = swiss)
summary(fit)


fit2 <- lm(Fertility ~ Examination*Catholic, data = swiss)
summary(fit2)


confint(fit2)


# categorical predictors

hist(swiss$Catholic, col = 'red')

swiss$religious <- ifelse(swiss$Catholic > 60, 'Lots', 'Few')
swiss$religious <- as.factor(swiss$religious)

fit3 <- lm(Fertility ~ Examination + religious, data = swiss)
summary(fit3)

fit4 <- lm(Fertility ~ religious*Examination, data = swiss)
summary(fit4)

# plots

ggplot(swiss, aes(x = Examination, y = Fertility)) + 
  geom_point() 

ggplot(swiss, aes(x = Examination, y = Fertility)) + 
  geom_point() + 
  geom_smooth()

ggplot(swiss, aes(x = Examination, y = Fertility)) + 
  geom_point() + 
  geom_smooth(method = 'lm')

ggplot(swiss, aes(x = Examination, y = Fertility, col = religious)) + 
  geom_point() 

ggplot(swiss, aes(x = Examination, y = Fertility, col = religious)) + 
  geom_point()  + 
  geom_smooth()

ggplot(swiss, aes(x = Examination, y = Fertility, col = religious)) + 
  geom_point()  + 
  geom_smooth(method = 'lm')


#

fit5 <- lm(Fertility ~ religious*Infant.Mortality*Examination, data = swiss)
summary(fit5)


# model comparison

rm(swiss)
swiss <- data.frame(swiss)

fit_full <- lm(Fertility ~ ., data = swiss)
summary(fit_full)

fit_reduced1 <- lm(Fertility ~ Infant.Mortality + Examination + Catholic + Education, data = swiss)
summary(fit_reduced1)

anova(fit_full, fit_reduced1)

fit_reduced2 <- lm(Fertility ~ Infant.Mortality + Education + Catholic + Agriculture, data = swiss)
summary(fit_reduced2)

anova(fit_full, fit_reduced2)


# model selection

optimal_fit <-  step(fit_full, direction = 'backward')
summary(optimal_fit)


#Памятка по интерпретации результатов регрессионного анализа с категориальными и непрерывными переменными
#Модель для примера: 
        
#DV ~ IV_numeric * IV_categorical

#IV_categorical - фактор с двумя уровнями (Level1 и Level2)

#Коэффициенты:
        
#Intercept — предсказанное значение DV для первого уровня IV_categorical с учётом того, что IV_numeric равна нулю.

#IV_numeric — насколько изменяется предсказанное значение DV при увеличении IV_numeric на одну единицу в группе, соответствующей первому уровню IV_categorical

#IV_categoricalLevel2 — насколько изменяется предсказанное значение DV при переходе от первого уровня IV_categorical ко второму уровню. С учётом того, что IV_numeric равна нулю. 

#IV_numeric:IV_categoricalLevel2 — насколько сильнее (или слабее) изменяется предсказанное значение DV при увеличении IV_numeric на одну единицу в группе, соответствующей второму уровню IV_categorical, по сравнению с первым уровнем. 

#Как предсказывать значения в новом датасете на основе полученных коэффициентов

#1). Предположим у нас есть новый объект, про который мы знаем, что он принадлежит к группе, соответствующей IV_categorical (Level1) и измеренный у него IV_numeric составил 10:
        
#Предсказанное значение DV = Intercept + 10 * IV_numeric
#        
#2). Предположим у нас есть новый объект, про который мы знаем, что он принадлежит к группе, соответствующей IV_categorical (Level2) и измеренный у него IV_numeric составил 6:
#        
#Предсказанное значение DV = Intercept + IV_categoricalLevel2 + 6 * (IV_numeric + IV_numeric:IV_categoricalLevel2)

