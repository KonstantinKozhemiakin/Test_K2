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


#������� �� ������������� ����������� �������������� ������� � ��������������� � ������������ �����������
#������ ��� �������: 
        
#DV ~ IV_numeric * IV_categorical

#IV_categorical - ������ � ����� �������� (Level1 � Level2)

#������������:
        
#Intercept � ������������� �������� DV ��� ������� ������ IV_categorical � ������ ����, ��� IV_numeric ����� ����.

#IV_numeric � ��������� ���������� ������������� �������� DV ��� ���������� IV_numeric �� ���� ������� � ������, ��������������� ������� ������ IV_categorical

#IV_categoricalLevel2 � ��������� ���������� ������������� �������� DV ��� �������� �� ������� ������ IV_categorical �� ������� ������. � ������ ����, ��� IV_numeric ����� ����. 

#IV_numeric:IV_categoricalLevel2 � ��������� ������� (��� ������) ���������� ������������� �������� DV ��� ���������� IV_numeric �� ���� ������� � ������, ��������������� ������� ������ IV_categorical, �� ��������� � ������ �������. 

#��� ������������� �������� � ����� �������� �� ������ ���������� �������������

#1). ����������� � ��� ���� ����� ������, ��� ������� �� �����, ��� �� ����������� � ������, ��������������� IV_categorical (Level1) � ���������� � ���� IV_numeric �������� 10:
        
#������������� �������� DV = Intercept + 10 * IV_numeric
#        
#2). ����������� � ��� ���� ����� ������, ��� ������� �� �����, ��� �� ����������� � ������, ��������������� IV_categorical (Level2) � ���������� � ���� IV_numeric �������� 6:
#        
#������������� �������� DV = Intercept + IV_categoricalLevel2 + 6 * (IV_numeric + IV_numeric:IV_categoricalLevel2)

