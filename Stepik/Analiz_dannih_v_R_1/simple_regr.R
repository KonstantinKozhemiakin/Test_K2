df  <- mtcars
df_numeric  <- df[,c(1,3:7)]

fit  <- lm(mpg ~ hp, df)
summary(fit)

ggplot(df, aes(hp, mpg))+
  geom_point(size = 5)+
  geom_smooth(method = "lm")+
  facet_grid(.~cyl)

ggplot(df, aes(hp, mpg))+
  geom_smooth(method = "lm", se = F)+
  facet_grid(.~cyl)

fitted_values_mpg  <- data.frame(mpg = df$mpg, fitted = fit$fitted.values )

new_hp <- data.frame(hp = c(100, 150, 129, 300))
new_hp$mpg  <- predict(fit, new_hp)

predict(fit, new_hp)


##################################

my_df  <- mtcars
my_df$cyl  <- factor(my_df$cyl, labels = c("four", "six", "eight"))
fit  <- lm(mpg ~ cyl, my_df)


#Памятка
cor.test(mtcars$mpg, mtcars$disp) # Расчет корреляции Пирсона 

cor.test(~ mpg + disp, mtcars) # запись через формулу

cor.test(mtcars$mpg, mtcars$disp, method = "spearman") # Расчет корреляции Спирмена 

cor.test(mtcars$mpg, mtcars$disp, method = "kendall") # Расчет корреляции Кендала 

cor(iris[, -5]) # построение корреляционной матрицы

fit <- lm(mpg ~ disp, mtcars) # построение линейной регрессии 

fit$coefficients # коэффициенты регрессии 

fit$fitted.values # предсказанные значения зависимой переменной 

#При наличии одинаковых значений в переменных расчет непараметрических корреляций будет сопровождаться предупреждением о невозможности рассчитать точное значение p - value.

#Если в ваших данных есть одинаковые наблюдения, но вы хотите рассчитать непараметрическую корреляцию, используйте функцию spearman_test  из пакета coin

library(coin)
spearman_test(~ mpg + disp, mtcars)

#Обратите внимание на различия в графиках. То что в первом aes() будет распространяться на все слои. А то, что в aes() конкретного geom - только на него.

ggplot(mtcars, aes(mpg, disp, col = factor(am)))+
        geom_point()+
        geom_smooth()

ggplot(mtcars, aes(mpg, disp))+
        geom_point(aes(col = factor(am)))+
        geom_smooth()

ggplot(mtcars, aes(mpg, disp))+
        geom_point()+
        geom_smooth(aes(col = factor(am)))
