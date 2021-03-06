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


#�������
cor.test(mtcars$mpg, mtcars$disp) # ������ ���������� ������� 

cor.test(~ mpg + disp, mtcars) # ������ ����� �������

cor.test(mtcars$mpg, mtcars$disp, method = "spearman") # ������ ���������� �������� 

cor.test(mtcars$mpg, mtcars$disp, method = "kendall") # ������ ���������� ������� 

cor(iris[, -5]) # ���������� �������������� �������

fit <- lm(mpg ~ disp, mtcars) # ���������� �������� ��������� 

fit$coefficients # ������������ ��������� 

fit$fitted.values # ������������� �������� ��������� ���������� 

#��� ������� ���������� �������� � ���������� ������ ����������������� ���������� ����� �������������� ��������������� � ������������� ���������� ������ �������� p - value.

#���� � ����� ������ ���� ���������� ����������, �� �� ������ ���������� ����������������� ����������, ����������� ������� spearman_test  �� ������ coin

library(coin)
spearman_test(~ mpg + disp, mtcars)

#�������� �������� �� �������� � ��������. �� ��� � ������ aes() ����� ���������������� �� ��� ����. � ��, ��� � aes() ����������� geom - ������ �� ����.

ggplot(mtcars, aes(mpg, disp, col = factor(am)))+
        geom_point()+
        geom_smooth()

ggplot(mtcars, aes(mpg, disp))+
        geom_point(aes(col = factor(am)))+
        geom_smooth()

ggplot(mtcars, aes(mpg, disp))+
        geom_point()+
        geom_smooth(aes(col = factor(am)))
