#' ---
#' title: "Final Project"
#' author: Shihan Qian
#' date: "`r format(Sys.time(), '%d %B, %Y')`"
#' ---
#' 
#' 
#' library loaded
library(car)
library(MASS)
library(glmnet)
library(splines)
library(ISLR)

#data loaded
data_tr <- read.table("data_tr.txt", header = TRUE, sep = "\t", dec = ".")[,-1]
head(data_tr)

####Section one: Data Exploration
summary(data_tr)
#We use the correlation test to check for the correlation between variables
round(cor(data_tr),digits = 2)
#We notice that hequity and hmort and hval have high correlation because hequity = hval-hmort, so we can decided to drop hequity.
#Also, we find that nohs+hs+smcol+col = 1, it's perfect collinearity, therefore, we decide to drop the col dummy. 
new_df = subset(data_tr, select = -c(hequity, col) )
#Shape of the data
y = new_df$tw
hist(y) #check the histogram for total wealth, we can see that the total wealth is heavily skewed to the right.

#check the histogram for individual retirement account (ira), non-401k financial assets (nifa), home mortgae(hmort),income(inc), given the high correlation between these variable and total wealth 
par(mfrow=c(2,3))
hist(new_df$ira)
hist(new_df$nifa)
hist(new_df$hmort)
hist(new_df$hval)
hist(new_df$inc)
hist(new_df$age)
#explore the basic relationship between total wealth and the individual independent vairable
plot(new_df$tw~new_df$ira)
abline(lm(new_df$tw~new_df$ira), col = 'red')
plot(new_df$tw~new_df$nifa)
abline(lm(new_df$tw~new_df$nifa), col = 'red')
plot(new_df$tw~new_df$hmort)
abline(lm(new_df$tw~new_df$hmort), col = 'red')
plot(new_df$tw~new_df$hval)
abline(lm(new_df$tw~new_df$hval), col = 'red')
plot(new_df$tw~new_df$inc)
abline(lm(new_df$tw~new_df$inc), col = 'red')
plot(new_df$tw~new_df$age)
abline(lm(new_df$tw~new_df$age), col = 'red')

# Multicollinearity check
all_var_model <- lm(tw ~ ., data=new_df)
summary(all_var_model)
vif_val = vif(all_var_model)
par(mfrow=c(1,1))
#visualize the vif value of the model
print(vif_val)
barplot(vif_val, las = 2, cex.names = 1, main = "VIF Values")
#multicollinearity problem happens when the covarients are highly correlated with each other, and a VIF value above 5 is high VIF value, it should be removed from the model.
df_1 = subset(new_df, select = -c(hs, nohs,educ))

#Run basic Regression on the data
#basic model based on the family situation
regression_1 <- lm(tw ~ age + marr + fsize + twoearn + hval, data = df_1)
par(mfrow=c(2,2))
plot(regression_1)
# a lot of outliers in the regression_1 given the huge deviation of dots in the qq plots
#normal qq plot is important, a perfect fit model, the qq plot will be a straight line
#model 2
#basic model based on retirement
regression_2 <- lm(tw~ira + e401 +age, data = df_1)
plot(regression_2)
#model 3
round(cor(df_1),digits = 2)
#model based on the variables that have high correlation with total wealth
regression_3 <- lm(tw~ira+nifa +hval +inc, data = df_1)
plot(regression_3)
#We can see that the outlier in the QQ plot decreases some compared to the first two model
#model 4
#model based on all the variables
regression_4 <- lm(tw~., data = df_1)
plot(regression_4)

####Section two: Model selection
#In Depth Analyses
#Comparing models with the dataset that with/without the high VIF variable
#using 10 folds corss validation to compare 4 model between forward stepwise,backward stepwise, lasso, and ridge.
n <- nrow(df_1)
k<-10
set.seed(1234)
ii <- sample(rep(1:k,floor(n/k)))
y <- df_1$tw

mspe1 <- mspe2 <- mspe3 <- mspe4 <- vector()
for (j in 1:k){
  
  hold <- (ii == j)
  train <- (ii != j)
  
  y_te <- y[hold]
  x_te <- df_1[hold,-1]
  y_tr <- y[train]
  x_tr <- df_1[train,-1]
  
  full.model <-lm(tw ~., data = df_1[train,]) 
  null.model <-lm(tw ~1, data = df_1[train,])
  reg_stepwise_forwards <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
  reg_stepwise_backwards <- stepAIC(full.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'backward')
  pr_forwards <- predict(reg_stepwise_forwards, newdata = df_1[hold,])
  pr_backwards <- predict(reg_stepwise_backwards, newdata = df_1[hold,])
  
  lambdas <- exp(seq(-4,1,by=0.1))
  LASSO.cv <- cv.glmnet(x = as.matrix(x_tr), y = y_tr, lambda = lambdas, alpha = 1)
  la_model <- glmnet(as.matrix(x_tr), y_tr, alpha = 0)
  best_lambda = LASSO.cv$lambda.min
  RR.cv <- cv.glmnet(x = as.matrix(x_tr), y = y_tr, lambda = lambdas, alpha = 0)
  best_lambda_rr = RR.cv$lambda.min
  rr_model <- glmnet(as.matrix(x_tr), y_tr, alpha = 0)
  pr_lasso <-  predict(la_model, s = best_lambda, newx = as.matrix(x_te))
  pr_ridge <- predict(rr_model, s = best_lambda_rr, newx = as.matrix(x_te))
  
  mspe1[j] <- mean((y_te - pr_forwards)^2)
  mspe2[j] <- mean((y_te - pr_backwards)^2)
  mspe3[j] <- mean((y_te - pr_lasso)^2)
  mspe4[j] <- mean((y_te - pr_ridge)^2)
}

## Return the mspe
print(mean(mspe1))
print(mean(mspe2))
print(mean(mspe3))
print(mean(mspe4))

#The first model forwards stepwise give the lowest mean, which is 1732037235
full.model <-lm(tw ~., data = df_1) 
null.model <-lm(tw ~1, data = df_1)
for_model <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
summary(for_model)
plot(for_model)

#Comparing models with the dataset that doesn't drop the high VIF variable
n <- nrow(data_tr)
k<-10
set.seed(12345)
ii <- sample(rep(1:k,floor(n/k)))
y <- data_tr$tw

mspe1 <- mspe2 <- mspe3 <- mspe4 <- vector()
for (j in 1:k){
  
  hold <- (ii == j)
  train <- (ii != j)
  
  y_te <- y[hold]
  x_te <- data_tr[hold,-1]
  y_tr <- y[train]
  x_tr <- data_tr[train,-1]
  
  full.model <-lm(tw ~., data = data_tr[train,]) 
  null.model <-lm(tw ~1, data = data_tr[train,])
  reg_stepwise_forwards <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
  reg_stepwise_backwards <- stepAIC(full.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'backward')
  pr_forwards <- predict(reg_stepwise_forwards, newdata = data_tr[hold,])
  pr_backwards <- predict(reg_stepwise_backwards, newdata = data_tr[hold,])
  
  lambdas <- exp(seq(-4,1,by=0.1))
  LASSO.cv <- cv.glmnet(x = as.matrix(x_tr), y = y_tr, lambda = lambdas, alpha = 1)
  la_model <- glmnet(as.matrix(x_tr), y_tr, alpha = 0)
  best_lambda = LASSO.cv$lambda.min
  RR.cv <- cv.glmnet(x = as.matrix(x_tr), y = y_tr, lambda = lambdas, alpha = 0)
  best_lambda_rr = RR.cv$lambda.min
  rr_model <- glmnet(as.matrix(x_tr), y_tr, alpha = 0)
  pr_lasso <-  predict(la_model, s = best_lambda, newx = as.matrix(x_te))
  pr_ridge <- predict(rr_model, s = best_lambda_rr, newx = as.matrix(x_te))
  
  mspe1[j] <- mean((y_te - pr_forwards)^2)
  mspe2[j] <- mean((y_te - pr_backwards)^2)
  mspe3[j] <- mean((y_te - pr_lasso)^2)
  mspe4[j] <- mean((y_te - pr_ridge)^2)
}

## Return the mspe
print(mean(mspe1))
print(mean(mspe2))
print(mean(mspe3))
print(mean(mspe4))

#The first model forwards stepwise give the lowest mean, which is 1733283481,which is a little bit higher than the model without high VIF,
#Therefore, we would prefer the first dataset that dropped the high VIF vairables.
full.model <-lm(tw ~., data = data_tr) 
null.model <-lm(tw ~1, data = data_tr)
for_model <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
summary(for_model)
plot(for_model)


####Transformation of the predictors
#Splines
#We can see that home mortgage, home value, and age has a flat relationship shown in the first section:data exploration.Therefore, we decide to to make splines of these variables.
plot(new_df$tw~new_df$hmort)
abline(lm(new_df$tw~new_df$hmort), col = 'red')
plot(new_df$tw~new_df$hval)
abline(lm(new_df$tw~new_df$hval), col = 'red')
plot(new_df$tw~new_df$age)
abline(lm(new_df$tw~new_df$age), col = 'red')
#Using for loop to choose the optimal degree of freedom
n <- 5
m <- 20
i <- 1
mse <- vector()
for (j in n:m){
  spline_hmort <- bs(df_1$hmort, df = j)
  spline_hmort <- as.data.frame(spline_hmort)
  names(spline_hmort) <- paste0("spline_hmort", 1:j)
  
  df_1_transform <- cbind(df_1, spline_hmort)
  df_1_transform <- subset(df_1_transform, select = -hmort)

  full.model <-lm(tw ~., data = df_1_transform[1:7139,]) 
  null.model <-lm(tw ~1, data = df_1_transform[1:7139,])
  for_model_trans <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
  
  pre <- predict(for_model_trans, newdata = df_1_transform[7139:7930,])
  mse[i] <- mean((df_1_transform$tw[7139:7930] - pre)^2)
  i <- i+ 1
}
min(mse)
print(mse)
#For the home mortage, the smallest mse is with degree of freedom of 20


plot(df_1$tw ~ df_1$hmort)
lines(df_1$hmort, predict(lm(tw~bs(hmort, df = 20), data = df_1)),col = 'red')
spline_hmort <- bs(df_1$hmort, df = 20)
spline_hmort <- as.data.frame(spline_hmort)
names(spline_hmort) <- paste0("spline_hmort", 1:20)

df_1_transform <- cbind(df_1, spline_hmort)
df_1_transform <- subset(df_1_transform, select = -hmort)
dim(df_1_transform)


n <- 5
m <- 20
i <- 1
mse <- vector()
for (j in n:m){
  spline_hval <- bs(df_1$hval, df = j)
  spline_hval <- as.data.frame(spline_hval)
  names(spline_hval) <- paste0("spline_hval", 1:j)
  
  df_1_transform_1 <- cbind(df_1_transform, spline_hval)
  df_1_transform_1 <- subset(df_1_transform_1, select = -hval)
  
  full.model <-lm(tw ~., data = df_1_transform_1[1:7139,]) 
  null.model <-lm(tw ~1, data = df_1_transform_1[1:7139,])
  for_model_trans <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
  
  pre <- predict(for_model_trans, newdata = df_1_transform_1[7139:7930,])
  mse[i] <- mean((df_1_transform_1$tw[7139:7930] - pre)^2)
  i <- i+ 1
}
min(mse)
print(mse)
#For the home value, the smallest mse is with degree of freedom of 17

plot(df_1$tw ~ df_1$hval)
lines(df_1$hval, predict(lm(tw~bs(hval, df = 17), data = df_1)),col = 'red')
spline_hval <- bs(df_1$hval, df = 17)
spline_hval <- as.data.frame(spline_hval)
names(spline_hval) <- paste0("spline_hval", 1:17)

df_1_trans_1 <- cbind(df_1_transform, spline_hval)
df_1_trans_1 <- subset(df_1_trans_1, select = -hval)
dim(df_1_trans_1)

n <- 5
m <- 10
i <- 1
mse <- vector()
for (j in n:m){
  spline_age <- bs(df_1$age, df = j)
  spline_age <- as.data.frame(spline_age)
  names(spline_age) <- paste0("spline_age", 1:j)
  
  df_1_transform_2 <- cbind(df_1_trans_1, spline_age)
  df_1_transform_2 <- subset(df_1_transform_2, select = -age)
  
  full.model <-lm(tw ~., data = df_1_transform_2[1:7139,]) 
  null.model <-lm(tw ~1, data = df_1_transform_2[1:7139,])
  for_model_trans <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
  
  pre <- predict(for_model_trans, newdata = df_1_transform_2[7139:7930,])
  mse[i] <- mean((df_1_transform_2$tw[7139:7930] - pre)^2)
  i <- i+ 1
}
min(mse)
print(mse)

plot(df_1$tw ~ df_1$age)
lines(df_1$age, predict(lm(tw~bs(age, df = 7), data = df_1)),col = 'red')
spline_age <- bs(df_1$age, df = 7)
spline_age <- as.data.frame(spline_age)
names(spline_age) <- paste0("spline_age", 1:7)

df_1_trans_2 <- cbind(df_1_trans_1, spline_age)
df_1_trans_2 <- subset(df_1_trans_2, select = -age)
dim(df_1_trans_2)
#test for the model selection
full.model <-lm(tw ~., data = df_1_trans_2) 
null.model <-lm(tw ~1, data = df_1_trans_2)
for_model_trans <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
summary(for_model_trans)
plot(for_model_trans)
summary(for_model_trans)$r.squared



#Interaction
#interaction between ira and age because ira is individual retirement account and it will vary very differently at different age since it's a variable related to retirement.
#Thus, we try to include this interaction
interaction <- df_1_trans_1$ira*df_1_trans_1$age
names(interaction) <- paste0('inter_', names(interaction))
#doesn't add on to the previous spline transformation
df_1_int <- cbind(df_1, interaction)
full.model <-lm(tw ~., data = df_1_int) 
null.model <-lm(tw ~1, data = df_1_int)
for_model_int <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
summary(for_model_int)
plot(for_model_int)
summary(for_model_int)$r.squared

#add on to the previous spline transformation
df_1_inter <- cbind(df_1_trans_2, interaction)
dim(df_1_inter)
full.model <-lm(tw ~., data = df_1_inter) 
null.model <-lm(tw ~1, data = df_1_inter)
for_model_inter <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
summary(for_model_inter)
plot(for_model_inter)
summary(for_model_inter)$r.squared



#interaction between male and marry status because "Married couples hold substantially higher wealth than single persons" and married man will have even higher avergae income,
#therefore, we decide to include the interaction between male and married status
interaction_1 <- df_1_inter$male*df_1_inter$marr
names(interaction_1) <- paste0('inter_', names(interaction_1))
#add on to the previous model with ira and age interaction
df_1_int_1 <- cbind(df_1_inter, interaction_1)
full.model <-lm(tw ~., data = df_1_int_1) 
null.model <-lm(tw ~1, data = df_1_int_1)
for_model_int_1 <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
summary(for_model_int_1)
plot(for_model_int_1)
summary(for_model_int_1)$r.squared



  
#using polynomial transformation
#We only include variable that has correlation around 0.5 or above
X <- subset(df_1, select = -c(tw,e401,hmort,male,twoearn,smcol,age,fsize,marr))


poly_df_1 <- as.data.frame(cbind(poly(as.matrix(X),degree = 3)))
full.model <-lm(df_1$tw ~., data = poly_df_1) 
null.model <-lm(df_1$tw ~1, data = poly_df_1)
for_model_poly <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
back_model_poly<- stepAIC(full.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'backward')
summary(for_model_poly)
summary(back_model_poly)
plot(for_model_poly)
summary(for_model_poly)$r.squared

#using polynomial transformation to select the optimal degree for only one variable
#explore the variable family size and we can see on the graph, there is actually a curve trend at family size of two, so linear regression may not be a best fit
#Therefore, we decide to do some polynomial transformation on family size
plot(df_1$tw~df_1$fsize)
abline(lm(df_1$tw~df_1$fsize), col = 'red')
n <- nrow(df_1)
k<-3
set.seed(123456)
ii <- sample(rep(1:k,floor(n/k)))
y <- df_1$tw
mspe<- mspe2<- vector()
for (j in 1:k){
  hold <- (ii == j)
  train <- (ii != j)
  
  y_te <- y[hold]
  
  poly_2 <- poly(df_1_int_1$fsize,2)
  poly_2 <- as.data.frame(poly_2)
  names(poly_2) <- paste0("poly_2", 1:2)
  
  df_poly_2 <- cbind(df_1_int_1, poly_2)
  df_poly_2 <- subset(df_poly_2, select = -fsize)
  
  poly_5 <- poly(df_1_int_1$fsize,5)
  poly_5 <- as.data.frame(poly_5)
  names(poly_5) <- paste0("poly_5", 1:5)
  
  df_poly_5 <- cbind(df_1_int_1, poly_5)
  df_poly_5 <- subset(df_poly_5, select = -fsize)
  
  full.model_2 <-lm(tw ~., data = df_poly_2[train,]) 
  null.model_2 <-lm(tw ~1, data = df_poly_2[train,])
  for_model_poly_2 <- stepAIC(null.model_2,scope = list(lower = null.model_2,upper = full.model_2), trace = FALSE, direction = 'forward')
  
  full.model_5 <-lm(tw ~., data = df_poly_5[train,]) 
  null.model_5 <-lm(tw ~1, data = df_poly_5[train,])
  for_model_poly_5 <- stepAIC(null.model_5,scope = list(lower = null.model_5,upper = full.model_5), trace = FALSE, direction = 'forward')
  
  pr2 <- predict(for_model_poly_2, newdata = df_poly_2[hold,])
  pr5 <- predict(for_model_poly_5, newdata = df_poly_5[hold,])
  
  mspe[j] <- mean((y_te - pr2)^2)
  mspe2[j] <- mean((y_te - pr5)^2)
  
}
mean(mspe)
mean(mspe2)
#It gives the same mse, therefore, we check the model
summary(for_model_poly_2)
#It doesn't include the polynomial transformation of fsize by using forward step wise selection, thus we know that it may not be a wise choice to include the polynomial transformation of variable fsize


#using GAM(Generalized Additive Model)
df_1
#I check for nifa ,ira and inc and find that they all have data concentrated in one part, therefore, i decide to use the ns() to generate generate Natural Cubic Splines on income and nifa.
plot(df_1$tw~df_1$nifa)
abline(lm(df_1$tw~df_1$nifa), col = 'red')
plot(df_1$tw~df_1$inc)
abline(lm(df_1$tw~df_1$inc), col = 'red')
plot(df_1$tw~df_1$ira)
abline(lm(df_1$tw~df_1$ira), col = 'red')
m<- 1
mspe<- vector()
for(i in 3:10){
  gam_df_1 <- lm(tw ~ ira + ns(age ,i)+e401 +nifa+ inc + hmort+ hval + male+ twoearn + smcol + fsize + marr ,data=df_1) 
  pr_gam <- predict(gam_df_1)
  mspe[m] <- mean((df_1$tw - pr_gam)^2)
  m <- m+1
}
print(mspe)
print(min(mspe))
#We choose the age of freedom of 10
m<- 1
mspe<- vector()
for(i in 3:5){
  gam_df_1 <- lm(tw ~ ira + ns(age ,10)+e401 +ns(nifa, i)+ inc + hmort+ hval + male+ twoearn + smcol + fsize + marr ,data=df_1) 
  pr_gam <- predict(gam_df_1)
  mspe[m] <- mean((df_1$tw - pr_gam)^2)
  m <- m+1
}
print(mspe)
print(min(mspe))
#We choose the degree of freedom of 5
m<- 1
mspe<- vector()
for(i in 3:5){
  gam_df_1 <- lm(tw ~ ira + ns(age ,10)+e401 +ns(nifa, 5)+ ns(inc, i) + hmort+ hval + male+ twoearn + smcol + fsize + marr ,data=df_1) 
  pr_gam <- predict(gam_df_1)
  mspe[m] <- mean((df_1$tw - pr_gam)^2)
  m <- m+1
}
print(mspe)
print(min(mspe))
#We choose the degree of freedom of 5

gam_df_1 <- lm(tw ~ ira + ns(age ,10)+e401 +ns(nifa, 5)+ ns(inc, 5) +hmort + hval + male+ twoearn + smcol + fsize + marr ,data=df_1) 
pr_gam <- predict(gam_df_1)
mspe <- mean((df_1$tw - pr_gam)^2)
print(mspe)
summary(gam_df_1)

#Try to add the interaction terms of marr and male
interaction_1 <- df_1$male*df_1$marr
names(interaction_1) <- paste0('inter_', names(interaction_1))
#add on to the previous model with ira and age interaction
df_2_int_1 <- cbind(df_1, interaction_1)
interaction <- df_1$ira*df_1$age
names(interaction) <- paste0('inter_', names(interaction))
#doesn't add on to the previous spline transformation
df_2_int_2 <- cbind(df_2_int_1, interaction)
gam_df_1 <- lm(tw ~  ira +ns(age ,10)+e401 +ns(nifa, 5)+ ns(inc, 5) + hmort + hval + twoearn + smcol + fsize +interaction_1 +interaction ,data=df_2_int_2) 
pr_gam <- predict(gam_df_1)
mspe <- mean((df_1$tw - pr_gam)^2)
print(mspe)
summary(gam_df_1)



#Using crossvalidation to test for the transformation model
n <- nrow(df_1)
k<-5
set.seed(2345)
ii <- sample(rep(1:k,floor(n/k)))
y <- df_1$tw
mspe1 <- mspe2 <- mspe3 <- mspe4 <- vector()
for (j in 1:k){
  hold <- (ii == j)
  train <- (ii != j)
  
  y_te <- y[hold]
  x_te <- df_1[hold,-1]
  y_tr <- y[train]
  x_tr <- df_1[train,-1]
  
  #the transformation only include using splines of three variable: age, hval, and hmort
  full.model <-lm(tw ~., data = df_1_trans_2[train,]) 
  null.model <-lm(tw ~1, data = df_1_trans_2[train,])
  for_model_trans <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
  #the transformation include splines and interactions: age*ira, male*marr 
  full.model <-lm(tw ~., data = df_1_int_1[train,]) 
  null.model <-lm(tw ~1, data = df_1_int_1[train,])
  for_model_int_1 <- stepAIC(null.model,scope = list(lower = null.model,upper = full.model), trace = FALSE, direction = 'forward')
  #the polynomial model
  full.model_2 <-lm(tw ~., data = df_poly_2[train,]) 
  null.model_2 <-lm(tw ~1, data = df_poly_2[train,])
  for_model_poly_2 <- stepAIC(null.model_2,scope = list(lower = null.model_2,upper = full.model_2), trace = FALSE, direction = 'forward')
  #GAM
  gam_df_1 <- lm(tw ~  ira +ns(age ,10)+e401 +ns(nifa, 5)+ ns(inc, 5) + hmort + hval + twoearn + smcol + fsize +interaction_1 +interaction ,data=df_2_int_2[train,]) 
  
  pr_spl <- predict(for_model_trans, newdata = df_1_trans_2[hold,])
  pr_spl_int <- predict(for_model_int_1, newdata = df_1_int_1[hold,])
  pr_poly <- predict(for_model_poly_2, newdata = df_poly_2[hold,])
  pr_gam <- predict(gam_df_1, newdata = df_2_int_2[hold,])
  
  mspe1[j] <- mean((y_te - pr_spl)^2)
  mspe2[j] <- mean((y_te - pr_spl_int)^2)
  mspe3[j] <- mean((y_te - pr_poly)^2)
  mspe4[j] <- mean((y_te - pr_gam)^2)
  
}
## Return the mspe
mean(mspe1)
mean(mspe2)
mean(mspe3)
mean(mspe4)

#The fourth model with gam and interaction give the smallest mspe4. Therefore, we decide to use this model for prediction.
gam_df_1 <- lm(tw ~  ira +ns(age ,10)+e401 +ns(nifa, 5)+ ns(inc, 5) + hmort + hval + twoearn + smcol + fsize +interaction_1 +interaction ,data=df_2_int_2)
summary(gam_df_1)


#Prediction
data_te <- read.table("data_for_prediction.txt", header = TRUE, sep = "\t", dec = ".")[,-1]
#Data transformation of test data
interaction_1 <- data_te$male*data_te$marr
names(interaction_1) <- paste0('inter_', names(interaction_1))
df_2_int_1 <- cbind(data_te, interaction_1)

interaction <- data_te$ira*data_te$age
names(interaction) <- paste0('inter_', names(interaction))
df_2_int_2 <- cbind(df_2_int_1, interaction)

my_predictions <- predict(gam_df_1, newdata = df_2_int_2)
length(my_predictions)
write.table(my_predictions, file = 'my_predictions.txt')





                                                                                                                                  