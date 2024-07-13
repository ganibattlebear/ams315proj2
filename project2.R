install.packages('mice')
library(mice)
#Read in data

setwd('/Users/gania/Documents/dataset')
getwd();
Dat <- read.csv('P2_095768.csv', header=TRUE)

#Form multiple regression model and residual plot of Y vs Ei
M_E <- lm(Y ~ E1+E2+E3+E4, data=Dat)
summary(M_E)
M_raw <- lm(Y ~ (E1+E2+E3+E4+G1+G2+G3+G4+G5+G6+G7+G8+G9+G10+G11+G12+G13+G14+G15+G16+G17+G18+G19+G20)^2, data=Dat)
plot(resid(M_raw) ~ fitted(M_raw), main='Residual Plot')

#Analyze regression with box-cox function
powah = 94/99
library(MASS)
b <- boxcox(M_raw)
lambda <- b$x[which.max(b$y)] 
# Transform entire equation with lambda = .9494949 = 94/99 according to Box-Cox
# Reference on extracting lambda: https://r-coder.com/box-cox-transformation-r/#Extracting_the_exact_lambda
M_trans <- lm( I(Y^(powah)) ~ (E1+E2+E3+E4+G1+G2+G3+G4+G5+G6+G7+G8+G9+G10+G11+G12+G13+G14+G15+G16+G17+G18+G19+G20)^2, data=Dat )
summary(M_raw)$adj.r.square
summary(M_trans)$adj.r.square

# New transformed residual plot
plot(resid(M_trans) ~ fitted(M_trans), main='New Residual Plot')

# Stepwise Regression of our model. We set nvmax to 24 when 
# wanting to see the full list of possible interactions.
install.packages("leaps")
library(leaps)
M <- regsubsets( model.matrix(M_trans)[,-1], I((Dat$Y)^powah),
                 nbest = 1 , nvmax=5, 
                 method = 'forward', intercept = TRUE )
temp <- summary(M)
# Now produce the model. We see the 3rd model has the highest jump in adj r^2.
# includes E2 + E1:E3 + E4:G11
install.packages("knitr")
library(knitr)
Var <- colnames(model.matrix(M_trans))
M_select <- apply(temp$which, 1, 
                  function(x) paste0(Var[x], collapse='+'))
kable(data.frame(cbind( model = M_select, adjR2 = temp$adjr2, BIC = temp$bic)),
      caption='Model Summary')

#Account for main contributions in model: E2 + E1:E3 + E4:G11 
# Split environment to get E1 + E2 + E3 + E4:G11
M_main <- lm( I(Y^powah) ~ E1+E2+E3+E4+G1+G2+G3+G4+G5+G6+G7+G8+G9+G10+G11+G12+G13+G14+G15+G16+G17+G18+G19+G20, data=Dat)
temp1 <- summary(M_main)
kable(temp1$coefficients[ abs(temp1$coefficients[,4]) <= 0.001, ], caption='Sig Coefficients')
# 2nd Stage = Checking whether main variables are significant and how so
M_2stage <- lm( I(Y^(powah)) ~ E1+E2+E3+E4:G11+G8:G19+G7:G15, data=Dat)
temp2 <- summary(M_2stage)
kable(temp2$coefficients[ abs(temp2$coefficients[,3]) >= 4, ])
#High adj R^2 of .6516 when using ~.95 BoxCox transformation with E1+E2+E3+E4:G11+G8:G19+G7:G15
#Also trying other transformations: Note use of no transformation (^1) does not 
#change which variables and interactions end up being significant. However 
#using a slightly lower exponent (.93 range) gives a lower p-value for the y-intercept.
exp <- .925
M_experiment <- lm(I(Y^exp) ~ E1+E2+E3+E4:G11+G8:G19+G7:G15, data=Dat)
summary(M_experiment)
#(4) If results given by different transformations are similar (they are) 
#you could try using all of the variables and combine terms 
#selected from all of these models. Highest adj R^2 at .6635 achieved below.
#However the Bonferonni inequality states the sum of alphas creates the total alpha
#so if we multiply each p value by 24 we would be at alpha = 4.0!!! So this
#tells us less is better and this model is useless.
M_exp <- lm(I(Y^powah) ~ E1+E2+E3+E1:E3+E2:E4+E2:G11+E4:G11+G1:G7+G1:G16+G1:G17+G2:G14+G2:G20+G4:G10+G4:G14+G5:G16+G7:G13+G7:G15+G8:G19+G10:G11+G11:G13, data=Dat)
summary(M_exp)
#Final model, sticking to what worked before.
M_att <- lm(I(Y^powah) ~ E1+E2+E3+E4:G11+G8:G19+G7:G15, data=Dat)
summary(M_att)
# Final model. Adj R^2 = .6516, F = 429.6, p-value < 2.2e-16. Significant at
# alpha = 24 * sum of p values = .015. There is a 1.5% chance we made a Type I error.
#
M_final <- lm(I(Y^powah) ~ E1 + E2 + E3 + E4:G11 + G8:G19 + G7:G15, data = Dat)
summary(M_final)
#Self-note put risk ratio AKA width of confidence interval
predict(M_final, Dat, interval="confidence") 
 