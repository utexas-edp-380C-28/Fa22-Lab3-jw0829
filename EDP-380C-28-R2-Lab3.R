#---------------------------------------------------#
# EDP 380C.26: Simulation in R
# Lab 3: An Overview of Linear Regression 
#        and Data Generation
#
# Name:  JIWON KIM
#---------------------------------------------------#

## SETUP -------------------------------------------
# Source rmvnorm()
source("scripts/rmvnorm.R")

## Begin Code --------------------------------------

# 1.a 
## (1) set param
p_x <- list(mu = -2, sd = 3)
p_y <- list(b0 = 12, b1 = 4, sd = 7)

## (2) generate x and y  
gen_xy <- function(n, p_x, p_y){

  x <- with(p_x, rnorm(n, mu, sd))
  y <- with(p_y, rnorm(n, (b0 + b1*x), sd)) 
  
  xy <- cbind(x = x, y = y)
  
  return(as.data.frame(xy))
}

## (3) analyze() 
analyze <- function(dataset){
  
  result <- summary(lm(y ~ x, data = dataset))

  return(c(m_x = mean(dataset$x) ,
           m_y = mean(dataset$y),
           s_x = sd(dataset$x),
           s_y = sd(dataset$y),
           cor = cor(dataset$x, dataset$y),
           b0 = result$coefficients[1],
           b1 = result$coefficients[2],
           s_e = result$sigma))
}

## (4) set seed
set.seed(17290)

## (5) run simul using replicate()
sim_data <- replicate(500, gen_xy(100, p_x, p_y), simplify = "array")

## (6) use analyze() 
sim_result <- apply(sim_data, MARGIN = 2, analyze)
out.1a <- list(mean = rowMeans(sim_result),
               sd = apply(sim_result, MARGIN = 1, sd))

#out.1a
# $mean
#       m_x        m_y        s_x        s_y        cor         b0         b1        s_e 
# -2.0061317  3.9810216  3.0034995 13.8596003  0.8631812 11.9780298  3.9873856  6.9749258 
# 
# $sd
#       m_x        m_y        s_x        s_y        cor         b0         b1        s_e 
# 0.31115637 1.41924855 0.20911389 1.03767609 0.02432744 0.86119359 0.24140973 0.49070350 



# 2.a 
## ols_reg()
ols_reg <- function(y, X){
  
  p <- ncol(X)
  n <- nrow(X)
  
  result <- matrix(NA, nrow = p + 2, ncol = 4)
  colnames(result) <- c("Estimate", "SE", "t value", "Pr(>|t|)")
  rownames(result) <- c(paste0('b', rep(0:(p-1))), "SD(e)", "R2")
  
  # Estimates 
  B <- solve(t(X) %*% X) %*% t(X) %*% y
  # solve(crossprod(X))
  
  # Hat Matrix and Residual 
  # H <- X %*% solve(t(X) %*% X)%*% t(X)
  # e <- (diag(1, n) - H) %*% y
  
  yhat <- X %*% B 
  e <- y - yhat 
  s2 <- (t(e) %*% e) / n
  
  cov <- (as.integer(s2)) * solve(t(X) %*% X)
  
  # SE 
  SE <- sqrt(diag(cov))
  
  # T values
  t_value <- B / SE
  
  # Pr
  dfs_residual <- n - p
  Pr <- 2 * pt(abs(t_value), dfs_residual, lower.tail = FALSE)
  
  # SD(e) and R2
  SD_e <- sqrt(s2)
  R2 <- (t(B) %*% cov(X) %*% B) / ((t(B) %*% cov(X) %*% B) + as.integer(s2))
  
  result[, 1] <- c(B, SD_e, R2)
  result[, 2] <- c(SE, NA, NA)
  result[, 3] <- c(t_value, NA, NA)
  result[, 4] <- c(Pr, NA, NA)
  
  return(result)
  
}

## using mtcar, fit reg model 
y <- as.matrix(mtcars$mpg)
X <- as.matrix(cbind(1, mtcars[c("wt", "cyl", "gear")]
  )
)

ols_reg(y, X)

#         Estimate        SE    t value     Pr(>|t|)
# b1    42.3863641 3.7774047 11.2210282 7.134111e-12
# b2    -3.3920819 0.7080399 -4.7908064 4.913555e-05
# b3    -1.5280010 0.3620872 -4.2199814 2.323616e-04
# b4    -0.5228629 0.6718684 -0.7782223 4.429649e-01
# SD(e)  2.4247668        NA         NA           NA
# R2     0.8581759        NA         NA           NA

summary(lm(mpg ~ wt + cyl + gear, data = mtcars))

# Residuals:
#     Min      1Q  Median      3Q     Max 
# -4.8443 -1.5455 -0.3932  1.4220  5.9416 
# 
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    42.3864     4.3790   9.679 1.97e-10 ***
#   wt           -3.3921     0.8208  -4.133 0.000294 ***
#   cyl          -1.5280     0.4198  -3.640 0.001093 ** 
#   gear         -0.5229     0.7789  -0.671 0.507524    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 2.592 on 28 degrees of freedom
# Multiple R-squared:  0.8329,	Adjusted R-squared:  0.815 
# F-statistic: 46.53 on 3 and 28 DF,  p-value: 5.262e-11


# 3.a
## (2) set seed and par
set.seed(21389)
n <- 100000
p_X <- list(mu = c(5, 10), sd = c(1, 2), rho = 0.3)

## transform par_list 
trans_p <- function(par_list){
  
  p <- length(par_list$mu)
  
  cor.mat <- matrix(par_list$rho, p, p)
  diag(cor.mat) <- 1
  Sigma <- diag(par_list$sd) %*% cor.mat %*% diag(par_list$sd)
  mu <- matrix(par_list$mu, nrow = p)
  
  return(list(mu = mu, Sigma = Sigma, cor = cor.mat))
  
}

p_Xt <- trans_p(p_X)

## (3) gen X 
X <- with(p_Xt, rmvnorm(n, mu, Sigma))



# 3.b
## set seed and param
set.seed(23921)
p_Y <- list(b = c(1, 1), R2 = 0.6, mu = 10, sd = 5)

## reg. coef
B <- with(p_Y, matrix(b, nrow = 2))

## var explained
V_exp <- with(p_Y, t(B) %*% p_Xt$Sigma %*% B * ((1/R2) - 1))
s_exp <- sqrt(V_exp)

## intercept 
b0 <- with(p_Y, mu - (t(p_Xt$mu) %*% B))

## gen Y 
Y <- with(p_Y, ((matrix(1, n) %*% b0) + (X %*% B)) + rnorm(n, 0, (s_exp)))

## run ols_reg()
ans_3b <- ols_reg(Y, cbind(1, X))

#         Estimate          SE   t value Pr(>|t|)
# b0    -4.9370600 0.039744587 -124.2197        0
# b1     0.9891621 0.006651245  148.7183        0
# b2     0.9997359 0.003313754  301.6929        0
# SD(e)  2.0361425          NA        NA       NA
# R2     0.6064754          NA        NA       NA

dif_3b <- list(b1 = ans_3b[2, 1] - p_Y$b[1],
               b2 = ans_3b[3, 1] - p_Y$b[2],
               R2 = ans_3b[5, 1] - p_Y$R2,
               mu = mean(Y) - p_Y$mu,
               sd = sd(Y) - p_Y$sd)
# $b1
# [1] -0.01083786
# 
# $b2
# [1] -0.0002641061
# 
# $R2
# [1] 0.006475361
# 
# $mu
# [1] 0.00164875
# 
# $sd
# [1] -1.789009



# 3.c 
## set seed and param
set.seed(123782)
p_c <- list(r = c(0.3, -0.4), mu = 10, sd = 5)

## get combined Sigma from combined cor. mat 
cor.mat <- with(p_c, rbind(c(1, r), cbind(r, p_Xt$cor)))
Sigma <- diag(c(p_c$sd, p_X$sd)) %*% cor.mat %*% diag(c(p_c$sd, p_X$sd))

## reg coef 
B_c <- solve(p_Xt$Sigma) %*% Sigma[2:3, 1]

## var explained 
v_exp_c <- with(p_c, sd^2 - t(B_c) %*% Sigma[2:3, 1])
s_exp_c <- sqrt(v_exp_c)

## R2 
R2 <- cor.mat[1, 2:3] %*% solve(p_Xt$cor) %*% cor.mat[2:3, 1]

## intercept 
b0_c <- with(p_c, mu - (t(p_Xt$mu) %*% B_c))
  
## gen Y 
Y_c <- with(p_c, (matrix(1, n) %*% b0_c) + X %*% B_c + rnorm(n, 0, s_exp_c))

ans_3c <- ols_reg(Y_c, cbind(1, X))

#         Estimate          SE   t value Pr(>|t|)
# b0    11.7786196 0.079489174  148.1789        0
# b1     2.3182149 0.013302490  174.2692        0
# b2    -1.3368715 0.006627508 -201.7156        0
# SD(e)  4.0162657          NA        NA       NA
# R2     0.3542879          NA        NA       NA

dif_3c <- list(r1 = cor(Y_c, X[,1]) - p_c$r[1],
               r2 = cor(Y_c, X[,2]) - p_c$r[2], 
               mu = mean(Y_c) - p_c$mu,
               sd = sd(Y) - p_c$sd)

# $r1
# [1,] 0.001778966
# 
# $r2
# [1,] 0.003319256
# 
# $mu
# [1] 0.000225353
# 
# $sd
# [1] -1.789009


# 3.d
## Method 1 
gen_m1 <- function(n, p_x, p_y){
  
  p <- length(p_x$mu)

  p_xt <- trans_p(p_x)
  X <- with(p_xt, rmvnorm(n, mu, Sigma))
  
  # Method 1 
  
  B <- with(p_y, matrix(b, nrow = p))

  V_exp <- with(p_y, (t(B) %*% p_xt$Sigma %*% B) * ((1/R2) - 1))
  s_exp <- sqrt(V_exp)
  
  b0 <- with(p_y, mu - (t(p_xt$mu) %*% B))
  
  Y <- with(p_y, ((matrix(1, n) %*% b0) + (X %*% B)) + rnorm(n, 0, (s_exp)))
 
  output <- cbind(Y,X)
  attr(output, 'par_xy') <- list(p_x, p_y)
  return(output) 
  

}

## set seed and param 
set.seed(6972)
p_x <- list(mu = rep(0, 5),
            sd = c(1, sqrt(2), sqrt(3), 2, sqrt(5)), 
            rho = 0.15)
p_y <- list(b = rep(1, 5), mu = 25, R2 = 0.5)

## gen data using method 1 
data_m1 <- gen_m1(100000, p_x, p_y)

ans_3d1 <- ols_reg(data_m1[,1], cbind(1, data_m1[,2:6]))
#         Estimate          SE    t value Pr(>|t|)
# b0    25.0219373 0.015166158 1649.85338        0
# b1     0.9812441 0.015593790   62.92531        0
# b2     1.0199428 0.011040332   92.38334        0
# b3     0.9881404 0.009026352  109.47284        0
# b4     0.9852023 0.007813859  126.08396        0
# b5     1.0050628 0.007029346  142.98097        0
# SD(e)  4.8183372          NA         NA       NA
# R2     0.5009573          NA         NA       NA

dif_3d1 <- list(b1 = ans_3d1[2, 1] - p_y$b[1],
                b2 = ans_3d1[3, 1] - p_y$b[2],
                b3 = ans_3d1[4, 1] - p_y$b[3],
                b4 = ans_3d1[5, 1] - p_y$b[4],
                b5 = ans_3d1[6, 1] - p_y$b[5],
                R2 = ans_3d1[8, 1] - p_y$R2,
                mu = mean(data_m1[,1]) - p_y$mu
)
# 
# $b1
# [1] -0.01875586
# 
# $b2
# [1] 0.01994276
# 
# $b3
# [1] -0.0118596
# 
# $b4
# [1] -0.01479772
# 
# $b5
# [1] 0.005062766
# 
# $R2
# [1] 0.0009573376
# 
# $mu
# [1] 0.03349281



## Method 2 
gen_m2 <- function(n, p_x, p_y){
  
  p <- length(p_x$mu)
  
  p_xt <- trans_p(p_x)
  X <- with(p_xt, rmvnorm(n, mu, Sigma))
  
  ## get combined Sigma from combined cor. mat 
  R <- with(p_y, rbind(c(1, rho), cbind(rho, p_xt$cor)))
  Sigma <- diag(c(p_y$sd, p_x$sd)) %*% R %*% diag(c(p_y$sd, p_x$sd))
  
  ## reg coef 
  B <- solve(p_xt$Sigma) %*% Sigma[2:(p + 1), 1]
  
  ## var explained 
  v_exp <- with(p_y, sd^2 - t(B) %*% Sigma[2:(p + 1), 1])
  s_exp <- sqrt(v_exp)
  
  ## R2 
  R2 <- R[1, 2:(p + 1)] %*% solve(p_xt$cor) %*% R[2:(p + 1), 1]
  
  ## intercept 
  b0 <- with(p_y, mu - (t(p_xt$mu) %*% B))
  
  ## gen Y 
  Y <- with(p_y, (matrix(1, n) %*% b0) + X %*% B + rnorm(n, 0, s_exp))

  output <- cbind(Y,X)
  attr(output, 'par_xy') <- list(p_x, p_y)
  return(output) 

}


## set seed and param
set.seed(1237)
p_x <- list(mu = rep(0, 5),
            sd = c(1, sqrt(2), sqrt(3), 2, sqrt(5)), 
            rho = 0.15)

p_y <- list(rho = c(-0.15, -0.5, 0.15, 0.3, 0.20), mu = 10, sd = 4)

## gen data using method 2 
data_m2 <- gen_m2(10000, p_x, p_y)
ans_3d2 <- ols_reg(data_m2[,1], cbind(1, data_m2[,2:6]))

#         Estimate         SE   t value      Pr(>|t|)
# b0     9.9958695 0.02645955 377.77921  0.000000e+00
# b1    -0.7106237 0.02723187 -26.09529 2.815207e-145
# b2    -1.6813177 0.01939977 -86.66691  0.000000e+00
# b3     0.3964946 0.01579660  25.09999 7.062746e-135
# b4     0.7044516 0.01358040  51.87266  0.000000e+00
# b5     0.4194456 0.01212141  34.60370 6.397081e-248
# SD(e)  2.8022547         NA        NA            NA
# R2     0.5318758         NA        NA            NA


dif_3d2 <- list(r1 = cor(data_m2[,1], data_m2[, 2]) - p_y$rho[1],
                r2 = cor(data_m2[,1], data_m2[, 3]) - p_y$rho[2],
                r3 = cor(data_m2[,1], data_m2[, 4]) - p_y$rho[3],
                r4 = cor(data_m2[,1], data_m2[, 5]) - p_y$rho[4],
                r5 = cor(data_m2[,1], data_m2[, 6]) - p_y$rho[5],
                mu = mean(data_m2[,1]) - p_y$mu,
                sd = sd(data_m2[,1]) - p_y$sd
)

# $r1
# [1] -0.009500006
# 
# $r2
# [1] -0.004060765
# 
# $r3
# [1] 0.002822168
# 
# $r4
# [1] -0.0114716
# 
# $r5
# [1] -0.005992614
# 
# $mu
# [1] -0.02163837
# 
# $sd
# [1] -0.02423421
