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
  s2 <- (t(e) %*% e) / (n - p) 
  
  cov <- (as.numeric(s2)) * solve(t(X) %*% X)
  
  # SE 
  SE <- sqrt(diag(cov))
  
  # T values
  t_value <- B / SE
  
  # Pr
  dfs_residual <- n - p
  Pr <- 2 * pt(abs(t_value), dfs_residual, lower.tail = FALSE)
  
  # SD(e) and R2
  SD_e <- sqrt(s2)
  R2 <- (t(B) %*% cov(X) %*% B) / ((t(B) %*% cov(X) %*% B) + as.numeric(s2))
  
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
# b0    42.3863641 4.3789952  9.6794726 1.965929e-10
# b1    -3.3920819 0.8208025 -4.1326409 2.941570e-04
# b2    -1.5280010 0.4197533 -3.6402363 1.092609e-03
# b3    -0.5228629 0.7788703 -0.6713094 5.075244e-01
# SD(e)  2.5921848        NA         NA           NA
# R2     0.8182681        NA         NA           NA

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
  
  cor_mat <- matrix(par_list$rho, p, p)
  diag(cor_mat) <- 1
  Sigma <- diag(par_list$sd) %*% cor_mat %*% diag(par_list$sd)
  mu <- matrix(par_list$mu, nrow = p)
  
  return(list(mu = mu, Sigma = Sigma, cor = cor_mat))
  
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
# b0    -4.9370600 0.040463429 -122.0129        0
# b1     0.9891621 0.006771543  146.0763        0
# b2     0.9997359 0.003373688  296.3332        0
# SD(e)  2.0361731          NA        NA       NA
# R2     0.5978875          NA        NA       NA

dif_3b <- rbind(b0 = cbind(b0, ans_3b[1, 1], b0 - ans_3b[1, 1]), 
                b1 = cbind(p_Y$b[1],  ans_3b[2, 1], p_Y$b[1] - ans_3b[2, 1]),
                b2 = cbind(p_Y$b[2], ans_3b[3, 1], p_Y$b[2] - ans_3b[3, 1]),
                R2 = cbind(p_Y$R2, ans_3b[5, 1], p_Y$R2 - ans_3b[5, 1]),
                mu = cbind(p_Y$mu, mean(Y), p_Y$mu - mean(Y)),
                sd = cbind(p_Y$sd, sd(Y), p_Y$sd - sd(Y)))

colnames(dif_3b) <- c("Pop", "Est", "Dif")
rownames(dif_3b) <- c("b0", "b1", "b2", "R2", "mu", "sd")

#     Pop        Est           Dif
# b0 -5.0 -4.9370600 -0.0629399882
# b1  1.0  0.9891621  0.0108378625
# b2  1.0  0.9997359  0.0002641061
# R2  0.6  0.5978875  0.0021125328
# mu 10.0 10.0016488 -0.0016487505
# sd  5.0  3.2109914  1.7890086150



# 3.c 
## set seed and param
set.seed(123782)
p_c <- list(r = c(0.3, -0.4), mu = 10, sd = 5)

## get combined Sigma from combined cor. mat 
cor_mat <- with(p_c, rbind(c(1, r), cbind(r, p_Xt$cor)))
Sigma <- diag(c(p_c$sd, p_X$sd)) %*% cor_mat %*% diag(c(p_c$sd, p_X$sd))

## reg coef 
B_c <- solve(p_Xt$Sigma) %*% Sigma[-1, 1]

## var explained 
v_exp_c <- with(p_c, sd^2 - t(B_c) %*% Sigma[-1, 1])
s_exp_c <- sqrt(v_exp_c)

## R2 
R2 <- cor_mat[1, -1] %*% solve(p_Xt$cor) %*% cor_mat[-1, 1]

## intercept 
b0_c <- with(p_c, mu - (t(p_Xt$mu) %*% B_c))
  
## gen Y 
Y_c <- with(p_c, (matrix(1, n) %*% b0_c) + X %*% B_c + rnorm(n, 0, s_exp_c))

ans_3c <- ols_reg(Y_c, cbind(1, X))

#         Estimate          SE   t value Pr(>|t|)
# b0    11.7786196 0.079813607  147.5766        0
# b1     2.3182149 0.013356784  173.5609        0
# b2    -1.3368715 0.006654558 -200.8956        0
# SD(e)  4.0163259          NA        NA       NA
# R2     0.3524265          NA        NA       NA

dif_3c <- rbind(b0 = cbind(b0_c, ans_3c[1], b0_c - ans_3c[1]), 
                b1 = cbind(B_c[1], ans_3c[2], B_c[1] - ans_3c[2]),
                b2 = cbind(B_c[2], ans_3c[3], B_c[2] - ans_3c[3]),
                sde = cbind(s_exp_c, ans_3c[4], s_exp_c - ans_3c[4]),
                R2 = cbind(R2, ans_3c[5], R2 - ans_3c[5])
)

colnames(dif_3c) <- c("Pop", "Est", "Dif")
rownames(dif_3c) <- c("b0", "b1", "b2", "SD(e)", "R2")

#             Pop        Est          Dif
# b0    11.9230769 11.7786196  0.144457282
# b1     2.3076923  2.3182149 -0.010522596
# b2    -1.3461538 -1.3368715 -0.009282314
# SD(e)  4.0191848  4.0163259  0.002858840
# R2     0.3538462  0.3524265  0.001419650


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

## get pop SD(e)
p <- length(p_x$mu)
p_xt <- trans_p(p_x)
X <- with(p_xt, rmvnorm(n, mu, Sigma))
B <- with(p_y, matrix(b, nrow = p))
V_exp <- with(p_y, (t(B) %*% p_xt$Sigma %*% B) * ((1/R2) - 1))
s_exp <- sqrt(V_exp)

## gen data using method 1 
data_m1 <- gen_m1(100000, p_x, p_y)

ans_3d1 <- ols_reg(data_m1[,1], cbind(1, data_m1[,2:6]))

#         Estimate          SE    t value Pr(>|t|)
# b0    25.0219373 0.015237786 1642.09794        0
# b1     0.9812441 0.015667438   62.62952        0
# b2     1.0199428 0.011092474   91.94908        0
# b3     0.9881404 0.009068983  108.95824        0
# b4     0.9852023 0.007850763  125.49128        0
# b5     1.0050628 0.007062545  142.30886        0
# SD(e)  4.8184817          NA         NA       NA
# R2     0.4986015          NA         NA       NA

dif_3d1 <- rbind(
  b1 = cbind(p_y$b[1], ans_3d1[2, 1], p_y$b[1] - ans_3d1[2, 1]),
  b2 = cbind(p_y$b[2], ans_3d1[3, 1], p_y$b[2] - ans_3d1[3, 1]),
  b3 = cbind(p_y$b[3], ans_3d1[4, 1], p_y$b[3] - ans_3d1[4, 1]),
  b4 = cbind(p_y$b[4], ans_3d1[5, 1], p_y$b[4] - ans_3d1[5, 1]),
  b5 = cbind(p_y$b[5], ans_3d1[6, 1], p_y$b[5] - ans_3d1[6, 1]),
  sde = cbind(s_exp, ans_3d1[7, 1], s_exp - ans_3d1[7, 1]),
  R2 = cbind(p_y$R2, ans_3d1[8, 1], p_y$R2 - ans_3d1[8, 1]),
  mu = cbind(p_y$mu, mean(data_m1[,1]), p_y$mu - mean(data_m1[,1]))
)

colnames(dif_3d1) <- c("Pop", "Est", "Dif")
rownames(dif_3d1) <- c("b1", "b2", "b3", "b4", "b5", "SD(e)", "R2", "mu")
 
#             Pop        Est          Dif
# b1     1.000000  0.9812441  0.018755865
# b2     1.000000  1.0199428 -0.019942757
# b3     1.000000  0.9881404  0.011859598
# b4     1.000000  0.9852023  0.014797719
# b5     1.000000  1.0050628 -0.005062766
# SD(e)  4.825922  4.8184817  0.007440416
# R2     0.500000  0.4986015  0.001398544
# mu    25.000000 25.0334928 -0.033492815

## Method 2 
gen_m2 <- function(n, p_x, p_y){
  
  p <- length(p_x$mu)
  
  p_xt <- trans_p(p_x)
  X <- with(p_xt, rmvnorm(n, mu, Sigma))
  
  ## get combined Sigma from combined cor_mat 
  
  R <- with(p_y, rbind(c(1, rho), cbind(rho, p_xt$cor)))
  Sigma <- diag(c(p_y$sd, p_x$sd)) %*% R %*% diag(c(p_y$sd, p_x$sd))
  
  ## reg coef 
  B <- solve(p_xt$Sigma) %*% Sigma[-1, 1]
  
  ## var explained 
  v_exp <- with(p_y, sd^2 - t(B) %*% Sigma[-1, 1])
  s_exp <- sqrt(v_exp)
  
  ## R2 
  R2 <- R[1, -1] %*% solve(p_xt$cor) %*% R[-1, 1]
  
  ## intercept 
  b0 <- with(p_y, mu - (t(p_xt$mu) %*% B))
  
  ## gen Y 
  Y <- with(p_y, (matrix(1, n) %*% b0) + X %*% B + rnorm(n, 0, s_exp))

  output <- cbind(Y,X)
  attr(output, 'par_xy') <- list(p_x, p_y)
  return(list(data = output,
              b0 = b0, B = B, s_exp = s_exp, R2 = R2)
  )

}


## set seed and param
set.seed(1237)
p_x <- list(mu = rep(0, 5),
            sd = c(1, sqrt(2), sqrt(3), 2, sqrt(5)), 
            rho = 0.15)

p_y <- list(rho = c(-0.15, -0.5, 0.15, 0.3, 0.20), mu = 10, sd = 4)

## gen data using method 2 
m2 <- gen_m2(10000, p_x, p_y)
data_m2 <- m2[[1]]
ans_3d2 <- ols_reg(data_m2[,1], cbind(1, data_m2[,2:6]))

#         Estimate         SE   t value      Pr(>|t|)
# b0     9.9958695 0.02803312 356.57357  0.000000e+00
# b1    -0.7106237 0.02885137 -24.63050 4.277579e-130
# b2    -1.6813177 0.02055348 -81.80209  0.000000e+00
# b3     0.3964946 0.01673604  23.69107 9.126511e-121
# b4     0.7044516 0.01438804  48.96093  0.000000e+00
# b5     0.4194456 0.01284228  32.66131 2.057206e-222(
# SD(e)  2.8030957         NA        NA            NA
# R2     0.5030344         NA        NA            NA


dif_3d2 <- rbind(b0 = cbind(m2$b0, ans_3d2[1], m2$b0 - ans_3d2[1]),
                 b1 = cbind(m2$B[1], ans_3d2[2], m2$B[1] - ans_3d2[2]),
                 b2 = cbind(m2$B[2], ans_3d2[3], m2$B[2] - ans_3d2[3]),
                 b3 = cbind(m2$B[3], ans_3d2[4], m2$B[3] - ans_3d2[4]),
                 b4 = cbind(m2$B[4], ans_3d2[5], m2$B[4] - ans_3d2[5]),
                 b5 = cbind(m2$B[5], ans_3d2[6], m2$B[5] - ans_3d2[6]),
                 sde = cbind(m2$s_exp, ans_3d2[7], m2$s_exp - ans_3d2[7]),
                 R2 = cbind(m2$R2, ans_3d2[8], m2$R2 - ans_3d2[8])
)


colnames(dif_3d2) <- c("Pop", "Est", "Dif")
rownames(dif_3d2) <- c("b0", "b1", "b2", "b3", "b4", "b5", "SD(e)", "R2")


#             Pop        Est          Dif
# b0    10.0000000  9.9958695  0.004130497
# b1    -0.7058824 -0.7106237  0.004741309
# b2    -1.6637807 -1.6813177  0.017536996
# b3     0.4075414  0.3964946  0.011046756
# b4     0.7058824  0.7044516  0.001430782
# b5     0.4209069  0.4194456  0.001461357
# SD(e)  2.8284271  2.8030957  0.025331412
# R2     0.5000000  0.5030344 -0.003034398