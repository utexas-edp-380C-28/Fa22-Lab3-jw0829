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

p_x <- list(mu = -2, sd = 3)
p_y <- list(b0 = 12, b1 = 4, sd = 7)


gen_xy <- function(n, p_x, p_y){

  x <- with(p_x, rnorm(n, mu, sd))
  y <- with(p_y, rnorm(n, (b0 + b1*x), sd)) 
  
  xy <- cbind(x = x, y = y)
  
  return(as.data.frame(xy))
}


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

set.seed(17290)

sim_data <- replicate(500, gen_xy(100, p_x, p_y), simplify = "array")

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

ols_reg <- function(y, X){
  
  p <- ncol(X)
  n <- nrow(X)
  
  result <- matrix(NA, nrow = p + 2, ncol = 4)
  colnames(result) <- c("Estimate", "SE", "t value", "Pr(>|t|)")
  rownames(result) <- c(paste0('b', rep(1:p)), "SD(e)", "R2")
  
  # Estimates 
  B <- solve(t(X) %*% X) %*% t(X) %*% y
  # solve(crossprod(X))
  
  # Hat Matrix and Residual 
  H <- X %*% solve(t(X) %*% X)%*% t(X)
  e <- (diag(1, n) - H) %*% y
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


y <- as.matrix(mtcars$mpg)
X <- as.matrix(cbind(rep(1, nrow(mtcars)), 
                     mtcars[c("wt", "cyl", "gear")]
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

set.seed(21389)

n <- 100000

p_X <- list(mu = c(5, 10), sd = c(1, 2), rho = 0.3)

trans_p <- function(par_list){
  
  p <- length(par_list$mu)
  
  cor.mat <- matrix(par_list$rho, p, p)
  diag(cor.mat) <- 1
  Sigma <- diag(par_list$sd) %*% cor.mat %*% diag(par_list$sd)
  mu <- matrix(par_list$mu, nrow = p)
  
  return(list(mu = mu, Sigma = Sigma))
  
}

p_Xt <- trans_p(p_X)

X <- with(p_Xt, rmvnorm(n, mu, Sigma))


# 3.b

set.seed(23921)

p_Y <- list(b1 = 1, b2 = 1, R2 = 0.6, mu = 10, sd = 5)

B <- with(p_Y, matrix(c(b1, b2), nrow = 2))

V_e <- with(p_Y, (t(B) %*% p_Xt$Sigma %*% B) * ((1/R2) - 1))
s_e <- sqrt(V_e)

Y <- with(p_Y, X %*% B + rnorm(n, 0, (s_e)))

ols_reg(Y, X)

# Error: cannot allocate vector of size 74.5 Gb




