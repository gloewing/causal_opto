library(data.table)
library(geepack)
library(rootSolve)
library(optimx)

dat <- data.table::fread("/Users/loewingergc/Desktop/NIMH Research/Causal/spontaneous_behavior_opto/example_data.csv")
source("/Users/loewingergc/Desktop/NIMH Research/Causal/msm_hr/Final Sims/Code/sandwich_normal.R")
n <- length(unique(dat$id))

set.seed(634)
data <- dat[sample(1:nrow(dat), 10000, replace = F)]
data$id <- as.numeric(as.factor(data$id))
data <- data[order(data$id),]

mod2 <- glm(YY ~ d_sum * preCount,
            family = "binomial",
            weights = ww,
            data = data)
beta_hat <- coef(mod2)

data$Y.cont <- rnorm(nrow(data), mean = predict(mod2))

mod2 <- lm(Y.cont ~ d_sum * preCount,
           weights = ww,
           data = data,x=TRUE, y = TRUE)
beta_hat <- coef(mod2)

# mod2.star <- geeglm(YY ~ d_sum * preCount,
#                     family = binomial,
#                     weights = ww,
#                     id = id,
#                     corstr = "independence",
#                     data = data)
# beta_hat.star <- coef(mod2.star)

mod3.star <- geese(Y.cont ~ d_sum * preCount,
                   mean.link = "identity",
                   variance = "gaussian",
                   # scale.fix = T,
                   # scale.value = 1,
                   weights = ww,
                   id = id,
                   corstr = "independence",
                   data = data)
beta_hat3 <- mod3.star$beta

# x_mat <- cbind(1, data$d_sum, data$preCount, data$d_sum * data$preCount)
x_mat <- model.matrix(mod2)
ww <- data$ww
YY <- data$Y.cont

## fitted value function
m <- function(x, beta) {
  (x %*% beta)[,1]
}

### point estimation
# approach 1: minimizing weighted distance

## takes in p-vector beta
## outputs value of objective function at that beta
obj_function <- function(beta) {
  sum(ww*(YY - m(x_mat, beta))^2)
}

## takes in p-vector beta
## outputs value of estimating function at that beta
est_function <- function(beta) {
  fit.vals <- m(x_mat, beta)
  colSums(x_mat * (ww * (YY - fit.vals)))
}

obj_gr <- function(beta) {
  - 2 * est_function(beta)
}

beta_hat.HC <- optim(beta_hat, fn = obj_function,
                     method = "L-BFGS-B", gr = obj_gr)$par

# beta_hat.HC <- optimr(rep(0, ncol(x_mat)), fn = obj_function, gr = obj_gr,
#                       method = "nlnm")$par

# approach 2: solving estimating equation directly

beta_hat.HC <- multiroot(est_function, rep(0, ncol(x_mat)), atol = 1e-9)$root


### standard error calculations
phi <- function(id, beta) {
  x_sub <- x_mat[data$id == id,]
  ww_sub <- ww[data$id == id]
  YY_sub <- YY[data$id == id]
  
  m_res <- m(x_sub, beta)
  return ( colSums(x_sub * ww_sub * (YY_sub - m_res)) )
}

phi_grad <- function(id, beta) {
  x_sub <- x_mat[data$id == id,]
  ww_sub <- ww[data$id == id]
  YY_sub <- YY[data$id == id]
  
  m_res <- m(x_sub, beta)
  out_prods <- lapply(1:nrow(x_sub), function(i) {
    x_sub[i,] %*% t(x_sub[i,])
  })
  scls <- - ww_sub
  Reduce('+', lapply(1:nrow(x_sub), function(i) {
    out_prods[[i]] * scls[i]
  }))
}

## estimating functions themselves
phi_mat <- t(sapply(unique(data$id), FUN = function(id) {
  phi(id, beta_hat3)
},simplify = 0))

## Jacobians of estimating functions
phi_grad_list <- lapply(unique(data$id), FUN = function(id) {
  phi_grad(id, beta_hat3)
})

## colSums(phi_mat) ## should be approximately zeroes!

A_hat <- ( t(phi_mat) %*% phi_mat ) / n
B_hat <- Reduce("+", phi_grad_list) / n
  
B_inv <- solve(B_hat)
V_hat <- B_inv %*% A_hat %*% B_inv

sqrt(diag(V_hat) / n) ## hand-coded standard errors
sqrt(diag(mod3.star$vbeta)) ## geese standard errors

sandwich_normal(data = (data), model = mod2, # requires either GLM model (warm-start) or X, Y, w
                beta_hat = as.numeric(beta_hat3), # beta estimate
                X = NULL,
                Y = NULL,
                w = NULL,
                id_name = "id")$SEs
                            
# cbind(beta_hat - qnorm(0.975) * sqrt(diag(V_hat) / n),
#       beta_hat + qnorm(0.975) * sqrt(diag(V_hat) / n)) ## hand-coded CIs
# 
# hc0 <- sandwich::vcovCL(mod2, cluster = ~id, type = "HC")
# sqrt(diag(hc0))
# hc1 <- sandwich::vcovCL(mod2, type = "HC1", cluster = ~id)
# sqrt(diag(hc1))
# hc2 <- sandwich::vcovCL(mod2, type = "HC2", cluster = ~id)
# sqrt(diag(hc2))
