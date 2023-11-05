library(invgamma)

ar1_simulation = function(N, rho){
  yt = rep(0,N)
  for(t in 2:N){
    yt[t] = rnorm(1, rho*yt[t-1],1)
  }
  return(yt)
}

N=500
yt = ar1_simulation(N,0.75)

burn_in <- 2000
thin <- 2
n_iter <- 6000

s2 <- 1
rho0 <- 0.5
rho <- rho0
lambda <- 2
a <- 2
b <- 4

ar1_gibbs_sampler = function(yt, lambda, rho0, rho, s2){
  samples <- matrix(0, n_iter, 2)
  for (i in 1:n_iter)
  {
    y_sum_squared = sum(yt^2)
    A = y_sum_squared - yt[N]^2 + lambda
    B = sum(yt[1:(N-1)] * yt[2:N])
    C = y_sum_squared - yt[1]^2 + 2*b + lambda * rho0^2
    rho = rnorm(1,mean = B/A, sd = sqrt(s2/A))
    s2 = rinvgamma(1, shape = a + N/2, scale = (C - B^2/A + A*(rho-B/A)^2)/2)
    samples[i,] <- c(rho, s2)
  }
  return(samples)
}
  
samples <- ar1_gibbs_sampler(yt, lambda, rho0, rho, s2)
burned_in <- samples[burn_in:n_iter,]
thinned <- burned_in[c(FALSE, TRUE),]

thinned <- data.frame(thinned)
colnames(thinned) = c("rho","sigma2")
par(mfrow = c(1,2))

plot(thinned$rho, xlab = "iterations",ylab= "sampled value", main = "Traceplot of rho", type='l', col = 'cyan')
plot(thinned$sigma2, xlab = "iterations", ylab = "sampled value", main = "Traceplot of sigma^2", type='l', col = 'red')

sampled_rho = c(thinned[1])
sampled_vars = c(thinned[2])

mean_rho = lapply(sampled_rho, mean)$rho
variance_rho = lapply(sampled_rho, var)$rho
nr1 = qt(0.95, df=2000-1)
mean_margin <- nr1*variance_rho/sqrt(2000)
rho_lower_bound = mean_rho - mean_margin
rho_upper_bound = mean_rho + mean_margin

mean_sigma = lapply(sampled_vars, mean)$sigma2
variance_sigma = lapply(sampled_vars, var)$sigma2
nr2 = qt(0.95, df=2000-1)
var_margin <- nr2*variance_sigma/sqrt(2000)
vars_lower_bound = mean_sigma - var_margin
vars_upper_bound = mean_sigma + var_margin