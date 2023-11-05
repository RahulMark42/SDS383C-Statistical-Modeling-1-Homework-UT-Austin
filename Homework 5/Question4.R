library(invgamma)
simulate_ar1 = function(n, rho){
  yt = rep(0,n)
  for(t in 2:n){
    yt[t] = rnorm(1, rho*yt[t-1],1)
  }
  return(yt)
}
N=500
yt = simulate_ar1(N,0.75)

burn_in <- 2000
thin <- 5
n_iter <- 10000

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
burned_in_samples <- samples[burn_in:n_iter,]
thinned_samples <- burned_in_samples[c(FALSE, FALSE, FALSE, TRUE),]

thinned_samples <- data.frame(thinned_samples)
colnames(thinned_samples) = c("rho","sigma2")
par(mfrow = c(1,2))

plot(thinned_samples$rho, xlab = "t value",ylab= "rho sampled value", main = "Traceplot of rho", type='l', col = 'green')
plot(thinned_samples$sigma2, xlab = "t value", ylab = "sigma2 sampled value", main = "Traceplot of sigma^2", type='l', col = 'blue')

sampled_means = c(thinned_samples[1])
sampled_vars = c(thinned_samples[2])

m_means = lapply(sampled_means, mean)$rho
s_means = lapply(sampled_means, var)$rho
nr1 = qt(0.95, df=2000-1)
mean_margin <- nr1*s_means/sqrt(2000)
mean_lower_bound = m_means - mean_margin
mean_upper_bound = m_means + mean_margin

m_vars = lapply(sampled_vars, mean)$sigma2
s_vars = lapply(sampled_vars, var)$sigma2
nr2 = qt(0.95, df=2000-1)
variance_margin <- nr2*s_vars/sqrt(2000)
variance_lower_bound = m_vars - variance_margin
variance_upper_bound = m_vars + variance_margin

variance_lower_bound
variance_upper_bound