library('MASS')
library(ggplot2)
library(mclust, quietly=TRUE)
velocities <- galaxies/1000
normalize = function(x){return(x/sum(x))}

z_sample = function(x,pi,mu){
  dmat = outer(mu,x,"-") 
  p.z.given.x = as.vector(pi) * dnorm(dmat,0,1) 
  p.z.given.x = apply(p.z.given.x,2,normalize)
  z = rep(0, length(x))
  for(i in 1:length(z)){
    z[i] = sample(1:length(pi), size=1,prob=p.z.given.x[,i],replace=TRUE)
  }
  return(z)
}


pi_sample = function(z,k){
  counts = colSums(outer(z,1:k,FUN="=="))
  pi = gtools::rdirichlet(1,counts+1)
  return(pi)
}

mu_sample = function(x, z, k, prior){
  df = data.frame(x=x,z=z)
  mu = rep(0,k)
  for(i in 1:k){
    sample.size = sum(z==i)
    sample.mean = ifelse(sample.size==0,0,mean(x[z==i]))
    
    post.prec = sample.size+prior$prec
    post.mean = (prior$mean * prior$prec + sample.mean * sample.size)/post.prec
    mu[i] = rnorm(1,post.mean,sqrt(1/post.prec))
  }
  return(mu)
}

gibbs_sampling <- function(x,k,niter =1000,muprior = list(mean=0,prec=0.1)){
  pi = rep(1/k,k) # initialize
  mu = rnorm(k,0,10)
  z = z_sample(x,pi,mu)
  res = list(mu=matrix(nrow=niter, ncol=k), pi = matrix(nrow=niter,ncol=k), z = matrix(nrow=niter, ncol=length(x)))
  res$mu[1,]=mu
  res$pi[1,]=pi
  res$z[1,]=z 
  for(i in 2:niter){
    pi = pi_sample(z,k)
    mu = mu_sample(x,z,k,muprior)
    z = z_sample(x,pi,mu)
    res$mu[i,] = mu
    res$pi[i,] = pi
    res$z[i,] = z
  }
  return(res)
}
k = 8
fit = gibbs_sampling (velocities, k)
fit = Mclust(velocities, G=k, model="V")
mu = as.numeric(fit$parameters$mean) 
pi = fit$parameters$pro
sigma = fit$parameters$variance$sigmasq
sigma = sigma^0.5

x_sample <- seq(7,40,length=1000)
pdf_calculator = function(x_value, mus, pis, sigma, k){
  sum <- 0
  for (k in 1:k){
    sum = sum + pis[k]*dnorm(x_sample, mus[k], sigma[k])
  }
  return(sum)
}
y_sample <- pdf_calculator(x_sample, mu, pi, sigma, k)

mean_sample = mean(y_sample)
variance_sample = var(y_sample)
nr1 = qt(0.95, df=1000-1)
mean_margin <- nr1*variance_sample/sqrt(1000)
sample_lower_bound = y_sample - mean_margin
sample_upper_bound = y_sample + mean_margin
  
plot(x_sample, y_sample, type='l', xlab = "Velocities", ylab='Densities', col='brown')
lines(x_sample, sample_upper_bound, type = 'l', lty=3, col='orange')
lines(x_sample, sample_lower_bound, type = 'l', lty=3, col='orange')
legend(x = "topright", legend=c("Actual Density", "90% CI Lower Bound", "90% CI Upper Bound", "Histogram of velocities"), 
       fill = c("brown","orange", "orange", "cyan"), cex = 0.6
)
hist(velocities, prob = TRUE, col = 'cyan', breaks = 50, add=TRUE)
title(main = paste(c("Density Plot for K = "),k),
      cex.main = 1,   font.main= 1, col.main= "orange",
)
