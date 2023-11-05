library('MASS')
library(ggplot2)
library(mclust, quietly=TRUE)
velocities <- galaxies/1000
normalize = function(x){return(x/sum(x))}

sample_z = function(x,pi,mu){
  dmat = outer(mu,x,"-") 
  p.z.given.x = as.vector(pi) * dnorm(dmat,0,1) 
  p.z.given.x = apply(p.z.given.x,2,normalize) 
  z = rep(0, length(x))
  for(i in 1:length(z)){
    z[i] = sample(1:length(pi), size=1,prob=p.z.given.x[,i],replace=TRUE)
  }
  return(z)
}


sample_pi = function(z,k){
  counts = colSums(outer(z,1:k,FUN="=="))
  pi = gtools::rdirichlet(1,counts+1)
  return(pi)
}

sample_mu = function(x, z, k, prior){
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

gibbs = function(x,k,niter =1000,muprior = list(mean=0,prec=0.1)){
  pi = rep(1/k,k) 
  mu = rnorm(k,0,10)
  z = sample_z(x,pi,mu)
  res = list(mu=matrix(nrow=niter, ncol=k), pi = matrix(nrow=niter,ncol=k), z = matrix(nrow=niter, ncol=length(x)))
  res$mu[1,]=mu
  res$pi[1,]=pi
  res$z[1,]=z 
  for(i in 2:niter){
    pi = sample_pi(z,k)
    mu = sample_mu(x,z,k,muprior)
    z = sample_z(x,pi,mu)
    res$mu[i,] = mu
    res$pi[i,] = pi
    res$z[i,] = z
  }
  return(res)
}
k = 8
fit = gibbs(velocities, k)
fit = Mclust(velocities, G=k, model="V")
mus = as.numeric(fit$parameters$mean) 
pis = fit$parameters$pro
sigma = fit$parameters$variance$sigmasq
sigma = sigma^0.5

x_sample <- seq(8,35,length=1000)
density_value = function(x_value, mus, pis, sigma, k){
  sum <- 0
  for (k in 1:k){
    sum = sum + pis[k]*dnorm(x_sample, mus[k], sigma[k])
  }
  return(sum)
}
y_sample <- density_value(x_sample, mus, pis, sigma, k)

means_sample = mean(y_sample)
variance_sample = var(y_sample)
nr1 = qt(0.95, df=1000-1)
mean_margin <- nr1*means_sample/sqrt(1000)
sample_lower_bound = y_sample - mean_margin
sample_upper_bound = y_sample + mean_margin

  
plot(x_sample, y_sample, type='l', xlab = "Velocities", ylab='Densities')
lines(x_sample, sample_upper_bound, type = 'l', lty=3, col='red')
lines(x_sample, sample_lower_bound, type = 'l', lty=3, col='red')
legend(x = "topright", legend=c("Fitted Density", "90% CI Lower Bound", "90% CI Upper Bound", "Histogram of Actual Values"), 
       fill = c("black","red", "red", "green"), cex = 0.6
)
hist(velocities, prob = TRUE, col = 'green', breaks = 25, add=TRUE)
title(main = paste(c("Fitted Density for K = "),k),
      cex.main = 1,   font.main= 1, col.main= "red",
)
