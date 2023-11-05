library(mvtnorm)
library(MASS)
library(MVNBayesian)
library(invgamma)
library(plotly)
library(mixtools)
dataset = faithful
dataset$waiting = dataset$waiting/100
dataset$eruptions = dataset$eruptions/10

n = 5000 

b = seq(1:n)
z = seq(1:n)
mu1 = 0.25
mu2 = 1.25
sigma1 = 0.2
sigma2 = 1.2
rho= 0.5

gibbs_sampling_multivariate_normal <- function(n,b,z,mu1,mu2,sigma1,sigma2,rho){
  for (i in 1:n)
  {
    z1=rbern(1,rho)
    if (z1==1)
    {b[i]=rnorm(1,mu1,sigma1)}
    else
    { b[i]=rnorm(1,mu2, sigma2)}
  }
  par(mar=c(2,2,2,2))
  par(mfrow = c(2, 1))
  hist(b, breaks=50,col="blue", xlab=" ",ylab=" ",yaxt="n",xaxt='n'  ,main=""  ) 
  title(main="Histogram with density curve for a mixture of two normal distributions", 
        font = 2, cex.main = 2, col.main="blue")
  par(new=TRUE)
  d <- density(b) 
  
  y=b
  d0 <- data.frame(y)       
  write.table(d0, "simulated-data-mix-of-normal.txt", row.names = FALSE,
              col.names = FALSE)
  
  
  d1 <- data.frame(rho, mu1, mu2, sigma1, sigma2)                
  write.table(d1, "fit-mixture-of-normal.txt", row.names = FALSE)
  
  m=5000 
  z = seq(1:n)
  for (i in 1:m)
  {
    print(i)
   
    for (j in 1:n)
    {
      temp1 = rho / sigma1 * exp(-(y[j] - mu1)^2 / 2 / sigma1 ^ 2)
      temp2 = (1 - rho) / sigma2 * exp(-(y[j] - mu2) ^ 2 / 2 / sigma2 ^ 2)
      temp = temp1 / (temp1 + temp2)
      z[j] =  rbinom(1, 1, temp)
    }
    

    t1=sum(z==1)
    t2=sum(z==0)
    rho=rbeta(1, t1+1, t2+1)
    
    a1= sum(z==1)/sigma1^2
    b1=0
    
    for (j in 1:n)
    {
      if (z[j] == 1)
      {
        b1 = b1 +  y[j] /sigma1 ^ 2
      }
    }
    mu1 = rnorm(1, b1/ a1, sqrt(1 / a1))
    

    a1= sum(z==0)/sigma2^2
    b1=0
    for (j in 1:n)
    {
      if (z[j] == 0)
      {
        b1 = b1 +  y[j] / sigma2 ^ 2
      }
    }
    mu2 = rnorm(1, b1/ a1, sqrt(1 / a1))
    

    a1= sum(z==1)
    b1=0
    for (j in 1:n)
    {
      if (z[j] == 1)
      {
        b1 = b1 +  (y[j] -mu1)^2 / 2
      }
    }
    sigma1= sqrt(rinvwishart(1, a1/2, b1))
    
    a1= sum(z==0)
    b1=0
    for (j in 1:n)
    {
      if (z[j] == 0)
      {
        b1 = b1 +  (y[j] - mu2)^2 / 2
      }
    }
    sigma2= sqrt(rinvgamma(1, a1/2, b1))
    
    d1 <- data.frame(rho, mu1, mu2, sigma1, sigma2)
    write.table(d1, "fit-mixture-of-normal.txt", row.names = FALSE,
                col.names = FALSE, append = TRUE)
  }
}

k=2
fit_density = mvnormalmixEM(dataset, k = k)
pi = fit_density$lambda
mu = fit_density$mu
cov = fit_density$sigma
x_vals <- seq(0.1, 0.6, length.out=100)
y_vals <- seq(0.3, 1, length.out=100)
density_calculator <- function(x_vals,y_vals,pi,mu,cov,k){
  count <- 1
  z_vals <- rep(0, length.out=100*100)
  for(i in 1:100){
    sum <- 0
    for(j in 1:100){
      for(K in 1:k){
        sum <- sum + pi[K]*dmvnorm(c(x_vals[i],y_vals[j]),mu[[k]],cov[[k]])
      }
    }
    z_vals[count] <- sum
    count = count + 1
  }
  return(z_vals)
}
z_vals = density_calculator(x_vals, y_vals, pi, mu, cov, k)
fig <-  plot_ly()
fig <- add_trace(
  fig,
  type = 'scatter',
  mode = 'markers',
  opacity=0.75,
  x = dataset$eruptions,
  y = dataset$waiting,
  name = 'Scatter plot'
)
fig <- add_trace(
  fig,
  type = 'contour',
  x = x_vals,
  y = y_vals,
  z= z_vals,
  name = 'Contours Plot', 
  colorscale='Jet'
)
print(fig)