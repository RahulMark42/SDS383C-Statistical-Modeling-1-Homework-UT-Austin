library(plotly)
n <- 5000
burn_in <- 2000
thin <- 3
mu1 <- 0
mu2 <- 2
rho <- 0.75
sigma <- (1-rho^2)
covars <- matrix(c(1, .75, .75, 1), 2)

gibbs_sampler <- function(obs){
  sample <- matrix(obs, nrow = n+1, ncol = 2, byrow = T)
  for(i in 1:n){
    obs[1] <- rnorm(1, mu1+rho*(obs[2] -mu2), sigma)
    obs[2] <- rnorm(1, mu2+rho*(obs[1] -mu1), sigma)
    sample[i+1,] <- obs
  }
  return(sample)
}

sampled_values = gibbs_sampler(c(0,0))
burned_in_sampled_values = sampled_values[burn_in:n,]
thinned_sampled_values = burned_in_sampled_values[c(FALSE, FALSE, TRUE),]


x <- seq(mu1-3, mu1+3, length= 500)
y <- seq(mu2-3, mu2+3, length= 500)

sigma21 <- covars[2,1]
sigma12 <- covars[1,2]
sigma11 <- covars[1,1]
sigma22 <- covars[2,2]
z <- function(x,y){ 
  z <- exp(-(sigma22*(x-mu1)^2+sigma11*(y-mu2)^2-2*sigma12*(x-mu1)*(y-mu2))/(2*(sigma11*sigma22-sigma12^2)))/(2*pi*sqrt(sigma11*sigma22-sigma12^2)) }
f <- t(outer(x,y,z))
df <- data.frame (u <- thinned_sampled_values[,1],
                  w <- thinned_sampled_values[,2])
fig <-  plot_ly()
fig <- add_trace(
  fig,
  type = 'contour',
  x = x,
  y = y,
  z=f,
  name = 'Contours of actual density' 
)
fig <- add_trace(
  fig,
  type = 'scatter',
  mode = 'markers',
  opacity=0.75,
  x = df$u,
  y = df$w,
  name = 'Scatter plot of sampled values'
)
fig %>%
  layout(title = 'Gibbs Sampling values scatter plot overlayed over contours of true density', plot_bgcolor = "#e5ecf6")
print(fig)
plot(thinned_sampled_values[,1], type = 'l', col='green', xlab='t values', ylab='Sampled Values of X')
plot(thinned_sampled_values[,2], type = 'l', col='blue', xlab='t values', ylab='Sampled Values of Y')
autocorrelation_1<-acf(thinned_sampled_values[,1],plot=TRUE)
autocorrelation_1<-acf(thinned_sampled_values[,2],plot=TRUE)