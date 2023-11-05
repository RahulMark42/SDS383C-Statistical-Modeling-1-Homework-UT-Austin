# Author: Rahul Nandakumar (E-ID: rn9355)

# Part 1: Checking the distributions.

# Making necessary imports
import scipy.stats
from scipy.stats import cauchy
import matplotlib.pyplot as plt
from random import sample
import numpy as np
import math

# Generating a cauchy distribution
x = np.linspace(-4, 4, num = 5000).tolist()
y = cauchy.pdf(x, loc = 0, scale = 1).tolist()
plt.plot(x, y, label = 'Cauchy Distribution')
plt.legend(loc="upper right")
plt.xlabel('x')
plt.ylabel('y')


# Calculating the means for each sample. We consider the sampling process 1000 times, and take 20 elements in our sample.
means = []
for i in range(1000):
    samples = sample(y, 20)
    mean = np.mean(samples)
    means.append(mean)

# Visualizing the distribution of means, and the cauchy distribution

fig,ax = plt.subplots()
ax.axis([-2, 2, 0, 1])
ax.hist(means, weights=np.ones_like(means)/np.sum(means), bins = 25, label='Histogram of means')
ax.plot(x, y, label = 'Cauchy Distribution')
ax.legend(loc="upper right")
ax.set_xlabel('x')
ax.set_ylabel('y')

# Part 2: Implementing the methods for calculating MLE for the cauchy (theta, 1) distribution

# Writing functions to calculate the gradient of the likelihood function, as derived in the problem. 

# First Gradient
def gradient_likelihood_function(theta, y):
    n = len(y)
    gradient = 0

    for i in range(n):
        gradient  = gradient + (2*(y[i] - theta))/(1 + (y[i] -theta)*(y[i] -theta))

    return gradient

# Second Gradient
def double_gradient_likelihood_function(theta, y):
    n = len(y)
    double_gradient = 0

    for i in range(n):
        double_gradient = double_gradient + 2*(-1 - 3*(y[i] -theta)*(y[i] -theta))/((1 + (y[i] -theta)*(y[i] -theta))*(1 + (y[i] -theta)*(y[i] -theta)))

    return double_gradient

# Function for stepwise gradient ascent
def stepwise_gradient_ascent(y):
    
    theta_old = 5
    c = 0

    while(c <= 10000):

        theta_temp = theta_old

        theta_old = theta_old + 0.5*(gradient_likelihood_function(theta_old, y))

        if(c%1000 == 0): print("c = {}, theta = {}".format(c, theta_old))
        c += 1

    return theta_old

# Function for newton-raphson method
def newton_raphson(y):
    
    theta_old = 5
    c = 0

    while(c <= 10000):

        theta_old = theta_old - (gradient_likelihood_function(theta_old, y))/((double_gradient_likelihood_function(theta_old, y)))

        if(c%1000 == 0): print("c = {}, theta = {}".format(c, theta_old))
        c += 1

    return theta_old

# Function for stochastic gradient ascent
def stochastic_gradient_ascent(y):
    
    theta_old = 5
    c = 0

    # We consider a sample of size 12 from the origial population to calculate the gradients. 
    y_sample = sample(y, 12)

    while(c <= 10000):

        theta_old= theta_old + 0.5*(gradient_likelihood_function(theta_old, y_sample))

        if(c%1000 == 0): print("c = {}, theta = {}".format(c, theta_old))
        c += 1

    return theta_old


# Considering our dataset, and calling the functions to get our answer. 
random_samples = [7.52, 9.92, 9.52, 21.97, 8.39, 8.09, 9.22, 9.37, 7.33, 15.32, 1.08, 8.51, 17.73, 11.20, 8.33, 10.83, 12.40, 14.49, 9.44, 3.67]

print("Newton Raphson Method Iterations: " + str(newton_raphson(random_samples)))
print("Stochastic Gradient Method Iterations: " + str(stochastic_gradient_ascent(random_samples)))
print("Stepwise Gradient Method Iterations: " + str(stepwise_gradient_ascent(random_samples)))
theta_fitted = scipy.stats.cauchy.fit(random_samples, fscale = 1)
print("Actual Theta Fitted using .fit() function: "+ str(theta_fitted))