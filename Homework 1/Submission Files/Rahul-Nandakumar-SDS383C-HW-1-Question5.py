# Author: Rahul Nandakumar (E-id: rn9355)

# Making necessary imports

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gamma
from scipy.stats import norm
from random import sample
import math


# Function to estimate the values of alpha and beta of gamma distribution with a population size of 100000
def estimators(alpha, beta, sample_size):

    # Generating a random variate of size 100000, and drawing a sample of size sample_size from it. 
    numbers = gamma.rvs(a = alpha, scale = beta, size = 100000)
    samples = sample(numbers.tolist(), sample_size)
    

    #Estimating alpha_hat and beta_hat
    mean = np.mean(samples)
    variance = np.var(samples)
    alpha_hat,_, beta_hat = gamma.fit(samples)
    x = np.linspace(np.min(numbers), np.max(numbers))

    return x, alpha_hat, alpha, beta_hat, beta, samples
    
# Function to plot the estimated density, true density and samples generated from the population
def plot(x, alpha_hat, alpha, beta_hat, beta, samples):

    plt.hist(samples, bins = 50, density=True, label='Sample Histogram')
    plt.plot(x, gamma.pdf(x, a = alpha_hat, scale = beta_hat), label = 'Estimated Density')
    plt.plot(x, gamma.pdf(x, a = alpha, scale = beta), label = 'True Density')
    plt.legend(loc="upper right")
    plt.title('50 samples')
    plt.show()


# Function to plot the values for samples when the convergence is required to be shown. 
def mvd_plotting(alpha, beta, B, sample_size):
    r = []
    for i in range(B):
        x, alpha_hat, alpha, beta_hat, beta, samples = estimators(alpha, beta, sample_size)
        r.append(math.sqrt(sample_size)*(alpha_hat - alpha)/math.sqrt(2*alpha_hat*(alpha_hat + 1)))
    x = np.linspace(np.min(r), np.max(r))
    plt.hist(r, bins = 50, density=True, label='Histogram calculated')    
    plt.plot(x, norm.pdf(x), label = 'Normal Density')
    plt.legend(loc="upper right")
    plt.title('50 samples')
    plt.show()

#Initializing the parameters. 
alpha = 4
beta = 0.6
sample_size = 50

# Estimating Gamma Parameters using method of moments. 
x, alpha_hat, alpha, beta_hat, beta, samples = estimators(alpha, beta, sample_size)

# Uncomment and run the remaining results.

#plot(x, alpha, alpha_hat, beta, beta_hat, samples)
#mvd_plotting(alpha, beta, 50, 50)
#mvd_plotting(alpha, beta, 500, 50)
mvd_plotting(alpha, beta, 1000, 50)