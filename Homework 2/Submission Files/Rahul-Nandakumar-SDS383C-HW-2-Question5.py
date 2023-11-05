# Author: Rahul Nandakumar (E-ID: rn9355)

# Making necessary imports
from scipy.stats import expon
from scipy.stats import uniform
import numpy as np
import matplotlib.pyplot as plt
from random import sample
import math


# We first generate our uniform distribution with theta = 10, and population size = 10000

xs = np.random.uniform(low=0.0, high=10, size=10000)
xs = np.sort(xs)

# We sample 1000 values from this population, for a total of 100000 iterations. Using this, we calculate the MLE value and zn. 
# Note, running time is large because of large population size (approx 3 mins). 
zn = []
for i in range(100000):
    samples = sample(xs.tolist(), 1000)
    zn.append(1000*(10 - max(samples))/10)

x = np.linspace(0, 10, 500)
y = [np.exp(-val) for val in x]

# We plot the histogram of zn's vs the exponential distribution where lambda = 1
plt.hist(np.sort(zn), bins = 50, density=True, label="Histogram of zn")
plt.plot(x, y, label = 'Exponential Distribution')
plt.legend(loc='upper right')