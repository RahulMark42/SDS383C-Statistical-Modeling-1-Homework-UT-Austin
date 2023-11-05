# Author: Rahul Nandakumar (E-id: rn9355)

# Making necessary imports

import matplotlib.pyplot as plt
from scipy.stats import binom
from scipy.stats import poisson

def visualize_distributions(m, p):

    # Specifying initial values. We specify a large value of m here.
    # Here r = range, m = number of trials or samples.
    # p = probability of success in each trial.

    r = list(range(1,m+1))
    l = m*p

    # Creating two distributions. 
    # First one is a Binomial Distribution with p = 0.4 for all values of r from 1 to m
    # Second one is a Poisson Distribution with mean = m*p for all values of r from 1 to m

    dist = [binom.pmf(r, m, p) for r in r]
    dist2 = [poisson.pmf(k = r, mu = l, loc = 0) for r in r]

    # We are taking a range of r = 3000-5000 to better view the details in the plot. 
    # Otherwise, the graphs would be too small to view
    # We plot the Binomial Distribution on a bar graph and the Poisson Distribution using the plot function.

    plt.bar(r[3000:5000], dist[3000:5000], label="Binomial")
    plt.plot(r[3000:5000], dist2[3000:5000], 'r', label = "Poisson")
    plt.legend(loc="upper left")
    plt.show()

visualize_distributions(10000, 0.4)