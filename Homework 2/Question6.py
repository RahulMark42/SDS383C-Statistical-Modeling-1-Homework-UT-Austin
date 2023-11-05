# Author: Rahul Nandakumar (E-ID: rn9355)

# Making necessary imports

import numpy as np
from scipy.stats import weibull_min
import math

# Functions to help us solve for k using newton raphson method
def find_sum_1(data, k):
    sum1 = 0
    for i in range(len(data)):
        sum1 = sum1 + (np.power(data[i], k)*math.log(data[i]))
    return sum1

def find_sum_2(data, k):
    sum2 = 0
    for i in range(len(data)):
        sum2 =  sum2 + np.power(data[i], k)
    return sum2

def find_sum_3(data):
    sum3 = 0
    for i in range(len(data)):
        sum3 = sum3 + np.log(data[i])
    return sum3

def function_value(data, k):
    function_value = (1/k) - (find_sum_1(data, k)/find_sum_2(data, k)) + (find_sum_3(data)/len(data))
    return function_value

def derivative_value(data, k):
    a = find_sum_1(data, k)/(find_sum_2(data, k))
    derivative_value = np.power(a, 2) - 2*(a) - 1/(np.power(k, 2))
    return derivative_value


# Given Data, and we take an initial guess for k
data = [225,171,198,189,189,135,162,135,117,162]
k = 3
c = 0

# Newton Raphson Iterations, we update K until the derivative of the function drops below 0.0001
while(function_value(data, k) >= 0.0001):
    k = k + function_value(data, k)/derivative_value(data, k)
    c += 1

# Estimating lambda once we find k using the closed form solution we formulated. 
sum = 0
for i in range(len(data)):
    sum = sum + math.pow(data[i], k)
lam = math.pow(sum/len(data), 1/k)

print("k,lambda = (" + str(k) + ", " + str(lam) + ")")