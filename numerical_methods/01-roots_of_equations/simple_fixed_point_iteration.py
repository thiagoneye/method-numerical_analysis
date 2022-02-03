"""
Simple Fixed-Point Iteration

The fixed-point iteration method is an open method, although open methods
usually show faster convergence than interval methods, they can show a
divergence of the result.
"""

# Imports

from math import exp
import matplotlib.pyplot as plt

# Inputs

x0 = 0       # Initial estimate
n = 100      # Number of iterations
ea = 0.001   # Absolute error

def f(x):
    "Curve Equation"
    y = exp(-x) - x
    return y

# Calculations

iterations = [0]
xlast = list()
xlast.append(float(x0))

def g(x,f):
    "Recursive Equation"
    y = f(x) + x
    return y

for i in range(n):
    xr = g(xlast[-1],f)

    iterations.append(i+1)
    xlast.append(xr)

    if (len(xlast) > 1):
        if (abs(xlast[-1] - xlast[-2]) < ea):
            break

# Outputs

plt.figure()
plt.plot(iterations, xlast)
plt.axvline(x=0,color='k')
plt.axhline(y=0,color='k')
plt.title('Convergence curve')
plt.xlabel('Number of iterations')
plt.ylabel('Root value')

print('The root of the equation is ' + str(round(xr, 3)))
