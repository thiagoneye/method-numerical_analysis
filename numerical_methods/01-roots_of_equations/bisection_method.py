"""
The Bisection Method

The bisection method is an interval method, so for its full operation, it is
necessary to initially have two points in such a way that in the interval
between them there is only one root.

Its use is recommended in conjunction with a graphical method for determining
intervals.
"""

# Imports

import matplotlib.pyplot as plt

# Inputs

xl = -100    # Initial estimate for lower value
xu = 100     # Initial estimate for higher value
n = 100      # Number of iterations
ea = 0.001   # Absolute error

def f(x):
    "Curve Equation"
    y = -0.5*x**2 + 2.5*x + 4.5
    return y

# Calculations

iterations = list()
xlast = list()
xl = float(xl)
xu = float(xu)
n = int(n)

for i in range(n):
    xr = (xl + xu)/2
    test = f(xl)*f(xr) # Signal change test at interval

    if (test < 0):
        xu = xr
    elif (test > 0):
        xl = xr
    else:
        if (f(xl) == 0):
            xr = xl
        else:
            xr = xu
        break

    iterations.append(i)
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
