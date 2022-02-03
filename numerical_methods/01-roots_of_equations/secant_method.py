"""
The Secant Method

The secant method is an open method, although open methods generally show
faster convergence than interval methods, they may show a divergence of the
result.
"""

# Import

from math import exp
import matplotlib.pyplot as plt

# Inputs

x0 = 0       # Initial estimate
x1 = 1       # Initial estimate
n = 100      # Number of iterations
ea = 0.001   # Absolute error

def f(x):
    "Curve equation"
    y = exp(-x) - x
    return y

# Calculations

iterations = [0, 1]
xlast = list()
xlast.append(float(x0))
xlast.append(float(x1))

for i in range(n):
    if ((f(xlast[-2]) - f(xlast[-1])) == 0):
        print("It is impossible to proceed, there is a division by zero.")
        break
    else:
        xr = xlast[-1] - (f(xlast[-1])*(xlast[-2] - xlast[-1]))/(f(xlast[-2]) - f(xlast[-1]))

        iterations.append(i+2)
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
