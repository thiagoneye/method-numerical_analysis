"""
The Müller's Method

The Müller's method is an open method for obtaining roots of polynomial
equations, although open methods generally show faster convergence than
interval methods, they may show a divergence of the result.
"""

# Imports

import matplotlib.pyplot as plt

# Inputs

x0 = 4.5     # Initial estimate
n = 100      # Number of iterations
ea = 0.001   # Absolute error

def f(x):
    """ Curve Equation """
    y = (x**3) - (13*x) - 12
    return y

# Calculations

x1 = x0 + 0.5
x2 = x0 + 1
iterations = [0, 1, 2]
xlast = [x0, x1, x2]

for i in range(n):
    f0 = f(x0)
    f1 = f(x1)
    f2 = f(x2)

    h0 = x1 - x0
    h1 = x2 - x1
    d0 = (f1 - f0)/h0
    d1 = (f2 - f1)/h1

    a = (d1 - d0)/(h1 + h0)
    b = a*h1 + d1
    c = f2

    disc = ((b**2) - 4*a*c)**(1/2)

    if (b >= 0):
        xr = x2 + (-2*c)/(b + disc)
    else:
        xr = x2 + (-2*c)/(b - disc)

    iterations.append(i+2)
    xlast.append(xr)

    if (abs(xlast[-1] - xlast[-2]) < ea):
        break

    x0 = x1
    x1 = x2
    x2 = xr

# Outputs

plt.figure()
plt.plot(iterations, xlast)
plt.axvline(x=0,color='k')
plt.axhline(y=0,color='k')
plt.title('Convergence curve')
plt.xlabel('Number of iterations')
plt.ylabel('Root value')

print('The root of the equation is ' + str(round(xr, 3)))
