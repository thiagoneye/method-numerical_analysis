"""
The Newton-Raphson Method

The Newton-Raphson method is an open method, although open methods generally
show faster convergence than interval methods, they may show a divergence of
the result.
"""

# Imports

import matplotlib.pyplot as plt

# Inputs

x0 = 10      # Initial estimate
n = 100      # Number of iterations
ea = 0.001   # Absolute error

def f(x):
    "Curve equation"
    y = -0.5*x**2 + 2.5*x + 4.5
    return y

def df(x):
    "Derivative of the curve equation"
    y = -x + 2.5
    return y

# Calculations

xlast = list()
iterations = [0]
xlast.append(float(x0))

for i in range(n):
    if (df(xlast[-1]) == 0):
        print("It is impossible to proceed, there is a division by zero.")
        break
    else:
        xr = xlast[-1] - f(xlast[-1])/df(xlast[-1])

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
