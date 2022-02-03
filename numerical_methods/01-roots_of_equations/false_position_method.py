"""
The Modified False Position Method

The false position method is an interval method, therefore, for its full
functioning, it is necessary initially to have two points so that in the
interval between them there is only one root.

Its use is recommended in conjunction with a graphical method for determining
intervals.
"""

# Imports

import matplotlib.pyplot as plt

# Inputs

xl = -4     # Initial estimate for lower value
xu = 4      # Initial estimate for higher value
n = 100     # Number of iterations
ea = 0.001  # Absolute error

def f(x):
    "Curve Equation"
    y = -0.5*x**2 + 2.5*x + 4.5
    return y

# Calculations

n = int(n)
xl = float(xl)
xu = float(xu)
fl = f(xl)
fu = f(xu)
il = 0
iu = 0
iterations = list()
xlast = list()

for i in range(n):
    xr = xu - fu*(xl - xu)/(fl - fu)
    fr = f(xr)
    test = fl*fr # Signal change test at interval

    if (test < 0):
        xu = xr
        fu = f(xu)
        iu = 0
        il += 1
        if (il >= 2):
            fl = fl/2
    elif (test > 0):
        xl = xr
        fl = f(xl)
        il = 0
        iu += 1
        if (iu >= 2):
            fu = fu/2
    else:
        if (fl == 0):
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
