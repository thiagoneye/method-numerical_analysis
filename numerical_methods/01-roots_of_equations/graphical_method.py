"""
The Graphical Method

Although it doesn't need input data and gives us the results quickly and
pleasantly, the graphical method is not as accurate as the others, since it is
obtained visually, therefore being subject to human reading errors.
"""

# Imports

import matplotlib.pyplot as plt
import numpy as np

# Inputs

def f(x):
    "Curve Equation"
    y = -0.5*x**2 + 2.5*x + 4.5
    return y

start = -5
stop = 5

# Calculations

x = np.linspace(start, stop)
y = f(x)

# Outputs

plt.figure()
plt.plot(x,y)
plt.axvline(x=0,color='k')
plt.axhline(y=0,color='k')