"""
Gauss Elimination with Partial Pivot

This method consists of applying successive elementary operations to a linear
system, to transform it into an easier-to-resolve system, which presents
exactly the same solutions as the original.
"""

# Imports

import numpy as np

# Inputs

AB = np.array([[3, -0.1, -0.2, 7.85],
               [0.1, 7, -0.3, -19.3],
               [0.3, -0.2, 10, 71.4]]) # Augmented matrix
precision = 3                          # Significant algharisms

# Calculations

# Forward elimination

 for i in range(len(AB)):
    if (i < len(AB)):            # Partial pivot
        AB_new = np.flip(AB[i:, i:], 0)
        AB[i:, i:] = AB_new
    for j in range(len(AB)):
        if (j > i):
            m = AB[j,i]/AB[i,i]  # Multiplication factor
            AB[j,:] -= m*AB[i,:] # Variable elimination

# Back substitution
A = AB[:,0:-1]
B = AB[:,-1]
X = np.zeros((len(AB),1))

for i in range(len(A)-1, -1, -1):
    for j in range(len(A)):
        X[i] -= A[i,j]*X[j]
    X[i] += B[i]
    X[i] /= A[i,i]

X = np.around(X, decimals=precision)

# Outputs

print("The problem solving vector is: \n")
print(X)
