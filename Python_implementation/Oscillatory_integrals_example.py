"""
Example file illustrating the use of Gauss-Legendre integration for oscillatory integrals
"""

import numpy as np
from Gauss_Legendre_integration import GaussLegendre_discretization, legpoly_int


"""
Example#1: f(x) = 1 + 2x^2 - x^3, exact values:
    
    int(f(x)*cos(x),-1, 1) = 2.6395
    int(f(x)*sin(x),-1, 1) = -0.3542

"""
def f(x):
    return 1 + 2*x**2 - x**3
# no. of sub-intervals
N = 2; 
# nodes
x_nodes = GaussLegendre_discretization(-1.0, 1.0, N)
f_vals = f(x_nodes)
w1_freq = 1.
f_int = legpoly_int(x_nodes, f_vals, w1_freq)

print(f"Oscillatory Integral of f(x) = 1 + 2*x**2 - x**3 with frequency w = {w1_freq}:")
print(f"cos term: {np.real(f_int):.4f} (exact: 2.6395),\
      sin term: {np.imag(f_int):.4f} (exact: -0.3542)")
      
      
"""
Example#2: g(x) = (1 + 2x^2 - x^3)*tanh(x), exact values:
    
    int(f(x)*cos(5x),-1, 1) = 0.1250
    int(f(x)*sin(5x),-1, 1) = -0.5831

"""
def g(x):
    return (1 + 2*x**2 - x**3) * np.tanh(x)

g_vals = g(x_nodes)
w2_freq = 5.
g_int = legpoly_int(x_nodes, g_vals, w2_freq)

print(f"Oscillatory Integral of g(x) = (1 + 2*x**2 - x**3)*tanh(x) with frequency w = {w2_freq}:")
print(f"cos term: {np.real(g_int):.4f} (exact: 0.1250),\
      sin term: {np.imag(g_int):.4f} (exact: -0.5831)")   
