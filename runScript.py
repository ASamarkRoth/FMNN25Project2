'''Set FMNN25Project2 as consoles working directory!!!'''

from rosenbrock import *
from optimizationNewton import *
from linesearch import *
import numpy as np

print("what is happening")

'''Test rosenbrock'''
f = rosenbrock
g = rosenbrock_grad


def a(x):
    return x[0]**2 + 2*x[1]**4
def b(x):
    return np.array([2*x[0], 8*x[1]**3])
def c(alpha):
	x0 = np.array([3,3])
	s = np.array([-1,-1])
	return a(x0 + alpha*s)

print(type(a), type(b))
print(line_search(c, 0))


prob = OptimizationProblem(a, b, np.array([7,11]))
solver = OriginalNewton(prob)
print(solver.newton_procedure())
>>>>>>> fb0ced4d9c40b6601f161ac1bbf34472d6d5166b
