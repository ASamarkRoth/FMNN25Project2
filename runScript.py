'''Set FMNN25Project2 as consoles working directory!!!'''

from optimizationNewton import *
import numpy as np

'''Test rosenbrock'''
from rosenbrock import *

f = rosenbrock
g = rosenbrock_grad
g(np.array([1,2]))

#Minimum is at (1,1) and f(1,1) = 0
guess = np.array([1.5,1.5])
OP = OptimizationProblem(f, x0 = guess)


