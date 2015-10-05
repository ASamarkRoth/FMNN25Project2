'''Set FMNN25Project2 as consoles working directory!!!'''

from rosenbrock import *
from optimizationNewton import *
from linesearch import *
import numpy as np

print("what is happening")

#%%
'''Test with easy function'''
def a(x):
    return x[0]**2 + 2*x[1]**4
def b(x):
    return np.array([2*x[0], 8*x[1]**3])
def c(alpha):
	x0 = np.array([3,3])
	s = np.array([-1,-1])
	return a(x0 + alpha*s)

#print(type(a), type(b))
#print(line_search(c, 0))
#print(exact_line_search(c, 0))



prob = OptimizationProblem(a, g = b, x0 = np.array([7,11]))

solver = OriginalNewton(prob)
#print("Linus testfunc")
#print(solver.newton_procedure(par_line_search="None"))
#print("test on a,b thing", solver.newton_procedure())

'''No gradient inserted'''
#prob = OptimizationProblem(a, x0 = np.array([7,11]))

#%%
'''Test rosenbrock'''
f = rosenbrock
g = rosenbrock_grad

prob = OptimizationProblem(f, x0 = np.array([0,0]))

solver = OriginalNewton(prob)

#print("Rosenbrock test")
#print(solver.newton_procedure(par_line_search= "exact"))

solver = OptimizationMethodsBroyden(prob)
print(solver.newton_procedure(par_line_search = "exact"))

