import numpy as np

def rosenbrock(x):
    return 100*(x[1] - x[0]**2)**2 + (1 - x[0])**2
    
def rosenbrock_function():
    def f(x):
        return rosenbrock(x)
    return f

def rosenbrock_grad(x):
    return np.array([-2*x[0]*(1 + 200*(x[1] - x[0]**2)), 200*(x[1] - x[0]**2)])
    
def rosenbrock_gradient():
    def g(x):
        return rosenbrock_grad(x)
    return g