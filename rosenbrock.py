import numpy as np

def rosenbrock(x):
    return 100*(x[1] - x[0]**2)**2 + (1 - x[0])**2
    
def rosenbrock_function(x):
    def f(x):
        return rosenbrock(x)
    return f

def rosenbrock_grad(x):
    return np.array([-2*x[0]*(1 + 200*(x[1] - x[0]**2)), 200*(x[1] - x[0]**2)])
    
def rosenbrock_gradient(x):
    def g(x):
        return rosenbrock_grad(x)
    return g

if __name__ == "__main__":
    print(rosenbrock((1,2)))
    print(rosenbrock_grad((1,2)))
