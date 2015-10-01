

import numpy as np
import scipy.linalg as sl
import abc #abstract base classes

'''Suggestion of how to compute grid for computation of g, perhaps also H?
Options: Create grid from the start which we use throughout whole optimisation or 
everytime we compute gk and Hk we compute new grid? Isn't it necessary to introduce a new grid everytime to assure
that xk+1 exists in the grid? '''

xk = np.array([3,4])
step =  0.001
bounds = [((xk[j]-step*(nbr_points//2)),xk[j]+1.1*step*(nbr_points//2)) for j in range(len(xk))]
gridk = np.mgrid[[slice(bounds[j][0], bounds[j][1], step) for j in range(len(xk))]]

f = gridk[0]** + gridk[1]*2
g,dx = np.gradient(f)

#%% 
'''Design suggestion'''

class OptimizationProblem:
    '''A class which generates the necessary components to handle and solve an
    optimization problem defined by an input function f'''
    def __init__(self, f, *args, **kwargs):
        self.f = f
        self.x0 = kwargs.get('x0')
        self.g = kwargs.get('g')
        if(self.x0 == None):
            self.x0 = self.guess_x(f)
        if(self.g == None):
            self.g = self.get_gradient(f,x0)
    
    @staticmethod
    def guess_x(f):
        return 1 #Can we do this?
        
    @classmethod
    def get_gradient(f,x0):
        return 2
    
'''This class is intended to be inherited'''
class OptimizationMethods:
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, OptimizationProblem):
        self.f = OptimizationProblem.f
        self.x0 = OptimizationProblem.x0
        self.g = OptimizationProblem.g
    
    def newton_procedure():
        xk = self.x0.copy()
        fk = self.f.copy()
        gk = self.g.copy()
        while True:
            gk = OptimizationProblem.get_gradient(fk,x0)
            Hk = self.hessian(gk,fk,x0)
            sk = np.dot(inv(Hk),gk)
            #Perhaps the two rows above should be one method, newton_direction
            alphak = self.line_search(fk, xk, tol = 1e-5)
            self.xk = self.xk + alphak
            if sl.norm(np.dot(alphak,sk)) < self.tol:
                x = xk
                fmin = f(xk)
                break
        return x, fmin
    
    @abc.abstractmethod
    def hessian():
        '''computes and returns hessian or hessian approx.'''
        
    @abc.abstractmethod
    def newton_direction():
        '''Computes newton_direction'''
    
    @abc.abstractmethod
    def line_search():
        '''Performs a line_search, finds the alpha which minimizes f(xk+alpha*sk)'''


#%%Testing abstract stuff
class Newton:
    
    def __init__(self, f, x0, tol, *args, **kwargs):    
        self.f = f
        self.x0 = x0
        self.tol = tol
        self.g = kwargs.get('g')   
    
    @abc.abstractmethod
    def get_new_x():
        '''Defines new x'''
    
class Newton2(Newton):
    
    def get_new_x(self,k):
        return self.x0+k