"""
@author: Anton Roth, Linus Jangland and Samuel Wiqvist 
"""
import numpy as np
import scipy.linalg as sl
import abc #abstract base classes
from linesearh import *

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
        return 1 #Can we do this? # dont we need to return a n-dim vector? 
        
    @classmethod 
    def get_gradient(cls, f,x0):
        return 2 # shoulden we compute the gradient numericaly now? 
    
'''This class is intended to be inherited'''
class OptimizationMethods:
    __metaclass__ = abc.ABCMeta
    
    def __init__(self, OptimizationProblem, par_line_search = "exact"):
        self.f = OptimizationProblem.f
        self.x0 = OptimizationProblem.x0
        self.g = OptimizationProblem.g
        self.par_line_search = par_line_search 
    
    def newton_procedure():
        xk = self.x0.copy()
        fk = self.f.copy()
        gk = self.g.copy()
        Gk = _initial_hessian(gk, xk)
        line_search = _get_line_search(fk, xk, tol, self.par_line_search)       
        while True:
            sk = _newton_direction(Gk, gk)  
            alphak = line_search(fk, xk, gk, tol = 1e-5)
            xk = xk + alphak*sk
            Gk = _update_hessian(Gk, xk)
            if sl.norm(alphak*sk) < self.tol:
                x = xk
                fmin = f(xk)
                break
        return [x, fmin]
        
    def _newton_direction(Gk, gk):
        '''Computes sk'''
        return (-1)*np.linalg.solve(Gk, gk)
    
    def _get_line_search(fk, xk, tol = 1e-5, par_line_search):
        '''Performs a line_search, finds the alpha which minimizes f(xk+alpha*sk)'''
        if par_line_search == "exact":
            def line_search(fk, xk, gk, tol = 1e-5):
                return _exact_line_search(fk, xk, gk, tol = 1e-5)
            return line_search
        elif par_line_search == "inexact":
            def line_search(fk, xk, tol = 1e-5):
                return _inexact_line_search(fk, xk, gk, tol = 1e-5)
            return line_search
        else:
            def line_search(fk, xk, tol = 1e-5):
                return 1
            return line_search
        
    def _exact_line_search(fk, xk, gk, tol = 1e-5):
        return line_search(fk, xk, tol)
    
    def _inexact_line_search(fk, xk, gk, tol = 1e-5):
        fp = np.dot(xk, gk) # detta borde vara ok 
        alpha_0 = 1 # detta är nog inte rätt
        return inexact_line_search(f, fp, alpha_0, tol)
        
    @abc.abstractmethod
    def _initial_hessian():
        '''Returns the inital hessian'''
    
    @abc.abstractmethod
    def _update_hessian(): # is this also the update method 
        '''Computes and returns hessian or hessian approx. 
            Updates the hessian
        '''
    

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
<<<<<<< HEAD
        
#%%
'''Cholesky on Hinv*gk = sk? or positive definiteness of H'''
A = np.array([[1,2,3],[4,5,6],[7,8,9]])
#L = sl.cho_factor(A)
try:
    L = sl.cho_factor(A)
except sl.LinAlgError:
        raise "The computed Hessian is not positive definite"
gk = np.array([3,3,3])
sk = sl.cho_solve(L,gk)

'''Cholesky decomposition, A = LL*, if unique --> A positive definite
H = J(g)'''

#%%

class OriginalNewton(OptimizationMethods):
    
    def _newton_direction(gk,Gk):
        Gk = 0.5*(G+G.T)
        try:
            L = sl.cho_factor(Gk)
        except sl.LinAlgError:
            print("The computed Hessian was not positive definite!")
        sk = sl.cho_solve(L,gk)
    
    def _update_hessian(G, xk, *args, **kwargs)
        return G(xk)
        
    
=======
#%%
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
>>>>>>> ba9abded135726b9452d6048504c698b6869fb7f
