"""
@author: Anton Roth, Linus Jangland and Samuel Wiqvist 
"""
import numpy as np
import scipy.linalg as sl
import abc #abstract base classes
<<<<<<< HEAD

'''Suggestion of how to compute grid for computation of g, perhaps also H?
Options: Create grid from the start which we use throughout whole optimisation or 
everytime we compute gk and Hk we compute new grid? Isn't it necessary to introduce a new grid everytime to assure
that xk+1 exists in the grid? '''

xk = np.array([3,4])
step =  0.001
nbr_points = 5
bounds = [((xk[j]-step*(nbr_points//2)),xk[j]+1.1*step*(nbr_points//2)) for j in range(len(xk))]
gridk = np.mgrid[[slice(bounds[j][0], bounds[j][1], step) for j in range(len(xk))]]

f = gridk[0]** + gridk[1]*2
g,dx = np.gradient(f)

#%% 
'''Design suggestion'''

=======
 
>>>>>>> ba9abded135726b9452d6048504c698b6869fb7f
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
    
    def __init__(self, OptimizationProblem, par_line = "exact"):
        self.f = OptimizationProblem.f
        self.x0 = OptimizationProblem.x0
        self.g = OptimizationProblem.g
        self.par_line = par_line 
    
    def newton_procedure():
        xk = self.x0.copy()
        fk = self.f.copy()
        gk = self.g.copy()
        Gk = _initial_hessian(gk, xk)
        while True:
<<<<<<< HEAD
            gk = OptimizationProblem.get_gradient(fk,x0)
            Hk = self.hessian(gk,fk,x0)
            sk = np.dot(inv(Hk),gk)
            #Perhaps the two rows above should be one method, newton_direction
            alphak = self.line_search(fk, xk, tol = 1e-5)
            self.xk = self.xk + alphak
            Hk = self.hessian()
            if sl.norm(np.dot(alphak,sk)) < self.tol:
=======
            sk = _newton_direction(Gk, gk)   
            if par_line == "exact": 
                alphak = self._line_search(fk, xk, tol = 1e-5, "exact")
            elif par_line == "inexact":
                alphak = self._line_search(fk, xk, tol = 1e-5, "inexact")
            else:
                alphak = 1
            xk = xk + alphak*sk
            Gk = _update_hessian(Gk, self.xk)
            if sl.norm(alphak*sk) < self.tol:
>>>>>>> ba9abded135726b9452d6048504c698b6869fb7f
                x = xk
                fmin = f(xk)
                break
        return [x, fmin]
        
    def _newton_direction(Gk, gk):
        '''Computes sk'''
        return (-1)*np.np.linalg.solve(Gk, gk)
        
    @abc.abstractmethod
    def _initial_hessian():
        '''Returns the inital hessian'''
    
    @abc.abstractmethod
    def _update_hessian(): # is this also the update method 
        '''Computes and returns hessian or hessian approx. 
            Updates the hessian
        '''
    
    @abc.abstractmethod
    def _line_search():
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
