"""
@author: Anton Roth, Linus Jangland and Samuel Wiqvist 
"""
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
            sk = _newton_direction(Gk, gk, xk)  
            alphak = line_search(fk, xk, gk, tol = 1e-5)
            xk = xk + alphak*sk
            Gk = _update_hessian(Gk, xk)
            if sl.norm(alphak*sk) < self.tol:
                x = xk
                fmin = f(xk)
                break
        return [x, fmin]
        
    def _newton_direction(Gk, gk, xk):
        '''Computes sk'''
        return (-1)*np.linalg.solve(Gk, gk(xk))
    
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
        fp = np.dot(xk, gk(xk)) # detta borde vara ok 
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
        
    def _initial_hessian(gk,xk):
        G = get_hessian(gk)
        return(G(xk))
    

 
