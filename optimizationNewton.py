import numpy as np
import scipy.linalg as sl
import gradhess
import linesearch
from abc import ABCMeta
from abc import abstractmethod

"""
@author: Anton Roth, Linus Jangland and Samuel Wiqvist 
"""
'''Design suggestion'''

class OptimizationProblem:
    '''A class which generates the necessary components to handle and solve an
    optimization problem defined by an input function f'''
    def __init__(self, f, g = None, x0 = None):
        self.f = f
        if (g == None):
            self.g = gradhess.get_gradient(f)
        else:
            self.g = g
        if (x0 == None):
            self.x0 = 3
        else:
            self.x0 = x0
    
    @staticmethod
    def guess_x(f):
        return 1 #Can we do this? # dont we need to return a n-dim vector? 
    
'''This class is intended to be inherited'''
class OptimizationMethods(metaclass=ABCMeta):
    
    def __init__(self, OptimizationProblem, par_line_search = "exact"):
        self.f = OptimizationProblem.f
        self.x0 = OptimizationProblem.x0
        self.g = OptimizationProblem.g
        self.par_line_search = par_line_search 
    
    def newton_procedure(self):
        xk = self.x0
        fk = self.f
        fp = fk # well ok..
        gk = self.g
        Gk = self._initial_hessian(xk, gk)
        line_search = self._get_line_search()
        while True:
            sk = self._newton_direction(xk, gk, Gk)
            def f_linear(alpha):
                return fk(xk + alpha*sk)
            alphak = line_search(f_linear, fp, 0)
            print("xk =", xk, ", sk = ", sk, ", alpha_k = ", alphak)
            xk = xk + alphak*sk
            Gk = self._update_hessian(xk, gk)
            if np.linalg.norm(alphak*sk) < 1e-5: # self.tol:
                x = xk
                fmin = fk(xk)
                break
        return [x, fmin]
        
    @abstractmethod
    def _newton_direction(self, xk, gk, Gk):
        '''Computes sk'''
    
    def _get_line_search(self):
        '''Performs a line_search, finds the alpha which minimizes f(xk+alpha*sk)'''
        if self.par_line_search == "exact":
            def line_search(f, fp, alpha_0):
                return linesearch.line_search(f, alpha_0)
            return line_search
        elif self.par_line_search == "inexact":
            def line_search(f, fp, alpha_0):
                return linesearch.inexact_line_search(f, fp, alpha_0)
            return line_search
        else:
            def line_search(f, fp, alpha_0):
                return 1
            return line_search

    def _exact_line_search(fk, xk, gk, tol = 1e-5):
        return line_search(fk, xk, tol)
    
    def _inexact_line_search(fk, xk, gk, tol = 1e-5):
        fp = np.dot(xk, gk(xk)) # detta borde vara ok 
        alpha_0 = 1 # detta är nog inte rätt
        return inexact_line_search(f, fp, alpha_0, tol)
        
    @abstractmethod
    def _initial_hessian(self, xk, g):
        pass
        '''Returns the inital hessian'''
    
    @abstractmethod
    def _update_hessian(self): # is this also the update method 
        pass
        '''Computes and returns hessian or hessian approx. 
            Updates the hessian
        '''
    

class OriginalNewton(OptimizationMethods):
    
    def _newton_direction(self, xk, g, G):
        Gk = gradhess.calc_hessian(g, xk)
        Gk = 0.5*(Gk + Gk.T)
        sk = -sl.solve(Gk, g(xk)) # just trying to see if it works
        return sk
        '''try:
            L = sl.cho_factor(Gk)
        except sl.LinAlgError:
            print("The computed Hessian was not positive definite!")
        sk = -sl.cho_solve(L, g(xk))
        return sk'''
    
    def _update_hessian(self, xk, g):
        return gradhess.calc_hessian(g, xk)
        
    def _initial_hessian(self, xk, g):
        return gradhess.calc_hessian(g, xk)
