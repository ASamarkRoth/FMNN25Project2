"""
@author: Anton Roth, Linus Jangland and Samuel Wiqvist 
"""

import numpy as np
import scipy.linalg as sl
from gradhess import *
import linesearch
from hessupdate import *
from abc import ABCMeta
from abc import abstractmethod


class OptimizationProblem:
    '''A class which generates the necessary components to handle and solve an
    optimization problem defined by an input function f'''
        
    def __init__(self, f, g = None, x0 = None):
        self.f = f
        if x0 == None:
            self.x0 = np.array([4,4])
        else:
            self.x0 = np.asarray(x0)
        if g == None:
            self.g = get_gradient(f)
        else:
            self.g = g
    
    @staticmethod
    def guess_x(f):
        return 1 #Can we do this? # dont we need to return a n-dim vector? 
    
class OptimizationMethods(metaclass=ABCMeta):
    '''This class is intended to be inherited'''
    
    def __init__(self, OptimizationProblem):
        self.f = OptimizationProblem.f
        self.x0 = OptimizationProblem.x0
        self.g = OptimizationProblem.g
    
    def newton_procedure(self, par_line_search = "exact"):
        xk = self.x0
        Gk = self._initial_hessian(xk, self.g)
        line_search = self._get_line_search(par_line_search)
        while True:
            sk = self._newton_direction(xk, self.g, Gk)
            def f_linear(alpha):
                return self.f(xk + alpha*sk)
            def f_linear_derivative(alpha):
                return self.g(xk + alpha*sk).dot(sk)
            alphak = line_search(f_linear, f_linear_derivative, 0)
            print("xk =", xk, ", fk = ", self.f(xk), ", sk = ", sk, ", alpha_k = ", alphak)
            xnew = xk + alphak*sk
            Gk = self._update_hessian(xk, xnew, self.g, Gk)
            xk = xnew
            if sl.norm(alphak*sk) < 1e-5: # self.tol:
                x = xk
                fmin = self.f(xk)
                break
        return [x, fmin]
        
    
    def _newton_direction(self, xk, g, G):
        '''Computes sk'''
        sk = -G.dot(g(xk))
        return sk
    

    def _get_line_search(self, par_line_search):
        '''Assign a line search algorithm to line_search'''
        if par_line_search == "exact":
            def line_search(f, fp, alpha_0):
                return linesearch.exact_line_search(f, alpha_0)
            return line_search
        elif par_line_search == "inexact":
            def line_search(f, fp, alpha_0):
                return linesearch.inexact_line_search(f, fp, alpha_0)
            return line_search
        elif par_line_search == "None":
            def line_search(f, fp, alpha_0):
                return 1
            return line_search

    def _exact_line_search(fk, xk, gk, tol = 1e-5):
        '''Performs a exact line search'''
        return line_search(fk, xk, tol)
    
    def _inexact_line_search(fk, xk, gk, tol = 1e-5):
        '''Performs a inexact line search, finds the alpha which minimizes f(xk+alpha*sk)'''
        fp = np.dot(xk, gk(xk)) # detta borde vara ok 
        alpha_0 = 1 # detta är nog inte rätt
        return inexact_line_search(f, fp, alpha_0, tol)
        
    @abstractmethod
    def _initial_hessian():
        pass
        '''Returns the inital hessian'''
    
    @abstractmethod
    def _update_hessian(): # is this also the update method 
        pass
        '''Computes and returns hessian or hessian approx. 
            Updates the hessian
        '''
    

class OriginalNewton(OptimizationMethods):
    
    def _newton_direction(self, xk, g, G):
        Gk = calc_hessian(g, xk)
        Gk = 0.5*(Gk + Gk.T)
        #sk = -sl.solve(Gk, g(xk)) # just trying to see if it works
        try:
            L = sl.cho_factor(Gk)
        except sl.LinAlgError:
            print("The computed Hessian was not positive definite!")
        sk = -sl.cho_solve(L, g(xk))
        return sk
        
    
    def _update_hessian(self, xk, xnew, g, G):
        return calc_hessian(g, xnew)
        
    def _initial_hessian(self, xk, g):
        return calc_hessian(g, xk)

class OptimizationMethodsBroyden(OptimizationMethods):

    def _update_hessian(self, xk, xk_old, gk, G): 
        delta = xk - xk_old
        gamma = gk(xk) - gk(sk)
        return broyden(H, delta, gamma)
        
    def _initial_hessian(self, xk, gk):
        return np.identity(len(self.x0))