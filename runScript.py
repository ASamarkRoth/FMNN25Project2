'''Set FMNN25Project2 as consoles working directory!!!'''
import numpy as np
import scipy.linalg as sl
import abc #abstract base classes
from gradhess import *
from linesearch import *
from rosenbrock import *

'''Test rosenbrock'''
f = rosenbrock_function(x)
g = rosenbrock_gradient(x)
