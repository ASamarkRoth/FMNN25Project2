''' different updating schemes for the inverse of the hessian '''
import numpy as np

# don't know if this is "good" or "bad" broyden, not sure about the difference....
def broyden(H, delta, gamma):
    ''' H is the current inverse of the hessian, returns the inverse for the next iteration '''
    u = delta - H.dot(gamma)
    a = 1/u.dot(gamma)
    return H + a*np.outer(u, u)

def DFP(H, delta, gamma):
    part2 = np.outer(delta, delta)/delta.dot(gamma)
    part3 = np.outer(H.dot(gamma), gamma.dot(H))/gamma.dot(H).dot(gamma)
    return H + part2 - part3

def BFGS(H, delta, gamma):
    part2 = (1 + gamma.dot(H).dot(gamma)/delta.dot(gamma))*np.outer(delta, delta)/delta.dot(gamma)
    part3 = (np.outer(delta, gamma.dot(H)) + np.outer(H.dot(gamma), delta))/delta.dot(gamma)
    return H + part2 - part3
