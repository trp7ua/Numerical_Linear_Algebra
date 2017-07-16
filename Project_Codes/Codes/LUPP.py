# LU Decomposition with partial pivoting

import numpy as np
from Operations import MMultiply
import time

def LUPP(A, b):
    t1 = time.time()
    LU = np.copy(A)
    n = len(b)
    piv = {}
    for i in xrange(n-1):
        mu = np.argmax(LU[i:, i])+i
        piv[i] = mu
        LU[[i, mu]] = LU[[mu, i]]
        if LU[i, i] != 0:
            LU[i+1:, i] = LU[i+1:, i]/LU[i, i]
            LU[i+1:, i+1:] = LU[i+1:, i+1:] - MMultiply(LU[i+1:, i], LU[i, i+1:], True)
            
    Pb = np.copy(b)
    for i in xrange(n):
        if i in piv:
            Pb[i], Pb[piv[i]] = Pb[piv[i]], Pb[i]
            
    ForwardSub(LU, Pb, True)
    BackwardSub(LU, Pb) 
    t2 = time.time() - t1
    return Pb, t2

def ForwardSub(L, b, unitdiag=False):
    n = len(b)
    if not unitdiag:
        b[0] = b[0]/L[0,0]
    for i in xrange(1, n):
        b[i] = (b[i] - MMultiply(L[i, :i], b[:i]))
        if not unitdiag:
            b[i] /= L[i, i]
            
    return b

def BackwardSub(U, b):
    n = len(b)
    b[n-1] = b[n-1]/U[n-1, n-1]
    for i in xrange(n-2, -1, -1):
        b[i] = (b[i] - MMultiply(U[i, i+1:], b[i+1:]))/U[i, i]
    
# Test cases	
if __name__ == "__main__":
    
    A = np.array([[10., -1., 2., 0.],
              [-1., 11., -1., 3.],
              [2., -1., 10., -1.],
              [0.0, 3., -1., 8.]], dtype=np.float64)
    """
    A = np.array([[3, 17, 10],
              [2, 4, -2],
              [6, 18, -12]], dtype=float)
    """
    b = np.array([6., 25., -11., 15.], dtype=np.float64)
    x = LUPP(A, b)
    print x    
    
            
