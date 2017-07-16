# The main test script for the project
# GenSparseMatrix takes an input density and generates an nxn sparse matrix using a uniform distribution.
# It ensures that the resulting sparse matrix is diagonally dominant

import numpy as np
from scipy import sparse
from SparseStorage import *
from Jacobi import *
from LUPP import *

def GenSparseMatrix(n, density):
    dok = sparse.rand(n, n, density, format='dok', dtype=np.float64, random_state=1)
    M = np.zeros((n, n))
    #print dok
    for item in dok.viewitems():
        M[item[0][0], item[0][1]] = item[1]
        
    for i in xrange(n):
        rowsum = sum(M[i,:])
        rowsum -= M[i, i]
        if rowsum > M[i, i]:
            M[i, i] = rowsum+0.01
    
    return M
    
if __name__ == "__main__":
    A = GenSparseMatrix(1000, 0.8)
    b = np.ones(1000, dtype=np.float64)
    Acsr = MatrixtoCSR(A)
    #x, t = LUPP(A, b)
    #print t
    """
    x, er, it, t = Jacobi(A, b)
    print er, it, t
    """
    x, er, it, t = SparseJacobi(Acsr, b, 1e-20)
    print er
    
