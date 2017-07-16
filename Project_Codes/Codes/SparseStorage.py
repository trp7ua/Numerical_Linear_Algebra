# Storage functions to exploit sparse structure
# MatrixtoCSR converts the matrix to CSR storage format
# SparseGet accesses the (i,j) element of the sparse matrix

import numpy as np
from pprint import pprint

def MatrixtoCSR(M):
    A = []
    IA = [0]
    JA = []
    (m, n) = np.shape(M)
    for i in xrange(m):
        nz_prev = 0
        for j in xrange(n):
            if M[i][j] != 0:
                A.append(M[i][j])
                JA.append(j)
                nz_prev += 1
                
        IA.append(IA[i]+nz_prev)
                
    csr = {'A':np.array(A), 'IA':np.array(IA, dtype=int), 'JA':np.array(JA, dtype=int), 'm':m, 'n':n}
    return csr

def SparseGet(Acsr, i, j):
    lo = Acsr['IA'][i]
    hi = Acsr['IA'][i+1]-1
    if Acsr['JA'][lo] > j or Acsr['JA'][hi] < j:
        return 0
    
    mid = (lo+hi)/2
        
    while hi-lo > 1:
        if j > Acsr['JA'][mid]:
            lo = mid
        else:
            hi = mid
        
        mid = (lo+hi)/2
               
    if Acsr['JA'][hi] == j:
        return Acsr['A'][hi]
    if Acsr['JA'][lo] == j:
        return Acsr['A'][lo]
        
    return 0
    
# Test cases	
if __name__ == "__main__":
    """
    M = np.array([[0, 0, 0, 0],
                 [5, 8, 0, 0],
                 [0, 0, 3, 0],
                 [0, 6, 0, 0]])
    """
    M = np.array([[10, 20, 0, 0, 0, 0],
                  [0, 30, 0, 40, 0, 0],
                  [0, 0, 50, 60, 70, 0],
                  [0, 0, 0, 0, 0, 80]])
    csr = MatrixtoCSR(M)
    pprint(csr)
    print SparseGet(csr, 3, 2)
                
