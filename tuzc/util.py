import numpy as np
import genericmatrix as gmat
from ffield import *

F = FField(8)

def numpy_to_generic(npmat):
    mat = gmat.GenericMatrix(npmat.shape, 0, F.unity, F.Add, F.Subtract,
                         F.Multiply, F.Divide)
    for i in range(npmat.shape[0]):
        mat.SetRow(i, np.array(npmat[i], dtype=int))
    return mat

def generic_to_numpy(mat):
    return np.array(mat.data, dtype=int)

def generate_rc_matrix(in_matrix, out_column_num):
    in_matrix_generic = numpy_to_generic(in_matrix)
    rand_matrix = gmat.GenericMatrix((in_matrix.shape[1], out_column_num), 0,
                                 F.unity, F.Add, F.Subtract, F.Multiply, F.Divide)
    for i in range(rand_matrix.rows):
        rand_matrix.SetRow(i, np.array([F.GetRandomElement() for x in range(out_column_num)], dtype=int))

    #return in_matrix.dot(rand_matrix)
    result = in_matrix_generic * rand_matrix
    result = generic_to_numpy(result)
    return result

def check_column_span(a, b):
    if not b.size:
        return False

    # check if all columns of a are in the column span of b
    ab_generic = numpy_to_generic(np.hstack((a,b)))
    b_generic = numpy_to_generic(b)
    #print "check span"
    #print ab_generic.Rank(), b_generic.Rank()

    #if np.linalg.matrix_rank(np.hstack((a, b)), tol=eps) > np.linalg.matrix_rank(b, tol=eps):
    if ab_generic.Rank() > b_generic.Rank():
        return False
    else:
        return True

def get_null_space(A):
    u, s, vh = np.linalg.svd(A)
    n = A.shape[1]   # the number of columns of A
    if len(s)<n:
        expanded_s = np.zeros(n, dtype = s.dtype)
        expanded_s[:len(s)] = s
        s = expanded_s
    null_mask = (s <= eps)
    null_space = np.compress(null_mask, vh, axis=0)
    return np.transpose(null_space)
