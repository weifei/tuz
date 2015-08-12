from util import *

def grank(code_matrix, S, T):
    if T[0]:
        H = numpy_to_generic(np.atleast_2d(code_matrix[T[0]].values))
        H2 = numpy_to_generic(np.atleast_2d(code_matrix[T[0]].ix[S[1]].values))
        rk_H = H.Rank()
        rk_H2 = H2.Rank()
    else:
        rk_H = 0
        rk_H2 = 0

    if T[0] or T[1]:
        H2G2 = np.hstack((code_matrix[T[0]].ix[S[1]].values, code_matrix[T[1]].ix[S[1]].values))
        H2G2 = numpy_to_generic(H2G2)
        rk_H2G2 = H2G2.Rank()
    else:
        rk_H2G2 = 0

    return rk_H + rk_H2G2 - rk_H2


def check_alignment_criterion(code_matrix, U1, U2, I1, I2, B, S):
    T0 = [U1, U2 + B]
    T1 = [U1 + I1, U2 + I2 + B]
    grank0 = grank(code_matrix, S, T0)
    grank1 = grank(code_matrix, S, T1)

    if grank1 > grank0:
        pi2 = code_matrix[I1].ix[S[1]].values
        H2U = code_matrix[U1].ix[S[1]].values

        H2UG2UB = np.hstack((code_matrix[U1].ix[S[1]].values, code_matrix[U2 + B].ix[S[1]].values))

        if (not check_column_span(pi2, H2U)) and (check_column_span(pi2, H2UG2UB)):
            return True
    else:
        return False

def get_alignment_solution(code_matrix, S, U1, I1, num):
    # need to identify basis
    T1 = sorted(set(U1).union(set(I1)))
    if U1:
        H2U = code_matrix[U1].ix[S[1]].values
    else:
        H2U = np.zeros((len(S[1]), 0))

    #print "I1", I1
    HI1 = code_matrix[I1].values
    HI1_generic = numpy_to_generic(HI1)
    H2I1 = code_matrix[I1].ix[S[1]].values
    H2 = np.hstack((H2U, H2I1))
    H2_generic = numpy_to_generic(H2)
    H2NULL = H2_generic.NullSpace()
    #print "NULLSPACE", H2NULL
    Nullsize = H2NULL.Size()

    # truncate the null basis
    # the last len(I1) rows
    FBase = H2NULL.SubMatrix(len(U1), len(U1)+len(I1)-1)
    #print "FBase", FBase

    rand_matrix = GenericMatrix((Nullsize[1], num), 0, F.unity, F.Add, F.Subtract, F.Multiply, F.Divide)
    for i in range(rand_matrix.rows):
        rand_matrix.SetRow(i, [F.GetRandomElement() for x in range(num)])

    #print "HI1_generic", HI1_generic
    result = HI1_generic * FBase * rand_matrix
    result = generic_to_numpy(result)
    return result