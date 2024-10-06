''' module which defines simple vector operations '''
import logging
import math

logger = logging.getLogger("piou.util.vector_op")

def cross_product(a,b,factor=1.0):
    ''' obtain the cross product of vector a and b. 
        if len(a)<3 or len(b)<3, will return None '''

    if a is None or b is None:
        return None

    if len(a)<3 or len(b)<3:
        return None

    a_cross_b=[]
    a_cross_b.append((a[1]*b[2]-a[2]*b[1])*factor)
    a_cross_b.append((a[2]*b[0]-a[0]*b[2])*factor)
    a_cross_b.append((a[0]*b[1]-a[1]*b[0])*factor)
    return a_cross_b

def diff(a,b):
    ''' calculate the difference of the two vector.
        mind that len(a) must be equal to len(b)'''
    if len(a)!=len(b):
        return None

    ret=[]
    n=len(a)
    for i in range(n):
        ret.append(0.0)
    for i in range(n):
        ret[i] = a[i]-b[i]
    return ret

def dot_product(a,b,factor=1.0):
    ''' obtain the dot product of vector a and b. 
        will return None if len(a)!=len(b) '''

    if a is None or b is None:
        return None

    if len(a)!=len(b):
        return None

    ret = 0
    for i in range(len(a)):
        ret += a[i]*b[i]*factor
    return ret

def norm(vec):
    ''' get the norm of the specified vector '''
    l = 0.0
    for v in vec:
        l += v*v
    return math.sqrt(l)

def get_inverse_matrix(mat):
    ''' get the inverse matrix of a 3x3 matrix. 
        beaware that this function will work only for a 3x3 matrix. '''
    a11=mat[0][0]
    a12=mat[0][1]
    a13=mat[0][2]
    a21=mat[1][0]
    a22=mat[1][1]
    a23=mat[1][2]
    a31=mat[2][0]
    a32=mat[2][1]
    a33=mat[2][2]
    deter = a11*a22*a33+a21*a32*a13+a31*a12*a23-\
            a11*a32*a23-a31*a22*a13-a21*a12*a33
    if math.fabs(deter) < 1e-14:
        logger.error("invalid cell; the determinant is 0.")
        return None
    ret = [] 
    for i in range(3):
        row = []
        for j in range(3):
            row.append([])
        ret.append(row)

    invdet = 1.0/deter
    ret[0][0] = invdet*(a22*a33-a23*a32)
    ret[0][1] = invdet*(a13*a32-a12*a33)
    ret[0][2] = invdet*(a12*a23-a13*a22)
    ret[1][0] = invdet*(a23*a31-a21*a33)
    ret[1][1] = invdet*(a11*a33-a13*a31)
    ret[1][2] = invdet*(a13*a21-a11*a23)
    ret[2][0] = invdet*(a21*a32-a22*a31)
    ret[2][1] = invdet*(a12*a31-a11*a32)
    ret[2][2] = invdet*(a11*a22-a12*a21)

    logger.debug("checking A x A^-1...")
    for i in range(3):
        for j in range(3):
            element=0
            for k in range(3):
                element += mat[i][k]*ret[k][j]
            logger.debug("matrix element "+str(i)+", "+str(j)+" : "+str(element))

    return ret

