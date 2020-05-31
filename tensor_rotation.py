""" Symbolic Tensor (Quaternion) Rotation

    The following script will perform symbolic tensor rotation using quaternions.
"""
# Author: Ken Sible
# Email:  ksible *at* outlook *dot* com

from sympy import Quaternion as quat
from sympy import Matrix

def rotate(tensor, axis, angle):
    """ Rotate symbolic vector or tensor about an arbitrary axis

        :arg:    3-vector or (3x3)-matrix
        :arg:    rotation axis (normal 3-vector)
        :arg:    rotation angle (in radians)
        :return: rotated tensor (of original type)
    """
    # Quaternion-Matrix Multiplication
    def mul(*args):
        if isinstance(args[0], list):
            q, M = args[1], args[0]
            for i, col in enumerate(M):
                M[i] = col * q
        else:
            q, M = args[0], args[1]
            for i, col in enumerate(M):
                M[i] = q * col
        return M
    # Rotation Quaternion (Axis, Angle)
    q = quat.from_axis_angle(axis, angle)
    if isinstance(tensor[0], list):
        tensor = Matrix(tensor)
        if tensor.shape != (3, 3):
            raise Exception('Invalid Matrix Size')
        # Rotation Formula: M' = (q.(q.M.q*)^T.q*)^T
        M = [quat(0, *tensor[:, i]) for i in range(tensor.shape[1])]
        M = mul(q, mul(M, q.conjugate()))
        for i in range(tensor.shape[1]):
            tensor[:, i] = [M[i].b, M[i].c, M[i].d]
        M = [quat(0, *tensor[i, :]) for i in range(tensor.shape[0])]
        M = mul(q, mul(M, q.conjugate()))
        for i in range(tensor.shape[0]):
            tensor[i, :] = [[M[i].b, M[i].c, M[i].d]]
        return tensor.tolist()
    if isinstance(tensor, list):
        if len(tensor) != 3:
            raise Exception('Invalid Vector Length')
        # Rotation Formula: v' = q.v.q*
        v = q * quat(0, *tensor) * q.conjugate()
        return [v.b, v.c, v.d]
    raise Exception('Unsupported Tensor Type')
