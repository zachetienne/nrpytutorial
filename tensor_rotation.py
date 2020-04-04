""" Symbolic Tensor (Quaternion) Rotation

The following script will perform symbolic tensor rotation using quaternions.
"""
# Author: Ken Sible
# Email:  ksible@outlook.com

from sympy import Quaternion as quat
from sympy import Matrix
from sympy.functions import transpose

# Input:  tensor = 3-vector or (3x3)-matrix
#         axis   = rotation axis (normal 3-vector)
#         angle  = rotation angle (in radians)
# Output: rotated tensor (of original type)
def rotate(tensor, axis, angle):
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
            raise Exception('Invalid Matrix Dimension')
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
    else:
        if len(tensor) != 3:
            raise Exception('Invalid Vector Length')
        # Rotation Formula: v' = q.v.q*
        tensor = q * quat(0, *tensor) * q.conjugate()
        return [tensor.b, tensor.c, tensor.d]
    raise Exception('Invalid Tensor Type: Matrix or Vector')
