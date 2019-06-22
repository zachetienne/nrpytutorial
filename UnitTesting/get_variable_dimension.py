# Takes in a tensor [tensor] and returns the rank of that tensor, along with the length of the tensor
# scalar -> rank 0 tensor with length 1 -> 0, 1
# vector with 5 elements -> rank 1 tensor with length 5 -> 1, 5
# tensor with NxN elements -> rank 2 tensor with length N -> 2, N
# ...
# Raises [IndexError] if the first argument of the tensor of any rank is the empty list.
# Assumes that the tensor being passed in is consistent in dimension, i.e. the 0'th index of the list has the same
# dimension as any other index of the list.

# Called by module_dict_to_list


def get_variable_dimension(tensor):
    
    dim = 0
    length = 1

    while isinstance(tensor, list):
        if dim == 0:
            length = len(tensor)
        dim += 1
        tensor = tensor[0]

    return dim, length
