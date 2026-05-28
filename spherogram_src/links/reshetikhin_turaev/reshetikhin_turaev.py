from .sparse_array import SparseTensor

import R_matrices

# tangle to graph where nodes are SparseTensors and edges are paris of
# tuples (SparseTensor, index in tensor)
# also, create a dictionary of labels of arcs to the edges

# For curls, do an honest implementation as SparseTensors, so that rot_num
# is only used when creating the tensor network. 
# The said dictionary should then use honest lables of edges instead of labels of arcs...

# The contraction sequence can then still be presented as a list of lists of labels

# When doing contractions, need to update the info in adjacent edges accordingly...ff

