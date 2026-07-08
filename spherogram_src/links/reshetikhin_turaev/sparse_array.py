from itertools import product as cartesian_product

class SparseArray:
    """
    A sparse array supporting arbitrary-dimensional indexing via tuples.

    Internally stores only non-default entries in a dict keyed by index tuples.
    Indices can be integers or tuples of integers of any length.
    """
    __slots__ = ('_data', '_default', '_rank', '_shape')

    def __init__(self, shape, data=None, default=0):
        self._data = {}
        self._default = default
        self._shape = tuple(shape)
        self._rank = len(self._shape)

        if data is not None:
            if isinstance(data, dict):
                data = data.items()
            for k, v in data:
                self[k] = v

    def __eq__(self, other):
        if not isinstance(other, SparseArray):
            return False
        return (self._shape == other._shape and
                self._default == other._default and
                self._data == other._data)

    def _key(self, index):
        if isinstance(index, tuple):
            return index
        return (index,)

    def __getitem__(self, index):
        return self._data.get(self._key(index), self._default)

    def __setitem__(self, index, value):
        key = self._key(index)
        if value == self._default:
            self._data.pop(key, None)
        else:
            if len(key) != self._rank:
                raise ValueError(
                    f'key length {len(key)} does not match rank {self._rank}')
            for i, (idx, dim) in enumerate(zip(key, self._shape)):
                if not (0 <= idx < dim):
                    raise ValueError(
                        f'index {idx} at axis {i} is out of range [0, {dim})')
            self._data[key] = value

    def __delitem__(self, index):
        key = self._key(index)
        if key not in self._data:
            raise KeyError(index)
        del self._data[key]

    def __contains__(self, index):
        return self._key(index) in self._data

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def __repr__(self):
        return f'SparseArray({self._data!r}, default={self._default!r})'

    def keys(self):
        return self._data.keys()

    def values(self):
        return self._data.values()

    def items(self):
        return self._data.items()

    def get(self, index, default=None):
        if default is None:
            default = self._default
        return self._data.get(self._key(index), default)

    def clear(self):
        self._data.clear()

    def copy(self):
        return SparseArray(self._shape, data=self._data.copy(), default=self._default)

    @property
    def rank(self):
        """Dimension of the indices."""
        return self._rank

    @property
    def shape(self):
        return self._shape
    
    @property
    def default(self):
        return self._default

    def nonzero_indices(self):
        """Return list of all indices with non-default values."""
        return list(self._data.keys())

    def to_dict(self):
        return dict(self._data)

    @classmethod
    def from_dict(cls, shape, d, default=0):
        result = cls(shape, default=default)
        for k, v in d.items():
            result[k] = v
        return result


class SparseTensor(SparseArray):
    """
    A sparse tensor supporting pairwise and multi-tensor contraction.

    Indices are tuples of integers. contract(other, pairs) contracts two
    tensors over specified index pairs (summing the product over those axes).
    multi_contract handles networks of arbitrarily many tensors at once.
    """

    def __repr__(self):
        return f'SparseTensor({self._data!r}, default={self._default!r})'

    def copy(self):
        return SparseTensor(self._shape, data=self._data.copy(), default=self._default)

    def _set(self, key, value):
        """Write self[key] = value, dropping the entry if it equals default."""
        if value == self._default:
            self._data.pop(key, None)
        else:
            self._data[key] = value

    def _accumulate(self, key, value):
        """Add value into self[key], maintaining sparsity."""
        if key in self._data:
            result = self._data[key] + value
            if not result:
                del self._data[key]
            else:
                self._data[key] = result
        else:
            self._data[key] = value

    def contract(self, other: 'SparseTensor', pairs):
        """
        Contract self with other over the specified index pairs, returning a
        new SparseTensor whose axes are the free axes of self followed by the
        free axes of other.

        pairs: iterable of (self_axis, other_axis) to be summed over.
               Only entries where self[..., v, ...] == other[..., v, ...] on
               every contracted axis contribute; their products are accumulated.

        Result axes: free axes of self (in order) + free axes of other (in order).

        Examples:
            # matrix multiply A[i,j] * B[j,k] -> C[i,k]
            A.contract(B, [(1, 0)])

            # double contraction A[i,j,k] * B[j,k,l] -> C[i,l]
            A.contract(B, [(1, 0), (2, 1)])
        """
        pairs = list(pairs)
        for ai, bj in pairs:
            if self._shape[ai] != other._shape[bj]:
                raise ValueError(
                    f'axis {ai} of self (size {self._shape[ai]}) is incompatible '
                    f'with axis {bj} of other (size {other._shape[bj]})')

        self_contracted  = {ai for ai, _  in pairs}
        other_contracted = {bj for _,  bj in pairs}
        self_free  = [i for i in range(self.rank)  if i not in self_contracted]
        other_free = [i for i in range(other.rank) if i not in other_contracted]
        result_shape = [self._shape[i] for i in self_free] + \
                       [other._shape[i] for i in other_free]

        if not self._data or not other._data:
            return SparseTensor(result_shape, default=self._default)

        # Group self entries by their values at the contracted axes (in pairs order).
        # For each contraction key, we only need to visit other entries that match.
        self_groups = {}
        for key, val in self.items():
            c_key = tuple(key[ai] for ai, _ in pairs)
            f_key = tuple(key[i]  for i in self_free)
            self_groups.setdefault(c_key, []).append((f_key, val))

        result = SparseTensor(result_shape, default=self._default)
        for key_b, val_b in other.items():
            c_key = tuple(key_b[bj] for _, bj in pairs)
            group = self_groups.get(c_key)
            if group is None:
                continue
            f_key_b = tuple(key_b[i] for i in other_free)
            for f_key_a, val_a in group:
                result._accumulate(f_key_a + f_key_b, val_a * val_b)
        return result
    
    def decorated_contract(self, other: 'SparseTensor', pairs):
        """
        An enhanced version of contract, where:

        pairs: a dictionary whose keys are (self_axis, other_axis) to be summed over,
                and values are tuples of the form (0/1, h[i,j]), where h[i,j] is a tensor of rank 2
                the 0 or 1 indicates whether i or j is contracted with self_axis.

        For example, let A[i,j] and B[k,l] be two tensors, then
            A.contract(B, {(1,0): (0, h)}) gives C[i,l]:= sum_{j,k} A[i,j] * h[j,k] * B[k,l]
            and
            A.contract(B, {(1,0): (1, h)}) gives C[i,l]:= sum_{j,k} A[i,j] * h[k,j] * B[k,l]
        """
        if self is other:
            return self._decorated_trace_pairs(pairs)

        for (ai, bj), (side, h) in pairs.items():
            h_self, h_other = (0, 1) if side == 0 else (1, 0)
            if self._shape[ai] != h._shape[h_self]:
                raise ValueError(
                    f'self axis {ai} (size {self._shape[ai]}) is incompatible '
                    f'with h axis {h_self} (size {h._shape[h_self]})')
            if other._shape[bj] != h._shape[h_other]:
                raise ValueError(
                    f'other axis {bj} (size {other._shape[bj]}) is incompatible '
                    f'with h axis {h_other} (size {h._shape[h_other]})')

        pair_list        = list(pairs.keys())
        self_contracted  = {ai for ai, _ in pair_list}
        other_contracted = {bj for _, bj in pair_list}
        self_free  = [i for i in range(self.rank)  if i not in self_contracted]
        other_free = [i for i in range(other.rank) if i not in other_contracted]
        result_shape = [self._shape[i] for i in self_free] + \
                       [other._shape[i] for i in other_free]

        if not self._data or not other._data:
            return SparseTensor(result_shape, default=self._default)

        # Group self entries by their contracted-axis values (in pair_list order).
        self_groups = {}
        for key, val in self._data.items():
            c_key = tuple(key[ai] for ai, _ in pair_list)
            f_key = tuple(key[i]  for i in self_free)
            self_groups.setdefault(c_key, []).append((f_key, val))

        # For each pair m, precompute: given k (other's contracted value),
        # which j values in self are reachable and with what h weight?
        # h_lookup[m][k] = [(j, h_val), ...]
        h_lookups = []
        for (ai, bj) in pair_list:
            side, h = pairs[(ai, bj)]
            lookup = {}
            for hkey, hval in h._data.items():
                hi, hj = hkey
                j, k = (hi, hj) if side == 0 else (hj, hi)
                lookup.setdefault(k, []).append((j, hval))
            h_lookups.append(lookup)

        result = SparseTensor(result_shape, default=self._default)

        for key_b, val_b in other._data.items():
            f_key_b = tuple(key_b[i] for i in other_free)

            # For each pair m, collect reachable (j_m, h_val_m) from h_m.
            per_pair = []
            for m, (ai, bj) in enumerate(pair_list):
                k = key_b[bj]
                js = h_lookups[m].get(k)
                if not js:
                    break
                per_pair.append(js)
            else:
                # Cartesian product: try every combination of j values across pairs.
                for combo in cartesian_product(*per_pair):
                    c_key = tuple(j for j, _ in combo)
                    group = self_groups.get(c_key)
                    if group is None:
                        continue
                    h_weight = 1
                    for _, hval in combo:
                        h_weight *= hval
                    hval_b = h_weight * val_b
                    for f_key_a, val_a in group:
                        result._accumulate(f_key_a + f_key_b, val_a * hval_b)

        return result

    def trace(self, i, j):
        """
        Single-tensor contraction: set indices i and j equal and sum,
        returning a SparseTensor of rank reduced by 2.

        Example: T[i,j,k].trace(0, 2) -> result[j] = sum_k T[k, j, k]
        """
        assert self.rank >= 2
        n = self.rank
        i, j = i % n, j % n
        if i == j:
            raise ValueError("trace indices must be distinct")

        result_shape = [v for idx, v in enumerate(self._shape) if idx != i and idx != j]
        if not self._data:
            return SparseTensor(result_shape, default=self._default)

        result = SparseTensor(result_shape, default=self._default)
        for key, value in self.items():
            if key[i] != key[j]:
                continue
            free_key = tuple(v for idx, v in enumerate(key) if idx != i and idx != j)
            result._accumulate(free_key, value)
        return result
    
    def decorated_trace(self, i, j, decoration):
        """
        Contract axes i and j of self with an edge tensor h inserted between
        them, returning a SparseTensor of rank reduced by 2.

        decoration: (side, h) where h is a rank-2 SparseTensor.
            side=0: result[free] = sum_{a,b} self[...,a at i,...,b at j,...] * h[a,b]
            side=1: result[free] = sum_{a,b} self[...,a at i,...,b at j,...] * h[b,a]

        Example: A[i,j,k].decorated_trace(0, 2, (0, h))
                 -> result[j] = sum_{i,k} A[i,j,k] * h[i,k]
        """
        return self._decorated_trace_pairs({(i, j): decoration})

    def _decorated_trace_pairs(self, pairs):
        """
        Multi-pair decorated trace. All axis indices refer to self's original axes.

        pairs: dict {(i, j): (side, h)} — each entry contracts axes i and j of
               self via the edge tensor h, simultaneously in a single pass.
        """
        contracted = set()
        for i, j in pairs:
            n = self.rank
            i, j = i % n, j % n
            if i == j:
                raise ValueError("trace indices must be distinct")
            for ax in (i, j):
                if ax in contracted:
                    raise ValueError(f"axis {ax} appears in more than one pair")
                contracted.add(ax)

        for (i, j), (side, h) in pairs.items():
            h_i, h_j = (0, 1) if side == 0 else (1, 0)
            if self._shape[i] != h._shape[h_i]:
                raise ValueError(
                    f'axis {i} of self (size {self._shape[i]}) is incompatible '
                    f'with h axis {h_i} (size {h._shape[h_i]})')
            if self._shape[j] != h._shape[h_j]:
                raise ValueError(
                    f'axis {j} of self (size {self._shape[j]}) is incompatible '
                    f'with h axis {h_j} (size {h._shape[h_j]})')

        free = [idx for idx in range(self.rank) if idx not in contracted]
        result_shape = [self._shape[idx] for idx in free]
        result = SparseTensor(result_shape, default=self._default)

        for key, val in self._data.items():
            weight = 1
            for (i, j), (side, h) in pairs.items():
                a, b = key[i], key[j]
                h_val = h[a, b] if side == 0 else h[b, a]
                if h_val == h._default:
                    weight = 0
                    break
                weight *= h_val
            if weight == 0:
                continue
            free_key = tuple(key[idx] for idx in free)
            result._accumulate(free_key, val * weight)

        return result

    def fixate(self, i, value):
        """
        Set index i to a fixed value, obtaining a new tensor with one less rank.

        T[i,j,k,l].fixate(i, 0) -> T[0,j,k,l]
        """
        n = self.rank
        i = i % n
        result_shape = [v for idx, v in enumerate(self._shape) if idx != i]
        result = SparseTensor(result_shape, default=self._default)
        for key, val in self._data.items():
            if key[i] != value:
                continue
            free_key = tuple(v for idx, v in enumerate(key) if idx != i)
            result._data[free_key] = val
        return result
    
    def permute(self, indices):
        """
        Reorder axes using pull-style indices: indices[i] is the axis of self
        that becomes axis i of the result.

        result[i0, i1, ...] = self[i_{indices[0]}, i_{indices[1]}, ...]

        Example: A.permute([2, 0, 1]) produces B where B[a,b,c] = A[b,c,a].
        """
        result_shape = [self._shape[i] for i in indices]
        result = SparseTensor(result_shape, default=self._default)
        for key, val in self._data.items():
            new_key = tuple(key[i] for i in indices)
            result._data[new_key] = val
        return result
