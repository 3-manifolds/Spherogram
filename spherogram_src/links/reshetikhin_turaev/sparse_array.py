class SparseArray:
    """
    A sparse array supporting arbitrary-dimensional indexing via tuples.

    Internally stores only non-default entries in a dict keyed by index tuples.
    Indices can be integers or tuples of integers of any length.
    """
    __slots__ = ('_data', '_default', '_rank')

    def __init__(self, data= None, default=0, rank = None):
        self._data = {}
        self._default = default
        self._rank = rank
        if data is not None:
            if isinstance(data, dict):
                data = data.items()
            for k, v in data:
                self[k] = v

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
            if self._rank is None:
                self._rank = len(key)
            elif len(key) != self._rank:
                raise ValueError(
                    f'key length {len(key)} does not match rank {self._rank}')
            self._data[key] = value

    def __delitem__(self, index):
        key = self._key(index)
        if key not in self._data:
            raise KeyError(index)
        del self._data[key]
        if not self._data:
            self._rank = None

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
        self._rank = None

    def copy(self):
        return SparseArray(data=self._data.copy(),
                            default=self._default,
                            rank = self.rank)

    @property
    def rank(self):
        """Dimension of the indices."""
        return self._rank

    def nonzero_indices(self):
        """Return list of all indices with non-default values."""
        return list(self._data.keys())

    def to_dict(self):
        return dict(self._data)

    @classmethod
    def from_dict(cls, d, default=0):
        result = cls(default=default)
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
        return SparseTensor(data=self._data.copy(), 
                            default=self._default, 
                            rank=self.rank)

    @property
    def rank(self):
        """Number of indices (tensor rank)"""
        return self._rank

    def _set(self, key, value):
        """Write self[key] = value, dropping the entry if it equals default."""
        assert self.rank == len(key)

        if value == self._default:
            self._data.pop(key, None)
            if not self._data:
                self._rank = None
        else:
            self._data[key] = value

    def _accumulate(self, key, value):
        """Add value into self[key], maintaining sparsity."""
        self._set(key, self._data.get(key, self._default) + value)

    def contract(self, other: SparseTensor, pairs):
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
        if not self._data or not other._data:
            return SparseTensor(default=self._default)

        self_contracted  = {ai for ai, _  in pairs}
        other_contracted = {bj for _,  bj in pairs}

        self_free  = [i for i in range(self.rank)  if i not in self_contracted]
        other_free = [i for i in range(other.rank) if i not in other_contracted]

        # Group self entries by their values at the contracted axes (in pairs order).
        # For each contraction key, we only need to visit other entries that match.
        self_groups = {}
        for key, val in self.items():
            c_key = tuple(key[ai] for ai, _ in pairs)
            f_key = tuple(key[i]  for i in self_free)
            self_groups.setdefault(c_key, []).append((f_key, val))

        result = SparseTensor(default=self._default, 
                              rank = len(self_free) + len(other_free))
        for key_b, val_b in other.items():
            c_key = tuple(key_b[bj] for _, bj in pairs)
            group = self_groups.get(c_key)
            if group is None:
                continue
            assert len(group) > 0

            f_key_b = tuple(key_b[i] for i in other_free)
            for f_key_a, val_a in group:
                result._accumulate(f_key_a + f_key_b, val_a * val_b)
        return result

    def trace(self, i, j):
        """
        Single-tensor contraction: set indices i and j equal and sum,
        returning a SparseTensor of rank reduced by 2.

        Example: T[i,j,k].trace(0, 2) -> result[j] = sum_k T[k, j, k]
        """
        assert self.rank - 2 >= 0
        if not self._data:
            return SparseTensor(default=self._default,
                                rank = self.rank - 2)

        n = self.rank
        i, j = i % n, j % n
        if i == j:
            raise ValueError("trace indices must be distinct")

        result = SparseTensor(default=self._default,
                              rank = self.rank - 2)
        for key, value in self.items():
            if key[i] != key[j]:
                continue
            free_key = tuple(v for idx, v in enumerate(key) if idx != i and idx != j)
            result._accumulate(free_key, value)
        return result