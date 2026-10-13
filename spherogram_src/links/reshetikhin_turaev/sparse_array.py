from itertools import product as cartesian_product


class SparseArray:
    """
    A sparse array supporting arbitrary-dimensional indexing via tuples.

    Internally stores only non-default entries in a dict keyed by index tuples.
    Indices can be integers or tuples of integers of any length.
    """

    __slots__ = ("_data", "_default", "_rank", "_shape")

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
        return (
            self._shape == other._shape
            and self._default == other._default
            and self._data == other._data
        )

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
                    f"key length {len(key)} does not match rank {self._rank}"
                )
            for i, (idx, dim) in enumerate(zip(key, self._shape)):
                if not (0 <= idx < dim):
                    raise ValueError(
                        f"index {idx} at axis {i} is out of range [0, {dim})"
                    )
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
        return f"SparseArray({self._data!r}, default={self._default!r})"

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
        return f"SparseTensor({self._data!r}, default={self._default!r})"

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

    def _accumulate_owned(self, key, value):
        """
        Add value into self[key], absorbing value's coefficients in place.

        Only for values this tensor owns: value must have been freshly created
        by the caller and must not be shared with any other tensor, since the
        entry it is merged into is mutated rather than rebuilt.  This avoids
        the copy of the whole accumulated polynomial that `x = x + value`
        performs on every hit, which is quadratic in the number of hits.

        Values that are not dict-based polynomials, and polynomials whose
        variables differ, fall back to ordinary addition.
        """
        data = self._data
        current = data.get(key)
        if current is None:
            data[key] = value
            return

        # Ask the type, not the instance: a missing attribute on a Sage
        # element is resolved through its category machinery and costs more
        # than the addition being guarded, while a lookup on the class is
        # about twenty times cheaper.
        cls = type(current)
        if (
            cls is not type(value)
            or not hasattr(cls, "poly_dict")
            or current.vars is not value.vars
        ):
            result = current + value
            if result:
                data[key] = result
            else:
                del data[key]
            return

        current_dict = current.poly_dict
        for exponent, coefficient in value.poly_dict.items():
            previous = current_dict.get(exponent)
            if previous is None:
                current_dict[exponent] = coefficient
            else:
                total = previous + coefficient
                if total:
                    current_dict[exponent] = total
                else:
                    del current_dict[exponent]
        if not current_dict:
            del data[key]

    def contract(self, other: "SparseTensor", pairs):
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
                    f"axis {ai} of self (size {self._shape[ai]}) is incompatible "
                    f"with axis {bj} of other (size {other._shape[bj]})"
                )

        self_contracted = {ai for ai, _ in pairs}
        other_contracted = {bj for _, bj in pairs}
        self_free = [i for i in range(self.rank) if i not in self_contracted]
        other_free = [i for i in range(other.rank) if i not in other_contracted]
        result_shape = [self._shape[i] for i in self_free] + [
            other._shape[i] for i in other_free
        ]

        if not self._data or not other._data:
            return SparseTensor(result_shape, default=self._default)

        # Group self entries by their values at the contracted axes (in pairs order).
        # For each contraction key, we only need to visit other entries that match.
        self_groups = {}
        for key, val in self.items():
            c_key = tuple(key[ai] for ai, _ in pairs)
            f_key = tuple(key[i] for i in self_free)
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

    def decorated_contract(
        self, other: "SparseTensor", pairs, consume=False, diagonal=None
    ):
        """
        An enhanced version of contract, where:

        pairs: a dictionary whose keys are (self_axis, other_axis) to be summed over,
                and values are tuples of the form (0/1, h[i,j]), where h[i,j] is a tensor of rank 2
                the 0 or 1 indicates whether i or j is contracted with self_axis.

        For example, let A[i,j] and B[k,l] be two tensors, then
            A.contract(B, {(1,0): (0, h)}) gives C[i,l]:= sum_{j,k} A[i,j] * h[j,k] * B[k,l]
            and
            A.contract(B, {(1,0): (1, h)}) gives C[i,l]:= sum_{j,k} A[i,j] * h[k,j] * B[k,l]

        consume: when True, the entries of whichever operand is streamed are
                freed as they are used, leaving that operand empty on return.
                Peak memory is then the resident operand plus the result,
                instead of both operands plus the result.  Only pass True when
                both operands are about to be discarded, and never for a
                tensor shared with anything else (an RMatrix's h tensors, say).

        diagonal: whether every h in pairs is diagonal, which allows a faster
                contraction.  Left as None it is worked out from the tensors;
                pass it when the caller already knows (see
                RMatrix.diagonal_h), since the answer is a property of the
                R matrices and does not change between contractions.  Passing
                True for a non-diagonal h silently gives wrong answers, as the
                off-diagonal terms are then never visited.
        """
        if self is other:
            return self._decorated_trace_pairs(pairs)

        for (ai, bj), (side, h) in pairs.items():
            h_self, h_other = (0, 1) if side == 0 else (1, 0)
            if self._shape[ai] != h._shape[h_self]:
                raise ValueError(
                    f"self axis {ai} (size {self._shape[ai]}) is incompatible "
                    f"with h axis {h_self} (size {h._shape[h_self]})"
                )
            if other._shape[bj] != h._shape[h_other]:
                raise ValueError(
                    f"other axis {bj} (size {other._shape[bj]}) is incompatible "
                    f"with h axis {h_other} (size {h._shape[h_other]})"
                )

        pair_list = list(pairs.keys())
        self_contracted = {ai for ai, _ in pair_list}
        other_contracted = {bj for _, bj in pair_list}
        self_free = [i for i in range(self.rank) if i not in self_contracted]
        other_free = [i for i in range(other.rank) if i not in other_contracted]
        result_shape = [self._shape[i] for i in self_free] + [
            other._shape[i] for i in other_free
        ]

        if not self._data or not other._data:
            return SparseTensor(result_shape, default=self._default)

        # Every entry of the grouped operand is revisited once per matching
        # entry of the streamed one, so the grouped operand has to stay
        # resident for the whole pass while the streamed one is touched a
        # single entry at a time.  Group whichever side is smaller.
        stream_self = len(self._data) > len(other._data)
        if stream_self:
            grouped, streamed = other, self
            grouped_axes = [bj for _, bj in pair_list]
            streamed_axes = [ai for ai, _ in pair_list]
            grouped_free, streamed_free = other_free, self_free
            # `side` says which axis of h meets self; the roles just swapped.
            grouped_sides = [1 - pairs[p][0] for p in pair_list]
        else:
            grouped, streamed = self, other
            grouped_axes = [ai for ai, _ in pair_list]
            streamed_axes = [bj for _, bj in pair_list]
            grouped_free, streamed_free = self_free, other_free
            grouped_sides = [pairs[p][0] for p in pair_list]

        # With diagonal h the two ends of a contracted edge always take the
        # same value, so the edge weight depends on that one value and can be
        # folded into the resident operand: one product per entry there,
        # instead of one per entry of the streamed operand, which is the
        # larger of the two.
        h_tensors = [pairs[p][1] for p in pair_list]
        if diagonal is None:
            diagonal = all(all(i == j for i, j in h._data) for h in h_tensors)

        # Group entries of the resident operand by their contracted-axis values.
        groups = {}
        if diagonal:
            for key, val in grouped._data.items():
                c_key = tuple(key[ax] for ax in grouped_axes)
                weight = 1
                unit_weight = True
                for m, h in enumerate(h_tensors):
                    hval = h._data.get((c_key[m], c_key[m]), h._default)
                    if type(hval) is int and hval == 1:
                        continue
                    if not hval:
                        weight = None
                        break
                    weight = hval if unit_weight else weight * hval
                    unit_weight = False
                if weight is None:
                    continue
                f_key = tuple(key[i] for i in grouped_free)
                groups.setdefault(c_key, []).append(
                    (f_key, val if unit_weight else val * weight)
                )
        else:
            for key, val in grouped._data.items():
                c_key = tuple(key[ax] for ax in grouped_axes)
                f_key = tuple(key[i] for i in grouped_free)
                groups.setdefault(c_key, []).append((f_key, val))

        # For each pair m, precompute: given k (the streamed operand's
        # contracted value), which values j of the grouped operand are
        # reachable and with what h weight?  h_lookup[m][k] = [(j, h_val), ...]
        # Not needed when h is diagonal: j is then k, and the weight has
        # already been folded into the grouped values.
        h_lookups = []
        if not diagonal:
            for m, h in enumerate(h_tensors):
                lookup = {}
                for hkey, hval in h._data.items():
                    hi, hj = hkey
                    j, k = (hi, hj) if grouped_sides[m] == 0 else (hj, hi)
                    lookup.setdefault(k, []).append((j, hval))
                h_lookups.append(lookup)

        result = SparseTensor(result_shape, default=self._default)
        if not groups:
            return result
        # Absorbing values in place only pays for dict-based polynomials; for
        # anything else (Sage polynomials, plain numbers) ordinary addition is
        # what happens anyway, so settle it once here rather than per entry.
        sample = next(iter(groups.values()))[0][1]
        accumulate = (
            result._accumulate_owned
            if hasattr(type(sample), "poly_dict")
            else result._accumulate
        )
        grouped_first = not stream_self
        streamed_data = streamed._data

        if consume:

            def streamed_entries():
                while streamed_data:
                    yield streamed_data.popitem()

        else:

            def streamed_entries():
                return iter(streamed_data.items())

        if diagonal:
            # The weight is already carried by the grouped values, so a
            # streamed entry simply meets the group with matching indices.
            for key_b, val_b in streamed_entries():
                c_key = tuple(key_b[ax] for ax in streamed_axes)
                group = groups.get(c_key)
                if group is not None:
                    f_key_b = tuple(key_b[i] for i in streamed_free)
                    if grouped_first:
                        for f_key_a, val_a in group:
                            accumulate(f_key_a + f_key_b, val_a * val_b)
                    else:
                        for f_key_a, val_a in group:
                            accumulate(f_key_b + f_key_a, val_a * val_b)
                del key_b, val_b
            return result

        for key_b, val_b in streamed_entries():
            f_key_b = tuple(key_b[i] for i in streamed_free)

            # For each pair m, collect reachable (j_m, h_val_m) from h_m.
            per_pair = []
            for m, ax in enumerate(streamed_axes):
                js = h_lookups[m].get(key_b[ax])
                if not js:
                    break
                per_pair.append(js)
            else:
                # Cartesian product: try every combination of j values across pairs.
                for combo in cartesian_product(*per_pair):
                    c_key = tuple(j for j, _ in combo)
                    group = groups.get(c_key)
                    if group is None:
                        continue
                    # Edges with no rotation carry the identity, whose entries
                    # are the plain integer 1.  Multiplying by those is not
                    # free: against a Sage polynomial it runs the full
                    # coercion machinery, and there is one such product per
                    # entry of the streamed operand.
                    h_weight = 1
                    unit_weight = True
                    for _, hval in combo:
                        if type(hval) is int and hval == 1:
                            continue
                        h_weight = hval if unit_weight else h_weight * hval
                        unit_weight = False
                    hval_b = val_b if unit_weight else h_weight * val_b
                    if grouped_first:
                        for f_key_a, val_a in group:
                            accumulate(f_key_a + f_key_b, val_a * hval_b)
                    else:
                        for f_key_a, val_a in group:
                            accumulate(f_key_b + f_key_a, val_a * hval_b)
            del key_b, val_b

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
                    f"axis {i} of self (size {self._shape[i]}) is incompatible "
                    f"with h axis {h_i} (size {h._shape[h_i]})"
                )
            if self._shape[j] != h._shape[h_j]:
                raise ValueError(
                    f"axis {j} of self (size {self._shape[j]}) is incompatible "
                    f"with h axis {h_j} (size {h._shape[h_j]})"
                )

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

    def permute(self, indices, consume=False):
        """
        Reorder axes using pull-style indices: indices[i] is the axis of self
        that becomes axis i of the result.

        result[i0, i1, ...] = self[i_{indices[0]}, i_{indices[1]}, ...]

        Example: A.permute([2, 0, 1]) produces B where B[a,b,c] = A[b,c,a].

        consume: when True, entries are moved out of self as the result is
                built, so the two never both hold the full tensor.  Self is
                left empty.
        """
        result_shape = [self._shape[i] for i in indices]
        result = SparseTensor(result_shape, default=self._default)
        if consume:
            data = self._data
            while data:
                key, val = data.popitem()
                result._data[tuple(key[i] for i in indices)] = val
        else:
            for key, val in self._data.items():
                new_key = tuple(key[i] for i in indices)
                result._data[new_key] = val
        return result
