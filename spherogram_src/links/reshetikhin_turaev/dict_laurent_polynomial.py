from ...sage_helper import _within_sage, sage_method

if _within_sage:
    from sage.all import PuiseuxSeriesRing, LaurentPolynomialRing, ZZ


@sage_method
def laurent_poly_from_dict(dict, vars, F):
    L = LaurentPolynomialRing(F, vars)
    return L(dict)


@sage_method
def puiseux_series_from_dict(poly_dict, var, F):
    """
    Build a Sage Puiseux series from a poly_dict and a single LaurentVariable.
    Key k represents var^(k / var.denominator).
    """
    P = PuiseuxSeriesRing(F, var.name)
    t = P.gen()
    result = P.zero()
    den = ZZ(var.denominator)
    for key, coef in poly_dict.items():
        (k,) = key
        result += coef * t ** (ZZ(k) / den)
    return result


class LaurentVariable:
    """
    A named variable with an optional denominator.
    Exponent key k represents var^(k / denominator).

    Examples:
        LaurentVariable('q')        # q^k for integer k
        LaurentVariable('q', 2)     # q^(k/2), so key 1 means q^(1/2)
    """

    __slots__ = ["name", "denominator"]

    def __init__(self, name, denominator=1):
        self.name = name
        self.denominator = denominator

    def __eq__(self, other):
        if isinstance(other, LaurentVariable):
            return self.name == other.name and self.denominator == other.denominator
        return NotImplemented

    def __hash__(self):
        return hash((self.name, self.denominator))

    def __repr__(self):
        if self.denominator == 1:
            return self.name
        return f"{self.name}[1/{self.denominator}]"

    def fmt_exp(self, k):
        """Format exponent k as a string for display."""
        if k == 0:
            return None
        num, den = k, self.denominator
        g = _gcd(abs(num), den)
        num, den = num // g, den // g
        if den == 1:
            if num == 1:
                return ""
            return f"^{num}"
        return f"^({num}/{den})"

    def latex_exp(self, k):
        """Format exponent k as a LaTeX superscript, or None if k == 0.

        >>> q = LaurentVariable('q', 2)
        >>> [q.latex_exp(k) for k in (0, 2, 4, -2, 1, -3)]
        [None, '', '^{2}', '^{-1}', '^{\\\\frac{1}{2}}', '^{-\\\\frac{3}{2}}']
        """
        if k == 0:
            return None
        num, den = k, self.denominator
        g = _gcd(abs(num), den)
        num, den = num // g, den // g
        if den == 1:
            if num == 1:
                return ""
            return f"^{{{num}}}"
        sign = "-" if num < 0 else ""
        return f"^{{{sign}\\frac{{{abs(num)}}}{{{den}}}}}"


def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a


def _lcm(a, b):
    return a * b // _gcd(a, b)


def _union_vars(vars1, vars2):
    """Return a variable tuple covering both vars1 and vars2.

    Variables are identified by name, so the two tuples may list a shared
    variable at different positions; a name occurring in both gets the LCM of
    the two denominators.  The result is symmetric in its arguments: the shared
    order is kept when both tuples list the same names in the same order, and
    the union is sorted by name otherwise.  Symmetry is what makes ``a + b``
    and ``b + a`` (or ``a * b`` and ``b * a``) land on the same variable tuple.

    >>> q2, t = LaurentVariable('q', 2), LaurentVariable('t')
    >>> _union_vars((t, LaurentVariable('q')), (q2, t))
    (q[1/2], t)
    """
    by_name1 = {v.name: v for v in vars1}
    by_name2 = {v.name: v for v in vars2}
    if [v.name for v in vars1] == [v.name for v in vars2]:
        names = [v.name for v in vars1]
    else:
        names = sorted(by_name1.keys() | by_name2.keys())
    union = []
    for name in names:
        v1, v2 = by_name1.get(name), by_name2.get(name)
        if v1 is None:
            union.append(v2)
        elif v2 is None:
            union.append(v1)
        else:
            union.append(LaurentVariable(name, _lcm(v1.denominator, v2.denominator)))
    return _intern_vars(tuple(union))


def _expand_dict_to(poly, target_vars):
    """Return poly's poly_dict re-expressed over target_vars.

    target_vars must contain every variable name used by poly, in any order,
    each with a denominator that is a multiple of poly's.  Exponent keys are
    permuted into the target order and scaled by new_denom / old_denom;
    variables of target_vars that poly does not use get exponent 0.
    """
    old_vars = poly.vars
    name_to_old = {v.name: i for i, v in enumerate(old_vars)}
    # slots[j] = (index in old_vars, exponent scale) or None for a new variable.
    slots = []
    for tv in target_vars:
        oi = name_to_old.get(tv.name)
        if oi is None:
            slots.append(None)
        else:
            slots.append((oi, tv.denominator // old_vars[oi].denominator))
    return {
        _intern_exponent(
            tuple(0 if s is None else key[s[0]] * s[1] for s in slots)
        ): coeff
        for key, coeff in poly.poly_dict.items()
    }


_vars_cache = {}


def _intern_vars(vars_tuple):
    return _vars_cache.setdefault(vars_tuple, vars_tuple)


# Exponent tuples are created afresh by every multiplication but take very few
# distinct values: contracting a Reshetikhin--Turaev network can hold millions
# of terms spread over a few hundred distinct exponent vectors.  Interning them
# keeps one tuple object per value, which is the largest single component of
# the memory a contraction uses.  The cache is bounded by the range of degrees
# a computation reaches; the limit is a backstop for long-lived sessions.
_exp_cache = {}

_EXP_CACHE_LIMIT = 1 << 20


def _intern_exponent(key):
    """Return the canonical tuple object equal to the exponent tuple `key`."""
    if len(_exp_cache) > _EXP_CACHE_LIMIT:
        _exp_cache.clear()
    return _exp_cache.setdefault(key, key)


class FastDictLaurentPolynomial:
    """
    A sparse Laurent polynomial in arbitrarily many variables.

    Represented as a dict mapping exponent tuples (of integers) to nonzero
    coefficients.  Each variable's denominator is encoded in its LaurentVariable,
    so exponent key k for variable v means v^(k / v.denominator).

    Warning: arithmetic operations (+, -, *, /, **) do NOT check that the two
    operands have compatible variables.  It is the caller's responsibility to
    ensure both polynomials share the same vars tuple (same names, same order,
    same denominators) before combining them. The denoiminators can be normalized
    via refactor_variables(new_vars) method.

    >>> p = DictLaurentPolynomial.from_str('q^2 - q^-1 + 3', ['q'])
    >>> p
    -q^-1 + 3 + q^2
    >>> r = DictLaurentPolynomial.from_str('q^2*t - 1', ['q', 't'])
    >>> r
    q^2*t - 1
    """

    __slots__ = ["vars", "poly_dict"]

    def __init__(self, vars, poly_dict):
        self.vars = _intern_vars(
            tuple(
                (
                    v
                    if isinstance(v, LaurentVariable)
                    else LaurentVariable(*((v,) if isinstance(v, str) else v))
                )
                for v in vars
            )
        )
        self.poly_dict = {k: v for k, v in poly_dict.items() if v != 0}

    def to_checked(self):
        """
        Return a DictLaurentPolynomial with the same variables and coefficients.

        >>> p = FastDictLaurentPolynomial.from_str('q^2 - 1', ['q'])
        >>> type(p).__name__
        'FastDictLaurentPolynomial'
        >>> c = p.to_checked()
        >>> type(c).__name__
        'DictLaurentPolynomial'
        >>> c == p
        True
        """
        return DictLaurentPolynomial._make(self.vars, self.poly_dict, _interned=True)

    @sage_method
    def to_sage(self):
        if len(self.vars) == 1:
            return puiseux_series_from_dict(self.poly_dict, self.vars[0], F=ZZ)
        elif all(var.denominator == 1 for var in self.vars):
            return laurent_poly_from_dict(
                self.poly_dict, [var.name for var in self.vars], F=ZZ
            )
        else:
            raise NotImplementedError(
                "Multi-variable Puiseux conversion to Sage is not supported."
            )

    @classmethod
    @sage_method
    def from_sage(cls, p, var_names=None):
        """
        Convert a Sage PuiseuxSeries or LaurentPolynomial to a DictLaurentPolynomial.

        var_names: optional list of variable name strings; defaults to the
        variable names from p's parent ring.

        For PuiseuxSeries the conversion relies on the internal _l (Laurent
        series) and _e (ramification index) attributes of Sage's implementation.

        sage: p = DictLaurentPolynomial.from_str('q^2 + 3 - q^-1', ['q'])
        sage: DictLaurentPolynomial.from_sage(p.to_sage()) == p
        True
        sage: r = DictLaurentPolynomial.from_str('q^2*t - 1', ['q', 't'])
        sage: DictLaurentPolynomial.from_sage(r.to_sage()) == r
        True
        sage: s = DictLaurentPolynomial.from_str('q^(1/2) + q^(-1/4)', ['q'])
        sage: DictLaurentPolynomial.from_sage(s.to_sage()) == s
        True
        """
        from sage.rings.puiseux_series_ring_element import PuiseuxSeries

        if isinstance(p, PuiseuxSeries):
            e = int(p.ramification_index())
            name = var_names[0] if var_names else str(p.variable())
            var = LaurentVariable(name, e)
            l = p.laurent_part()
            poly_dict = {
                (int(k),): int(v)
                for k, v in zip(l.exponents(), l.coefficients())
                if v != 0
            }
            return cls._make((var,), poly_dict)

        # LaurentPolynomial (univariate or multivariate, all denominators 1).
        parent = p.parent()
        names = var_names or [str(v) for v in parent.gens()]
        vars_tuple = tuple(LaurentVariable(name) for name in names)
        nvars = len(names)
        poly_dict = {}
        for exp, v in p.dict().items():
            if v == 0:
                continue
            key = (int(exp),) if nvars == 1 else tuple(int(k) for k in exp)
            poly_dict[key] = int(v)
        return cls._make(vars_tuple, poly_dict)

    @classmethod
    def _make(cls, vars, poly_dict, _interned=False):
        """Construct without cleaning — caller guarantees no zero values.
        Pass _interned=True when vars is already a canonical interned tuple."""
        obj = object.__new__(cls)
        obj.vars = (
            vars
            if _interned
            else _intern_vars(vars if isinstance(vars, tuple) else tuple(vars))
        )
        obj.poly_dict = poly_dict
        return obj

    @classmethod
    def generator(cls, vars, index=0):
        """
        Return the unit monomial for the variable at `index`.
        Exponent key 1 represents var^(1/denominator).

        >>> q = LaurentVariable('q', 2)
        >>> g = DictLaurentPolynomial.generator([q])
        >>> g
        q^(1/2)
        >>> g ** 3
        q^(3/2)
        """
        exp = tuple(1 if i == index else 0 for i in range(len(vars)))
        return cls._make(vars, {exp: 1})

    def __bool__(self):
        """
        >>> bool(DictLaurentPolynomial.from_str('q + 1', ['q']))
        True
        >>> bool(DictLaurentPolynomial.from_str('q - q', ['q']))
        False
        """
        return bool(self.poly_dict)

    def __eq__(self, other):
        """
        >>> p = DictLaurentPolynomial.from_str('q + 1', ['q'])
        >>> p == DictLaurentPolynomial.from_str('1 + q', ['q'])
        True
        >>> p == 0
        False
        >>> DictLaurentPolynomial.from_str('q - q', ['q']) == 0
        True
        >>> DictLaurentPolynomial.from_str('1', ['q']) == 1
        True
        >>> v4 = LaurentVariable('q', 4)
        >>> p4 = DictLaurentPolynomial._make((v4,), {(4,): 1})
        >>> p4 == DictLaurentPolynomial.from_str('q', ['q'])
        True
        >>> qt = DictLaurentPolynomial.from_str('q^2*t - 1', ['q', 't'])
        >>> qt == DictLaurentPolynomial.from_str('t*q^2 - 1', ['t', 'q'])
        True
        >>> qt == DictLaurentPolynomial.from_str('q^2*s - 1', ['q', 's'])
        False
        """
        if isinstance(other, FastDictLaurentPolynomial):
            if self.vars is other.vars:
                return self.poly_dict == other.poly_dict
            if {v.name for v in self.vars} != {v.name for v in other.vars}:
                return False
            common_vars = _union_vars(self.vars, other.vars)
            return _expand_dict_to(self, common_vars) == _expand_dict_to(
                other, common_vars
            )
        if other == 0:
            return not self.poly_dict
        if other == 1:
            return (
                len(self.poly_dict) == 1
                and self.poly_dict.get((0,) * len(self.vars), 0) == 1
            )

        if isinstance(other, int):
            return (
                len(self.poly_dict) == 1
                and self.poly_dict.get((0,) * len(self.vars), 0) == other
            )
        
        return NotImplemented

    __hash__ = None

    def __neg__(self):
        """
        >>> p = DictLaurentPolynomial.from_str('q + 2', ['q'])
        >>> -p
        -2 - q
        """
        return type(self)._make(
            self.vars, {k: -v for k, v in self.poly_dict.items()}, _interned=True
        )

    def __add__(self, other):
        """
        >>> p = DictLaurentPolynomial.from_str('q + 1', ['q'])
        >>> r = DictLaurentPolynomial.from_str('q^-1 - 1', ['q'])
        >>> p + r
        q^-1 + q
        >>> p + 3
        4 + q
        """
        if isinstance(other, FastDictLaurentPolynomial):
            sp, op = self.poly_dict, other.poly_dict
            # Copy the larger dict; iterate over the smaller to minimise lookups.
            if len(sp) >= len(op):
                result, iterate = dict(sp), op
            else:
                result, iterate = dict(op), sp
            for k, v in iterate.items():
                if k in result:
                    s = result[k] + v
                    if s:
                        result[k] = s
                    else:
                        del result[k]
                else:
                    result[k] = v
            return type(self)._make(self.vars, result, _interned=True)

        if other == 0:
            return type(self)._make(self.vars, dict(self.poly_dict), _interned=True)

        zero_key = (0,) * len(self.vars)
        result = dict(self.poly_dict)
        s = result.get(zero_key, 0) + other
        if s:
            result[zero_key] = s
        else:
            result.pop(zero_key, None)
        return type(self)._make(self.vars, result, _interned=True)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        """
        >>> p = DictLaurentPolynomial.from_str('q^2 + q', ['q'])
        >>> r = DictLaurentPolynomial.from_str('q', ['q'])
        >>> p - r
        q^2
        >>> p - 1
        -1 + q + q^2
        """
        if isinstance(other, FastDictLaurentPolynomial):
            result = dict(self.poly_dict)
            for k, v in other.poly_dict.items():
                if k in result:
                    s = result[k] - v
                    if s:
                        result[k] = s
                    else:
                        del result[k]
                else:
                    result[k] = -v
            return type(self)._make(self.vars, result, _interned=True)
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        """
        >>> p = DictLaurentPolynomial.from_str('q + 1', ['q'])
        >>> r = DictLaurentPolynomial.from_str('q - 1', ['q'])
        >>> p * r
        -1 + q^2
        >>> p * 3
        3 + 3*q
        """
        if isinstance(other, FastDictLaurentPolynomial):
            result = {}
            # Share one tuple object per distinct exponent vector: the tuples
            # built below are overwhelmingly duplicates of ones already alive.
            if len(_exp_cache) > _EXP_CACHE_LIMIT:
                _exp_cache.clear()
            intern = _exp_cache.setdefault
            if len(self.poly_dict) == 1:
                # monomial * polynomial: no cancellation possible, assign directly
                ((k1, v1),) = self.poly_dict.items()
                n = len(k1)
                if n == 1:
                    for k2, v2 in other.poly_dict.items():
                        k = (k1[0] + k2[0],)
                        result[intern(k, k)] = v1 * v2
                elif n == 2:
                    for k2, v2 in other.poly_dict.items():
                        k = (k1[0] + k2[0], k1[1] + k2[1])
                        result[intern(k, k)] = v1 * v2
                else:
                    for k2, v2 in other.poly_dict.items():
                        k = tuple(a + b for a, b in zip(k1, k2))
                        result[intern(k, k)] = v1 * v2
            elif len(other.poly_dict) == 1:
                # polynomial * monomial: no cancellation possible, assign directly
                ((k2, v2),) = other.poly_dict.items()
                n = len(k2)
                if n == 1:
                    for k1, v1 in self.poly_dict.items():
                        k = (k1[0] + k2[0],)
                        result[intern(k, k)] = v1 * v2
                elif n == 2:
                    for k1, v1 in self.poly_dict.items():
                        k = (k1[0] + k2[0], k1[1] + k2[1])
                        result[intern(k, k)] = v1 * v2
                else:
                    for k1, v1 in self.poly_dict.items():
                        k = tuple(a + b for a, b in zip(k1, k2))
                        result[intern(k, k)] = v1 * v2
            else:
                n = len(self.vars)
                if n == 1:
                    for k1, v1 in self.poly_dict.items():
                        for k2, v2 in other.poly_dict.items():
                            k = (k1[0] + k2[0],)
                            prod = v1 * v2
                            old = result.get(k)
                            if old is None:
                                result[intern(k, k)] = prod
                            else:
                                s = old + prod
                                if s:
                                    result[k] = s
                                else:
                                    del result[k]
                elif n == 2:
                    for k1, v1 in self.poly_dict.items():
                        for k2, v2 in other.poly_dict.items():
                            k = (k1[0] + k2[0], k1[1] + k2[1])
                            prod = v1 * v2
                            old = result.get(k)
                            if old is None:
                                result[intern(k, k)] = prod
                            else:
                                s = old + prod
                                if s:
                                    result[k] = s
                                else:
                                    del result[k]
                else:
                    for k1, v1 in self.poly_dict.items():
                        for k2, v2 in other.poly_dict.items():
                            k = tuple(a + b for a, b in zip(k1, k2))
                            prod = v1 * v2
                            old = result.get(k)
                            if old is None:
                                result[intern(k, k)] = prod
                            else:
                                s = old + prod
                                if s:
                                    result[k] = s
                                else:
                                    del result[k]
            return type(self)._make(self.vars, result, _interned=True)

        if self == 0 or other == 0:
            return type(self)._make(self.vars, {}, _interned=True)

        if other == 1:
            return type(self)._make(self.vars, dict(self.poly_dict), _interned=True)

        return type(self)._make(
            self.vars, {k: v * other for k, v in self.poly_dict.items()}, _interned=True
        )

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        """
        Divide by a monomial DictLaurentPolynomial, or by a scalar that evenly
        divides all coefficients.  Raises ValueError otherwise.

        >>> p = DictLaurentPolynomial.from_str('q^2 - 1', ['q'])
        >>> q = DictLaurentPolynomial.from_str('q', ['q'])
        >>> p / q
        -q^-1 + q
        >>> r = DictLaurentPolynomial.from_str('2*q^2 + 4*q', ['q'])
        >>> r / 2
        2*q + q^2
        >>> r / 3
        Traceback (most recent call last):
            ...
        ValueError: DictLaurentPolynomial division: scalar 3 does not divide all coefficients
        >>> p / DictLaurentPolynomial.from_str('q + 1', ['q'])
        Traceback (most recent call last):
            ...
        ValueError: DictLaurentPolynomial division: divisor must be a monomial
        """
        if isinstance(other, FastDictLaurentPolynomial):
            if len(other.poly_dict) != 1:
                raise ValueError(
                    "DictLaurentPolynomial division: divisor must be a monomial"
                )
            return self * other**-1
        if other == 0:
            raise ZeroDivisionError("DictLaurentPolynomial division by zero")
        if not all(c % other == 0 for c in self.poly_dict.values()):
            raise ValueError(
                f"DictLaurentPolynomial division: scalar {other!r} does not "
                f"divide all coefficients"
            )
        return type(self)._make(
            self.vars,
            {k: c // other for k, c in self.poly_dict.items()},
            _interned=True,
        )

    def simplify_denominator(self):
        """
        Return a simplified copy where each variable's denominator is divided
        by the GCD it shares with all exponents for that variable.

        >>> v = LaurentVariable('q', 4)
        >>> p = DictLaurentPolynomial._make((v,), {(2,): 1, (-2,): -1})
        >>> p.vars[0].denominator
        4
        >>> s = p.simplify_denominator()
        >>> s.vars[0].denominator
        2
        >>> p == s
        True
        """
        reductions = []
        for i, var in enumerate(self.vars):
            g = 0
            for key in self.poly_dict:
                g = _gcd(abs(key[i]), g)
                if g == 1:
                    break
            reductions.append(_gcd(g, var.denominator))

        new_vars = tuple(
            LaurentVariable(var.name, var.denominator // r)
            for var, r in zip(self.vars, reductions)
        )
        new_poly_dict = {
            tuple(k // r for k, r in zip(key, reductions)): coef
            for key, coef in self.poly_dict.items()
        }
        return type(self)._make(new_vars, new_poly_dict)

    def refactor_variables(self, new_vars):
        """
        Return a copy rescaled to use new_vars.

        new_vars must list the same variable names in the same order, and each
        new denominator must be a multiple of the current denominator.  Exponent
        keys are multiplied by the ratio new_denom / old_denom so the rational
        exponents (key / denom) remain unchanged.

        >>> v1 = LaurentVariable('q', 1)
        >>> p = DictLaurentPolynomial._make((v1,), {(1,): 1, (-1,): -1})
        >>> p
        -q^-1 + q
        >>> v4 = LaurentVariable('q', 4)
        >>> p4 = p.refactor_variables([v4])
        >>> p4.vars[0].denominator
        4
        >>> p4
        -q^-1 + q
        """
        new_vars = tuple(
            v if isinstance(v, LaurentVariable) else LaurentVariable(v)
            for v in new_vars
        )
        scales = []
        for old_var, new_var in zip(self.vars, new_vars):
            if old_var.name != new_var.name:
                raise ValueError(
                    f"refactor_variables: variable name mismatch: "
                    f"{old_var.name!r} vs {new_var.name!r}"
                )
            if new_var.denominator % old_var.denominator != 0:
                raise ValueError(
                    f"refactor_variables: new denominator {new_var.denominator} "
                    f"is not a multiple of old denominator {old_var.denominator} "
                    f"for variable {old_var.name!r}"
                )
            scales.append(new_var.denominator // old_var.denominator)
        new_poly_dict = {
            _intern_exponent(tuple(k * s for k, s in zip(key, scales))): coef
            for key, coef in self.poly_dict.items()
        }
        return type(self)._make(new_vars, new_poly_dict)

    def change_vars(self, rules, new_var_names=None):
        """
        Return a new DictLaurentPolynomial with variables substituted by rules.

        rules: dict mapping LaurentVariable (or variable name string) ->
        DictLaurentPolynomial (or expression string), or a list of
        DictLaurentPolynomials (one per variable, in order).  Variables absent
        from a dict-style rules are kept as themselves (identity substitution).

        String values are parsed using the source variable names as the ring,
        or new_var_names if provided.  A string value describes the image of
        the full variable (e.g. q^1); for source variables with denominator
        d > 1 the image must be a monomial.

        new_var_names: optional list of variable name strings for the target
        ring.  Source variables absent from rules are auto-mapped to themselves
        (by name); if their name is not already in new_var_names it is appended.

        >>> q_var = LaurentVariable('q')
        >>> t_var = LaurentVariable('t')
        >>> p = DictLaurentPolynomial.from_str('q^2 + q', ['q'])
        >>> t = DictLaurentPolynomial.generator([t_var])
        >>> p.change_vars({q_var: t})
        t + t^2
        >>> p.change_vars([DictLaurentPolynomial.from_str('t^-1', ['t'])])
        t^-2 + t^-1
        >>> p.change_vars({'q': 'q^-1'})
        q^-2 + q^-1
        >>> r = DictLaurentPolynomial.from_str('q + t', ['q', 't'])
        >>> r.change_vars({'q': 'q*t^2', 't': 't^-1'})
        q*t^2 + t^-1
        >>> s = DictLaurentPolynomial.from_str('q^(1/2)*t^(1/3)', ['q', 't'])
        >>> s.change_vars({'q': 'q*t^2', 't': 't^-1'})
        q^(1/2)*t^(2/3)
        >>> p.change_vars({'q': 'a^2'}, new_var_names=['a'])
        a^2 + a^4
        >>> r.change_vars({'q': 'a^2'}, new_var_names=['a'])
        a^2 + t
        >>> r.change_vars({'q': 'a^2', 't': 'b^-1'}, new_var_names=['a', 'b'])
        a^2 + b^-1
        >>> half_q = DictLaurentPolynomial.from_str('q^(1/2)', ['q'])
        >>> half_q.change_vars({'q': 'a^4'}, new_var_names=['a'])
        a^2
        """
        if isinstance(rules, list):
            rules = dict(zip(self.vars, rules))

        # Normalize string keys to LaurentVariable via name lookup.
        if any(isinstance(k, str) for k in rules):
            name_to_var = {v.name: v for v in self.vars}
            rules = {
                (name_to_var[k] if isinstance(k, str) else k): v
                for k, v in rules.items()
            }

        # When new_var_names is provided, auto-fill unspecified source vars
        # with string identity rules, extending new_var_names as needed.
        if new_var_names is not None:
            rules = dict(rules)
            new_var_names = list(new_var_names)
            for var in self.vars:
                if var not in rules:
                    if var.name not in new_var_names:
                        new_var_names.append(var.name)
                    rules[var] = var.name

        # Parse string values.
        var_names = (
            new_var_names if new_var_names is not None else [v.name for v in self.vars]
        )
        has_str_values = any(isinstance(v, str) for v in rules.values())
        if has_str_values:
            parsed_rules = {}
            for src_var, img in rules.items():
                if isinstance(img, str):
                    parsed = type(self).from_str(img, var_names)
                    d = src_var.denominator
                    if d > 1:
                        # String describes image of the full variable (src_var^1).
                        # We need the image of the unit generator (src_var^(1/d)).
                        # Only valid when image is a monomial (can be raised to 1/d power).
                        if len(parsed.poly_dict) > 1:
                            raise ValueError(
                                f"change_vars: image of {src_var!r} has denominator "
                                f"{d} > 1 but the image is not a monomial"
                            )
                        # Scale image variable denominators by d; keys unchanged.
                        new_img_vars = tuple(
                            LaurentVariable(v.name, v.denominator * d)
                            for v in parsed.vars
                        )
                        parsed = type(self)._make(new_img_vars, parsed.poly_dict)
                    parsed_rules[src_var] = parsed
                else:
                    parsed_rules[src_var] = img
            rules = parsed_rules

        # Build identity images for variables not in rules.
        full_rules = {}
        for i, var in enumerate(self.vars):
            if var in rules:
                full_rules[var] = rules[var]
            else:
                full_rules[var] = type(self).generator(self.vars, index=i)

        # Unify all image DLPs to a common vars tuple when string values were used.
        if has_str_values and full_rules:
            all_imgs = list(full_rules.values())
            ref_vars = all_imgs[0].vars

            def _lcm2(a, b):
                return a * b // _gcd(a, b)

            common_denoms = [v.denominator for v in ref_vars]
            for img in all_imgs[1:]:
                for i, v in enumerate(img.vars):
                    common_denoms[i] = _lcm2(common_denoms[i], v.denominator)
            common_vars = tuple(
                LaurentVariable(v.name, d) for v, d in zip(ref_vars, common_denoms)
            )
            full_rules = {
                src_var: img.refactor_variables(common_vars)
                for src_var, img in full_rules.items()
            }

        result = None
        for key, coef in self.poly_dict.items():
            term = None
            for var, k in zip(self.vars, key):
                if k == 0:
                    continue
                factor = full_rules[var] ** k
                term = factor if term is None else term * factor

            if term is None:
                zero_key = (0,) * len(next(iter(full_rules.values())).vars)
                term = type(self)._make(
                    next(iter(full_rules.values())).vars, {zero_key: coef}
                )
            else:
                term = term * coef

            result = term if result is None else result + term

        if result is None:
            img = next(iter(full_rules.values()))
            return type(self)._make(img.vars, {})
        return result

    @classmethod
    def from_str(cls, s, vars):
        """
        Parse a string into a DictLaurentPolynomial.

        vars: list of variable name strings, e.g. ['q'] or ['q', 't'].
        The denominator for each LaurentVariable is the LCM of all
        denominators appearing in its exponents in the string.

        >>> DictLaurentPolynomial.from_str('q^2 - q^-1 + 3', ['q'])
        -q^-1 + 3 + q^2
        >>> DictLaurentPolynomial.from_str('q^(1/2) + q^(-1/2)', ['q'])
        q^(-1/2) + q^(1/2)
        >>> DictLaurentPolynomial.from_str('1 - 1/(q*t)', ['q', 't'])
        1 - q^-1*t^-1
        >>> DictLaurentPolynomial.from_str('(q^2 - 1) / q', ['q'])
        -q^-1 + q
        >>> DictLaurentPolynomial.from_str('t1^2 + t2^-1', ['t1', 't2'])
        t1^2 + t2^-1
        >>> DictLaurentPolynomial.from_str('t^2 + t1', ['t', 't1'])
        t^2 + t1
        >>> DictLaurentPolynomial.from_str('t1^2*t10 + t1 - t10^-1', ['t1', 't10'])
        t1^2*t10 + t1 - t10^-1
        >>> DictLaurentPolynomial.from_str('(t1^2 - 1) / t1', ['t1', 't2'])
        t1 - t1^-1
        >>> DictLaurentPolynomial.from_str('3', [])
        3
        """
        if len(vars) != len(set(vars)):
            raise ValueError(f"from_str: duplicate variable names in {vars!r}")

        s = s.replace(" ", "")

        # --- tokeniser state (mutable via single-element list) ---
        pos = [0]

        def expect(ch):
            if pos[0] >= len(s) or s[pos[0]] != ch:
                got = repr(s[pos[0]]) if pos[0] < len(s) else "end of string"
                raise ValueError(
                    f"from_str: expected {ch!r}, got {got} "
                    f"at position {pos[0]} in {s!r}"
                )
            pos[0] += 1

        def parse_pos_int():
            start = pos[0]
            while pos[0] < len(s) and s[pos[0]].isdigit():
                pos[0] += 1
            if pos[0] == start:
                got = repr(s[pos[0]]) if pos[0] < len(s) else "end of string"
                raise ValueError(
                    f"from_str: expected integer at position {pos[0]} "
                    f"in {s!r}, got {got}"
                )
            return int(s[start : pos[0]])

        def parse_signed_int():
            sign = 1
            if pos[0] < len(s) and s[pos[0]] in "+-":
                if s[pos[0]] == "-":
                    sign = -1
                pos[0] += 1
            return sign * parse_pos_int()

        def parse_exponent():
            # Called after '^' has been consumed.
            if pos[0] < len(s) and s[pos[0]] == "(":
                pos[0] += 1
                num = parse_signed_int()
                den = 1
                if pos[0] < len(s) and s[pos[0]] == "/":
                    pos[0] += 1
                    den = parse_pos_int()
                    if den == 0:
                        raise ValueError(
                            f"from_str: zero denominator in exponent in {s!r}"
                        )
                expect(")")
                return num, den
            return parse_signed_int(), 1

        # --- grammar (builds AST as nested tuples) ---
        # Nodes: ('v', name, num, den)  variable^(num/den)
        #        ('c', n)               integer constant
        #        ('u-', expr)           unary negation
        #        (op, left, right)      op in '+', '-', '*', '/'

        sorted_vars = sorted(vars, key=len, reverse=True)

        def parse_expr():
            result = parse_term()
            while pos[0] < len(s) and s[pos[0]] in "+-":
                op = s[pos[0]]
                pos[0] += 1
                result = (op, result, parse_term())
            return result

        def parse_term():
            result = parse_factor()
            while pos[0] < len(s) and s[pos[0]] in "*/":
                op = s[pos[0]]
                pos[0] += 1
                result = (op, result, parse_factor())
            return result

        def parse_factor():
            sign = 1
            while pos[0] < len(s) and s[pos[0]] in "+-":
                if s[pos[0]] == "-":
                    sign = -sign
                pos[0] += 1
            result = parse_atom()
            return ("u-", result) if sign == -1 else result

        def parse_atom():
            if pos[0] >= len(s):
                raise ValueError(f"from_str: unexpected end of expression in {s!r}")
            c = s[pos[0]]

            if c == "(":
                pos[0] += 1
                result = parse_expr()
                expect(")")
                return result

            # Match a variable name (longest first to handle ambiguous prefixes).
            for var in sorted_vars:
                end = pos[0] + len(var)
                if s[pos[0] : end] == var:
                    # Require a non-identifier character to follow (avoid prefix match).
                    if end < len(s) and (s[end].isalnum() or s[end] == "_"):
                        continue
                    pos[0] = end
                    num, den = 1, 1
                    if pos[0] < len(s) and s[pos[0]] == "^":
                        pos[0] += 1
                        num, den = parse_exponent()
                    return ("v", var, num, den)

            # Must be an integer coefficient.
            if c.isdigit():
                n = parse_pos_int()
                if pos[0] < len(s) and s[pos[0]] == ".":
                    raise ValueError(
                        f"from_str: decimal numbers not supported "
                        f"at position {pos[0]} in {s!r}"
                    )
                return ("c", n)

            raise ValueError(
                f"from_str: unexpected character {c!r} "
                f"at position {pos[0]} in {s!r}; known variables: {vars!r}"
            )

        # --- parse ---

        ast = parse_expr()

        if pos[0] != len(s):
            raise ValueError(
                f"from_str: unexpected content {s[pos[0]:]!r} "
                f"at position {pos[0]} in {s!r}"
            )

        # --- scan AST for per-variable denominator LCMs ---

        def lcm(a, b):
            return a * b // _gcd(a, b)

        var_lcms = {v: 1 for v in vars}

        def collect_denoms(node):
            tag = node[0]
            if tag == "v":
                _, name, num, den = node
                if num != 0:
                    var_lcms[name] = lcm(var_lcms[name], den)
            elif tag == "u-":
                collect_denoms(node[1])
            elif tag != "c":
                collect_denoms(node[1])
                collect_denoms(node[2])

        collect_denoms(ast)

        # --- build generators in the joint variable space ---

        vars_list = tuple(LaurentVariable(v, var_lcms[v]) for v in vars)
        n = len(vars)
        generators = {
            v: cls._make(vars_list, {tuple(1 if j == i else 0 for j in range(n)): 1})
            for i, v in enumerate(vars)
        }

        # --- evaluate AST using DictLaurentPolynomial arithmetic ---

        def evaluate(node):
            tag = node[0]
            if tag == "c":
                return node[1]
            if tag == "v":
                _, name, num, den = node
                if num == 0:
                    return 1
                return generators[name] ** (num * var_lcms[name] // den)
            if tag == "u-":
                return -evaluate(node[1])
            if tag == "+":
                return evaluate(node[1]) + evaluate(node[2])
            if tag == "-":
                return evaluate(node[1]) - evaluate(node[2])
            if tag == "*":
                return evaluate(node[1]) * evaluate(node[2])
            # tag == '/'
            right_val = evaluate(node[2])
            if isinstance(right_val, int):
                if abs(right_val) != 1:
                    raise ValueError(
                        f"from_str: division by non-unit coefficient {right_val} "
                        f"in {s!r}; write the coefficient in the numerator"
                    )
                return evaluate(node[1]) * right_val
            if len(right_val.poly_dict) != 1:
                raise ValueError(f"from_str: can only divide by a monomial in {s!r}")
            return evaluate(node[1]) * (right_val**-1)

        result = evaluate(ast)

        if isinstance(result, int):
            zero_key = (0,) * n
            return cls._make(vars_list, {zero_key: result} if result != 0 else {})
        return result

    def __pow__(self, n):
        """
        >>> p = DictLaurentPolynomial.from_str('q', ['q'])
        >>> p ** 3
        q^3
        >>> p ** -2
        q^-2
        >>> p ** 0
        1
        """
        try:
            n = n.__index__()
        except (AttributeError, TypeError):
            raise ValueError(f"exponent must be an integer, got {n!r}")
        if n == 0:
            zero_key = (0,) * len(self.vars)
            return type(self)._make(self.vars, {zero_key: 1}, _interned=True)
        if n < 0:
            if len(self.poly_dict) != 1:
                raise ValueError("negative powers only supported for monomials")
            ((key, coef),) = self.poly_dict.items()
            if coef != 1 and coef != -1:
                raise ValueError(
                    f"negative powers require a ±1 leading coefficient, got {coef!r}"
                )
            inv_key = tuple(k * n for k in key)
            inv_coef = coef ** (-n)  # -n > 0, so int**int stays int; (±1)^k = (±1)^{-k}
            return type(self)._make(self.vars, {inv_key: inv_coef}, _interned=True)
        result = self
        base = self
        n -= 1
        while n:
            if n & 1:
                result = result * base
            base = base * base
            n >>= 1
        return result

    def _sorted_terms(self):
        """Return (exponent, coefficient) pairs in display order.

        Univariate polynomials are sorted by ascending degree (matching
        PuiseuxSeries); multivariate ones by descending total degree, then
        lex descending.  Degrees are compared exactly as rationals.
        """
        L = 1
        for var in self.vars:
            L = L * var.denominator // _gcd(L, var.denominator)
        scales = tuple(L // var.denominator for var in self.vars)

        def _sort_key(item):
            exp = item[0]
            total = sum(k * s for k, s in zip(exp, scales))
            if len(self.vars) == 1:
                return total  # ascending for univariate (matches PuiseuxSeries)
            return (-total, tuple(-k for k in exp))  # descending for multivariate

        return sorted(self.poly_dict.items(), key=_sort_key)

    def __repr__(self):
        if not self.poly_dict:
            return "0"
        terms = []
        for exp, coef in self._sorted_terms():
            parts = []
            for var, k in zip(self.vars, exp):
                fmt = var.fmt_exp(k)
                if fmt is not None:
                    parts.append(var.name + fmt)
            monomial = "*".join(parts)
            if not monomial:
                terms.append(str(coef))
            elif coef == 1:
                terms.append(monomial)
            elif coef == -1:
                terms.append(f"-{monomial}")
            else:
                terms.append(f"{coef}*{monomial}")
        return " + ".join(terms).replace("+ -", "- ")


class DictLaurentPolynomial(FastDictLaurentPolynomial):
    """
    A checked variant of FastDictLaurentPolynomial that normalises denominators
    to their LCM and expands to the union of variable sets before each binary
    arithmetic operation.  Variables are matched by name, so operands may list
    the same variables in different orders.

    >>> a = DictLaurentPolynomial.from_str('q^2', ['q'])
    >>> b = DictLaurentPolynomial.from_str('q^(1/2)', ['q'])
    >>> a + b
    q^(1/2) + q^2
    >>> t = DictLaurentPolynomial.from_str('t', ['t'])
    >>> a + t
    q^2 + t
    >>> a + DictLaurentPolynomial.from_str('t*q^-2', ['t', 'q'])
    q^2 + q^-2*t
    """

    def _match(self, other):
        """Return (self, other) refactored to share a common vars tuple.

        Variables are matched by name, so the two operands may list a shared
        variable at different positions.  If the variable sets differ, both are
        expanded to their union; denominators for shared variables are
        normalised to their LCM.  The common tuple is chosen symmetrically (see
        :func:`_union_vars`), so ``a op b`` and ``b op a`` produce results over
        the same variables rather than two transposed representations.

        >>> a = DictLaurentPolynomial.from_str('q^2*t + 1', ['q', 't'])
        >>> b = DictLaurentPolynomial.from_str('t^3*q - 2', ['t', 'q'])
        >>> a * b
        q^3*t^4 + q*t^3 - 2*q^2*t - 2
        >>> b * a
        q^3*t^4 + q*t^3 - 2*q^2*t - 2
        >>> a + b
        q*t^3 + q^2*t - 1
        >>> a * b == b * a
        True
        >>> a + b == b + a
        True
        """
        if self.vars is other.vars:
            return self, other

        union_vars = _union_vars(self.vars, other.vars)
        return (
            FastDictLaurentPolynomial._make(
                union_vars, _expand_dict_to(self, union_vars), _interned=True
            ),
            FastDictLaurentPolynomial._make(
                union_vars, _expand_dict_to(other, union_vars), _interned=True
            ),
        )

    def __add__(self, other):
        lhs, rhs = (
            self._match(other)
            if isinstance(other, FastDictLaurentPolynomial)
            else (self, other)
        )
        result = FastDictLaurentPolynomial.__add__(lhs, rhs)
        return DictLaurentPolynomial._make(
            result.vars, result.poly_dict, _interned=True
        )

    def __sub__(self, other):
        lhs, rhs = (
            self._match(other)
            if isinstance(other, FastDictLaurentPolynomial)
            else (self, other)
        )
        result = FastDictLaurentPolynomial.__sub__(lhs, rhs)
        return DictLaurentPolynomial._make(
            result.vars, result.poly_dict, _interned=True
        )

    def __mul__(self, other):
        lhs, rhs = (
            self._match(other)
            if isinstance(other, FastDictLaurentPolynomial)
            else (self, other)
        )
        result = FastDictLaurentPolynomial.__mul__(lhs, rhs)
        return DictLaurentPolynomial._make(
            result.vars, result.poly_dict, _interned=True
        )

    def simplify_variables(self):
        """
        Return a copy with any variable whose exponent is 0 in every term removed.

        >>> p = DictLaurentPolynomial.from_str('q^2 + 1', ['q', 't'])
        >>> p.simplify_variables()
        1 + q^2
        >>> p.simplify_variables().vars
        (q,)
        >>> DictLaurentPolynomial.from_str('3', ['q', 't']).simplify_variables()
        3
        >>> DictLaurentPolynomial.from_str('q*t + q', ['q', 't']).simplify_variables()
        q*t + q
        """
        if not self.poly_dict or not self.vars:
            return self
        keep = [
            any(key[i] != 0 for key in self.poly_dict) for i in range(len(self.vars))
        ]
        if all(keep):
            return self
        new_vars = _intern_vars(tuple(v for v, k in zip(self.vars, keep) if k))
        new_dict = {
            tuple(exp for exp, k in zip(key, keep) if k): coeff
            for key, coeff in self.poly_dict.items()
        }
        return DictLaurentPolynomial._make(new_vars, new_dict, _interned=True)

    def coefficients(self, var):
        """
        Return a dict mapping each exponent key of `var` to its coefficient,
        a DictLaurentPolynomial in the remaining variables (in their original
        order).  As elsewhere, key k means var^(k / var.denominator).  Keys
        are sorted in increasing order and only nonzero coefficients appear.

        >>> p = DictLaurentPolynomial.from_str('q^2 + 1', ['q'])
        >>> p.coefficients('q')
        {0: 1, 2: 1}
        >>> p = DictLaurentPolynomial.from_str('q^2*t + q*t^2 + 3', ['q', 't'])
        >>> p.coefficients('q')
        {0: 3, 1: t^2, 2: t}
        >>> p.coefficients('t')
        {0: 3, 1: q^2, 2: q}
        >>> p.coefficients('t')[1].vars
        (q,)
        >>> r = DictLaurentPolynomial.from_str('q^(1/2)*t - q^(-1/2) + t', ['q', 't'])
        >>> r.coefficients('q')
        {-1: -1, 0: t, 1: t}
        >>> r.coefficients(LaurentVariable('t'))
        {0: -q^(-1/2), 1: 1 + q^(1/2)}
        >>> p.coefficients('s')
        Traceback (most recent call last):
            ...
        ValueError: coefficients: variable 's' not found
        """
        if isinstance(var, str):
            names = [v.name for v in self.vars]
            if var not in names:
                raise ValueError(f"coefficients: variable {var!r} not found")
            idx = names.index(var)
        elif var in self.vars:
            idx = self.vars.index(var)
        else:
            raise ValueError(f"coefficients: variable {var!r} not found")
        rest_vars = _intern_vars(self.vars[:idx] + self.vars[idx + 1 :])
        grouped = {}
        for key, coeff in self.poly_dict.items():
            rest_key = _intern_exponent(key[:idx] + key[idx + 1 :])
            grouped.setdefault(key[idx], {})[rest_key] = coeff
        return {
            exp: DictLaurentPolynomial._make(rest_vars, grouped[exp], _interned=True)
            for exp in sorted(grouped)
        }

    def latex(self):
        """
        Return a LaTeX string for the polynomial, with terms in the same
        order as repr().

        >>> p = DictLaurentPolynomial.from_str('q^2 + 1', ['q'])
        >>> p.latex()
        '1 + q^{2}'
        >>> r = DictLaurentPolynomial.from_str('q^(1/2)*t - q^(-1/2) + t', ['q', 't'])
        >>> r.latex()
        'q^{\\\\frac{1}{2}} t + t - q^{-\\\\frac{1}{2}}'
        >>> DictLaurentPolynomial.from_str('-3*q^-2*t^3 + 2 - q', ['q', 't']).latex()
        '-q - 3 q^{-2} t^{3} + 2'
        >>> DictLaurentPolynomial.from_str('q - q', ['q']).latex()
        '0'
        """
        if not self.poly_dict:
            return "0"
        terms = []
        for exp, coef in self._sorted_terms():
            parts = []
            for var, k in zip(self.vars, exp):
                fmt = var.latex_exp(k)
                if fmt is not None:
                    parts.append(f"{var.name}{fmt}")
            monomial = " ".join(parts)
            if not monomial:
                terms.append(str(coef))
            elif coef == 1:
                terms.append(monomial)
            elif coef == -1:
                terms.append(f"-{monomial}")
            else:
                terms.append(f"{coef} {monomial}")
        return " + ".join(terms).replace("+ -", "- ")
