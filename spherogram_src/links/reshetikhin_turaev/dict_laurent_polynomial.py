from ...sage_helper import _within_sage, sage_method

if _within_sage:
    from sage.all import PuiseuxSeriesRing, LaurentPolynomialRing, ZZ

@sage_method
def laurent_poly_from_dict(dict, vars, F = ZZ):
    L = LaurentPolynomialRing(F, vars)
    return L(dict)

@sage_method
def puiseux_series_from_dict(poly_dict, var, F=ZZ):
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
    __slots__ = ['name', 'denominator']

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
        return f'{self.name}[1/{self.denominator}]'

    def fmt_exp(self, k):
        """Format exponent k as a string for display."""
        if k == 0:
            return None
        num, den = k, self.denominator
        g = _gcd(abs(num), den)
        num, den = num // g, den // g
        if den == 1:
            if num == 1:  return ''
            return f'^{num}'
        return f'^({num}/{den})'

def _gcd(a, b):
    while b:
        a, b = b, a % b
    return a

_vars_cache = {}

def _intern_vars(vars_tuple):
    return _vars_cache.setdefault(vars_tuple, vars_tuple)

class DictLaurentPolynomial:
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
    q^2 + 3 - q^-1
    >>> r = DictLaurentPolynomial.from_str('q^2*t - 1', ['q', 't'])
    >>> r
    q^2*t - 1
    """
    __slots__ = ['vars', 'poly_dict']

    def __init__(self, vars, poly_dict):
        self.vars = _intern_vars(tuple(
            v if isinstance(v, LaurentVariable) else LaurentVariable(*((v,) if isinstance(v, str) else v))
            for v in vars
        ))
        self.poly_dict = {k: v for k, v in poly_dict.items() if v != 0}

    @sage_method
    def to_sage(self):
        if len(self.vars) == 1:
            return puiseux_series_from_dict(self.poly_dict, self.vars[0])
        elif all(var.denominator == 1 for var in self.vars):
            return laurent_poly_from_dict(self.poly_dict, [var.name for var in self.vars])
        else:
            raise NotImplementedError('Multi-variable Puiseux conversion to Sage is not supported.')

    @classmethod
    @sage_method
    def from_sage(cls, p, var_names=None):
        """
        Convert a Sage PuiseuxSeries or LaurentPolynomial to a DictLaurentPolynomial.

        var_names: optional list of variable name strings; defaults to the
        variable names from p's parent ring.

        For PuiseuxSeries the conversion relies on the internal _l (Laurent
        series) and _e (ramification index) attributes of Sage's implementation.

        >>> p = DictLaurentPolynomial.from_str('q^2 + 3 - q^-1', ['q'])
        >>> DictLaurentPolynomial.from_sage(p.to_sage()) == p
        True
        >>> r = DictLaurentPolynomial.from_str('q^2*t - 1', ['q', 't'])
        >>> DictLaurentPolynomial.from_sage(r.to_sage()) == r
        True
        >>> s = DictLaurentPolynomial.from_str('q^(1/2) + q^(-1/4)', ['q'])
        >>> DictLaurentPolynomial.from_sage(s.to_sage()) == s
        True
        """
        from sage.rings.puiseux_series_ring_element import PuiseuxSeries
        
        if isinstance(p, PuiseuxSeries):
            e = int(p.ramification_index())
            name = var_names[0] if var_names else str(p.variable())
            var = LaurentVariable(name, e)
            l = p.laurent_part()
            poly_dict = {(int(k),): int(v) for k, v in zip(l.exponents(), l.coefficients()) if v != 0}
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
        obj.vars = vars if _interned else _intern_vars(vars if isinstance(vars, tuple) else tuple(vars))
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
        >>> p4 = DictLaurentPolynomial._make((v4,), {(4,): 1})   # q stored with denom 4
        >>> p4 == DictLaurentPolynomial.from_str('q', ['q'])      # same value, different denom
        True
        """
        if isinstance(other, DictLaurentPolynomial):
            if self.vars is other.vars:
                return self.poly_dict == other.poly_dict
            if len(self.vars) != len(other.vars):
                return False
            if any(v1.name != v2.name for v1, v2 in zip(self.vars, other.vars)):
                return False
            def lcm(a, b): return a * b // _gcd(a, b)
            common_vars = tuple(
                LaurentVariable(v1.name, lcm(v1.denominator, v2.denominator))
                for v1, v2 in zip(self.vars, other.vars)
            )
            return (self.refactor_variables(common_vars).poly_dict ==
                    other.refactor_variables(common_vars).poly_dict)
        if other == 0:
            return not self.poly_dict
        if other == 1:
            return len(self.poly_dict) == 1 and \
                    self.poly_dict.get((0,) * len(self.vars), 0) == 1
        return NotImplemented

    __hash__ = None

    def __neg__(self):
        """
        >>> p = DictLaurentPolynomial.from_str('q + 2', ['q'])
        >>> -p
        -q - 2
        """
        return DictLaurentPolynomial._make(
            self.vars, {k: -v for k, v in self.poly_dict.items()}, _interned=True)

    def __add__(self, other):
        """
        >>> p = DictLaurentPolynomial.from_str('q + 1', ['q'])
        >>> r = DictLaurentPolynomial.from_str('q^-1 - 1', ['q'])
        >>> p + r
        q + q^-1
        >>> p + 3
        q + 4
        """
        if isinstance(other, DictLaurentPolynomial):
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
            return DictLaurentPolynomial._make(self.vars, result, _interned=True)
    
        if other == 0:
            return DictLaurentPolynomial._make(self.vars, dict(self.poly_dict), _interned=True)

        zero_key = (0,) * len(self.vars)
        result = dict(self.poly_dict)
        s = result.get(zero_key, 0) + other
        if s:
            result[zero_key] = s
        else:
            result.pop(zero_key, None)
        return DictLaurentPolynomial._make(self.vars, result, _interned=True)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        """
        >>> p = DictLaurentPolynomial.from_str('q^2 + q', ['q'])
        >>> r = DictLaurentPolynomial.from_str('q', ['q'])
        >>> p - r
        q^2
        >>> p - 1
        q^2 + q - 1
        """
        if isinstance(other, DictLaurentPolynomial):
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
            return DictLaurentPolynomial._make(self.vars, result, _interned=True)
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        """
        >>> p = DictLaurentPolynomial.from_str('q + 1', ['q'])
        >>> r = DictLaurentPolynomial.from_str('q - 1', ['q'])
        >>> p * r
        q^2 - 1
        >>> p * 3
        3*q + 3
        """
        if isinstance(other, DictLaurentPolynomial):
            result = {}
            if len(self.poly_dict) == 1:
                # monomial * polynomial: no cancellation possible, assign directly
                (k1, v1), = self.poly_dict.items()
                n = len(k1)
                if n == 1:
                    for k2, v2 in other.poly_dict.items():
                        result[(k1[0] + k2[0],)] = v1 * v2
                elif n == 2:
                    for k2, v2 in other.poly_dict.items():
                        result[(k1[0] + k2[0], k1[1] + k2[1])] = v1 * v2
                else:
                    for k2, v2 in other.poly_dict.items():
                        result[tuple(a + b for a, b in zip(k1, k2))] = v1 * v2
            elif len(other.poly_dict) == 1:
                # polynomial * monomial: no cancellation possible, assign directly
                (k2, v2), = other.poly_dict.items()
                n = len(k2)
                if n == 1:
                    for k1, v1 in self.poly_dict.items():
                        result[(k1[0] + k2[0],)] = v1 * v2
                elif n == 2:
                    for k1, v1 in self.poly_dict.items():
                        result[(k1[0] + k2[0], k1[1] + k2[1])] = v1 * v2
                else:
                    for k1, v1 in self.poly_dict.items():
                        result[tuple(a + b for a, b in zip(k1, k2))] = v1 * v2
            else:
                n = len(self.vars)
                if n == 1:
                    for k1, v1 in self.poly_dict.items():
                        for k2, v2 in other.poly_dict.items():
                            k = (k1[0] + k2[0],)
                            prod = v1 * v2
                            if k in result:
                                s = result[k] + prod
                                if s:
                                    result[k] = s
                                else:
                                    del result[k]
                            else:
                                result[k] = prod
                elif n == 2:
                    for k1, v1 in self.poly_dict.items():
                        for k2, v2 in other.poly_dict.items():
                            k = (k1[0] + k2[0], k1[1] + k2[1])
                            prod = v1 * v2
                            if k in result:
                                s = result[k] + prod
                                if s:
                                    result[k] = s
                                else:
                                    del result[k]
                            else:
                                result[k] = prod
                else:
                    for k1, v1 in self.poly_dict.items():
                        for k2, v2 in other.poly_dict.items():
                            k = tuple(a + b for a, b in zip(k1, k2))
                            prod = v1 * v2
                            if k in result:
                                s = result[k] + prod
                                if s:
                                    result[k] = s
                                else:
                                    del result[k]
                            else:
                                result[k] = prod
            return DictLaurentPolynomial._make(self.vars, result, _interned=True)
        
        if self == 0 or other == 0:
            return DictLaurentPolynomial._make(self.vars, {}, _interned=True)

        if other == 1:
            return DictLaurentPolynomial._make(self.vars, dict(self.poly_dict), _interned=True)

        return DictLaurentPolynomial._make(
            self.vars, {k: v * other for k, v in self.poly_dict.items()}, _interned=True)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        """
        Divide by a monomial DictLaurentPolynomial, or by a scalar that evenly
        divides all coefficients.  Raises ValueError otherwise.

        >>> p = DictLaurentPolynomial.from_str('q^2 - 1', ['q'])
        >>> q = DictLaurentPolynomial.from_str('q', ['q'])
        >>> p / q
        q - q^-1
        >>> r = DictLaurentPolynomial.from_str('2*q^2 + 4*q', ['q'])
        >>> r / 2
        q^2 + 2*q
        >>> r / 3
        Traceback (most recent call last):
            ...
        ValueError: DictLaurentPolynomial division: scalar 3 does not divide all coefficients
        >>> p / DictLaurentPolynomial.from_str('q + 1', ['q'])
        Traceback (most recent call last):
            ...
        ValueError: DictLaurentPolynomial division: divisor must be a monomial
        """
        if isinstance(other, DictLaurentPolynomial):
            if len(other.poly_dict) != 1:
                raise ValueError('DictLaurentPolynomial division: divisor must be a monomial')
            return self * other ** -1
        if other == 0:
            raise ZeroDivisionError('DictLaurentPolynomial division by zero')
        if not all(c % other == 0 for c in self.poly_dict.values()):
            raise ValueError(
                f'DictLaurentPolynomial division: scalar {other!r} does not '
                f'divide all coefficients')
        return DictLaurentPolynomial._make(
            self.vars, {k: c // other for k, c in self.poly_dict.items()},
            _interned=True)

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
        return DictLaurentPolynomial._make(new_vars, new_poly_dict)

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
        q - q^-1
        >>> v4 = LaurentVariable('q', 4)
        >>> p4 = p.refactor_variables([v4])
        >>> p4.vars[0].denominator
        4
        >>> p4
        q - q^-1
        """
        new_vars = tuple(
            v if isinstance(v, LaurentVariable) else LaurentVariable(v)
            for v in new_vars
        )
        scales = []
        for old_var, new_var in zip(self.vars, new_vars):
            if old_var.name != new_var.name:
                raise ValueError(
                    f'refactor_variables: variable name mismatch: '
                    f'{old_var.name!r} vs {new_var.name!r}')
            if new_var.denominator % old_var.denominator != 0:
                raise ValueError(
                    f'refactor_variables: new denominator {new_var.denominator} '
                    f'is not a multiple of old denominator {old_var.denominator} '
                    f'for variable {old_var.name!r}')
            scales.append(new_var.denominator // old_var.denominator)
        new_poly_dict = {
            tuple(k * s for k, s in zip(key, scales)): coef
            for key, coef in self.poly_dict.items()
        }
        return DictLaurentPolynomial._make(new_vars, new_poly_dict)

    def change_vars(self, rules, new_var_names = None):
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
        t^2 + t
        >>> p.change_vars([DictLaurentPolynomial.from_str('t^-1', ['t'])])
        t^-1 + t^-2
        >>> p.change_vars({'q': 'q^-1'})
        q^-1 + q^-2
        >>> r = DictLaurentPolynomial.from_str('q + t', ['q', 't'])
        >>> r.change_vars({'q': 'q*t^2', 't': 't^-1'})
        q*t^2 + t^-1
        >>> s = DictLaurentPolynomial.from_str('q^(1/2)*t^(1/3)', ['q', 't'])
        >>> s.change_vars({'q': 'q*t^2', 't': 't^-1'})
        q^(1/2)*t^(2/3)
        >>> p.change_vars({'q': 'a^2'}, new_var_names=['a'])
        a^4 + a^2
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
        var_names = new_var_names if new_var_names is not None else [v.name for v in self.vars]
        has_str_values = any(isinstance(v, str) for v in rules.values())
        if has_str_values:
            parsed_rules = {}
            for src_var, img in rules.items():
                if isinstance(img, str):
                    parsed = DictLaurentPolynomial.from_str(img, var_names)
                    d = src_var.denominator
                    if d > 1:
                        # String describes image of the full variable (src_var^1).
                        # We need the image of the unit generator (src_var^(1/d)).
                        # Only valid when image is a monomial (can be raised to 1/d power).
                        if len(parsed.poly_dict) > 1:
                            raise ValueError(
                                f'change_vars: image of {src_var!r} has denominator '
                                f'{d} > 1 but the image is not a monomial')
                        # Scale image variable denominators by d; keys unchanged.
                        new_img_vars = tuple(
                            LaurentVariable(v.name, v.denominator * d)
                            for v in parsed.vars
                        )
                        parsed = DictLaurentPolynomial._make(new_img_vars, parsed.poly_dict)
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
                full_rules[var] = DictLaurentPolynomial.generator(self.vars, index=i)

        # Unify all image DLPs to a common vars tuple when string values were used.
        if has_str_values and full_rules:
            all_imgs = list(full_rules.values())
            ref_vars = all_imgs[0].vars
            def _lcm2(a, b): return a * b // _gcd(a, b)
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
                term = DictLaurentPolynomial._make(
                    next(iter(full_rules.values())).vars, {zero_key: coef})
            else:
                term = term * coef

            result = term if result is None else result + term

        if result is None:
            img = next(iter(full_rules.values()))
            return DictLaurentPolynomial._make(img.vars, {})
        return result

    @classmethod
    def from_str(cls, s, vars):
        """
        Parse a string into a DictLaurentPolynomial.

        vars: list of variable name strings, e.g. ['q'] or ['q', 't'].
        The denominator for each LaurentVariable is the LCM of all
        denominators appearing in its exponents in the string.

        >>> DictLaurentPolynomial.from_str('q^2 - q^-1 + 3', ['q'])
        q^2 + 3 - q^-1
        >>> DictLaurentPolynomial.from_str('q^(1/2) + q^(-1/2)', ['q'])
        q^(1/2) + q^(-1/2)
        >>> DictLaurentPolynomial.from_str('1 - 1/(q*t)', ['q', 't'])
        1 - q^-1*t^-1
        >>> DictLaurentPolynomial.from_str('(q^2 - 1) / q', ['q'])
        q - q^-1
        >>> DictLaurentPolynomial.from_str('t1^2 + t2^-1', ['t1', 't2'])
        t1^2 + t2^-1
        >>> DictLaurentPolynomial.from_str('t^2 + t1', ['t', 't1'])
        t^2 + t1
        >>> DictLaurentPolynomial.from_str('t1^2*t10 + t1 - t10^-1', ['t1', 't10'])
        t1^2*t10 + t1 - t10^-1
        >>> DictLaurentPolynomial.from_str('(t1^2 - 1) / t1', ['t1', 't2'])
        t1 - t1^-1
        """
        if len(vars) != len(set(vars)):
            raise ValueError(f'from_str: duplicate variable names in {vars!r}')

        s = s.replace(' ', '')

        # --- tokeniser state (mutable via single-element list) ---
        pos = [0]

        def expect(ch):
            if pos[0] >= len(s) or s[pos[0]] != ch:
                got = repr(s[pos[0]]) if pos[0] < len(s) else 'end of string'
                raise ValueError(
                    f'from_str: expected {ch!r}, got {got} '
                    f'at position {pos[0]} in {s!r}')
            pos[0] += 1

        def parse_pos_int():
            start = pos[0]
            while pos[0] < len(s) and s[pos[0]].isdigit():
                pos[0] += 1
            if pos[0] == start:
                got = repr(s[pos[0]]) if pos[0] < len(s) else 'end of string'
                raise ValueError(
                    f'from_str: expected integer at position {pos[0]} '
                    f'in {s!r}, got {got}')
            return int(s[start:pos[0]])

        def parse_signed_int():
            sign = 1
            if pos[0] < len(s) and s[pos[0]] in '+-':
                if s[pos[0]] == '-':
                    sign = -1
                pos[0] += 1
            return sign * parse_pos_int()

        def parse_exponent():
            # Called after '^' has been consumed.
            if pos[0] < len(s) and s[pos[0]] == '(':
                pos[0] += 1
                num = parse_signed_int()
                den = 1
                if pos[0] < len(s) and s[pos[0]] == '/':
                    pos[0] += 1
                    den = parse_pos_int()
                    if den == 0:
                        raise ValueError(
                            f'from_str: zero denominator in exponent in {s!r}')
                expect(')')
                return num, den
            return parse_signed_int(), 1

        # --- polynomial representation ---
        # Each poly is a list of (coef: int, exps: dict[var_name → (num, den)]).
        # Exponents use exact rational arithmetic; variable names are strings.

        sorted_vars = sorted(vars, key=len, reverse=True)

        def _rat_add(n1, d1, n2, d2):
            n = n1 * d2 + n2 * d1
            d = d1 * d2
            if n == 0:
                return 0, 1
            g = _gcd(abs(n), d)
            return n // g, d // g

        def _mono_key(exps):
            return tuple(sorted(exps.items()))

        def _combine(terms):
            acc = {}
            for c, e in terms:
                k = _mono_key(e)
                if k in acc:
                    acc[k] = (acc[k][0] + c, acc[k][1])
                else:
                    acc[k] = (c, e)
            return [(c, e) for c, e in acc.values() if c != 0]

        def poly_neg(p):
            return [(-c, e) for c, e in p]

        def poly_add(p1, p2):
            return _combine(p1 + p2)

        def poly_mul(p1, p2):
            result = []
            for c1, e1 in p1:
                for c2, e2 in p2:
                    c = c1 * c2
                    e = dict(e1)
                    for v, (n2, d2) in e2.items():
                        if v in e:
                            n1, d1 = e[v]
                            rn, rd = _rat_add(n1, d1, n2, d2)
                            if rn == 0:
                                del e[v]
                            else:
                                e[v] = (rn, rd)
                        else:
                            e[v] = (n2, d2)
                    result.append((c, e))
            return _combine(result)

        def poly_inv(p):
            if len(p) != 1:
                raise ValueError(
                    f'from_str: can only divide by a monomial in {s!r}')
            c, e = p[0]
            if c == 0:
                raise ValueError(f'from_str: division by zero in {s!r}')
            if abs(c) != 1:
                raise ValueError(
                    f'from_str: division by non-unit coefficient {c} '
                    f'in {s!r}; write the coefficient in the numerator')
            return [(c, {v: (-n, d) for v, (n, d) in e.items()})]

        # --- grammar ---

        def parse_expr():
            result = parse_term()
            while pos[0] < len(s) and s[pos[0]] in '+-':
                op = s[pos[0]]; pos[0] += 1
                right = parse_term()
                result = poly_add(result, poly_neg(right) if op == '-' else right)
            return result

        def parse_term():
            result = parse_factor()
            while pos[0] < len(s) and s[pos[0]] in '*/':
                op = s[pos[0]]; pos[0] += 1
                right = parse_factor()
                result = poly_mul(result, poly_inv(right) if op == '/' else right)
            return result

        def parse_factor():
            sign = 1
            while pos[0] < len(s) and s[pos[0]] in '+-':
                if s[pos[0]] == '-':
                    sign = -sign
                pos[0] += 1
            result = parse_atom()
            return poly_neg(result) if sign == -1 else result

        def parse_atom():
            if pos[0] >= len(s):
                raise ValueError(
                    f'from_str: unexpected end of expression in {s!r}')
            c = s[pos[0]]

            if c == '(':
                pos[0] += 1
                result = parse_expr()
                expect(')')
                return result

            # Match a variable name (longest first to handle ambiguous prefixes).
            for var in sorted_vars:
                end = pos[0] + len(var)
                if s[pos[0]:end] == var:
                    # Require a non-identifier character to follow (avoid prefix match).
                    if end < len(s) and (s[end].isalnum() or s[end] == '_'):
                        continue
                    pos[0] = end
                    num, den = 1, 1
                    if pos[0] < len(s) and s[pos[0]] == '^':
                        pos[0] += 1
                        num, den = parse_exponent()
                    return [(1, {var: (num, den)} if num != 0 else {})]

            # Must be an integer coefficient.
            if c.isdigit():
                n = parse_pos_int()
                if pos[0] < len(s) and s[pos[0]] == '.':
                    raise ValueError(
                        f'from_str: decimal numbers not supported '
                        f'at position {pos[0]} in {s!r}')
                return [(n, {})]

            raise ValueError(
                f'from_str: unexpected character {c!r} '
                f'at position {pos[0]} in {s!r}; known variables: {vars!r}')

        # --- parse and convert ---

        poly = parse_expr()

        if pos[0] != len(s):
            raise ValueError(
                f'from_str: unexpected content {s[pos[0]:]!r} '
                f'at position {pos[0]} in {s!r}')

        def lcm(a, b):
            return a * b // _gcd(a, b)

        var_lcms = {v: 1 for v in vars}
        for _, exps in poly:
            for v, (n, d) in exps.items():
                var_lcms[v] = lcm(var_lcms[v], d)

        vars_list = [LaurentVariable(v, var_lcms[v]) for v in vars]

        poly_dict = {}
        for coef, exps in poly:
            key = tuple(
                exps[v][0] * (var_lcms[v] // exps[v][1]) if v in exps else 0
                for v in vars
            )
            if key in poly_dict:
                new_v = poly_dict[key] + coef
                if new_v:
                    poly_dict[key] = new_v
                else:
                    del poly_dict[key]
            else:
                poly_dict[key] = coef

        return cls._make(vars_list, poly_dict)

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
            raise ValueError(f'exponent must be an integer, got {n!r}')
        if n == 0:
            zero_key = (0,) * len(self.vars)
            return DictLaurentPolynomial._make(self.vars, {zero_key: 1}, _interned=True)
        if n < 0:
            if len(self.poly_dict) != 1:
                raise ValueError('negative powers only supported for monomials')
            (key, coef), = self.poly_dict.items()
            if coef != 1 and coef != -1:
                raise ValueError(
                    f'negative powers require a ±1 leading coefficient, got {coef!r}')
            inv_key = tuple(k * n for k in key)
            inv_coef = coef ** (-n)  # -n > 0, so int**int stays int; (±1)^k = (±1)^{-k}
            return DictLaurentPolynomial._make(self.vars, {inv_key: inv_coef}, _interned=True)
        result = self
        base = self
        n -= 1
        while n:
            if n & 1:
                result = result * base
            base = base * base
            n >>= 1
        return result

    def __repr__(self):
        if not self.poly_dict:
            return '0'
        # Sort by total degree descending (exact rational arithmetic), then lex descending.
        L = 1
        for var in self.vars:
            L = L * var.denominator // _gcd(L, var.denominator)
        scales = tuple(L // var.denominator for var in self.vars)
        def _sort_key(item):
            exp = item[0]
            total = sum(k * s for k, s in zip(exp, scales))
            return (-total, tuple(-k for k in exp))
        terms = []
        for exp, coef in sorted(self.poly_dict.items(), key=_sort_key):
            parts = []
            for var, k in zip(self.vars, exp):
                fmt = var.fmt_exp(k)
                if fmt is not None:
                    parts.append(var.name + fmt)
            monomial = '*'.join(parts)
            if not monomial:
                terms.append(str(coef))
            elif coef == 1:
                terms.append(monomial)
            elif coef == -1:
                terms.append(f'-{monomial}')
            else:
                terms.append(f'{coef}*{monomial}')
        return ' + '.join(terms).replace('+ -', '- ')

