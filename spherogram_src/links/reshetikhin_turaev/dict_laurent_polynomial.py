from ...sage_helper import _within_sage, sage_method

if _within_sage:
    from sage.all import LaurentPolynomialRing, ZZ

@sage_method
def laurent_poly_from_dict(dict, vars, F = ZZ):
    L = LaurentPolynomialRing(F, vars)
    return L(dict)

import re

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
        if all(var.denominator == 1 for var in self.vars):
            return laurent_poly_from_dict(self.poly_dict, [var.name for var in self.vars])
        else:
            raise NotImplementedError('Can only convert DictLaurentPolynomial with integer exponentials to LaurentPolynomial in Sage.')

    @classmethod
    def _make(cls, vars, poly_dict):
        """Construct without cleaning — caller guarantees no zero values."""
        obj = object.__new__(cls)
        obj.vars = _intern_vars(vars if isinstance(vars, tuple) else tuple(vars))
        obj.poly_dict = poly_dict
        return obj

    @classmethod
    def generator(cls, vars, index=0):
        """
        Return the unit monomial for the variable at `index`.
        Exponent key 1 represents var^(1/denominator).

        Example:
            q = LaurentVariable('q', 2)
            gen = DictLaurentPolynomial.generator([q])
            # gen represents q^(1/2); gen**3 represents q^(3/2)
        """
        exp = tuple(1 if i == index else 0 for i in range(len(vars)))
        return cls._make(vars, {exp: 1})

    def __bool__(self):
        return bool(self.poly_dict)

    def __eq__(self, other):
        if isinstance(other, DictLaurentPolynomial):
            return self.poly_dict == other.poly_dict
        if other == 0:
            return not self.poly_dict
        return NotImplemented

    __hash__ = None

    def __neg__(self):
        return DictLaurentPolynomial._make(
            self.vars, {k: -v for k, v in self.poly_dict.items()})

    def __add__(self, other):
        if isinstance(other, DictLaurentPolynomial):
            result = dict(self.poly_dict)
            for k, v in other.poly_dict.items():
                if k in result:
                    s = result[k] + v
                    if s:
                        result[k] = s
                    else:
                        del result[k]
                else:
                    result[k] = v
            return DictLaurentPolynomial._make(self.vars, result)
        if other == 0:
            return DictLaurentPolynomial._make(self.vars, dict(self.poly_dict))
        zero_key = (0,) * len(self.vars)
        result = dict(self.poly_dict)
        s = result.get(zero_key, 0) + other
        if s:
            result[zero_key] = s
        else:
            result.pop(zero_key, None)
        return DictLaurentPolynomial._make(self.vars, result)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
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
            return DictLaurentPolynomial._make(self.vars, result)
        return self.__add__(-other)

    def __rsub__(self, other):
        return (-self).__add__(other)

    def __mul__(self, other):
        if isinstance(other, DictLaurentPolynomial):
            result = {}
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
            return DictLaurentPolynomial._make(self.vars, result)
        if other == 0:
            return DictLaurentPolynomial._make(self.vars, {})
        if other == 1:
            return DictLaurentPolynomial._make(self.vars, dict(self.poly_dict))
        return DictLaurentPolynomial._make(
            self.vars, {k: v * other for k, v in self.poly_dict.items()})

    def __rmul__(self, other):
        return self.__mul__(other)

    def change_vars(self, rules):
        """
        Return a new DictLaurentPolynomial with variables substituted by rules.

        rules must have the same length as self.vars.  Exponent keys are
        rescaled when denominators change: exponent k (meaning var^(k/old_denom))
        rules: dict mapping LaurentVariable -> DictLaurentPolynomial,
        or a list of DictLaurentPolynomials (one per variable, in order).

        Each variable is substituted by the corresponding polynomial.
        Variables absent from a dict-style rules are kept as themselves
        (identity substitution).

        Examples:
            q = LaurentVariable('q', 2)
            t = LaurentVariable('t', 1)
            t_half = DictLaurentPolynomial.generator([t])
            p.change_vars({q: t_half})   # substitute q^(1/2) -> t
        """
        if isinstance(rules, list):
            rules = dict(zip(self.vars, rules))

        # Build identity images for variables not in rules.
        full_rules = {}
        for i, var in enumerate(self.vars):
            if var in rules:
                full_rules[var] = rules[var]
            else:
                full_rules[var] = DictLaurentPolynomial.generator(self.vars, index=i)

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

    @staticmethod
    def _preprocess_division(s, sym_names):
        """Expand /var and /(var1*var2*...) into *var^(-1)*... so the main parser handles it."""
        sorted_syms = sorted(sym_names, key=len, reverse=True)
        sym_alt = '|'.join(re.escape(sym) for sym in sorted_syms)
        exp_pat = r'(?:\^(?:\([+-]?\d+(?:/\d+)?\)|[+-]?\d+))?'
        factor_pat = r'(?:' + sym_alt + r')' + exp_pat
        product_pat = factor_pat + r'(?:\*' + factor_pat + r')*'
        denom_re = re.compile(r'/\((' + product_pat + r')\)|/(' + factor_pat + r')')
        factor_re = re.compile(r'(' + sym_alt + r')(' + exp_pat + r')')

        def negate_exp(exp):
            if not exp:
                return '^(-1)'
            inner = exp[1:]  # strip '^'
            if inner.startswith('(') and inner.endswith(')'):
                inner = inner[1:-1]
            if inner.startswith('-'):
                inner = inner[1:]
                return f'^({inner})' if '/' in inner else f'^{inner}'
            return f'^(-{inner})'

        def replace_denom(m):
            content = m.group(1) if m.group(1) is not None else m.group(2)
            result = []
            for part in content.split('*'):
                fm = factor_re.fullmatch(part)
                if fm is None:
                    raise ValueError(f'Cannot parse denominator factor: {part!r}')
                result.append(fm.group(1) + negate_exp(fm.group(2)))
            return '*' + '*'.join(result)

        # Strip outer parens from pure-product numerators: (a*b)/c → a*b/c.
        # A "pure product" means no + or - at the top level inside the parens.
        expanded = []
        i, n = 0, len(s)
        while i < n:
            if s[i] != '(':
                expanded.append(s[i]); i += 1; continue
            depth, j, pure = 1, i + 1, True
            while j < n and depth > 0:
                if s[j] == '(':   depth += 1
                elif s[j] == ')': depth -= 1
                elif s[j] in '+-' and depth == 1: pure = False
                j += 1
            if pure and j < n and s[j] == '/':
                expanded.append(s[i+1:j-1])  # content without outer parens
            else:
                expanded.append(s[i:j])
            i = j
        s = ''.join(expanded)

        return denom_re.sub(replace_denom, s)

    @classmethod
    def from_str(cls, s, vars):
        """
        Parse a string into a DictLaurentPolynomial.

        vars: list of variable name strings, e.g. ['q'] or ['q', 't'].
        The denominator for each LaurentVariable is the LCM of all
        denominators appearing in its exponents in the string.

        Supported formats:
            'q^(1/2) + q^(-3/2) - 2'
            '2*q^(1/2) - q^(-1) + 3*q^2'
            'q^2*t^(-1/3) + 1'
            '1 - 1/(q*t)'
            '(q - 1) / q'
        """
        def lcm(a, b):
            return a * b // _gcd(a, b)
        def lcm_list(lst):
            r = 1
            for x in lst:
                r = lcm(r, x)
            return r

        s = s.replace(' ', '')
        s = cls._preprocess_division(s, vars)

        # Split at + or - that are not inside parentheses.
        term_strs = []
        depth = 0
        start = 0
        for i, c in enumerate(s):
            if c == '(':
                depth += 1
            elif c == ')':
                depth -= 1
            elif c in '+-' and depth == 0 and i > start and s[i-1] != '^':
                term_strs.append(s[start:i])
                start = i
        term_strs.append(s[start:])
        term_strs = [t for t in term_strs if t]

        # Pattern for each variable: sym optionally followed by ^(num/den) or ^num.
        sym_pats = {
            sym: re.compile(
                re.escape(sym) +
                r'(?:\^(?:\(([+-]?\d+)(?:/(\d+))?\)|([+-]?\d+)))?'
            )
            for sym in vars
        }

        parsed = []               # [(coef, {sym: (num, den)})]
        all_denoms = {sym: {1} for sym in vars}

        for ts in term_strs:
            # Extract leading coefficient (handles: 2*, -3*, -, +, 2, -2).
            m = re.match(r'^([+-]?\d*)\*?', ts)
            prefix = m.group(1)
            rest = ts[m.end():]
            coef = 1 if prefix in ('', '+') else (-1 if prefix == '-' else int(prefix))

            var_exps = {}
            for sym in vars:
                pm = sym_pats[sym].search(rest)
                if pm:
                    if pm.group(1) is not None:
                        num = int(pm.group(1))
                        den = int(pm.group(2)) if pm.group(2) else 1
                    elif pm.group(3) is not None:
                        num, den = int(pm.group(3)), 1
                    else:
                        num, den = 1, 1
                    var_exps[sym] = (num, den)
                    all_denoms[sym].add(den)

            parsed.append((coef, var_exps))

        var_lcms = {sym: lcm_list(all_denoms[sym]) for sym in vars}
        vars_list = [LaurentVariable(sym, var_lcms[sym]) for sym in vars]

        poly_dict = {}
        for coef, var_exps in parsed:
            key = tuple(
                var_exps[sym][0] * (var_lcms[sym] // var_exps[sym][1])
                if sym in var_exps else 0
                for sym in vars
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
        if not isinstance(n, int):
            raise ValueError(f'exponent must be an integer, got {n!r}')
        if n == 0:
            zero_key = (0,) * len(self.vars)
            return DictLaurentPolynomial._make(self.vars, {zero_key: 1})
        if n < 0:
            if len(self.poly_dict) != 1:
                raise ValueError('negative powers only supported for monomials')
            (key, coef), = self.poly_dict.items()
            inv_key = tuple(k * n for k in key)
            inv_coef = coef ** n  # works when coef is ±1 or a symbolic type
            return DictLaurentPolynomial._make(self.vars, {inv_key: inv_coef})
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
            monomial = ''.join(parts)
            if not monomial:
                terms.append(str(coef))
            elif coef == 1:
                terms.append(monomial)
            elif coef == -1:
                terms.append(f'-{monomial}')
            else:
                terms.append(f'{coef}*{monomial}')
        return ' + '.join(terms).replace('+ -', '- ')

