from .dict_laurent_polynomial import FastDictLaurentPolynomial, LaurentVariable

from .sparse_array import SparseTensor

import csv, ast, pathlib, os, math

dir_path = pathlib.Path(__file__).resolve().parent

_cache = dict()

def laurent_sparse_tensor_from_file(file, vars = ['t', 'q'], sage_polynomials = False):
    reader = csv.reader(file)
    header = next(reader)
    shape = ast.literal_eval(header[0])

    if header[1] != 'ZZ':
        raise NotImplementedError

    data = dict()
    for line in reader:
        key, value = line
        key = tuple(ast.literal_eval(key))
        assert key not in data.keys(), f'{key} appeared multiple times in {file.name}'
        data[key] = FastDictLaurentPolynomial.from_str(value, vars = vars)

    # Unify variable denominators: compute LCM across all loaded polynomials so
    # every value shares the same vars tuple (enabling interning and consistent arithmetic).
    if data:
        common_denoms = [1] * len(vars)
        for poly in data.values():
            for i, var in enumerate(poly.vars):
                d = common_denoms[i]
                common_denoms[i] = d * var.denominator // math.gcd(d, var.denominator)
        common_vars = tuple(LaurentVariable(v, d) for v, d in zip(vars, common_denoms))
        data = {key: poly.refactor_variables(common_vars) for key, poly in data.items()}

    if sage_polynomials:
        data = {key: poly.to_sage() for key, poly in data.items()}

    return SparseTensor(shape = shape, data = data)

def laurent_sparse_tensor_from_path(path, vars = ['t', 'q'], compressed = False, sage_polynomials = False):
    if compressed:
        import bz2
        with bz2.open(path, 'rt') as f:
            return laurent_sparse_tensor_from_file(f, vars = vars, sage_polynomials = sage_polynomials)
    else:
        with open(path, 'r') as f:
            return laurent_sparse_tensor_from_file(f, vars = vars, sage_polynomials = sage_polynomials)

class RMatrix:
    __slots__ = ['_R', '_h', '_id']

    def __init__(self, Rp, Rm, hp, hm):
        self._R = (Rp, Rm)
        self._h = (hp, hm)

        self._id = SparseTensor(hp.shape, data = {(i,i): 1 for i in range(hp.shape[0])})

    def R(self, sign):
        if sign == 1:
            return self._R[0].copy()
        else:
            assert sign == -1
            return self._R[1].copy()
        
    def h(self, sign):
        if sign == 1:
            return self._h[0].copy()
        elif sign == -1:
            return self._h[1].copy()
        else:
            assert sign == 0
            return self._id
    
    @staticmethod
    def from_directory(dir_path, vars = ['t', 'q'], compressed = False, sage_polynomials = False):
        names = [name + '.csv' + ('.bz2' if compressed else '') 
                 for name in ['Rp', 'Rn', 'hp', 'hn']]
        
        tensors = [laurent_sparse_tensor_from_path(os.path.join(dir_path, name),
                                                   vars = vars,
                                                   compressed = compressed,
                                                   sage_polynomials = sage_polynomials)
                   for name in names]
        
        return RMatrix(*tensors)

def colored_links_gould_R_matrices(n, sage_polynomials = False):
    if 0 < n <= 4:
        key = (f'V{n}', sage_polynomials)
        if key in _cache.keys():
            return _cache[key]
        else:
            _cache[key] = RMatrix.from_directory(dir_path = os.path.join(dir_path, f'R_matrices/V{n}/'),
                                                           vars = ['t', 'q'], 
                                                           compressed = False,
                                                           sage_polynomials = sage_polynomials)
            return _cache[key]
    else:
        raise NotImplementedError
    
def _q_binomial(n, k, q):
    if k < 0 or k > n:
        return 0
    # table[i][j] = q_binomial(i, j, q)
    table = [[0] * (k + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        table[i][0] = 1
    for i in range(1, n + 1):
        for j in range(1, min(i, k) + 1):
            table[i][j] = table[i-1][j-1] + q**j * table[i-1][j]
    return table[n][k]

def _q_pochhammer(a, q, n):
    result = 1
    for k in range(n):
        result = result * (1 - a * q**k)
    return result

def _q_pow(e):
    """
    DictLaurentPolynomial representing q^(e/2). e must be an integer.
    """
    return FastDictLaurentPolynomial._make((LaurentVariable('q', 4),), {(e,): 1})

def colored_jones_R_matrices(n, sage_polynomials=False):
    """
    The R matrices for the n-colored Jones polynomial.
    In particular, n = 1 gives the Jones polynomial.
    """
    if n < 0:
        raise NotImplementedError

    n = n + 1

    q_actual = _q_pow(4)   # q^1
    q_inv    = _q_pow(-4)  # q^(-1)

    shape  = (n, n, n, n)
    data_p = {}
    data_n = {}

    for i in range(n):
        for j in range(n):
            for k in range(n):
                l = i + j - k
                if l < 0 or l >= n:
                    continue

                # JRRp[i,j,k,l] = JRp[i, j, m_p, n]  with m_p = j - k
                # = q^(-(n-1)^2/4) * q^(-(i+j-k)*k) * q^((n-1)*(i+k)/2)
                #   * qbin[j, m_p] * qp[q^(n-1-i), q^-1, m_p]
                m_p  = j - k
                e4_p = -(n-1)**2 - 4*(i+j-k)*k + 2*(n-1)*(i+k)
                mono = _q_pow(e4_p)
                qb   = _q_binomial(j, m_p, q_actual)    # 0 when m_p < 0 or m_p > j
                qp   = _q_pochhammer(_q_pow(4*(n-1-i)), q_inv, m_p)
                val  = mono * qb * qp
                if val:
                    if not sage_polynomials:
                        data_p[(i, j, k, l)] = val
                    else:
                        data_p[(i, j, k, l)] = val.to_sage()

                # JRRn[i,j,k,l] = JRn[i, j, m_n, n]  with m_n = k - j
                # = q^((n-1)^2/4) * (-1)^m_n * q^(i*j + m_n*(m_n-1)/2) * q^(-(n-1)*(i+k)/2)
                #   * qbin[i, m_n] * qp[q^(n-1-j), q^-1, m_n]
                m_n  = k - j
                # m_n*(m_n-1) is always even (product of consecutive integers)
                e4_n = (n-1)**2 + 4*(i*j + m_n*(m_n-1)//2) - 2*(n-1)*(i+k)
                sign = (-1)**m_n
                mono = _q_pow(e4_n) * sign
                qb   = _q_binomial(i, m_n, q_actual)    # 0 when m_n < 0 or m_n > i
                qp   = _q_pochhammer(_q_pow(4*(n-1-j)), q_inv, m_n)
                val  = mono * qb * qp
                if val:
                    if not sage_polynomials:
                        data_n[(i, j, k, l)] = val
                    else:
                        data_n[(i, j, k, l)] = val.to_sage()

    Rp = SparseTensor(shape, data=data_p)
    Rn = SparseTensor(shape, data=data_n)

    # hp[i,i] = q^(i + (1-n)/2) = q^(i - (n-1)/2),  key e4 = 4*i - 2*(n-1)
    # hn[i,i] = 1 / hp[i,i]                        ,  key e4 = 2*(n-1) - 4*i
    if not sage_polynomials:
        hp = SparseTensor((n, n), data={(i, i): _q_pow(4*i - 2*(n-1)) for i in range(n)})
        hn = SparseTensor((n, n), data={(i, i): _q_pow(2*(n-1) - 4*i) for i in range(n)})
    else:
        hp = SparseTensor((n, n), data={(i, i): _q_pow(4*i - 2*(n-1)).to_sage() for i in range(n)})
        hn = SparseTensor((n, n), data={(i, i): _q_pow(2*(n-1) - 4*i).to_sage() for i in range(n)})

    return RMatrix(Rp, Rn, hp, hn)

def prefactor_colored_jones(n, writhe, sage_polynomial = False):
    n = n + 1

    if not sage_polynomial:
        return _q_pow(writhe * ((n**2) -1))
    else:
        return _q_pow(writhe * ((n**2) -1)).to_sage()
