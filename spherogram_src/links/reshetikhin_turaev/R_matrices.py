from .dict_laurent_polynomial import DictLaurentPolynomial

from .sparse_array import SparseTensor

import csv, ast, pathlib, os

dir_path = pathlib.Path(__file__).resolve().parent

_cache = dict()

def laurent_sparse_tensor_from_file(file, vars = ['t', 'q']):
    reader = csv.reader(file)
    header = next(reader)
    shape = ast.literal_eval(header[0])

    assert header[1] == 'LaurentPolynomial', f'Expected type LaurentPolynomial, got {header[1]}'

    data = dict()
    for line in reader:
        key, value = line
        key = tuple(ast.literal_eval(key))
        assert key not in data.keys(), f'{key} appeared multiple times in {file.name}'
        value = DictLaurentPolynomial.from_str(value, vars = vars)
        data[key] = value

    return SparseTensor(shape = shape, data = data)

def laurent_sparse_tensor_from_path(path, vars = ['t', 'q'], compressed = False):
    if compressed:
        import bz2
        with bz2.open(path, 'rt') as f:
            return laurent_sparse_tensor_from_file(f, vars = vars)
    else:
        with open(path, 'r') as f:
            return laurent_sparse_tensor_from_file(f, vars = vars)

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
    def laurent_R_from_directory(dir_path, vars = ['t', 'q'], compressed = False):
        names = [name + '.csv' + ('.bz2' if compressed else '') 
                 for name in ['Rp', 'Rn', 'hp', 'hn']]
        
        tensors = [laurent_sparse_tensor_from_path(os.path.join(dir_path, name),
                                                   vars = vars,
                                                   compressed = compressed)
                   for name in names]
        
        return RMatrix(*tensors)

def colored_links_gould_R_matrices(n):
    if 0 < n <= 4:
        key = f'V{n}'
        if key in _cache.keys():
            return _cache[key]
        else:
            _cache[key] = RMatrix.laurent_R_from_directory(dir_path = os.path.join(dir_path, f'R_matrices/V{n}/'))
            return _cache[key]
    else:
        raise NotImplementedError