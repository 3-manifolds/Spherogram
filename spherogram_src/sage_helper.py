"""
Helper code for dealing with additional functionality when Sage is
present.
"""
try:
    import sage.all
    _within_sage = True
except:
    _within_sage = False

import functools

class SageNotAvailable(Exception):
    pass

if _within_sage:
    def sage_method(function):
        @functools.wraps(function)
        def wrapper(*args, **kwds):
            return function(*args, **kwds)
        wrapper._sage_method = True
        return wrapper
else:
    def sage_method(function):
        @functools.wraps(function)
        def wrapper(*args, **kwds):
            raise SageNotAvailable('Sorry, this feature requires using SnapPy inside Sage.')
        wrapper._sage_method = True
        return wrapper

def sage_methods(obj):
    ans = []
    for attr in dir(obj):
        try:
            methods = getattr(obj, attr)
            if methods._sage_method == True:
                ans.append(methods)
        except AttributeError:
            pass
    return ans



