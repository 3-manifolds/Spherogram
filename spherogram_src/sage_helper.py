"""
Helper code for dealing with additional functionality when Sage is
present.

Any method which works only in Sage should be decorated with
"@sage_method" and any doctests (in Sage methods or not) which should
be run only in Sage should be styled with input prompt "sage:" rather
than the usual ">>>".

Similarly, doctests which require SnapPy should be styled in a block
where the first non-whitespace character is | followed by a space.
"""
try:
    import sage.all
    _within_sage = True
except Exception:
    _within_sage = False
    import decorator


import doctest


class SageNotAvailable(Exception):
    pass


if _within_sage:
    def sage_method(function):
        function._sage_method = True
        return function
    from packaging.version import parse
    from sage.version import version as _sage_version
    sage_pd_clockwise = (parse(_sage_version) <= parse('10.1.beta0'))
else:
    def _sage_method(function, *args, **kw):
        raise SageNotAvailable('Sorry, this feature requires using SnapPy inside Sage.')

    def sage_method(function):
        return decorator.decorator(_sage_method, function)

    sage_pd_clockwise = None
