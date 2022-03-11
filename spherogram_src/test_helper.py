"""
Helper code for dealing with additional functionality when Sage is
present.

Any method which works only in Sage should be decorated with
"@sage_method" and any doctests (in Sage methods or not) which should
be run only in Sage should be styled with input prompt "sage:" rather
than the usual ">>>".

Additionally, doctests which require SnapPy marked with the "+SNAPPY"
flag.
"""
try:
    import sage.all
    _within_sage = True
except Exception:
    _within_sage = False
    import decorator

try:
    import snappy
    _have_snappy = True
except ImportError:
    print('WARNING: Could not import snappy in test_helper.py.')
    _have_snappy = False


import doctest
import re
import types


class SageNotAvailable(Exception):
    pass


if _within_sage:
    def sage_method(function):
        function._sage_method = True
        return function
else:
    def _sage_method(function, *args, **kw):
        raise SageNotAvailable('Sorry, this feature requires using SnapPy inside Sage.')

    def sage_method(function):
        return decorator.decorator(_sage_method, function)


SNAPPY_FLAG = doctest.register_optionflag('SNAPPY')


def filter_out_snappy(pieces):
    ans = []
    for piece in pieces:
        if _have_snappy or not isinstance(piece, doctest.Example):
            ans.append(piece)
        elif not piece.options.get(SNAPPY_FLAG, False):
            ans.append(piece)
    return ans


if _within_sage:
    class DocTestParser(doctest.DocTestParser):
        def parse(self, string, name='<string>'):
            string = re.subn(r'(\n\s*)sage:|(\A\s*)sage:', r'\g<1>>>>', string)[0]
            pieces = doctest.DocTestParser.parse(self, string, name)
            return filter_out_snappy(pieces)

    globs = dict()
else:
    class DocTestParser(doctest.DocTestParser):
        def parse(self, string, name='<string>'):
            pieces = doctest.DocTestParser.parse(self, string, name)
            return filter_out_snappy(pieces)

    globs = dict()


def print_results(module, results):
    print(module.__name__ + ':')
    print('   %s failures out of %s tests.' % (results.failed, results.attempted))


def doctest_modules(modules, verbose=False, print_info=True, extraglobs=dict()):
    finder = doctest.DocTestFinder(parser=DocTestParser())
#    full_extraglobals = dict(globs.items() + extraglobs.items())
    full_extraglobals = globs.copy()
    full_extraglobals.update(extraglobs)
    failed, attempted = 0, 0
    for module in modules:
        if isinstance(module, types.ModuleType):
            runner = doctest.DocTestRunner(verbose=verbose)
            for test in finder.find(module, extraglobs=full_extraglobals):
                runner.run(test)
            result = runner.summarize()
        else:
            result = module(verbose=verbose)
        failed += result.failed
        attempted += result.attempted
        if print_info:
            print_results(module, result)

    if print_info:
        print('\nAll doctests:\n   %s failures out of %s tests.' % (failed, attempted))
    return doctest.TestResults(failed, attempted)
