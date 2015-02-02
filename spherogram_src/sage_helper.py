"""
Helper code for dealing with additional functionality when Sage is
present.

Any method which works only in Sage should be decorated with
"@sage_method" and any doctests (in Sage methods or not) which should
be run only in Sage should be styled with input prompt "sage:" rather
than the usual ">>>".
"""
try:
    import sage.all
    _within_sage = True
except:
    _within_sage = False

import functools, doctest, re

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


# Not currently used, but could be exploited by an interpeter to hide
# sage_methods when in plain Python.

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

# Used for doctesting

if _within_sage:
    class DocTestParser(doctest.DocTestParser):
        def parse(self, string, name='<string>'):
            string = re.subn('([\n\A]\s*)sage:', '\g<1>>>>', string)[0]
            return doctest.DocTestParser.parse(self, string, name)
else:
    DocTestParser = doctest.DocTestParser


def print_results(module, runner):
    print(module.__name__ + ':')
    print('   %s failures out of %s tests.' %  (runner.failures, runner.tries))
    
def doctest_modules(modules, verbose=False):
    finder = doctest.DocTestFinder(parser=DocTestParser())
    failures, tries = 0, 0
    for module in modules:
        runner = doctest.DocTestRunner(verbose=verbose)
        for test in finder.find(module):
            runner.run(test)
        failures += runner.failures
        tries += runner.tries
        print_results(module, runner)
    print('\nAll doctests:\n   %s failures out of %s tests.' % (failures, tries))
    

        
    

