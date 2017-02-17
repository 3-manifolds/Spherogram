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
except:
    _within_sage = False
    import decorator

try:
    import snappy
    _have_snappy = True
except ImportError:
    print('Could not import snappy in sage_helper.py.')
    _have_snappy = False
    
import doctest, re, types

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

if _within_sage:
    class DocTestParser(doctest.DocTestParser):
        def parse(self, string, name='<string>'):
            string = re.subn('([\n\A]\s*)sage:', '\g<1>>>>', string)[0]
            if _have_snappy:
                string = re.subn('(^\s*)\| ', '\g<1>', string, flags=re.M)[0]
            return doctest.DocTestParser.parse(self, string, name)

    globs = dict()
else:
    class DocTestParser(doctest.DocTestParser):
        def parse(self, string, name='<string>'):
            if _have_snappy:
                string = re.subn('(^\s*)\| ', '\g<1>', string, flags=re.M)[0]
            return doctest.DocTestParser.parse(self, string, name)
        
    globs = dict()

def print_results(module, results):
    print(module.__name__ + ':')
    print('   %s failures out of %s tests.' %  (results.failed, results.attempted))
    
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
    

        
    

