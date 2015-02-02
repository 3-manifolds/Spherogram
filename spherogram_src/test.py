import spherogram, spherogram.links.test, spherogram.links.simplify
from spherogram.sage_helper import _within_sage
import snappy, doctest
SAGE_METHOD = doctest.register_optionflag('SAGE_METHOD')

snappy.number.Number._accuracy_for_testing = 8
modules = [spherogram.codecs.DT, spherogram.graphs, spherogram.presentations,
           spherogram.links.links, spherogram.links.links_base,
           spherogram.links.random_links, spherogram.links.orthogonal,
           spherogram.links.simplify, spherogram.links.invariants]

spherogram.links.links_base.Link.exterior = snappy._link_exterior

if _within_sage:
    snappy.Manifold.use_field_conversion('snappy')
    snappy.ManifoldHP.use_field_conversion('snappy')
    import spherogram.links.morse, spherogram.links.invariants
    modules += [spherogram.links.morse]
else:
    from doctest import Example, SKIP
    class SageAwareExample(Example):
        """Subclass of doctest.Example which treats the SAGE_METHOD flag as
        being equivalent to SKIP
        
        """
        def __init__(self, source, want, exc_msg=None, lineno=0, indent=0,
                     options=None):
            Example.__init__(self, source, want, exc_msg, lineno, indent,
                             options)
            if self.options.get(SAGE_METHOD, False):
                self.options[SKIP] = True

    doctest.Example = SageAwareExample

for module in modules:
    results = doctest.testmod(module)
    print('%s: \n    %s failures out of %s tests\n' % (module.__name__, results[0], results[1]))

spherogram.links.test.run()
