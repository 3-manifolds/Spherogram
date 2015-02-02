import spherogram, spherogram.links, spherogram.links.test
import spherogram.links.simplify, spherogram.links.morse

from spherogram.sage_helper import _within_sage, doctest_modules
import snappy, re, getopt, sys

snappy.number.Number._accuracy_for_testing = 8
modules = [spherogram.codecs.DT, spherogram.graphs, spherogram.presentations,
           spherogram.links.links, spherogram.links.links_base,
           spherogram.links.random_links, spherogram.links.orthogonal,
           spherogram.links.simplify, spherogram.links.invariants, spherogram.links.morse]

spherogram.links.links_base.Link.exterior = snappy._link_exterior

if _within_sage:
    snappy.Manifold.use_field_conversion('snappy')
    snappy.ManifoldHP.use_field_conversion('snappy')


optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
verbose = len(optlist) > 0
doctest_modules(modules, verbose)

spherogram.links.test.run()
