from __future__ import print_function
import spherogram, spherogram.links, spherogram.links.test
import spherogram.links.simplify, spherogram.links.morse
import spherogram.links.seifert

import spherogram.test_helper as test_helper
import re, getopt, sys


modules = [spherogram.codecs.DT, spherogram.codecs.Base64LikeDT,
           spherogram.graphs, spherogram.presentations,
           spherogram.links.links, spherogram.links.links_base,
           spherogram.links.random_links, spherogram.links.orthogonal,
           spherogram.links.simplify, spherogram.links.invariants,
           spherogram.links.morse, spherogram.links.seifert]

# Apply the monkey-patches that snappy applies when it is imported.
if test_helper._have_snappy:
    import snappy
    spherogram.links.links_base.Link.exterior = snappy._link_exterior

def run_doctests(verbose=False, print_info=True):
    if test_helper._have_snappy:
        snappy.number.Number._accuracy_for_testing = 8
    if test_helper._have_snappy and test_helper._within_sage:
        snappy.Manifold.use_field_conversion('snappy')
        snappy.ManifoldHP.use_field_conversion('snappy')
    return test_helper.doctest_modules(modules, verbose, print_info)

def run_all_tests():
    optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
    verbose = len(optlist) > 0
    results = run_doctests(verbose)
    print()
    spherogram.links.test.run()
    return results.failed

if __name__ == '__main__':
    run_all_tests()
