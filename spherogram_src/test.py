from __future__ import print_function
import snappy
import spherogram, spherogram.links, spherogram.links.test
import spherogram.links.simplify, spherogram.links.morse
import spherogram.links.seifert

import spherogram.sage_helper as sage_helper
import re, getopt, sys


modules = [spherogram.codecs.DT, spherogram.codecs.Base64LikeDT,
           spherogram.graphs, spherogram.presentations,
           spherogram.links.links, spherogram.links.links_base,
           spherogram.links.random_links, spherogram.links.orthogonal,
           spherogram.links.simplify, spherogram.links.invariants,
           spherogram.links.morse, spherogram.links.seifert]

# Apply the monkey-patches that snappy applies when it is imported.
spherogram.links.links_base.Link.exterior = snappy._link_exterior
spherogram.links.links_base.Link._lookup_DT = snappy._link_lookup_DT

def run_doctests(verbose=False, print_info=True):
    snappy.number.Number._accuracy_for_testing = 8
    if sage_helper._within_sage:
        snappy.Manifold.use_field_conversion('snappy')
        snappy.ManifoldHP.use_field_conversion('snappy')
    return sage_helper.doctest_modules(modules, verbose, print_info)

if __name__ == '__main__':
    optlist, args = getopt.getopt(sys.argv[1:], 'v', ['verbose'])
    verbose = len(optlist) > 0
    run_doctests(verbose)
    print()
    spherogram.links.test.run()
