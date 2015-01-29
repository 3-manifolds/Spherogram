import spherogram, spherogram.links.test, spherogram.links.simplify
import snappy, doctest

snappy.number.Number._accuracy_for_testing = 8
modules = [spherogram.codecs.DT, spherogram.graphs, spherogram.presentations,
           spherogram.links.links, spherogram.links.links_base,
           spherogram.links.random_links, spherogram.links.orthogonal,
           spherogram.links.simplify]

if snappy.SnapPy._within_sage:
    snappy.Manifold.use_field_conversion('snappy')
    snappy.ManifoldHP.use_field_conversion('snappy')
    spherogram.links.links_base.Link.exterior = snappy._link_exterior
    import spherogram.links.morse, spherogram.links.invariants
    modules += [spherogram.links.invariants, spherogram.links.morse]

for module in modules:
    results = doctest.testmod(module)
    print('%s: \n    %s failures out of %s tests\n' % (module.__name__, results[0], results[1]))

spherogram.links.test.run()
