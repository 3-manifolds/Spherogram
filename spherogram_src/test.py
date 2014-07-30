import spherogram, spherogram.links.test
import snappy, doctest

snappy.number.Number._accuracy_for_testing = 8
modules = [spherogram.codecs.DT, spherogram.graphs, spherogram.presentations,
           spherogram.links.links, spherogram.links.orthogonal]

if snappy.SnapPy._within_sage:
    snappy.Manifold.use_field_conversion('snappy')
    snappy.ManifoldHP.use_field_conversion('snappy')
    modules.append(spherogram.links.invariants)

for module in modules:
    results = doctest.testmod(module)
    print('%s: \n    %s failures out of %s tests\n' % (module.__name__, results[0], results[1]))

spherogram.links.test.run()
