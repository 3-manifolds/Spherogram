"""

Only the test for K8n1 is run as part of the doctest suite.  To test
them all, do::

  python -m spherogram.links.bands.regsession

You should see::

  Starting K6a3
      OK
  Starting K8n1
      OK
  Starting K14n21673
      OK
  Starting 16n1008880
      OK

"""

from .core import banded_links
from ...test_helper import snappy


original = {'K6a3':['L7a4', 'L7a5', 'L7a6', 'L8a11', 'L8a3', 'L8a6', 'L8n2', 'L9n13'],
            'K8n1':['L10a41', 'L10a60', 'L10n18', 'L10n21', 'L10n22', 'L10n24',
                    'L10n27', 'L10n30', 'L10n36', 'L10n38', 'L10n40', 'L10n41',
                    'L10n44', 'L10n48', 'L10n59', 'L10n61', 'L10n64', 'L11n198',
                    'L11n37', 'L11n53', 'L11n57', 'L12n1043', 'L12n1093', 'L12n1109',
                    'L12n1183', 'L12n1325', 'L12n428', 'L12n802', 'L12n822', 'L12n851',
                    'L14n15889', 'L7a5', 'L7a6', 'L7n2', 'L8n2', 'L9a16', 'L9a27',
                    'L9a34', 'L9n11', 'L9n13', 'L9n14', 'L9n16', 'L9n9'],
            'K14n21673':['L10n14', 'L10n18', 'L10n24', 'L11a221', 'L11a232', 'L11a247',
                         'L11n137', 'L11n6', 'L12a473', 'L12a549', 'L12a671', 'L12n105',
                         'L12n1270', 'L12n1375', 'L12n242', 'L12n424', 'L12n551', 'L12n627',
                         'L12n634', 'L12n651', 'L12n684', 'L13a1628', 'L13a2812',
                         'L13a3563', 'L13a4143', 'L13a4280', 'L13n1400', 'L13n2070',
                         'L13n3058', 'L13n3059', 'L13n3188', 'L13n3204', 'L13n3355',
                         'L13n3492', 'L13n3843', 'L13n388', 'L13n4077', 'L13n4218',
                         'L13n425', 'L13n4393', 'L13n4476', 'L13n4802', 'L13n4891',
                         'L13n5149', 'L13n5277', 'L13n5682', 'L13n6457', 'L13n723',
                         'L13n773', 'L14a12973', 'L14a13404', 'L14a19559', 'L14a3177',
                         'L14a5723', 'L14n13143', 'L14n13800', 'L14n15423', 'L14n15748',
                         'L14n15964', 'L14n16063', 'L14n16313', 'L14n17248', 'L14n17411',
                         'L14n22250', 'L14n22561', 'L14n24030', 'L14n24832', 'L14n25099',
                         'L14n27706', 'L14n3180', 'L14n32161', 'L14n36351', 'L14n365',
                         'L14n39756', 'L14n8190', 'L7a5', 'L9a27', 'L9n11'],
            '16n1008880':['L11a321', 'L12a1299', 'L12a227', 'L13a3112', 'L13a4312',
                          'L13n3416', 'L13n5029', 'L14a12316', 'L14a13220',
                          'L14a14675', 'L14a16563', 'L14a17068', 'L14a3399',
                          'L14a4618', 'L14a6045', 'L14a7390', 'L14a9531', 'L14n1816',
                          'L14n25590', 'L14n2624', 'L14n2625', 'L14n27107', 'L14n28138',
                          'L14n33194', 'L14n3977', 'L7a5', 'L9a21']}



def hyperbolize(L):
    for i in range(2):
        E = L.exterior()
        for j in range(5):
            if E.solution_type(enum=True) in {1, 2}:
                return E
            E.randomize()


def regression_links_only(L):
    ans = dict()
    for data in banded_links(L, max_twists=2, max_band_len=None, paths='shortest'):
        E = hyperbolize(data[0])
        if E is not None:
            ident = snappy.HTLinkExteriors.identify(E, extends_to_link=True)
            if ident:
                ident = ident.name()
                if ident not in ans:
                    ans[ident] = data[1]
    return ans


def test_all():
    for knot, links in original.items():
        print('Starting ' + knot)
        L = snappy.Manifold(knot).link()
        new_links = sorted(regression_links_only(L))
        if new_links == links:
            print('    OK')
        else:
            diff = set(new_links).symmetric_difference(set(links))
            print('    Error: ', diff)

            
def test_one():
    """
    >>> test_one()   #doctest: +SNAPPY
    True
    """
    knot = 'K8n1'
    L = snappy.Manifold(knot).link()
    new_links = sorted(regression_links_only(L))
    return new_links == original[knot]


if __name__ == '__main__':
    test_all()
