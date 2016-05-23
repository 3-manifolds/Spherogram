import snappy, sys
sys.path.append('/Users/dunfield/s/dev')
import spherogram.graphs as graphs
from spherogram import DTcodec

"""
Morwen's code triangulates 1100 links/second
Marc's code triangulates 
"""

dtcodes = [line.strip() for line in open('DTcodes')][:1000]

def test_Morwen():
    for dt in dtcodes:
        snappy.Manifold('DT[' + dt + ']')

def test_Marc():
    for dt in dtcodes:
        L = DTcodec(dt)

DT = DTcodec('mamchgilkjbamfed')
G = DT.fat_graph

def test1():
    v = list(G.vertices)[0]
    D = G.incidence_dict
    for i in xrange(155000):
        D[v]
        
def test2():
    keys = G.incidence_dict.keys()
    v = keys[0]
    D = dict( (k, (i, range(10*i, 10*(i+10)))) for i,k in enumerate(keys))
    for i in xrange(155000):
        D[v] 

def test3():
    keys = G.incidence_dict.keys()
    v = keys[0]
    D = dict( (k, None) for i,k in enumerate(keys))
    for i in xrange(155000):
        D[v] 

def test4():
    e = list(G.edges)[0]
    v = e[0]
    w = e[1]
    for i in xrange(50000):
        v == w
    return v


class TestObject(object):
    pass

def test5():
    D = dict([(TestObject(),i) for i in range(10)])
    k = D.keys()[0]
    for i in xrange(10000000):
        D[k]

def test6():
    D = dict([(i, TestObject()) for i in range(10)])
    k = D.keys()[0]
    for i in xrange(10000000):
        D[k]

def test7():
    a, b = TestObject(), TestObject()
    for i in xrange(1000000):
        a == b


def test8():
    L = graphs.CyclicList(range(10))
    L = range(10)
    for i in xrange(298128):
        L[2], len(L)
