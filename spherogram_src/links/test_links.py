import spherogram
import unittest
from random import randrange

class TestLinkFunctions(unittest.TestCase):
    
    def setUp(self):
        # Trefoil
        a = spherogram.links.links.Crossing('a')
        b = spherogram.links.links.Crossing('b')
        c = spherogram.links.links.Crossing('c')
        a[2] = b[1]
        b[3] = c[0]
        c[2] = a[1]
        a[3] = b[0]
        b[2] = c[1]
        c[3] = a[0]
        self.K3a1 = spherogram.links.links.Link([a,b,c])
        self.K3a1_prime = spherogram.links.links.Link('K3a1')

        self.K7a2 = spherogram.links.links.Link('K7a2')
        self.K8a3 = spherogram.links.links.Link('K8a3')
        self.K8a13 = spherogram.links.links.Link('K8a13')

        # Hopf Link
        a = spherogram.links.links.Crossing('a') 
        b = spherogram.links.links.Crossing('b') 
        a[0]=b[1] 
        a[1]=b[0] 
        a[2]=b[3] 
        a[3]=b[2] 
        self.L2a1 = spherogram.links.links.Link([a,b])

        #Borromean Link (3) 
        a = spherogram.links.links.Crossing('a') 
        b = spherogram.links.links.Crossing('b') 
        c = spherogram.links.links.Crossing('c') 
        d = spherogram.links.links.Crossing('d') 
        e = spherogram.links.links.Crossing('e') 
        f = spherogram.links.links.Crossing('f') 
        a[2] = f[1]
        a[3] = e[0] 
        b[1] = a[0] 
        b[2] = e[3] 
        c[0] = a[1] 
        c[1] = b[0] 
        d[3] = c[2] 
        d[0] = b[3] 
        e[2] = d[1] 
        e[1] = f[0] 
        f[2] = c[3] 
        f[3] = d[2] 
        self.Borr = spherogram.links.links.Link([a,b,c,d,e,f])
        
        self.knots = [self.K3a1, self.K3a1_prime, self.K7a2, self.K8a3, self.K8a13]
        self.links = [self.L2a1, self.Borr]
        self.all_links = self.knots + self.links
               
    def testWrithe(self):
        self.assertEqual(abs(self.K3a1.writhe()),       3)
        self.assertEqual(abs(self.K3a1_prime.writhe()), 3)
        self.assertEqual(abs(self.K7a2.writhe()),       3)
        self.assertEqual(abs(self.K8a3.writhe()),       2)
        self.assertEqual(abs(self.K8a13.writhe()),      4)
        self.assertEqual(abs(self.L2a1.writhe()),       2)
        self.assertEqual(abs(self.Borr.writhe()),       0)

    def testLinkingNumber(self):
        for k in self.knots:
            self.assertEqual(k.linking_number(),           0)
        self.assertEqual(abs(self.L2a1.linking_number()),  1)
        self.assertEqual(abs(self.Borr.linking_number()),  0)

    def testLinkingMatrix(self):
        M = [[0]]
        N = [[0, -1], [-1, 0]]
        P = [[0,0,0],[0,0,0],[0,0,0]]
        for k in self.knots:
            self.assertEqual(k.linking_matrix(),      M)
        self.assert_(self.L2a1.linking_matrix() == N or self.L2a1.linking_matrix() == -1*N)
        self.assertEqual(self.Borr.linking_matrix(),       P)

    def testKnotGroup(self):
        # Correct number of generators? abelian invariants?
        for k in self.all_links:
            self.assertEqual(len(k.knot_group().generators()),             len(k.crossings))
            self.assertEqual(k.knot_group().abelian_invariants()[0],       0)
            self.assertEqual(len(k.knot_group().abelian_invariants()),     len(k.link_components))
       
    def testAlexanderPoly(self): 
        t = var('t')
        a = var('a')
        b = var('b')
        c = var('c')
        # method = 'wirt'
        self.assertEqual(self.K3a1.alexander_poly(),       t - 1 + t**-1)
        self.assertEqual(self.K3a1_prime.alexander_poly(), t - 1 + t**-1)
        # is this really creating the knot 7_2 ?
        self.assertEqual(self.K7a2.alexander_poly(),       -t**2 + 5*t + -7 + 5*t**-1 - t**-2)
        self.assertEqual(self.K8a3.alexander_poly(),       t**3 - 3*t**2 + 6*t - 7 + 6*t**-1 - 3*t**-2 + t**-3)
        self.assertEqual(self.K8a13.alexander_poly(),      -t**3 + 3*t**2 - 4*t + 5 - 4*t**-1 + 3*t**-2 - t**-3)
        self.assert_(self.L2a1.alexander_poly() == sqrt(t) - 1/sqrt(t) or self.L2a1.alexander_poly() == - sqrt(t) + 1/sqrt(t))
        self.assertEqual(self.Borr.alexander_poly(),       -t**2 + 4*t + 4*t**-1 - t**-2 - 6)
        # method = 'snappy'
        self.assertEqual(self.K3a1.alexander_poly(method='snappy'),       a**2 - a + 1)
        self.assertEqual(self.K3a1_prime.alexander_poly(method='snappy'), a**2 - a + 1)
        self.assertEqual(self.K7a2.alexander_poly(method='snappy'),       a**4 - 5*a**3 + 7*a**2 - 5*a + 1)
        self.assertEqual(self.K8a3.alexander_poly(method='snappy'),       a**6 - 3*a**5 + 6*a**4 - 7*a**3 + 6*a**2 - 3*a + 1)
        self.assertEqual(self.K8a13.alexander_poly(method='snappy'),      a**6 - 3*a**5 + 4*a**4 - 5*a**3 + 4*a**2 - 3*a + 1)
        
    def testConnectedSum(self):
        random_index = randrange(0,len(self.knots))
        k1 = self.knots[random_index]
        random_index = randrange(0, len(self.knots))
        k2 = self.knots[random_index]
        Sum = k1.connected_sum(k2)
        self.assertEqual(Sum.alexander_poly(), k1.alexander_poly()*k2.alexander_poly())
        self.assertEqual(Sum.signature(),      k1.signature()+k2.signature())
        self.assertEqual(Sum.determinant(),    k1.determinant()*k2.determinant())

    def testGoeritzMatrix(self):
        self.assert_(self.K3a1.goeritz_matrix() == Matrix([3])       or self.K3a1.goeritz_matrix() == -1*Matrix([3]))
        self.assert_(self.K3a1_prime.goeritz_matrix() == Matrix([3]) or self.K3a1_prime.goeritz_matrix() == -1*Matrix([3]))
        # Goeritz matrix's columns, etc. get rearranged ?
#        self.assert_(self.K7a2.goeritz_matrix() == Matrix([0])       or self.K7a2.goeritz_matrix() == -1*Matrix([0])
#        self.assert_(self.K8a3.goeritz_matrix() == Matrix([0])       or self.K8a3.goeritz_matrix() == -1*Matrix([0]))
#        self.assert_(self.K8a13.goeritz_matrix() == Matrix([0])      or self.K8a13.goeritz_matrix() == -1*Matrix([0]))
       
    def testSignature(self):
        self.assertEqual(abs(self.K3a1.signature()),       2)
        self.assertEqual(abs(self.K3a1_prime.signature()), 2)
        self.assertEqual(abs(self.K7a2.signature()),       2)
        self.assertEqual(abs(self.K8a3.signature()),       2)
        self.assertEqual(abs(self.K8a13.signature()),      4)
        self.assertEqual(abs(self.L2a1.signature()),       1)
        self.assertEqual(abs(self.Borr.signature()),       0)

    def testCopy(self):
        random_index = randrange(0,len(self.knots))
        k1 = self.knots[random_index]
        k1_prime = k1.copy()
        self.assertEqual(k1.determinant(),    k1_prime.determinant())
        self.assertEqual(k1.writhe(),         k1_prime.writhe())
        self.assertEqual(k1.signature(),      k1_prime.signature())
        self.assertEqual(k1.alexander_poly(), k1_prime.alexander_poly())
        #component of randomness in goeritz matrix ? rows/columns mixed up
        #self.assertEqual(k1.goeritz_matrix(), k1_prime.goeritz_matrix())

        random_index = randrange(0,len(self.links))
        k2 = self.links[random_index]
        k2_prime = k2.copy()
        self.assertEqual(k2.determinant(),    k2_prime.determinant())
        self.assertEqual(k2.writhe(),         k2_prime.writhe())
        self.assertEqual(k2.signature(),      k2_prime.signature())
        self.assertEqual(k2.alexander_poly(), k2_prime.alexander_poly())

    def testMirror(self):
        random_index = randrange(0,len(self.knots))
        k1 = self.knots[random_index]
        k1_prime = k1.mirror()
        self.assert_(k1.signature() == -1*k1_prime.signature())
        self.assert_(k1.writhe() == -1*k1_prime.writhe())

        # What does work for mirror of a link?
#        random_index = randrange(0,len(self.links))
#        k2 = self.links[random_index]
#        k2_prime = k2.mirror()
#        self.assert_(k2.signature() == -1*k2_prime.signature())
#        self.assert_(k2.writhe() == -1*k2_prime.writhe())

    def testDet(self):
        self.assertEqual(self.K3a1.determinant(),                 3)
        self.assertEqual(self.K3a1.determinant(method='color'),   3)
        self.assertEqual(self.K3a1.determinant(method='goeritz'), 3)
        
        self.assertEqual(self.K3a1_prime.determinant(),                 3)
        self.assertEqual(self.K3a1_prime.determinant(method='color'),   3)
        self.assertEqual(self.K3a1_prime.determinant(method='goeritz'), 3)

        self.assertEqual(self.K7a2.determinant(),                 19)
        self.assertEqual(self.K7a2.determinant(method='color'),   19)
        self.assertEqual(self.K7a2.determinant(method='goeritz'), 19)

        self.assertEqual(self.K8a3.determinant(),                 27)
        self.assertEqual(self.K8a3.determinant(method='color'),   27)
        self.assertEqual(self.K8a3.determinant(method='goeritz'), 27)

        self.assertEqual(self.K8a13.determinant(),                 21)
        self.assertEqual(self.K8a13.determinant(method='color'),   21)
        self.assertEqual(self.K8a13.determinant(method='goeritz'), 21)

        self.assertEqual(self.L2a1.determinant(),                 2)
        self.assertEqual(self.L2a1.determinant(method='color'),   2)
        self.assertEqual(self.L2a1.determinant(method='goeritz'), 2)

        self.assertEqual(self.Borr.determinant(),                 16)
        self.assertEqual(self.Borr.determinant(method='color'),   16)
        self.assertEqual(self.Borr.determinant(method='goeritz'), 16)

suite = unittest.TestLoader().loadTestsFromTestCase(TestLinkFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)
