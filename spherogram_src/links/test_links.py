import spherogram
import unittest
from random import randrange

class TestLinkFunctions(unittest.TestCase):
    
    def setUp(self):
        # Trefoil
        a = spherogram.Crossing('a')
        b = spherogram.Crossing('b')
        c = spherogram.Crossing('c')
        a[2] = b[1]
        b[3] = c[0]
        c[2] = a[1]
        a[3] = b[0]
        b[2] = c[1]
        c[3] = a[0]
        self.Tref = spherogram.Link([a,b,c])

        self.K3_1 = spherogram.Link('3_1')
        self.K7_2 = spherogram.Link('7_2')
        self.K8_3 = spherogram.Link('8_3')
        self.K8_13 = spherogram.Link('8_13')

        # Hopf Link
        a = spherogram.Crossing('a') 
        b = spherogram.Crossing('b') 
        a[0]=b[1] 
        a[1]=b[0] 
        a[2]=b[3] 
        a[3]=b[2] 
        self.L2a1 = spherogram.Link([a,b])

        #Borromean Link (3) 
        a = spherogram.Crossing('a') 
        b = spherogram.Crossing('b') 
        c = spherogram.Crossing('c') 
        d = spherogram.Crossing('d') 
        e = spherogram.Crossing('e') 
        f = spherogram.Crossing('f') 
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
        self.Borr = spherogram.Link([a,b,c,d,e,f])
        
        self.L6a2 = spherogram.Link('L6a2')
        self.L6a4 = spherogram.Link('L6a4')
        self.L7a3 = spherogram.Link('L7a3')

        self.knots = [self.K3_1, self.K7_2, self.K8_3, self.K8_13]
        self.links = [self.L2a1, self.L6a4, self.L6a2, self.L7a3]
        self.all_links = self.knots + self.links
               
    def random_knot(self):
        random_index = randrange(0,len(self.knots))
        return self.knots[random_index]

    def random_link(self):
        random_index = randrange(0,len(self.links))
        return self.links[random_index]

    def testWrithe(self):
        self.assertEqual(abs(self.Tref.writhe()),       3)
        self.assertEqual(abs(self.K3_1.writhe()),       3)
        self.assertEqual(abs(self.K7_2.writhe()),       7)
        self.assertEqual(abs(self.K8_3.writhe()),       0)
        self.assertEqual(abs(self.K8_13.writhe()),      2)
        self.assertEqual(abs(self.L2a1.writhe()),       2)
        self.assertEqual(abs(self.L6a2.writhe()),       6)
        self.assertEqual(abs(self.Borr.writhe()),       0)
        self.assertEqual(abs(self.L6a4.writhe()),       0)
        self.assertEqual(abs(self.L7a3.writhe()),       3)

    def testLinkingNumber(self):
        for k in self.knots:
            self.assertEqual(k.linking_number(),           0)
        self.assertEqual(abs(self.L2a1.linking_number()),  1)
        self.assertEqual(abs(self.L6a2.linking_number()),  3)
        self.assertEqual(abs(self.Borr.linking_number()),  0)
        self.assertEqual(abs(self.L6a4.linking_number()),  0)
        self.assertEqual(abs(self.L7a3.linking_number()),  0)

    def testLinkingMatrix(self):
        A = [[0]]
        B = [[0, -1], [-1, 0]]
        B2 = [[0, 1], [1, 0]]
        C = [[0, -3], [-3, 0]]
        C2 = [[0, 3], [3, 0]]
        D = [[0,0,0],[0,0,0],[0,0,0]]
        E = [[0, 0], [0, 0]]
        for k in self.knots:
            self.assertEqual(k.linking_matrix(),      A)
        self.assert_(self.L2a1.linking_matrix() == B or self.L2a1.linking_matrix() == B2)
        self.assert_(self.L6a2.linking_matrix() == C or self.L6a2.linking_matrix() == C2)
        self.assertEqual(self.Borr.linking_matrix(),       D)
        self.assertEqual(self.L6a4.linking_matrix(),       D)
        self.assertEqual(self.L7a3.linking_matrix(),       E)

    def testKnotGroup(self):
        # Correct number of generators? abelian invariants?
        for k in self.all_links:
            self.assertEqual(len(k.knot_group().generators()),             len(k.crossings))
            self.assertEqual(k.knot_group().abelian_invariants()[0],       0)
            self.assertEqual(len(k.knot_group().abelian_invariants()),     len(k.link_components))
       
    def testAlexanderPoly(self): 
        t = var('t')
        a = var('a')

        # method = 'wirt'
        self.assertEqual(self.Tref.alexander_poly(),       t - 1 + t**-1)
        self.assertEqual(self.K3_1.alexander_poly(),       t - 1 + t**-1)
        self.assertEqual(self.K7_2.alexander_poly(),       3*t - 5 + 3*t**-1)
        self.assertEqual(self.K8_3.alexander_poly(),       -4*t + 9 - 4*t**-1)
        self.assertEqual(self.K8_13.alexander_poly(),      2*t**2 - 7*t + 11 - 7*t**-1 + 2*t**-2)
        self.assert_(self.L2a1.alexander_poly() == sqrt(t) - 1/sqrt(t) or self.L2a1.alexander_poly() == - sqrt(t) + 1/sqrt(t))
        self.assertEqual(self.Borr.alexander_poly(),       -t**2 + 4*t + 4*t**-1 - t**-2 - 6)
        self.assertEqual(self.L6a4.alexander_poly(),       -t**2 + 4*t + 4*t**-1 - t**-2 - 6)

        # method = 'snappy'
        self.assertEqual(self.Tref.alexander_poly(method='snappy'),       a**2 - a + 1)
        self.assertEqual(self.K3_1.alexander_poly(method='snappy'), a**2 - a + 1)
        self.assertEqual(self.K7_2.alexander_poly(method='snappy'),       3*a**2 - 5*a + 3)
        self.assertEqual(self.K8_3.alexander_poly(method='snappy'),       4*a**2 - 9*a + 4)
        self.assertEqual(self.K8_13.alexander_poly(method='snappy'),      2*a**4 - 7*a**3 + 11*a**2 - 7*a + 2)
        
    def testConnectedSum(self):
        repeat = 3
        while repeat > 0:
            k1 = self.random_knot()
            k2 = self.random_knot()
            Sum = k1.connected_sum(k2)
            self.assertEqual(Sum.alexander_poly(), k1.alexander_poly()*k2.alexander_poly())
            self.assertEqual(Sum.signature(),      k1.signature()+k2.signature())
            self.assertEqual(Sum.determinant(),    k1.determinant()*k2.determinant())
            repeat -=1

    def testSignature(self):
        self.assertEqual(abs(self.Tref.signature()),       2)
        self.assertEqual(abs(self.K3_1.signature()),       2)
        self.assertEqual(abs(self.K7_2.signature()),       2)
        self.assertEqual(abs(self.K8_3.signature()),       0)
        self.assertEqual(abs(self.K8_13.signature()),      0)
        self.assertEqual(abs(self.L2a1.signature()),       1)
        self.assertEqual(abs(self.L6a2.signature()),       3)
        self.assertEqual(abs(self.Borr.signature()),       0)
        self.assertEqual(abs(self.L6a4.signature()),       0)
        self.assertEqual(abs(self.L7a3.signature()),       3)

    def testCopy(self):
        repeat = 3
        while repeat > 0:
            k1 = self.random_knot()
            k1_prime = k1.copy()
            self.assertEqual(k1.determinant(),    k1_prime.determinant())
            self.assertEqual(k1.writhe(),         k1_prime.writhe())
            self.assertEqual(k1.signature(),      k1_prime.signature())
            self.assertEqual(k1.alexander_poly(), k1_prime.alexander_poly())
            repeat-=1

        repeat = 3
        while repeat > 0:
            k2 = self.random_link()
            k2_prime = k2.copy()
            self.assertEqual(k2.determinant(),    k2_prime.determinant())
            self.assertEqual(k2.writhe(),         k2_prime.writhe())
            self.assertEqual(k2.signature(),      k2_prime.signature())
            self.assertEqual(k2.alexander_poly(), k2_prime.alexander_poly())
            repeat-=1

    def testMirror(self):
        repeat = 3
        while repeat > 0:
            k1 = self.random_knot()
            k1_prime = k1.mirror()
            self.assert_(k1.signature() == -1*k1_prime.signature(), msg="knot signature failed for "+ repr(k1))
            self.assert_(k1.writhe() == -1*k1_prime.writhe(), msg="knot writhe failed for "+repr(k1))
            repeat-=1
            
        repeat = 3
        while repeat > 0:
            k2 = self.random_link()
            k2_prime = k2.mirror()
            self.assert_(k2.signature() == -1*k2_prime.signature(), msg="link signature failed for " + repr(k2))
            self.assert_(k2.writhe() == -1*k2_prime.writhe(), msg="link writhe failed for " + repr(k2))
            repeat-=1

    def testDet(self):
        self.assertEqual(self.K3_1.determinant(),                 3)
        self.assertEqual(self.K3_1.determinant(method='color'),   3)
        self.assertEqual(self.K3_1.determinant(method='goeritz'), 3)
        
        self.assertEqual(self.Tref.determinant(),                 3)
        self.assertEqual(self.Tref.determinant(method='color'),   3)
        self.assertEqual(self.Tref.determinant(method='goeritz'), 3)

        self.assertEqual(self.K7_2.determinant(),                 11)
        self.assertEqual(self.K7_2.determinant(method='color'),   11)
        self.assertEqual(self.K7_2.determinant(method='goeritz'), 11)

        self.assertEqual(self.K8_3.determinant(),                 17)
        self.assertEqual(self.K8_3.determinant(method='color'),   17)
        self.assertEqual(self.K8_3.determinant(method='goeritz'), 17)

        self.assertEqual(self.K8_13.determinant(),                 29)
        self.assertEqual(self.K8_13.determinant(method='color'),   29)
        self.assertEqual(self.K8_13.determinant(method='goeritz'), 29)

        self.assertEqual(self.L2a1.determinant(),                 2)
        self.assertEqual(self.L2a1.determinant(method='color'),   2)
        self.assertEqual(self.L2a1.determinant(method='goeritz'), 2)

        self.assertEqual(self.L6a2.determinant(),                 10)
        self.assertEqual(self.L6a2.determinant(method='color'),   10)
        self.assertEqual(self.L6a2.determinant(method='goeritz'), 10)

        self.assertEqual(self.Borr.determinant(),                 16)
        self.assertEqual(self.Borr.determinant(method='color'),   16)
        self.assertEqual(self.Borr.determinant(method='goeritz'), 16)

        self.assertEqual(self.L6a4.determinant(),                 16)
        self.assertEqual(self.L6a4.determinant(method='color'),   16)
        self.assertEqual(self.L6a4.determinant(method='goeritz'), 16)

        self.assertEqual(self.L7a3.determinant(),                 16)
        self.assertEqual(self.L7a3.determinant(method='color'),   16)
        self.assertEqual(self.L7a3.determinant(method='goeritz'), 16)

    def testBlackGraph(self):
        repeat = 3
        while repeat > 0:
            k1 = self.random_knot()
            self.assert_(k1.black_graph().is_planar())
            repeat-=1

        repeat = 3
        while repeat > 0:
            k2 = self.random_link()
            self.assert_(k2.black_graph().is_planar())
            repeat-=1


suite = unittest.TestLoader().loadTestsFromTestCase(TestLinkFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)
