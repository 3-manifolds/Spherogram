from spherogram.links.alexander import Exhaustion
import sympy

def kauffman_bracket(self):
    n = len(self.link)
    PD = self.link.PD_code()

    A, B = sympy.symbols('A,B')
    web = 1
    for c in self.crossings:
        print('web before')
        print(web)

        a,b,c,d = c.strand_labels

        web *= A*(strand_symb(a,d)*strand_symb(b,c)) + (1/A)*(strand_symb(a,b)*strand_symb(c,d))
        web = web.expand()
        print('web about to simplify')
        print(web)

        web = simplify_jones_expression(web)

        print('web simplified')
        print(web)

    return web.expand()

def jones_polynomial(self):
    B = self.kauffman_bracket()
    w = self.link.writhe()
    A,q = sympy.symbols('A,q')
    return ((-A**3)**w*B/(-A**2-(1/A)**2)).subs(A,q**(1/4))

def test_kauffman():
    from spherogram import Link
    L = Link('L11a548')
    E = good_exhaustion(L)
    return E.kauffman_bracket()

def test_kauffman_big():
    from spherogram import random_link
    L = random_link(50,alternating=True)
    E = good_exhaustion(L)
    return E.kauffman_bracket()


def test_jones():
    from spherogram import Link
    L = Link('L11a548')
    E = good_exhaustion(L, max_failed_tries=100)
    print(E.jones_polynomial())

def strand_symb(a,b):
    a, b = sorted([a,b])
    return sympy.Symbol('P'+str(a)+'c'+str(b))

def test_simplify():
    A,B,P0c0,P0c1,P1c2,P3c5,P3c4,P20c100 = sympy.symbols('A,B,P0c0,P0c1,P1c2,P3c5,P3c4,P20c100')
    poly = P0c0+P0c0*P0c1+P0c1*P1c2+P3c5*P3c4+P3c4*P20c100
    print(poly)
    simplify_jones_expression(poly)

    P13c19,P14c18,P14c19,P15c20 = sympy.symbols('P13c19,P14c18,P14c19,P15c20')
    monomial = P13c19*P14c18*P14c19*P15c20
    print(simplify_monomial(monomial))

def simplify_jones_expression(poly):
    monomials = poly.args
    for monomial in monomials[:]:
        print('monomial:')
        print(monomial)

#        print('monomials:')
#        print(monomials)

        simplified_mon = simplify_monomial(monomial)
#        print('simplified:')
#        print(simplified_mon)
        poly = poly.subs(monomial,simplified_mon)
#        print(poly)
    return poly


def simplify_monomial(monomial):
    monomial = remove_squares(monomial)
    monomial = remove_loops(monomial)
    A = sympy.Symbol('A')
    strand_vars = monomial.free_symbols - set(sympy.symbols('A,B'))
    strand_labels = all_labels(strand_vars)

    while len(set(strand_labels)) < len(strand_labels):
#        collections.Counter
#        print(monomial)
        for l in strand_labels:
            if strand_labels.count(l) > 1:
                matching_vars = [v for v in strand_vars if l in var_to_strand_labels(v)]
                v1, v2 = matching_vars[0], matching_vars[1]
                new_v = combine_strands(v1,v2, l)
                monomial = monomial.subs(v1*v2,new_v)
                monomial = monomial.subs(new_v*new_v,-A**2-(1/A)**2)
                sl = var_to_strand_labels(new_v)
                if sl[0] == sl[1]:
                    monomial = monomial.subs(new_v,-A**2-(1/A)**2)
                break

#        monomial = remove_squares(monomial)
#        monomial = remove_loops(monomial)

        strand_vars = monomial.free_symbols - {A}
        strand_labels = all_labels(strand_vars)

    return monomial

def all_labels(strand_vars):
    strand_labels = []
    for v in strand_vars:
        strand_labels.extend(var_to_strand_labels(v))
    return strand_labels


def combine_strands(v1, v2, common_label):
    labels1 = var_to_strand_labels(v1)
    labels2 = var_to_strand_labels(v2)
    l1 = labels1[1 - labels1.index(common_label)]
    l2 = labels2[1 - labels2.index(common_label)]
    l1, l2 = sorted([l1,l2])
    return sympy.Symbol('P'+str(l1)+'c'+str(l2))


def remove_squares(monomial):
    A, B = sympy.symbols('A,B')
    strand_vars = monomial.free_symbols - {A,B}
    for v in strand_vars:
        monomial = monomial.subs(v*v,-A**2-(1/A)**2)
    return monomial


def remove_loops(monomial):
    A, B = sympy.symbols('A,B')
    strand_vars = monomial.free_symbols - {A,B}
    for v in strand_vars:
        l1, l2 = var_to_strand_labels(v)
        if l1 == l2:
            monomial = monomial.subs(v,-A**2-(1/A)**2)
    return monomial


def var_to_strand_labels(v):
    return map(int, str(v)[1:].split('c'))


def join_strands(v1, v2):
    pass


"""
def jones_ring(n):
    names_dict = {(i,j):'P'+str(i)+'c'+str(j) for i in range(2*n) for j in range(i,2*n)}
    names_list = list(names_dict.values())
    names_list.append('A')
    names_list.append('B')
    R = PolynomialRing(QQ,names_list)
    return R

def condense_polynomial(poly, A, B, R, gd):
    while True:
        no_simplifications = True
        for monomial in poly.monomials():
            print(poly)
            poly, fs = condense_term(monomial,poly,A,B, R, gd)
            if fs:
                no_simplifications = False
        if no_simplifications:
            break

    return poly

def condense_term(monomial, poly, A, B, R, gd):
    variables = list(monomial.variables())
    degrees = monomial.degrees()
    gens = R.gens()
    found_simplification = False

    if A in variables:
        variables.remove(A)
    if B in variables:
        variables.remove(B)
    var_labels = [map(int,str(v)[1:].split('c')) for v in variables]

    seen_var_indices = []
    for i,d in enumerate(degrees):
        if d == 2:
            v = gens[i]
            if v not in [A,B]:
                poly = poly.subs({v*v: -A*A-B*B})
                found_simplification = True
                seen_var_indices.append( variables.index(v) )

    for index in sorted(seen_var_indices,reverse=True):
        variables.pop(index)
        var_labels.pop(index)
    seen_var_indices = []


    for i,v in enumerate(var_labels):
        if v[0] == v[1]:
            poly = poly.subs({variables[i]:-A*A-B*B})
            found_simplification = True
            seen_var_indices.append(i)

    for index in sorted(seen_var_indices,reverse=True):
        variables.pop(index)
        var_labels.pop(index)

    all_strand_labels = []
    for vl in var_labels:
        all_strand_labels.extend(vl)

    found_matching_pair = False
    for label in all_strand_labels:
        if all_strand_labels.count(label) > 1:
            found_matching_pair = True
            found_simplification = True
            matching_label = label
            break
    if found_matching_pair:
        matched_indices = [i for (i, pair) in enumerate(var_labels) if matching_label in pair]
        p1, p2 = var_labels[matched_indices[0]], var_labels[matched_indices[1]]
        v1, v2 = variables[matched_indices[0]], variables[matched_indices[1]]
        l1 = p1[ 1 - p1.index(matching_label) ]
        l2 = p2[ 1 - p2.index(matching_label) ]
        l1, l2 = sorted([l1, l2])
        new_variable = gd['P'+str(l1)+'c'+str(l2)]
        print(v1,v2,new_variable)
        print('poly before')
        print(poly)
        poly = poly.subs({v1*v2:new_variable})
        print('poly after')
        print(poly)
    return poly, found_simplification

def strand_contract(v1, v2):
    s = set(v1).intersection(set(v2)).pop()
    if s:
        return s.pop()


def test_jones():
    R = jones_ring(5)
    gd = R.gens_dict()
    A = gd['A']
    B = gd['B']
    poly = 2*gd['P0c0']*gd['P1c1']+gd['P1c2']*gd['P2c3']-A*gd['P1c4']*gd['P2c3']+B*gd['P4c5']*gd['P3c4']+gd['P0c1']*gd['P0c1']
    poly = condense_polynomial(poly,A,B, R, gd)

"""
