from link import *
from tangles import *

def figure8():
    a, b, c, d = [Crossing(x) for x in 'abcd']
    a[0] = c[1]
    a[1] = d[0]
    a[2] = b[1]
    a[3] = b[0]
    b[2] = d[3]
    b[3] = c[2]
    c[3] = d[2]
    c[0] = d[1]
    return Link([a,b,c,d])

def punct_torus():
    a = Crossing('a')
    a[0] = a[2]
    a[1] = a[3]
    return Link([a], check_planarity=False)

def whitehead():
    a, b, c, d, e =  crossings = [Crossing(x) for x in 'abcde']
    a[0] = b[3]
    a[1] = d[0]
    a[2] = d[3]
    a[3] = c[0]
    b[0] = c[3]
    b[1] = e[2]
    b[2] = e[1]
    c[1] = d[2]
    c[2] = e[3]
    d[1] = e[0]
    return Link(crossings)

K, W, T = figure8(), whitehead(), punct_torus()
print K.is_planar(), W.is_planar(), punct_torus().is_planar()
print K.PD_code(True)
print W.PD_code(True)
print K.exterior().volume(), W.exterior().volume()



