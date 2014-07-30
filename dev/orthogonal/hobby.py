"""
Smoothing curves using Hobby's splines.  Mostly copied from code of
Dylan Thurston.
"""

from sage.all import vector, sqrt, RR, cos, sin
import math

def solve_tri(a, b, c, v):
    """
    Solve the tridiagonal matrix a[i-1]*w[i-1]+b[i]*w[i]+c[i]*w[i+1] = v[i]
    """
    if len(a) == 0:
        return [v[0]/b[0]]
    else:
        b[1] = b[1] - c[0]*a[0]/b[0]
        v[1] = v[1] - v[0]*a[0]/b[0]
        w = solve_tri(a[1:], b[1:], c[1:], v[1:])
        w.insert(0, (v[0]-c[0]*w[0])/b[0])
        return w
    
def solve_cyclic (a, b, c, v):
    if (len(a) == 2):
        a[0] = a[0] + c[1]
        c[0] = c[0] + a[1]
        return solve_tri(a[0:1], b, c[0:1], v)
    else:
        [a0, a1] = a[0:2]
        [c0, c1] = c[0:2]
        [b0, b1] = b[0:2]
        [v0, v1] = v[0:2]
        a[1] = -a0*a1/b1
        b[1] = b0 - c0*a0/b1
        b[2] = b[2] - c1*a1/b1
        c[1] = -c0*c1/b1
        v[1] = v0 - v1*c0/b1
        v[2] = v[2] - v1*a1/b1
        w = solve_cyclic(a[1:],b[1:],c[1:],v[1:])
        w.insert(1, (v1-a0*w[0] - c1*w[1])/b1)
        return w

def polar(r, phi):
    return vector(RR, [r*cos(phi), r*sin(phi)])

def to_polar((x,y)):
    r = sqrt(x**2+y**2)
    phi = 0 if r == 0 else RR(math.atan2(y,x))
    return r, phi
                  
    
def good_bezier(p1, angle1, angle2, p2, tension1=1.0, tension2=1.0):
    """
    Compute a nice curve from p1 to p2 with specified tangents
    """
    l, psi = to_polar(p2-p1)
    ctheta, stheta = polar(1.0, angle1-psi)
    cphi,   sphi   = polar(1.0, psi-angle2)
    a = sqrt(2.0)
    b = 1.0/16.0
    c = (3.0 - sqrt(5.0))/2.0
    alpha = a*(stheta - b*sphi) * (sphi - b*stheta) * (ctheta - cphi)
    rho = (2 + alpha) / (1 + (1-c)*ctheta + c*cphi) / tension1
    sigma = (2 - alpha) / (1 + (1-c)*cphi + c*ctheta) / tension2
    return (p1 + polar(l*rho/3, angle1), p2 - polar(l*sigma/3, angle2))

def mock_curvature (dist, angle1, angle2, tension1=1.0, tension2=1.0):
    """
    Compute an affine linear approximation to the curvature of good_bezier
    at endpoints
    """
    a = tension1 ** 2/dist * (3-1/tension2)
    b = -(tension1 ** 2/dist)*(1/tension2)
    c = -(tension2 ** 2/dist)*(1/tension1)
    d = (tension2 ** 2/dist)*(3-1/tension1)
    return ((a, -b, -c, d), (a*angle1 - b*angle2, - c*angle1 + d*angle2))

def accum_matrix (ml):
    l = []
    d = [0]
    u = []
    v = [0]
    for m in ml:
        d[-1] = d[-1] +  m[0][0]
        u.append(m[0][1])
        l.append(m[0][2])
        d.append(m[0][3])
        v[-1] = v[-1] + m[1][0]
        v.append(m[1][1])
    return (l,d,u,v)

def accum_cyclic_matrix(ml):
    (l,d,u,v) = accum_matrix(ml)
    d[0] = d[0] + d[-1]
    v[0] = v[0] + v[-1]
    return (l,d[:-1],u,v[:-1])

def normalize_angle(phi):
    pi = math.pi
    if phi < -pi or phi > pi:
        return float(phi + pi) % float(2*pi) - pi
    else:
        return phi

def hobby_dirs(points, tangent1 = None, tangent2 = None, curl1 = 1.0, curl2 = 1.0):
    """
    Compute the optimal angles for the tangents, using specified boundary
    conditions
    """
    if tangent1 is None:
        curl1 = 1.0
    if tangent2 is None:
        curl2 = 1.0
    matrices = []
    phi = 0
    for i in range(len(points)-1):
        (dist, theta) = to_polar(points[i+1] - points[i])
        phi = phi + normalize_angle(theta - phi)
        matrices.append(mock_curvature(dist, phi, phi))
    tangent2 = phi + normalize_angle(tangent2 - phi)
    ((a,b,c,d),(v,w)) = matrices[0]
    if tangent1 is not None:
        # Change equations to set initial tangent
        matrices[0] = ((1,0,c,d),(tangent1,w))
    else:
        # Curvature at first vertex is (1/curl1) * curvature at second
        matrices[0] = ((c+curl1*a, d+curl1*b, c, d), (w+curl1*v, w))
    ((a,b,c,d),(v,w)) = matrices[-1]
    if tangent2 is not None:
        # Change equations to set final tangent
        matrices[-1] = ((a,b,0,1),(v,tangent2))
    else:
        # Curvature at final vertex is (1/curl2) * curvature at previous
        matrices[-1] = ((a, b, a+curl2*c, b+curl2*d), (v, v+curl2*w))
    tridiag = accum_matrix(matrices)
    return apply(solve_tri,tridiag)

def hobby_cyclic_dirs(points):
    matrices = []
    phi = 0
    phis = []
    for i in range(len(points)-1):
        (d, theta) = to_polar(points[i+1] - points[i])
        phi = phi + normalize_angle(theta - phi)
        phis.append(phi)
        matrices.append(mock_curvature(d, phi, phi))
    (d,theta) = to_polar(points[0] - points[-1])
    phi = phi + normalize_angle(theta - phi)
    phi0 = phis[0] + normalize_angle(phi - phis[0])
    matrices.append(mock_curvature(d, phi, phi0))
    tridiag = accum_cyclic_matrix(matrices)
    return apply(solve_cyclic,tridiag)

