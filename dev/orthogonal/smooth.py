import snappy, orthogonal, plink

L = orthogonal.link_from_manifold(snappy.Manifold('16n12345'))
M = snappy.Manifold()
M.LE.load_from_spherogram(L)
polylines = [ [vector(v) for v in segment] for segment in  M.LE.polylines()[0][0]]

def smooth_segment(seg):
    if len(seg) <= 4:
        return bezier_path([seg])
    
    # If we have more points we break things by inserting
    # vertices into the *middle* of certain edges.

    a, b, c, d = seg[:4]
    n = (c+d)/2
    return bezier_path( [ [a, b, c, n]]) + smooth_segment([n] + seg[3:])

def smooth(polylines):
    return sum( [ smooth_segment(s) for s in polylines ], Graphics())

def smooth_segment2(seg):
    if len(seg) <= 4:
        return bezier_path([seg])
    
    # If we have more points we break things by inserting
    # vertices into the *middle* of certain edges.

    a, b, c, d = seg[:4]
    n = (c+d)/2
    return bezier_path( [ [a, b, c, n]]) + smooth_segment([n] + seg[3:])

def smooth(polylines):
    return sum( [ smooth_segment2(s) for s in polylines ], Graphics())


