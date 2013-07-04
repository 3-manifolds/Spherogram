import snappy, orthogonal, plink, hobby

L = orthogonal.link_from_manifold(snappy.Manifold('16n12345'))
M = snappy.Manifold()
M.LE.load_from_spherogram(L)


def my_polylines(LE):
    segments = M.LE.polylines(gapsize=0)[0][0]
    starts = {segment[-1]:segment for segment in segments}
    ans = starts.pop(segments[0][0])
    while len(starts):
        ans += starts.pop(ans[-1])
    return [vector(RR, k) for k in ans[:-1]]
    
#polylines = my_polylines(M.LE)

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


sq = matrix(RR, [(0,0), (1,0), (1,1), (0, 1)]).rows()

def smooth(polylines):
    return sum( [ smooth_segment2(s) for s in polylines ], Graphics())

def hobby_smooth(points):
    dirs = hobby.hobby_cyclic_dirs(points)
    print dirs
    path = []
    for k in range(len(points)-1):
        (c0, c1) = hobby.good_bezier(points[k],dirs[k],dirs[k+1],points[k+1])
        if k == 0:
            path.append([points[0], c0, c1, points[1]])
        else:
            path.append([c0, c1, points[k+1]])
        (c0, c1) = hobby.good_bezier(points[-1], dirs[-1], dirs[0], points[0])
        path.append([c0, c1, points[0]])
    return path
            

def test_gb(p1, a1, a2, p2, t1=1.0, t2=1.0):
    c0, c1 = hobby.good_bezier(p1, a1, a2, p2, t1, t2)

    G = bezier_path([[p1, c0, c1, p2]]) + point(c0) + point(c1)
    G.set_aspect_ratio(1.0)
    return G

def test_hd(points, t):
    #points = matrix(RR, [(0,0), (1, 1), (2, -2), (3, 0), (2, 3), (0, 2)]).rows()
    dirs = hobby.hobby_dirs(points, 0, RR.pi()/2)
    print dirs
    path = []
    for k in range(len(points)-1):
        (c0, c1) = hobby.good_bezier(points[k],dirs[k],dirs[k+1],points[k+1], t, t)
        if k == 0:
            path.append([points[0], c0, c1, points[1]])
        else:
            path.append([c0, c1, points[k+1]])
    
    #+ line(sum(path, []), color='red')
    #G = line(points)  + point2d(points) + bezier_path(path, color='green')
    G = point2d(points) + bezier_path(path, color='green')
    G.set_aspect_ratio(1.0)
    return G
    
