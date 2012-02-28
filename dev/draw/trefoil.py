import os, sys, re, tempfile

def test_peer_code(code):
    fname = tempfile.mktemp()
    file = open(fname, 'w')
    file.write(code + "\nexit")
    file.close()
    os.popen("draw " + fname).read()
    return len(open('draw.out').readlines()) > 20 

def test(): 
    orders = [[1, 3, 5], [1, 5, 3], [3, 1, 5], [3, 5, 1], [5, 1, 3], [5, 3, 1]]
    signs = [ (i, j, k) for i in [1, -1] for j in [1, -1] for k in [1, -1]]
    for a, b, c in orders:
        for i, j, k in signs:
            s = "[%d %d %d]/# # #" % (i*a, j*b, c*k)
            print s.split('/')[0], test_peer_code(s)
    
test_peer_code('[-5 7 -1 3] / + + + +')
