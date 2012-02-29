"""
Plotting pictures of links via Andrew Bartholomew's "Draw programme"
http://www.layer8.co.uk/maths/draw/index.htm
"""

import os, sys, re, tempfile, subprocess
from subprocess import PIPE
from tangles import RationalTangle


def knot(fractions):
    if len(fractions) == 1:
        return RationalTangle(*fractions[0]).denominator_closure()
    else:
        A, B, C = [RationalTangle(*f) for f in fractions]
        T = A + B + C
        return T.numerator_closure()


def create_file_with_contents(file_name, data):
    file = open(file_name, 'w')
    file.write(data)
    file.close()

def run_draw(peer_code, infile, outfile, *options):
    create_file_with_contents(infile, peer_code + '\nexit')
    p = subprocess.Popen( ('draw', infile, outfile) + options, stdout=PIPE, stderr=PIPE)
    return p.stderr.read()

def run_silent(*command):
    subprocess.Popen(command, stdout=PIPE, stderr=PIPE).wait()
    
def link_pdf(peer_code):
    curr_dir = os.path.abspath(os.path.curdir)
    tmp_dir = tempfile.mkdtemp()
    os.chdir(tmp_dir)
    ans = run_draw(peer_code, 'peer_code', 'link.mps', '--pen-size=4', '--disc-size=40')
    if len(ans):
        raise ValueError('draw failed: ' + ans + 'for ' + peer_code)
    run_silent('mpost', 'link.mps')
    run_silent('epstopdf', 'link.1')
    data = open('link.pdf', 'rb').read()
    os.system('rm -rf ' + tmp_dir)
    os.chdir(curr_dir)
    return data

def test():
    import ntools
    knots = ntools.DataInFile('montesinos_knots')
    for K, d in knots:
        L = knot(d)
        pc = L.peer_code()
        try:
            open('/tmp/knots/' + K + '.pdf', 'wb').write( link_pdf(pc) )
        except ValueError:
            print K, d, pc

def test2():
    for i in range(10):
        L = knot([(8,23)])
        pc = L.peer_code()
        try:
            open('/tmp/knots/' + '8_23' + '.pdf', 'wb').write( link_pdf(pc) )
        except ValueError:
            print pc

def test3():
    for i in range(10):
        L = knot([(6,29)])
        pc = L.peer_code()
        print pc
        open('/tmp/knots/' + '10_8' + '.pdf', 'wb').write( link_pdf(pc) )
        os.system('open /tmp/knots/10_8.pdf')
        raw_input('hit any key')


