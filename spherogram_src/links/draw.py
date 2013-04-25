from __future__ import print_function
"""
Plotting pictures of links via Andrew Bartholomew's "Draw programme"
http://www.layer8.co.uk/maths/draw/index.htm
"""

import os, sys, re, tempfile, subprocess
from subprocess import PIPE

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
        os.chdir(curr_dir)
        raise ValueError('draw failed: ' + ans + 'for ' + peer_code)
    run_silent('mpost', 'link.mps')
    run_silent('env', 'epstopdf', 'link.1')
    data = open('link.pdf', 'rb').read()
    os.chdir(curr_dir)
    os.system('rm -rf ' + tmp_dir)
    return data

def test():
    import spherogram.links.test
    knots = spherogram.links.test.some_knots()
    tmp_dir = tempfile.mkdtemp()
    files = []
    for name, K in knots:
        pc = K.peer_code()
        try:
            filename = tmp_dir + os.sep + name + '.pdf'
            open(filename, 'wb').write( link_pdf(pc) )
            files.append(filename)
        except ValueError:
            print("Problem with:",  name, pc)

    os.system("open -a Preview.app " + tmp_dir + "/*.pdf")

if __name__ == '__main__':
    test()

