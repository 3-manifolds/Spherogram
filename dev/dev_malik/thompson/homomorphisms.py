import snappy
from thompson import *
from graphknotter import *
import itertools

a = TreeMap( TreeSequence([DyadicRational(3,2)]) , TreeSequence([DyadicRational(1,2)]) )
b = TreeMap( TreeSequence([DyadicRational(7,3)]) , TreeSequence([DyadicRational(5,3)]) )
A = a.inverse()
B = b.inverse()

ball_size = 11

def AB_ball(max_size):
    spheres = [AB_sphere(size) for size in range(1,max_size)]
    return itertools.chain(*spheres)


def AB_sphere(size):
    letters = 'abAB'
    all_words_itertool = itertools.product(letters, repeat=size)
    all_words = [''.join(word) for word in all_words_itertool]
    return filter(is_reduced, all_words)


def is_reduced(word):
    for i in range(len(word) - 1):
        u, v = word[i], word[i + 1]
        if u == 'a' and v == 'A':
            return False
        if u == 'A' and v == 'a':
            return False
        if u == 'b' and v == 'B':
            return False
        if u == 'B' and v == 'b':
            return False
    return True


def word_image(word, image_a, image_b):
    images = []
    for letter in word:
        if letter == 'a':
            image = image_a
        if letter == 'b':
            image = image_b
        if letter == 'A':
            image = inverse(image_a)
        if letter == 'B':
            image = inverse(image_b)
        images.append(image)
    return ''.join(images)

def inverse(word):
    return word[::-1].swapcase()

def link(word):
    return link_diagram(ABWord(word,a,b).evaluate().planar_graph())
"""
f = open('knotted_words.txt','r')
split_lines = [line.split() for line in f]
f.close()
powers = []
n = 0
for word, cn in split_lines:
    print(n)
    n += 1
    cn_powers = []
    try:
        for i in range(10):
            K = link(word*(i+1))
            K.simplify(mode='global')
            cn_powers.append(len(K))
    except:
        continue
    powers.append( [word, cn_powers] )

g = open('powers.txt','w')
g.write(str(powers))
g.close()
"""
"""
knots = []
f = open('knotted_words.txt','w')
for word in AB_ball(ball_size):
    M = ABWord(word,a,b).evaluate()
    G = M.planar_graph()
    K = link_diagram(G)
    try:
        K.simplify(mode='global')
    except:
        continue
    if len(K)>0:
#        knots.append([K,word])
        f.write(word+' ' + str(len(K)) +'\n')
f.close()
"""
"""
M = snappy.Manifold('m004')
G = M.fundamental_group()
if G.num_generators() == 2:
    possible_word_pairs = itertools.product(AB_ball(ball_size),repeat=2)
    for image_a, image_b in possible_word_pairs:
        not_homomorphism = False
        for rel in G.relators():
            image = word_image(rel,image_a,image_b)
            M = ABWord(image,a,b).evaluate()
            M.reduce()
            if M.size() != 2:
                not_homomorphism = True
                break
        if not not_homomorphism:
            print('image_a: '+image_a)
            print('image_b: '+image_b)
            print('Manifold: '+M.name())
"""


"""
for i in range(10):
    print('i: '+str(i))
    M = snappy.OrientableClosedCensus[i]
    G = M.fundamental_group()
    if G.num_generators() == 2:
        possible_word_pairs = itertools.product(AB_ball(ball_size),repeat=2)
        for image_a, image_b in possible_word_pairs:
            not_homomorphism = False
            for rel in G.relators():
                image = word_image(rel,image_a,image_b)
                M = ABWord(image,a,b).evaluate()
                M.reduce()
                if M.size() != 2:
                    not_homomorphism = True
                    break
            if not not_homomorphism:
                print('image_a: '+image_a)
                print('image_b: '+image_b)
                print('Manifold: '+M.name())
"""
