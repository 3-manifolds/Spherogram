import spherogram as spherogram
from spherogram.graphs import FatGraph
from collections import defaultdict
from spherogram.dev.tangle_patch import *
from spherogram.links.simplify import reverse_type_I
import random

def faces(G):
    i_dict = G.incidence_dict

    #direction = (e,v) #on edge e looking towards v (i.e. left arm holding edge)
    #when you traverse the faces counterclockwise

    unseen_pairs = { (edge,vertex) for edge in G.edges for vertex in edge.incident_to()}
    face_list = []

    while len(unseen_pairs) > 0:

        direction = unseen_pairs.pop()
        face = [direction]
        next_edge = i_dict[direction[1]].succ(direction[0])
        vert_pair = next_edge.incident_to()
        next_direction = (next_edge, vert_pair[1 - vert_pair.index(direction[1])])
        while next_direction in unseen_pairs:
            unseen_pairs.remove(next_direction)
            face.append(next_direction)
            direction = next_direction

            next_edge = i_dict[direction[1]].succ(direction[0])
            vert_pair = next_edge.incident_to()
            next_direction = (next_edge, vert_pair[1 - vert_pair.index(direction[1])])

        face_list.append(face)

    return face_list


def link_diagram(G):
    face_list = faces(G)
    crossing_dict = {}

    for edge in G.edges:  # create one crossing per edge
        l = label(edge)
        crossing_dict[l] = spherogram.Crossing(l)

    for face in face_list:  # connect along faces
        for n, direction in enumerate(face):
            edge = direction[0]
            next_edge = face[(n + 1) % len(face)][0]
            c = crossing_dict[label(edge)]
            cnext = crossing_dict[label(next_edge)]
            o = open_overposition(c)
            u = open_underposition(cnext)
#            print('c,o: ')
#            print(c,o)
#            print('cnext,u: ')
#            print(cnext,u)
            c[o] = cnext[u]

    for edge in G.edges:  # switch twisted edges
        c = crossing_dict[label(edge)]
        if edge.twisted:
            c.rotate_by_90()

    return spherogram.Link(crossing_dict.values())


def label(edge):
    vs = edge.incident_to()
    s1 = vs[0]+'['+str(edge.slot(vs[0]))+']'
    s2 = vs[1]+'['+str(edge.slot(vs[1]))+']'
    return frozenset( (s1,s2) )


def open_underposition(crossing):
    if crossing.adjacent[0] is None:
        return 0
    elif crossing.adjacent[2] is None:
        return 2
    print('Both under positions occupied')
    return


def open_overposition(crossing):
    if crossing.adjacent[1] is None:
        return 1
    elif crossing.adjacent[3] is None:
        return 3
    print('No open over position')
    return


def link_to_fat_graph(link, alternate_around_crossing=False,
                      alternate_strand=False):
    G = FatGraph()
    B = link.black_graph()
    visited_slots = []
    for c in link.crossings:
        for i, adj in enumerate(c.adjacent):
            if (c.label, i) not in visited_slots:
                label = str(c.label)
                adjlabel = str(adj[0].label)
                if alternate_around_crossing:
                    G.add_edge((label, i), (adjlabel, adj[1]), i % 2)
                elif alternate_strand:
                    diff = (i * adj[1]) % 2
                    G.add_edge((label, i), (adjlabel, adj[1]), diff)
                else:
                    G.add_edge((label, i), (adjlabel, adj[1]), 0)
                visited_slots.append((c.label, i))
                visited_slots.append((adj[0].label, adj[1]))
    return G


def graph_knot_sequence(link, length, alternate_around_crossing=False,
                        alternate_strand=False):
    sequence = [link]
    link_copy = link.copy()
    for i in range(length):
        link_copy = link_diagram(link_to_fat_graph(link_copy,
                                                   alternate_around_crossing))
        sequence.append(link_copy)
    return sequence


def double(link):
    G = link_to_fat_graph(link)
    f = faces(G)
    colors = face_two_coloring(G, f)
    D = link_diagram(link_to_fat_graph(link.copy()))
    H = link_to_fat_graph(D)
    for edge in H.edges:
        edge.twisted = (edge.twisted + find_face_color(edge, colors)) % 2
    return link_diagram(H)


def group_by_quadruples(crossings):
    labels = defaultdict(list)
    for i,c in enumerate(crossings):
        label1, label2 = list(c.label)
        l1 = label1.split('\'')[1].split('[')[0], label1.split('\'')[3].split('[')[0]
        l2 = label2.split('\'')[1].split('[')[0], label2.split('\'')[3].split('[')[0]
        common = set(l1).intersection(set(l2)).pop()
        labels[common].append(i)
    return labels

def face_two_coloring(G,faces):
    colors = {}
    colors[tuple(faces[0])] = 0
    faces.pop(0)
    while faces:
        print(colors)
        print(faces)
        for f in faces:
            for other_face in colors:
                if share_edge(f,other_face):
                    colors[tuple(f)] = 1 - colors[other_face]
                    faces.remove(f)
                    break
    return colors

def share_edge(face1,face2):
    s1 = set([e for e,v in face1])
    s2 = set([e for e,v in face2])
    return len(s1.intersection(s2)) > 0


def edge_to_str_pair(edge):
    return set(map(lambda s: s.strip(), str(edge).split('---')))


def doubled_edge_to_str_pair(e):
    s1 = {e[0].split()[0].split('\'')[1], e[0].split()[1].split('\'')[1]}
    s2 = {e[1].split()[0].split('\'')[1], e[1].split()[1].split('\'')[1]}
    return [s1,s2]

def simple_cut(link, cs1, cs2):
    return Tangle(2,link.crossings,[cs1,cs2.opposite(),cs1.opposite(),cs2])

def four_connect(link1, css1, link2, css2):
    T1 = simple_cut(link1, css1[0], css1[1])
    T2 = simple_cut(link2, css2[0], css2[1])
    L = T1.circular_sum(T2,0)
    L._rebuild()
    return L

#find a pair of parallel strands in a doubled link
def identify_pair(link):
    f = link.faces()
    for face in f:
        if len(face) == 4:
            signs = [(cs.strand_index % 2) for cs in face]
            comps = [component_index(link,cs) for cs in face]
            if sorted(signs) == [0, 0, 1, 1] and len(set(comps)) > 1:
                second_zero = 3 - list(reversed(signs)).index(0)
                if comps[second_zero] != comps[(second_zero+2) % 4]:
                    return face[second_zero],face[(second_zero+2) % 4]

    raise Exception()


def build_face(cs):
    next_cs = cs.next_corner()
    f = [cs, next_cs]
    while next_cs != cs:
        next_cs = next_cs.next_corner()
        f.append(next_cs)
    f.pop()
    return f

def component_index(link,cs):
    cs = cs.oriented()
    for i,comp in enumerate(link.link_components):
        if cs in comp:
            return i


def whitehead_double(link):
    D = double(link)
    cs1, cs2 = identify_pair(D)
    T1 = simple_cut(D, cs1, cs2)
    T2 = clasp()
    clear_orientations(T1)
    clear_orientations(T2)
    L = T1.circular_sum(T2, 1)
    L._rebuild()
    return L


def make_writhe_zero(link):
    writhe = link.writhe()
    if writhe > 0:
        hand = 'left'
    else:
        hand = 'right'
    writhe = abs(writhe)
    for i in range(writhe):
        cs = random.choice(link.crossing_strands())
        reverse_type_I(link, cs, i, hand, False)
        link._rebuild()
        print(link.writhe())


def clear_orientations(link):
    for c in link.crossings:
        c._clear()


def clasp():
    c, d = spherogram.Crossing('c'),spherogram.Crossing('d')
    c[0] = d[1]
    c[1] = d[0]
    d2, d3 = d.crossing_strands()[2],d.crossing_strands()[3]
    c2, c3 = c.crossing_strands()[2],c.crossing_strands()[3]
    return Tangle(2,[c,d],[d2, d3, c3, c2])


def random_four_connect(link1, link2):
    from random import choice, sample
    css1 = sample(choice(link1.faces()), 2)
    css2 = sample(choice(link2.faces()), 2)
    return four_connect(link1, css1, link2, css2)


def find_face_color(e, colors):
    pair = doubled_edge_to_str_pair(e)
    for face in colors:
        edges = [edge_to_str_pair(edge) for edge, v in face]
        if pair[0] in edges and pair[1] in edges:
            return colors[face]


def relabel(link):
    for i, c in enumerate(link.crossings):
        c.label = i


def cycle(n, signs):
    G = FatGraph()
    for i in range(n):
        G.add_edge((str(i), 0), (str((i + 1) % n), 1), signs[i])
    return G
