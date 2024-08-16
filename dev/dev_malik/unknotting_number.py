import spherogram
from spherogram.links.random_links import map_to_link, random_map
from random import choice, randint


def doubly_connected_crossing(link):
    doubly_connected = None
    for c in link.crossings:
        num_neighbors = len({x[0] for x in c.adjacent})
        if num_neighbors == 2:
            return c
        elif num_neighbors == 3:
            doubly_connected = c
    return doubly_connected


def unknot_sequence_optimized(link, simp_mode='global'):
    link.simplify(mode=simp_mode)
    switching_count = 0
    while len(link):
        most_crossings_removed = 0
        best_crossing = None
        for i, c in enumerate(link.crossings):
            link_copy = link.copy()
            num_crossings_before = len(link_copy)
            crossing_copy = link_copy.crossings[i]
            switch_crossing(link_copy, crossing_copy)
            link_copy.simplify(mode=simp_mode)
            num_crossings_after = len(link_copy)
            crossings_removed = (num_crossings_before - num_crossings_after)
            if crossings_removed > most_crossings_removed:
                best_crossing = c
                most_crossings_removed = crossings_removed
        switch_crossing(link, best_crossing)
        link.simplify(mode=simp_mode)
        switching_count += 1
    return switching_count


def min_unknotting_optimized(link, simp_mode='global', num_attempts=10, backtrack_lower=10, backtrack_upper=50):
    known_lower_bound = abs(link.signature()) / 2
    #    print('Known lower bound: ' + str(known_lower_bound))
    unknotting_num = len(link) / 2
    for i in range(num_attempts):
        print(float(i) / num_attempts * 100)
        link_copy = link.copy()
        link_copy.backtrack(randint(backtrack_lower, backtrack_upper))
        link_copy.simplify(mode=simp_mode)
        x = unknot_sequence_optimized(link_copy)
        if x < unknotting_num:
            unknotting_num = x
            if unknotting_num == known_lower_bound:
                break
    return (known_lower_bound, unknotting_num)


def unknot_sequence(link, simp_mode='level'):
    link.simplify(mode=simp_mode)
    switching_count = 0
    while len(link):
        doubly_connected = doubly_connected_crossing(link)
        if doubly_connected:
            switch_crossing(link, doubly_connected)
        else:
            switch_crossing(link, choice(link.crossings))
            print('random')
        link.simplify(mode=simp_mode)
        switching_count += 1
    return switching_count


def unknot_sequence(link, simp_mode='level'):
    # duplicate function, why ?
    link.simplify(mode=simp_mode)
    switching_count = 0
    while len(link):
        switch_crossing(link, choice(link.crossings))
        link.simplify(mode=simp_mode)
        switching_count += 1
    return switching_count


def min_unknotting_sequence(link, simp_mode='level', num_attempts=100,
                            backtrack_lower=10, backtrack_upper=40):
    known_lower_bound = abs(link.signature()) / 2
    # print('Known lower bound: ' + str(known_lower_bound))
    unknotting_num = len(link) / 2
    for i in range(num_attempts):
        # print(float(i)/num_attempts*100)
        link_copy = link.copy()
        link_copy.backtrack(randint(backtrack_lower, backtrack_upper))
        link_copy.simplify(mode=simp_mode)
        x = unknot_sequence(link_copy)
        if x < unknotting_num:
            unknotting_num = x
            if unknotting_num == known_lower_bound:
                break
    return (known_lower_bound, unknotting_num)


def switch_crossing(link, crossing):
    crossing.rotate(1)
    link._rebuild()


def random_alternating_knot(n):
    while True:
        L = map_to_link(random_map(n, edge_conn_param=4))
        if len(L.link_components) == 1:
            break
    return L.alternating()


def mean_unknotting_num(cns, num_links_per_crossing,
                        simp_mode='level', num_attempts=100,
                        backtrack_lower=10, backtrack_upper=40):
    means = []
    for i in cns:
        print(i)
        sig, un = 0, 0
        for j in range(num_links_per_crossing):
            s, u = min_unknotting_sequence(random_alternating_knot(i),
                                           simp_mode=simp_mode,
                                           num_attempts=num_attempts,
                                           backtrack_lower=backtrack_lower,
                                           backtrack_upper=backtrack_upper)
            sig += s
            un += u
            print(s, u)
        means.append((float(sig) / num_links_per_crossing,
                      float(un) / num_links_per_crossing))
    return means
