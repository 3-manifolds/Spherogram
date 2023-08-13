import spherogram as sg
import random


def random_petaluma_diagram(size):
    if size % 2:
        return petaluma_knot(permutation(size))
    else:
        print("Size must be even")


def petaluma_knot(height_perm):
    size = len(height_perm)
    crossing_dict = {}
    visited_dict = {}
    for i in range(size-1):
        for j in range(i+1,size):
            crossing_dict[ (i,j) ] = sg.Crossing('c'+str(i)+str(j))
            visited_dict[ (i,j) ] = False
    end_position = 2
    if height_perm[1] < height_perm[0]:
        end_position = 3
    old_crossing = crossing_dict[0,1]
    visited_dict[0,1] = True

    for i in range(size):
        for j in strands_to_cross(size,i):
            if (i,j) == (0,1):
                continue
            a,b = sorted([i,j]) #must be ordered
            next_crossing = crossing_dict[a,b]
            ordered_or_backward = 0
            if i > j:
                ordered_or_backward = 1
            if height_perm[i] < height_perm[j]: # strand i passes under strand j
                if not visited_dict[a,b]: #if first time crossing has come up
                    next_crossing[0] = old_crossing[end_position]
                    end_position = 2
                else: #crossing has been seen before
                    if (b-a) % 2 == 1: #strand j should cross strand i right to left
                        next_crossing[2] = old_crossing[end_position]
                        end_position = 0
                    else: #left to right
                        next_crossing[0] = old_crossing[end_position]
                        end_position = 2
            else: #strand i passes over strand j
                if not visited_dict[a,b]:
                    next_crossing[1] = old_crossing[end_position]
                    end_position = 3
                else: #crossing has been seen before
                    if (b-a) % 2 == 1: #right to left
                        next_crossing[1] = old_crossing[end_position]
                        end_position = 3
                    else: #left to right
                        next_crossing[3] = old_crossing[end_position]
                        end_position = 1
            visited_dict[a,b] = True
            old_crossing = next_crossing

    first = crossing_dict[0, 1]
    last = crossing_dict[size-2, size-1]
    first_open = first.adjacent.index(None)  # last open spot
    last_open = last.adjacent.index(None)
    first[first_open] = last[last_open]
    return sg.Link(crossing_dict.values())


def strands_to_cross(size, i):
    mid = (size - 1) / 2
    L = []
    for j in range(mid):
        L.insert(0, (i-(2*(j+1))) % size)
        L.append((i+(2*(j+1))) % size)
    return L


def permutation(size):
    L = range(size)
    random.shuffle(L)
    return L
