import spherogram
import itertools


def test_link(link, sign):
    assert link.signature() == -sign
    assert all(c.sign == sign for c in link.crossings)


# From linkinfo, to make sure that we don't undo the progress of
# d852beb76a33798a6

L2a1_0 = spherogram.Link([[4, 1, 3, 2], [2, 3, 1, 4]])
L2a1_1 = spherogram.Link([[4, 2, 3, 1], [2, 4, 1, 3]])

test_link(L2a1_0, -1)
test_link(L2a1_1, 1)

base0 = [[3, 0, 2, 1], [1, 2, 0, 3]]
base1 = [[3, 1, 2, 0], [1, 3, 0, 2]]

for perm in itertools.permutations(range(4)):
    L = spherogram.Link([[perm[i] for i in PD] for PD in base0])
    test_link(L, -1)
    L = spherogram.Link([[perm[i] for i in PD] for PD in base0])
    test_link(L, -1)

for perm in itertools.permutations(range(4)):
    L = spherogram.Link([[perm[i] for i in PD] for PD in base1])
    test_link(L, 1)
    L = spherogram.Link([[perm[i] for i in PD] for PD in base1])
    test_link(L, 1)


# From issue 35

A = spherogram.Link([(1, 2, 0, 3), (2, 1, 3, 0)])
B = spherogram.Link([(1, 3, 0, 2), (3, 1, 2, 0)])
C = spherogram.Link([(1, 2, 0, 3), (3, 0, 2, 1)])

test_link(A, 1)
test_link(B, 1)
test_link(C, -1)


# A positive trefoil plus a 2-crossing linking component:
base2 = [(2, 0, 3, 7), (6, 2, 7, 1), (0, 6, 1, 5), (9, 4, 8, 3), (4, 9, 5, 8)]

# linking component reversed:
base3 = [(2, 0, 3, 7), (6, 2, 7, 1), (0, 6, 1, 5), (8, 3, 9, 4), (4, 9, 5, 8)]


for perm in itertools.permutations([3, 4, 5, 8, 9]):
    full_perm = list(range(10))
    for i, j in zip([3, 4, 5, 8, 9], perm):
        full_perm[i] = j
    assert len(set(full_perm)) == 10

    L = spherogram.Link([[full_perm[i] for i in PD] for PD in base2])
    assert L.signature() == -3
    assert all(c.sign == 1 for c in L.crossings)

    L = spherogram.Link([[full_perm[i] for i in PD] for PD in base3])
    assert L.signature() == -1
