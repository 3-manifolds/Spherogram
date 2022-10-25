import spherogram
import string

# Unknots from

unknot_28 = [1,-4,-3,6,5,-2,-7,8,4,-5,-9,10,2,-1,-11,7,12,-13,-6,3,14,-12,
             -10,9,13,-14,-8,11,15,-18,-17,20,19,-16,-21,17,22,-23,-20,21,
             24,-15,-25,26,16,-19,-27,28,18,-24,-26,27,23,-22,-28,25]


unknot_43 = [-1,2,-3,-4,-5,6,-7,8,-9,10,-11,12,-13,3,14,-15,16,-17,18,
             -19,20,1,-21,22,-23,24,-25,7,-26,27,-10,28,-29,13,4,30,31,
             -32,17,-33,34,35,-2,-36,-22,-37,-38,25,-8,39,-27,11,-40,
             29,41,23,37,-31,42,-16,33,-43,19,-20,-35,-14,-30,5,
             -6,38,-24,9,-39,26,-12,40,-28,-41,36,21,32,-42,15,-34,43,-18]

unknot_78 = [1,2,-3,-4,5,6,-7,-8,9,10,-2,-11,12,13,-14,-5,15,16,17,18,-19,
            -20,21,22,-6,-23,24,25,-26,-1,27,28,-29,-30,23,14,-31,-32,11,
            26,-33,-34,30,7,-22,-15,4,31,-13,-24,34,36,-28,-9,20,-17,37,
            38,35,39,-40,41,42,-43,-44,45,46,-47,-48,49,-39,-50,51,-52,
            -53,54,43,-55,-56,48,57,-58,-38,-61,-59,60,52,-62,-41,56,63,
            -46,-64,65,58,-66,-49,40,67,-51,-68,59,69,-70,-45,71,55,-42,-72,
            53,74,-69,-73,-37,-65,75,47,-63,-71,44,76,-74,-60,68,-77,-35,66,
            -57,-75,64,70,-76,-54,72,62,-67,50,77,61,73,-18,19,-10,-27,78,33,
            -25,-12,32,3,-16,-21,8,29,-36,-78]


def from_gauss_code(code):
    """
    This is the basic unsigned/unoriented variant.

    >>> L = from_gauss_code(unknot_28)
    >>> L.simplify('pickup')
    True
    >>> L
    <Link: 0 comp; 0 cross>
    >>> L.unlinked_unknot_components
    1
    """
    n = len(code)
    labels = list(range(1, n + 1))
    code_to_label = dict(zip(code, labels))
    label_to_code = dict(zip(labels, code))
    dt = []
    for i in range(1, n, 2):
        a = label_to_code[i]
        j = code_to_label[-a]
        if a < 0:
            j = -j
        dt.append(j)

    return spherogram.Link(f'DT: {dt}')


def regina_DT(code):
    """
    >>> L = regina_DT('flmnopHKRQGabcdeJI')
    >>> L
    <Link: 1 comp; 18 cross>
    >>> L = regina_DT('hPONrMjsLfiQGDCBKae')
    >>> L
    <Link: 1 comp; 19 cross>
    """
    cross = len(code)
    l = string.ascii_letters[cross - 1]
    dt_prefix = 'DT: ' + l + 'a' + l
    return spherogram.Link(dt_prefix + code)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
