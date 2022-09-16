class OrderedSet():
    def __init__(self, iterable=None):
        if iterable is not None:
            self.elts = {e: None for e in iterable}
        else:
            self.elts = dict()

    def __len__(self):
        return len(self.elts)

    def __contains__(self, key):
        return key in self.elts

    def add(self, key):
        self.elts[key] = None

    def discard(self, key):
        if key in self.elts:
            self.elts.pop(key)

    def __iter__(self):
        return iter(self.elts)

    def pop(self):
        return self.elts.popitem()[0]

    def update(self, sequence):
        for s in sequence:
            self.elts[s] = None

    def difference_update(self, sequence):
        for s in sequence:
            if s in self.elts:
                self.elts.pop(s)

    def remove(self, key):
        if key not in self.elts:
            raise KeyError('element not in set')
        self.elts.pop(key)

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and all(x == y
                                                   for x, y in zip(self, other))
        return set(self) == set(other)


if __name__ == '__main__':
    s = OrderedSet('abracadaba')
    t = OrderedSet('simsalabim')
    print(s == t, t == t)
