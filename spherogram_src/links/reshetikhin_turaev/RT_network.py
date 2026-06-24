class DirectedEdge:
    __slots__ = ['label', 'index', 'sign', 'reversed_edge']

    def __init__(self, label, reversed_edge = None):
        self.label = label
        self.index = max(label, ~label)
        self.sign = 1 if label == self.index else -1

        if reversed_edge is None: 
            reversed_edge = DirectedEdge(~self.label, self)    
        self.reversed_edge = reversed_edge

    def __str__(self):
        return ('' if self.sign == 1 else '~') + str(self.index)

    def __repr__(self):
        return str(self)
    
    def __hash__(self):
        return hash(self.label)
    
    def __eq__(self, other):
        return self.label == other.label
    
    def __invert__(self):
        return self.reversed_edge

class RTNetwork:
    
    def __init__(self, tensors, T = None, network = None, rot_num = None, boundary = None, boundary_labels = None):
        """
        Represent the tensor network obtained by applying 
        the Reshetikhin--Turaev functor determined by the 
        given RMatrix tensors to the tangle T.
    
        The network is represented as a list of pairs (tensor, legs)
        """
        self.tensors = tensors

        if T is not None:
            assert T.is_upward(), 'Tangle should be upward for the Reshetikhin--Turaev functor to apply'
            self.rot_num = T.rot_num()
            self.tangle = T.copy()

            self.boundary = T.boundary
            self.boundary_labels = T.strand_labels

            self.edge = edge = dict()

            self.network = network = []
            for c in T.crossings:
                labels = c.strand_labels
                for lab in labels:
                    if lab not in edge.keys():
                        edge[lab] = DirectedEdge(lab)
                        
                if c.sign == 1:
                    key = (~edge[labels[3]], 
                           ~edge[labels[0]], 
                            edge[labels[2]], 
                            edge[labels[1]])
                else:
                    assert c.sign == -1, f'Crossing {c} is not oriented'
                    key = (~edge[labels[0]], 
                           ~edge[labels[1]], 
                            edge[labels[3]], 
                            edge[labels[2]])

                network.append((tensors.R(c.sign), key))

            self.idle_labels = set(self.boundary_labels)

            for arc in self.idle_labels:
                if arc not in edge.keys():
                    edge[arc] = DirectedEdge(arc)
                    network.append((tensors.h(0), (~edge[arc], edge[arc])))
        else:
            assert all(item is not None for item in (network, rot_num, boundary, boundary_labels))
            self.network = network
            self.rot_num = rot_num
            self.boundary = boundary
            self.boundary_labels = boundary_labels
            self.idle_labels = set(boundary_labels)

            self.edge = edge = dict()
            for _, key in network:
                for e in key:
                    if e.index not in edge.keys():
                        edge[e.index] = e if e.sign == 1 else ~e

    def optimal_contraction_sequence(self):
        try:
            import numpy as np
            import opt_einsum as oe
        except ImportError:
            raise ModuleNotFoundError('Modules numpy and opt_einsum is required for computing the optimal contraction sequences')
        
        oe_network = []
        idle = set(self.idle_labels)
        for tensor, key in self.network:
            if len(key) == 2 and key[0].index == key[1].index:
                try:
                    idle.remove(key[0].index)
                except:
                    raise ValueError(f'key {key[0].index} not found in {idle}')
            else:
                oe_network.append(np.empty(tensor.shape))
                oe_network.append([edge.index for edge in key])

        return oe.contract_path(*oe_network, idle)[0]

    def contract_nodes(self, indices):
        idx1, idx2 = indices
        tensor1, key1 = self.network[idx1]
        tensor2, key2 = self.network[idx2]

        pairs = {}
        for pos_i, ei in enumerate(key1):
            for pos_j, ej in enumerate(key2):
                if ei.index == ej.index and ei.sign * ej.sign == -1:
                    side = 0 if ei.sign == 1 else 1
                    pairs[(pos_i, pos_j)] = (side, self.tensors.h(self.rot_num[ei.index]))

        result_tensor = tensor1.decorated_contract(tensor2, pairs)

        contracted1 = {pos_i for pos_i, _ in pairs}
        contracted2 = {pos_j for _, pos_j in pairs}

        if idx1 == idx2:
            contracted_all = contracted1 | contracted2
            new_key = tuple(e for pos, e in enumerate(key1) if pos not in contracted_all)
            self.network.pop(idx1)
        else:
            new_key = (tuple(e for pos, e in enumerate(key1) if pos not in contracted1) +
                       tuple(e for pos, e in enumerate(key2) if pos not in contracted2))
            hi, lo = max(idx1, idx2), min(idx1, idx2)
            self.network.pop(hi)
            self.network.pop(lo)

        del tensor1, tensor2
        
        self.network.append((result_tensor, new_key))
        self._resolve_self_loops(len(self.network) - 1)

    def _resolve_self_loops(self, idx):
        while True:
            _, key = self.network[idx]
            pairs = {}
            seen = {}
            for pos, e in enumerate(key):
                if e.index in seen:
                    other_pos, other_e = seen[e.index]
                    if e.sign * other_e.sign == -1:
                        side = 0 if other_e.sign == 1 else 1
                        pairs[(other_pos, pos)] = (side, self.tensors.h(self.rot_num[e.index]))
                else:
                    seen[e.index] = (pos, e)
            if not pairs:
                break
            self.contract_nodes((idx, idx))
            idx = len(self.network) - 1
 
    def contract_sequence(self, seq):
        for indices in seq:
            self.contract_nodes(indices)

    def contract_all(self):
        for i in range(len(self.network)):
            self._resolve_self_loops(i)
        self.contract_sequence(self.optimal_contraction_sequence())

    def evaluate(self):
        """
        Fixate all idle labels at value 0, obtaining a new RTNework with (0,0) boundary,
        contract all and return the product of all values of the resulting tensors.
        """
        assert self.boundary == (1,1)

        new_network = []
        prefactor = 1

        for tensor, key in self.network:
            idle_positions = sorted(
                [pos for pos, e in enumerate(key) if e.index in self.idle_labels],
                reverse=True
            )
            non_idle_key = tuple(e for e in key if e.index not in self.idle_labels)
            t = tensor
            for pos in idle_positions:
                t = t.fixate(pos, 0)
            if t.rank == 0:
                prefactor *= t[()]
            else:
                new_network.append((t, non_idle_key))

        reduced = RTNetwork(
            self.tensors,
            network=new_network,
            rot_num=self.rot_num,
            boundary=(0, 0),
            boundary_labels=[]
        )
        reduced.contract_all()

        result = prefactor
        for tensor, _ in reduced.network:
            result *= tensor[()]
        return result
