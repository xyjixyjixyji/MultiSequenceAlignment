"""
Algorithms implemented for Multiple Sequence Alignment (MSA)
* Dynamic Programming (in other file)
* A* Algorithm (in other file)
* Genetic Algorithm
    * Genetic Algorithm impl in this file
"""
import numpy as np
import json
import time
from MSA_DP import DP2d, DP3d

json_file_path = "hparam.json"
with open(json_file_path, 'r') as f:
    params = json.load(f)

scores = params['scores']


class Genetic2D():
    
    ''' 
    @param:
        query, seq: 2 sequences to be aligned
        nr_eachGen: # of pops in each gen
        nr_gen:: # of gens the algorithm shall inherit
        poss_gap: the possibility a match (or mismatch) will happen in the sequence
        scores: dictionary to count the cost
    '''   
    def __init__(self, query, seq, kwargs, scores=scores):
        self.query = query
        self.seq = seq
        # pops is just a bunch of seqs of [x, xy, y, ....]
        # the indices are corresponding
        self.pops = []
        self.costs = []
        
        self.nr_eachGen = kwargs["nr_eachGen"]
        self.nr_gen = kwargs["nr_gen"]
        
        # (0, poss_xy) -> xy (poss_xy + poss_x) -> x (poss_xy + poss_x + poss_y) -> y
        self.poss_xy = kwargs["poss_xy"]
        self.poss_x = (1 - self.poss_xy) / 2
        self.poss_y = self.poss_x
        
        self.mutation_rate = kwargs["poss_mut"]
        self.crossover_rate = kwargs["poss_cross"]
        
        self.scores = scores
        
    def gen_pop(self):
        #! call when the class is defined the first time
        # to generate the initial pops
        query, seq = self.query, self.seq
        m, n = len(query), len(seq)
        count = 0
        while count < self.nr_eachGen:
            pop = []
            ptr_x, ptr_y = 0, 0
            while (ptr_x != m or ptr_y != n):
                random_num = np.random.uniform(low=0.0, high=1.0)
                if random_num < self.poss_xy:
                    if ptr_x < m and ptr_y < n:
                        # move is XY
                        pop.append('xy')
                        ptr_x += 1
                        ptr_y += 1
                    else:   continue
                elif random_num < (self.poss_xy + self.poss_x):
                    if ptr_x < m:
                        pop.append('x')
                        ptr_x += 1
                    else:   continue
                else:
                    if ptr_y < n:
                        pop.append('y')
                        ptr_y += 1
                    else:   continue
            self.pops.append(pop)
            count += 1
        
    def get_cost(self):
        # called each time the generation is updated
        self.costs.clear()
        scores = self.scores
        s_match, s_sub, s_gap = scores.get('match', 0), scores.get('sub', 3), scores.get('gap', 2)
        for pop in self.pops:
            cost = 0
            ptr_x, ptr_y = 0, 0
            for move in pop:
                if move == 'x':
                    ptr_x += 1
                    cost += s_gap
                elif move == 'y':
                    ptr_y += 1
                    cost += s_gap
                else: # move == 'xy'
                    ptr_x += 1
                    ptr_y += 1
                    __cost_xy = s_match if (self.query[ptr_x - 1] == self.seq[ptr_y - 1]) else s_sub
                    cost += __cost_xy
            self.costs.append(cost)
                    
    def get_Fitness(self):
        # find the best cost individual
        # The fitness of each individual is the best minus their cost
        # lower cost -> higher fitness
        max = float('-inf')
        for cost in self.costs:
            if cost > max:
                max = cost
        fit = np.array(self.costs)
        fit = (max - fit) + 1e-3
        return fit
    
    def select(self):
        fit = self.get_Fitness()
        indices = np.random.choice(np.arange(self.nr_eachGen), size=self.nr_eachGen, replace=True,
                           p=(fit / fit.sum()))
        new_pops = []
        for index in indices:
            new_pops.append(self.pops[index])
        self.pops.clear()
        self.pops.extend(new_pops)
        self.get_cost()

    def mutation(self):
        mutation_rate = self.mutation_rate
        # best fit individual do not incoperate the mutation and crossover
        minindex, mincost = -1, float('inf')
        for index, cost in enumerate(self.costs):
            if cost < mincost:
                mincost = cost
                minindex = index

        for index, indiv in enumerate(self.pops):
            if index == minindex:
                continue

            if np.random.rand() < mutation_rate:
                len_indiv = len(indiv)
                mut_pos = np.random.randint(0, len_indiv)
                if indiv[mut_pos] == 'xy':
                    indiv[mut_pos] = 'x'
                    indiv.insert(mut_pos + 1, 'y')
                elif indiv[mut_pos] == 'x':
                    ptrl = mut_pos - 1
                    ptrr = mut_pos + 1
                    while(ptrl >= 0 and ptrr < len_indiv):
                        if ptrl >= 0 and indiv[ptrl] == 'y':
                            indiv[mut_pos] = 'xy'
                            indiv.pop(ptrl)
                            break
                        elif ptrr < len_indiv and indiv[ptrr] == 'y':
                            indiv[mut_pos] = 'xy'
                            indiv.pop(ptrr)
                            break
                        else:
                            ptrl -= 1
                            ptrr += 1
                elif indiv[mut_pos] == 'y':
                    ptrl = mut_pos - 1
                    ptrr = mut_pos + 1
                    while(ptrl >= 0 or ptrr < len_indiv):
                        if ptrl >= 0 and indiv[ptrl] == 'x':
                            indiv[mut_pos] = 'xy'
                            indiv.pop(ptrl)
                            break
                        elif ptrr < len_indiv and indiv[ptrr] == 'x':
                            indiv[mut_pos] = 'xy'
                            indiv.pop(ptrr)
                            break
                        else:
                            ptrl -= 1
                            ptrr += 1
                else:   pass
        
    def crossover(self):
        new_pops = []
        crossover_rate = self.crossover_rate
        mincost = float('inf')
        minindex = -1
        # best fit individual do not incoperate the mutation and crossover
        for index, cost in enumerate(self.costs):
            if cost < mincost:
                mincost = cost
                minindex = index
        
        for index, father in enumerate(self.pops):
            if index == minindex:
                new_pops.append(father)
                continue
            child = father
            mother = self.pops[np.random.randint(low=0, high=self.nr_eachGen)]
            # find the cross points
            cross_point = np.random.randint(0, len(father))
            # count the m and n
            m, n = 0, 0
            for move in father[:cross_point]:
                if move == 'x' or move == 'xy':
                    m += 1
                if move == 'y' or move == 'xy':
                    n += 1
            count_m, count_n = 0, 0
            for index, move in enumerate(mother):
                if move == 'x' or move == 'xy':
                    count_m += 1
                if move == 'y' or move == 'xy':
                    count_n += 1
                if count_m == m and count_n == n:
                    if np.random.rand() < crossover_rate:
                        child[cross_point:] = mother[index + 1:]
            new_pops.append(child)
        self.pops.clear()
        self.pops.extend(new_pops)

    # below are not only for genetic algorithm, not in general
    def getcost(self, moves):
        scores = self.scores
        s_match, s_sub, s_gap = scores.get('match', 0), scores.get('sub', 3), scores.get('gap', 2)
        cost = 0
        ptr_x, ptr_y = 0, 0
        for move in moves:
            if move == 'x':
                ptr_x += 1
                cost += s_gap
            elif move == 'y':
                ptr_y += 1
                cost += s_gap
            else: # move == 'xy'
                ptr_x += 1
                ptr_y += 1
                __cost_xy = s_match if (self.query[ptr_x - 1] == self.seq[ptr_y - 1]) else s_sub
                cost += __cost_xy
        return cost

    def findmin(self):
        mincost, minindex = float('inf'), -1
        for i, c in enumerate(self.costs):
            if c < mincost:
                mincost = c
                minindex = i
        return self.pops[minindex]

    def run(self):
        self.gen_pop()
        self.get_cost()
        iter = 0
        while (iter < self.nr_gen):
            self.crossover()
            self.mutation()
            self.get_cost()
            self.select()
            iter += 1

    def resolve2d(self, moves):
        query, data = self.query, self.seq
        seq1, seq2 = "", ""
        ptr1, ptr2 = 0, 0
        for move in moves:
            if move == 'x':
                seq1 += query[ptr1]
                seq2 += '-'
                ptr1 += 1
            elif move == 'y':
                seq2 += data[ptr2]
                seq1 += '-'
                ptr2 += 1
            else:
                seq1 += query[ptr1]
                seq2 += data[ptr2]
                ptr1 += 1
                ptr2 += 1
        return seq1, seq2

def genetic2d(seq1, seq2, kwargs, scores):
    gen = Genetic2D(seq1, seq2, kwargs, scores=scores)
    gen.run()
    minimal = gen.findmin()
    score = gen.getcost(minimal)
    query_a, seq_a = gen.resolve2d(minimal)
    return score, query_a, seq_a

class Genetic3D():
    
    def __init__(self, query, seq1, seq2, kwargs, scores=scores):
        self.query = query
        self.seq1 = seq1
        self.seq2 = seq2
        self.m, self.n, self.p = len(query), len(seq1), len(seq2)
        # pops is just a bunch of seqs of [x, xy, y, ....]
        # the indices are corresponding
        self.pops = []
        self.costs = []
        
        self.nr_eachGen = kwargs["nr_eachGen"]
        self.nr_gen = kwargs["nr_gen"]
        
        # possibilities
        self.poss_xyz = kwargs["poss_xyz"]
        self.poss_xy = kwargs["poss_xy"]
        self.poss_xz = kwargs["poss_xy"]
        self.poss_yz = kwargs["poss_xy"]
        self.poss_x = (1 - self.poss_xyz - 3 * self.poss_xy) / 3
        self.poss_y = self.poss_x
        self.poss_z = self.poss_x
        
        self.mutation_rate = kwargs["poss_mut"]
        self.crossover_rate = kwargs["poss_cross"]
        
        self.scores = scores
    
    def gen_pop(self):
        #! call when the class is defined the first time
        # to generate the initial pops
        m, n, p = self.m, self.n, self.p
        # print(self.poss_xyz, self.poss_xy, self.poss_xz, self.poss_yz, self.poss_y, self.poss_z, self.poss_x)
        count = 0
        while count < self.nr_eachGen:
            pop = []
            ptr_x, ptr_y, ptr_z = 0, 0, 0
            while (ptr_x != m or ptr_y != n or ptr_z != p):
                # print((ptr_x, m), (ptr_y, n), (ptr_z, p))
                random_num = np.random.uniform(low=0.0, high=1.0)
                if random_num < self.poss_xyz:
                    if ptr_x < m and ptr_y < n and ptr_z < p:
                        pop.append('xyz')
                        ptr_x, ptr_y, ptr_z = \
                        ptr_x + 1, ptr_y + 1, ptr_z + 1
                    else:   continue
                elif random_num < self.poss_xyz + self.poss_xy:
                    if ptr_x < m and ptr_y < n:
                        pop.append('xy')
                        ptr_x, ptr_y = ptr_x + 1, ptr_y + 1
                    else:    continue
                elif random_num < self.poss_xyz + self.poss_xy + self.poss_xz:
                    if ptr_x < m and ptr_z < p:
                        pop.append('xz')
                        ptr_x, ptr_z = ptr_x + 1, ptr_z + 1
                    else:    continue
                elif random_num < self.poss_xyz + self.poss_xy + self.poss_xz + self.poss_yz:
                    if ptr_y < n and ptr_z < p:
                        pop.append('yz')
                        ptr_y, ptr_z = ptr_y + 1, ptr_z + 1
                    else:    continue
                elif random_num < self.poss_xyz + self.poss_xy + self.poss_xz + self.poss_yz + \
                                  self.poss_x:
                    if ptr_x < m:
                        pop.append('x')
                        ptr_x += 1
                    else:   continue
                elif random_num < self.poss_xyz + self.poss_xy + self.poss_xz + self.poss_yz + \
                                  self.poss_x + self.poss_y:
                    if ptr_y < n:
                        pop.append('y')
                        ptr_y += 1
                    else:   continue
                else:
                    if ptr_z < p:
                        pop.append('z')
                        ptr_z += 1
                    else:   continue
                        
            self.pops.append(pop)
            # print(pop)
            count += 1
        
    def get_cost(self):
        self.costs.clear()
        scores = self.scores
        s_match, s_sub, s_gap = scores.get('match', 0), scores.get('sub', 3), scores.get('gap', 2)
        seq1, seq2, seq3 = self.query, self.seq1, self.seq2
        
        # for pop in self.pops:
        #     x = pop.count('x')
        #     y = pop.count('y')
        #     z = pop.count('z')
        #     xy = pop.count('xy')
        #     xz = pop.count('xz')
        #     yz = pop.count('yz')
        #     xyz = pop.count('xyz')
        #     print(x + xy + xz + xyz, y + yz + xy + xyz, z + xz + yz + xyz)
        #     print(self.m, self.n, self.p)
        
        for pop in self.pops:
            cost = 0
            ptr_x, ptr_y, ptr_z = 0, 0, 0
            for move in pop:
                # print(ptr_x, ptr_y, ptr_z)
                # print(move)
                if move == 'x':
                    ptr_x += 1
                    cost += 2 * s_gap
                elif move == 'y':
                    ptr_y += 1
                    cost += 2 * s_gap
                elif move == 'z':
                    ptr_z += 1
                    cost += 2 * s_gap
                elif move == 'xy':
                    ptr_x, ptr_y = ptr_x + 1, ptr_y + 1
                    __score_xy = s_match if (seq1[ptr_x - 1] == seq2[ptr_y - 1]) else s_sub
                    cost += 2 * s_gap + __score_xy
                elif move == 'xz':
                    ptr_x, ptr_z = ptr_x + 1, ptr_z + 1
                    __score_xz = s_match if (seq1[ptr_x - 1] == seq3[ptr_z - 1]) else s_sub
                    cost += 2 * s_gap + __score_xz
                elif move == 'yz':
                    ptr_y, ptr_z = ptr_y + 1, ptr_z + 1
                    __score_yz = s_match if (seq2[ptr_y - 1] == seq3[ptr_z - 1]) else s_sub
                    cost += 2 * s_gap + __score_yz
                elif move == 'xyz':
                    ptr_x, ptr_y, ptr_z = \
                    ptr_x + 1, ptr_y + 1, ptr_z + 1
                    __score_xy = s_match if (seq1[ptr_x - 1] == seq2[ptr_y - 1]) else s_sub
                    __score_xz = s_match if (seq1[ptr_x - 1] == seq3[ptr_z - 1]) else s_sub
                    __score_yz = s_match if (seq2[ptr_y - 1] == seq3[ptr_z - 1]) else s_sub
                    cost += __score_xy + __score_xz + __score_yz
            self.costs.append(cost)
            
    def get_Fitness(self):
        max = float('-inf')
        for cost in self.costs:
            if cost > max:
                max = cost
        fit = np.array(self.costs)
        fit = (max - fit) + 1e-3
        return fit
        
    def select(self):
        fit = self.get_Fitness()
        indices = np.random.choice(np.arange(self.nr_eachGen), size=self.nr_eachGen, replace=True,
                           p=(fit / fit.sum()))
        new_pops = []
        for index in indices:
            new_pops.append(self.pops[index])
        self.pops.clear()
        self.pops.extend(new_pops)
        self.get_cost()
    
    def mutation(self):
        '''
        7 cases:
            * x:   f(x) = xy, delete the closest y
            * y:   f(y) = yz, delete the closest z
            * z:   f(z) = xz, delete the closest x
            * xy:  f(xy) = xyz, delete the closest z
            * xz:  f(xz) = xyz, delete the closest y
            * yz:  f(yz) = xyz, delete the closest x
            * xyz: f(xyz) = xy, insert a z next to it
        '''
        mutation_rate = self.mutation_rate
        # best fit individual do not incorperate the mutation
        mincost = float('inf')
        minindex = -1
        for index, cost in enumerate(self.costs):
            if cost < mincost:
                mincost = cost
                minindex = index

        for index, indiv in enumerate(self.pops):
            # operation each individual except the best
            if index == minindex:
                continue

            if np.random.rand() < mutation_rate:
                len_indiv = len(indiv)
                mut_pos = np.random.randint(0, len_indiv)
                # start doing mutation op
                if indiv[mut_pos] == 'xyz':
                    indiv[mut_pos] = 'xy'
                    indiv.insert(mut_pos + 1, 'z')
                elif indiv[mut_pos] == 'xy':
                    ptrl = mut_pos - 1
                    ptrr = mut_pos + 1
                    while (ptrl >= 0 and ptrr < len_indiv):
                        if ptrl >= 0 and indiv[ptrl] == 'z':
                            indiv[mut_pos] = 'xyz'
                            indiv.pop(ptrl)
                            break
                        elif ptrr < len_indiv and indiv[ptrr] == 'z':
                            indiv[mut_pos] = 'xyz'
                            indiv.pop(ptrr)
                            break
                        else:
                            ptrl -= 1
                            ptrr += 1
                elif indiv[mut_pos] == 'xz':
                    ptrl = mut_pos - 1
                    ptrr = mut_pos + 1
                    while (ptrl >= 0 and ptrr < len_indiv):
                        if ptrl >= 0 and indiv[ptrl] == 'y':
                            indiv[mut_pos] = 'xyz'
                            indiv.pop(ptrl)
                            break
                        elif ptrr < len_indiv and indiv[ptrr] == 'y':
                            indiv[mut_pos] = 'xyz'
                            indiv.pop(ptrr)
                            break
                        else:
                            ptrl -= 1
                            ptrr += 1
                elif indiv[mut_pos] == 'yz':
                    ptrl = mut_pos - 1
                    ptrr = mut_pos + 1
                    while (ptrl >= 0 and ptrr < len_indiv):
                        if ptrl >= 0 and indiv[ptrl] == 'x':
                            indiv[mut_pos] = 'xyz'
                            indiv.pop(ptrl)
                            break
                        elif ptrr < len_indiv and indiv[ptrr] == 'x':
                            indiv[mut_pos] = 'xyz'
                            indiv.pop(ptrr)
                            break
                        else:
                            ptrl -= 1
                            ptrr += 1
                elif indiv[mut_pos] == 'x':
                    ptrl = mut_pos - 1
                    ptrr = mut_pos + 1
                    while (ptrl >= 0 and ptrr < len_indiv):
                        if ptrl >= 0 and indiv[ptrl] == 'y':
                            indiv[mut_pos] = 'xy'
                            indiv.pop(ptrl)
                            break
                        elif ptrr < len_indiv and indiv[ptrr] == 'y':
                            indiv[mut_pos] = 'xy'
                            indiv.pop(ptrr)
                            break
                        else:
                            ptrl -= 1
                            ptrr += 1
                elif indiv[mut_pos] == 'y':
                    ptrl = mut_pos - 1
                    ptrr = mut_pos + 1
                    while (ptrl >= 0 and ptrr < len_indiv):
                        if ptrl >= 0 and indiv[ptrl] == 'z':
                            indiv[mut_pos] = 'yz'
                            indiv.pop(ptrl)
                            break
                        elif ptrr < len_indiv and indiv[ptrr] == 'z':
                            indiv[mut_pos] = 'yz'
                            indiv.pop(ptrr)
                            break
                        else:
                            ptrl -= 1
                            ptrr += 1
                elif indiv[mut_pos] == 'z':
                    ptrl = mut_pos - 1
                    ptrr = mut_pos + 1
                    while (ptrl >= 0 and ptrr < len_indiv):
                        if ptrl >= 0 and indiv[ptrl] == 'x':
                            indiv[mut_pos] = 'xz'
                            indiv.pop(ptrl)
                            break
                        elif ptrr < len_indiv and indiv[ptrr] == 'x':
                            indiv[mut_pos] = 'xz'
                            indiv.pop(ptrr)
                            break
                        else:
                            ptrl -= 1
                            ptrr += 1
                else:   pass
                
    def crossover(self):
        new_pops = []
        crossover_rate = self.crossover_rate
        # best fit individual do not incorperate the crossover
        mincost = float('inf')
        minindex = -1
        for index, cost in enumerate(self.costs):
            if cost < mincost:
                mincost = cost
                minindex = index
        
        for index, father in enumerate(self.pops):
            if index == minindex:
                new_pops.append(father)
                continue
            
            child = father
            mother = self.pops[np.random.randint(low=0, high=self.nr_eachGen)]
            # find the crosspoints
            cross_point = np.random.randint(0, len(father))
            # count the m, n, p
            m, n, p = 0, 0, 0
            for move in father[:cross_point]:
                if move == 'x' or move == 'xy' or move == 'xz' or move == 'xyz':
                    m += 1
                if move == 'y' or move == 'xy' or move == 'yz' or move == 'xyz':
                    n += 1
                if move == 'z' or move == 'xz' or move == 'yz' or move == 'xyz':
                    p += 1
            # find valid crosspoint in mother individual
            count_m, count_n, count_p = 0, 0, 0
            for index, move in enumerate(mother):
                if move == 'x' or move == 'xy' or move == 'xz' or move == 'xyz':
                    count_m += 1
                if move == 'y' or move == 'xy' or move == 'yz' or move == 'xyz':
                    count_n += 1
                if move == 'z' or move == 'xz' or move == 'yz' or move == 'xyz':
                    count_p += 1
                if count_m == m and count_n == n and count_p == p:
                    if np.random.rand() < crossover_rate:
                        child[cross_point:] = mother[index + 1:]
            new_pops.append(child)
        self.pops.clear()
        self.pops.extend(new_pops)
    
    def run(self):
        self.gen_pop()
        self.get_cost()
        iter = 0
        while (iter < self.nr_gen):
            self.crossover()
            self.mutation()
            self.get_cost()
            self.select()
            iter += 1
    
    def findmin(self):
        mincost, minindex = float('inf'), -1
        for i, c in enumerate(self.costs):
            if c < mincost:
                mincost = c
                minindex = i
        return self.pops[minindex]
    
    def resolve3d(self, moves):
        query, seq1, seq2 = self.query, self.seq1, self.seq2
        query_o, seq1_o, seq2_o = "", "", ""
        ptr_x, ptr_y, ptr_z = 0, 0, 0
        for move in moves:
            if move == 'x':
                query_o += query[ptr_x]
                ptr_x += 1
            elif move == 'y':
                seq1_o += seq1[ptr_y]
                ptr_y += 1
            elif move == 'z':
                seq2_o += seq2[ptr_z]
                ptr_z += 1
            elif move == 'xy':
                query_o += query[ptr_x]
                seq1_o += seq1[ptr_y]
                ptr_x, ptr_y = ptr_x + 1, ptr_y + 1
            elif move == 'xz':
                query_o += query[ptr_x]
                seq2_o += seq2[ptr_z]
                ptr_x, ptr_z = ptr_x + 1, ptr_z + 1
            elif move == 'yz':
                seq1_o += seq1[ptr_y]
                seq2_o += seq2[ptr_z]
                ptr_y, ptr_z = ptr_y + 1, ptr_z + 1
            else:
                query_o += query[ptr_x]
                seq1_o += seq1[ptr_y]
                seq2_o += seq2[ptr_z]
                ptr_x, ptr_y, ptr_z = ptr_x + 1, ptr_y + 1, ptr_z + 1
        return query_o, seq1_o, seq2_o
    
    def getcost(self, moves):      
        scores = self.scores
        s_match, s_sub, s_gap = \
        scores.get('match', 0), scores.get('sub', 3), scores.get('gap', 2)
        
        seq1, seq2, seq3 = self.query, self.seq1, self.seq2
        
        ptr_x, ptr_y, ptr_z = 0, 0, 0
        cost = 0
        for move in moves:
            # print(ptr_x, ptr_y, ptr_z)
            # print(move)
            if move == 'x':
                ptr_x += 1
                cost += 2 * s_gap
            elif move == 'y':
                ptr_y += 1
                cost += 2 * s_gap
            elif move == 'z':
                ptr_z += 1
                cost += 2 * s_gap
            elif move == 'xy':
                ptr_x, ptr_y = ptr_x + 1, ptr_y + 1
                __score_xy = s_match if (seq1[ptr_x - 1] == seq2[ptr_y - 1]) else s_sub
                cost += 2 * s_gap + __score_xy
            elif move == 'xz':
                ptr_x, ptr_z = ptr_x + 1, ptr_z + 1
                __score_xz = s_match if (seq1[ptr_x - 1] == seq3[ptr_z - 1]) else s_sub
                cost += 2 * s_gap + __score_xz
            elif move == 'yz':
                ptr_y, ptr_z = ptr_y + 1, ptr_z + 1
                __score_yz = s_match if (seq2[ptr_y - 1] == seq3[ptr_z - 1]) else s_sub
                cost += 2 * s_gap + __score_yz
            elif move == 'xyz':
                ptr_x, ptr_y, ptr_z = \
                ptr_x + 1, ptr_y + 1, ptr_z + 1
                __score_xy = s_match if (seq1[ptr_x - 1] == seq2[ptr_y - 1]) else s_sub
                __score_xz = s_match if (seq1[ptr_x - 1] == seq3[ptr_z - 1]) else s_sub
                __score_yz = s_match if (seq2[ptr_y - 1] == seq3[ptr_z - 1]) else s_sub
                cost += __score_xy + __score_xz + __score_yz
        return cost

def genetic3d(seq1, seq2, seq3, kwargs, scores):
    gen3d = Genetic3D(seq1, seq2, seq3, kwargs, scores = scores)
    gen3d.run()
    minimal = gen3d.findmin()
    score = gen3d.getcost(minimal)
    query_a, seq1_a, seq2_a = gen3d.resolve3d(minimal)
    return score, query_a, seq1_a, seq2_a

def test2d():
    seq1 = 'KJXXJAJKPXKJJXJKPXKJXXJAJKPXKJJXJKPXKJXXJAJKPXKJXXJAJKHXKJXXJAJKPXKJXXJAJKHXKJXX'
    seq2 = 'VXTLKZOKMOKAPHXHMLOWZHTPPHKPKIAXPOXKSKSWJSTSGNSHIOTTLPLLMZKUJHXTPWOWHZGAHLWKKPKMPXOTMZJUOPJ'
    kwargs = params["genetic2d"]

    gen2d= Genetic2D(seq1, seq2, kwargs)
    
    gen2d.run()
    
    print(gen2d.getcost(gen2d.findmin()))
    
    _, score = DP2d(seq1, seq2, scores)
    
    print("optimal: %d" % score)

def test3d():
    seq1 = 'KJXXJAJKPXKJJXJKPXKJXXJAJKPXKJJXJKPXKJXXJAJKPXKJXXJAJKHXKJXXJAJKPXKJXXJAJKHXKJXX'
    seq2 = 'VXTLKZOKMOKAPHXHMLOWZHTPPHKPKIAXPOXKSKSWJSTSGNSHIOTTLPLLMZKUJHXTPWOWHZGAHLWKKPKMPXOTMZJUOPJ'
    seq3 = 'PJJAPJJPPJJPJJAPJJPPJJPJJAPJJPPJJPJJAPJJPPJJPJJAPJJPPJJPJJAPJJPJJKJJP'
    kwargs = params["genetic3d"]

    start = time.time()
    gen3d = Genetic3D(seq1, seq2, seq3, kwargs)
    gen3d.run()
    
    minimal = gen3d.findmin()
    end = time.time()
    # print(minimal)
    print(gen3d.getcost(minimal))
    
    _, score = DP3d(seq1, seq2, seq3, scores)
    print("optimal: %d" % score)
    print("TIME EXPIRED: %.6f" % (end - start))

if __name__ == '__main__':
    test2d() 
    test3d()