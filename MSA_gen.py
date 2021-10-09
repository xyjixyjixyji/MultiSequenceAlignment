"""
Algorithms implemented for Multiple Sequence Alignment (MSA)
* Dynamic Programming (in other file)
* A* Algorithm (in other file)
* Genetic Algorithm
    * Genetic Algorithm impl in this file
"""
import time
import numpy as np
import copy
from MSA_DP import DP2d
from numpy import random


scores = {
    'match': 0,
    'sub': 3,
    'gap': 2,
}

class Genetic2D():
    
    ''' 
    @param:
        query, seq: 2 sequences to be aligned
        nr_eachGen: # of pops in each gen
        nr_gen:: # of gens the algorithm shall inherit
        poss_gap: the possibility a match (or mismatch) will happen in the sequence
        scores: dictionary to count the cost
    '''   
    def __init__(self, query, seq, nr_eachGen=50, nr_gen=300, poss_xy=0.8, poss_mut = 0.001, poss_cross=0.8, scores=scores):
        self.query = query
        self.seq = seq
        # pops is just a bunch of seqs of [x, xy, y, ....]
        # the indices are corresponding
        self.pops = []
        self.costs = []
        
        self.nr_eachGen = nr_eachGen
        self.nr_gen = nr_gen
        
        # (0, poss_xy) -> xy (poss_xy + poss_x) -> x (poss_xy + poss_x + poss_y) -> y
        self.poss_xy = poss_xy
        self.poss_x = (1 - poss_xy) / 2
        self.poss_y = self.poss_x
        
        self.mutation_rate = poss_mut
        self.crossover_rate = poss_cross
        
        self.scores = scores
        
    def gen_pop(self):
        #! call when the class is defined the first time
        # to generate the initial pops
        query, seq = self.query, self.seq
        m, n = len(query), len(seq)
        # the # of xy should not be more than min_len
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
        for indiv in self.pops:
            if np.random.rand() < self.mutation_rate:
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

def genetic2d(seq1, seq2, scores):
    gen = Genetic2D(seq1, seq2, scores=scores)
    gen.run()
    minimal = gen.findmin()
    score = gen.getcost(minimal)
    query_a, seq_a = gen.resolve2d(minimal)
    return score, query_a, seq_a

class Genetic3D():
    
    def __init__(self, query, seq1, seq2, nr_eachGen=50, nr_gen=300, \
        poss_xyz=0.7, poss_xy = 0.06,\
        poss_mut = 0.001, poss_cross=0.8, \
        scores=scores):
        self.query = query
        self.seq1 = seq1
        self.seq2 = seq2
        self.x, self.y, self.z = len(query), len(seq1), len(seq2)
        # pops is just a bunch of seqs of [x, xy, y, ....]
        # the indices are corresponding
        self.pops = []
        self.costs = []
        
        self.nr_eachGen = nr_eachGen
        self.nr_gen = nr_gen
        
        # possibilities
        self.poss_xyz = poss_xyz
        self.poss_xy = poss_xy
        self.poss_xz = poss_xy
        self.poss_yz = poss_xy
        self.poss_x = 1 - poss_xyz - 3 * poss_xy
        self.poss_y = self.poss_x
        self.poss_z = self.poss_x
        
        self.mutation_rate = poss_mut
        self.crossover_rate = poss_cross
        
        self.scores = scores

def test():
    seq1 = 'KJXXJAJKPXKJJXJKPXKJXXJAJKPXKJJXJKPXKJXXJAJKPXKJXXJAJKHXKJXXJAJKPXKJXXJAJKHXKJXX'
    seq2 = 'VXTLKZOKMOKAPHXHMLOWZHTPPHKPKIAXPOXKSKSWJSTSGNSHIOTTLPLLMZKUJHXTPWOWHZGAHLWKKPKMPXOTMZJUOPJ'
    
    gen2d= Genetic2D(seq1, seq2)
    
    gen2d.run()
    
    print(gen2d.getcost(gen2d.findmin()))
    
    _, score = DP2d(seq1, seq2, scores)
    
    print("optimal: %d" % score)

if __name__ == '__main__':
    test()