"""
Algorithms implemented for Multiple Sequence Alignment (MSA)
* Dynamic Programming （in other file)
* A* Algorithm
    * AStar2d() and AStar3d() implemented in this file
* Genetic Algorithm (in other file)
"""
import heapq
import time
from MSA_DP import DP2d

from util.trace_path import resolve_path2d, resolve_path3d

scores = {
        'match': 0,
        'sub': 3,
        'gap': 2,
    }

class PriorityQueue():
    
    def __init__(self):
        self.q = []
        self.size = 0
        
    def push(self, node, prev_op, prev_point, cost, prio):
        heapq.heappush(self.q, (prio, (cost, prev_op, prev_point, node)))
        self.size += 1
        
    def pop(self):
        # we need the prio of node to add it to successors
        self.size -= 1
        return heapq.heappop(self.q)

    def empty(self):
        return self.q == []
    
def find_successors2d(point, m, n):
    '''
    grid: (m + 1) by (n + 1)
    In 2d case, a point (x, y) successor should be (x + 1, y), (x, y + 1) or (x + 1, y + 1)
    However, when x is m or y is n, some of the successors are out of bound
    '''
    x, y = point
    suc_x, suc_y, suc_xy = (x + 1, y), (x, y + 1), (x + 1, y + 1)
    if x == m:
        suc_x = None
        suc_xy = None
    if y == n:
        suc_y = None
        suc_xy = None
    return suc_x, suc_y, suc_xy

def heuristic2d(point, m, n, scores):
    x, y = point
    s_gap = scores.get('gap', 2)
    return ((m - x) - (n - y)) * s_gap

def AStar2d(seq1, seq2, scores=None):
    '''
    Using Astar algorithm in 2d case of MSA
    f(n) = g(n) + h(n)
    f() denotes cost, g() denotes cost of a path from start node
    heuristic: for node (x, y), h((x, y)) = [(g_x - x) - (g_y - y)] * s_gap
    '''
    
    m, n = len(seq1), len(seq2)
    path = {}
    close_list = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    # seq1: x
    # seq2: y
    s_match, s_sub, s_gap = scores.get('match', 0), scores.get('sub', 3), scores.get('gap', 2)
    # to start with, the first point is (0, 0) and the goal point is (m, n)
    # so the f = 0 + (m - n) * s_gap
    init_point = (0, 0)
    init_op = None
    init_cost = 0
    init_prev = None
    init_f = heuristic2d((0, 0), m, n, scores)
    # use a priority queue to maintain the order of traversing
    #! PUSH ORDER: (node, prev_op, prev_point, cost, prio):
    #! POP ORDER:  (prio, (cost, prev_op, prev_point, node))
    q = PriorityQueue()
    q.push(init_point, init_op, init_prev, init_cost, init_f)
    goal_point = (m, n)
    opt_res = -1
    traversed = 0
    while (q.empty() is False):
        cur_cost, prev_op, prev_point, cur_point = q.pop()[1]
        x, y = cur_point
        if close_list[x][y] != 0:
            continue
        # record the close list and path list
        path[cur_point] = (prev_op, prev_point)
        close_list[x][y] = 1
        
        traversed += 1
        # goal reached
        # print("curpoint:", cur_point, cur_cost, _ )
        if cur_point == goal_point:
            opt_res = cur_cost
            break
        #! WHEN PUSHING:
        #! next prio = cur_cost + thiscost + heuristic
        #! cost = cur_cost + thiscost
        suc_x, suc_y, suc_xy = find_successors2d(cur_point, m, n)
        #! PUSH ORDER: push(node, prev_op, prev_point, cost, prio):
        if suc_x is not None:
            # deal with x
            thiscost = s_gap
            h = heuristic2d(suc_x, m, n, scores)
            g = cur_cost + thiscost
            f = g + h
            q.push(suc_x, 'x', cur_point, g, f)
            # print("suc_x:", suc_x, cur_cost, thiscost, heuristic)
        if suc_y is not None:
            # deal with y
            thiscost = s_gap
            h = heuristic2d(suc_y, m, n, scores)
            g = cur_cost + thiscost
            f = g + h
            q.push(suc_y, 'y', cur_point, g, f)
            # print("suc_y:", suc_y, cur_cost, thiscost, heuristic)
        if suc_xy is not None:
            # deal with xy
            x, y = suc_xy
            thiscost = s_match if (seq1[x - 1] == seq2[y - 1]) else s_sub
            h = heuristic2d(suc_xy, m, n, scores)
            g = cur_cost + thiscost
            f = g + h
            q.push(suc_xy, 'xy', cur_point, g, f)
            # print("suc_xy:", suc_xy, cur_cost, thiscost, heuristic)
    score = opt_res
    return path, score 

def find_successors3d(point, m, n, p):
    '''
    Find all 7 successors in 3d cube
    if out of bound, the successors are set None
    '''
    x, y, z = point
    suc_x, suc_y, suc_z = (x + 1, y, z), (x, y + 1, z), (x, y, z + 1)
    suc_xy, suc_xz, suc_yz = (x + 1, y + 1, z), (x + 1, y, z + 1), (x, y + 1, z + 1)
    suc_xyz = (x + 1, y + 1, z + 1)
    if x == m:
        suc_x, suc_xy, suc_xz, suc_xyz = None, None, None, None
    if y == n:
        suc_y, suc_xy, suc_yz, suc_xyz = None, None, None, None
    if z == p:
        suc_z, suc_xz, suc_yz, suc_xyz = None, None, None, None
    return suc_x, suc_y, suc_z, suc_xy, suc_xz, suc_yz, suc_xyz


    m, n, p = len(seq1), len(seq2), len(seq3)
    tab_xy = [[float('inf') for _ in range(n + 1)] for _ in range(m + 1)]
    tab_xz = [[float('inf') for _ in range(p + 1)] for _ in range(m + 1)]
    tab_yz = [[float('inf') for _ in range(p + 1)] for _ in range(n + 1)]

    for x in range(m + 1):
        for y in range(n + 1):
            _, tab_xy[x][y] = DP2d(seq1[x: ], seq2[y: ], scores)
            
    for x in range(m + 1):
        for z in range(p + 1):
            _, tab_xz[x][z] = DP2d(seq1[x: ], seq3[z: ], scores)
            
    for y in range(n + 1):
        for z in range(p + 1):
            _, tab_yz[y][z] = DP2d(seq2[y: ], seq3[z: ], scores)

    return tab_xy, tab_xz, tab_yz

def heuristic3d(point, seq1, seq2, seq3, scores, option: str):
    if option == 'opt':
        x, y, z = point
        postfix1, postfix2, postfix3 = seq1[x: ], seq2[y: ], seq3[z: ]
        _, h1 = DP2d(postfix1, postfix2, scores)
        _, h2 = DP2d(postfix2, postfix3, scores)
        _, h3 = DP2d(postfix1, postfix3, scores)
        return h1 + h2 + h3
    
    elif option == 'L1':
        x, y, z = point
        m, n, p = len(seq1), len(seq2), len(seq3)
        L1norm = abs(m - x) + abs(n - y) + abs(p - z)
        return 2 * L1norm
    
    elif option == 'pw_L1':
        x, y, z = point
        m, n, p = len(seq1), len(seq2), len(seq3)
        s_gap = scores.get('gap', 2)
        return (abs(m - x - n + y) + abs(m - x - p + z) + abs(n - y - p + z)) * s_gap

def AStar3d(seq1, seq2, seq3, scores=None):
    '''
    Using Astar algorithm in 3d case of MSA
    f(n) = g(n) + h(n)
    f() denotes cost, g() denotes cost of a path from start node
    heuristic: for node (x, y, z) use the sum of pairwise 2D-MSA cost between
    seq1[x: ]    seq2[y: ]    seq3[z: ]
    
    #TODO: backtrack through path dict
    #TODO: remember to save a close list!
    # PUSH ORDER: (node, prev_op, prev_point, cost, prio)
    # POP ORDER:  (prio, (cost, prev_op, prev_point, node))
    '''
    
    m, n, p = len(seq1), len(seq2), len(seq3)
    path = {}
    close_list = [[[0 for _ in range(p + 1)] for _ in range(n + 1)] for _ in range(m + 1)]
    # open_list = [[[0 for _ in range(p + 1)] for _ in range(n + 1)] for _ in range(m + 1)]
    
    # seq1: x    seq2: y    seq3: z
    s_match, s_sub, s_gap = scores.get('match', 0), scores.get('sub', 3), scores.get('gap', 2)
    q = PriorityQueue()
    # PUSH ORDER: (node, prev_op, prev_point, cost, prio)
    # POP ORDER:  (prio, (cost, prev_op, prev_point, node))
    
    option = 'pw_L1'
    
    init_point = (0, 0, 0)
    init_op = None
    init_prev = None
    init_cost = 0
    init_prio = heuristic3d(init_point, seq1, seq2, seq3, scores, option)
    
    # init_prio = heuristic3d(init_point, tab_xy, tab_xz, tab_yz)
    q.push(init_point, init_op, init_prev, init_cost, init_prio)
    #! WHEN PUSHING:
    #! next prio = cur_cost + thiscost + heuristic
    #! cost = cur_cost + thiscost
    opt_res = 0
    goal_point = (m, n, p)
    traverse = 0
    
    while (q.empty() is False):
        
        cur_prio, (cur_cost, prev_op, prev_point, cur_point) = q.pop()
        x, y, z = cur_point
        
        if close_list[x][y][z] != 0:
            continue

        traverse += 1

        # record the close list and path
        
        close_list[x][y][z] = 1
        
        path[cur_point] = (prev_op, prev_point)
        
        # DEBUG: print(cur_point, cur_cost)
        
        if cur_point == goal_point:
            opt_res = cur_cost
            break
        
        suc_x, suc_y, suc_z, \
        suc_xy, suc_xz, suc_yz, \
        suc_xyz = find_successors3d(cur_point, m, n, p)
        
        if suc_x is not None:
            x, y, z = suc_x
            # if open_list[x][y][z] == 0:
                # x step, y z gap
            thiscost = 2 * s_gap
            h = heuristic3d(suc_x, seq1, seq2, seq3, scores, option)
            g = cur_cost + thiscost
            f = g + h
            q.push(suc_x, 'x', cur_point, g, f)
            # open_list[x][y][z] = 'pass'
        
        if suc_y is not None:
            x, y, z = suc_y
            # if open_list[x][y][z] == 0:
            thiscost = 2 * s_gap
            h = heuristic3d(suc_y, seq1, seq2, seq3, scores, option)
            g = cur_cost + thiscost
            f = g + h
            q.push(suc_y, 'y', cur_point, g, f)
            # open_list[x][y][z] = 'pass'
        
        if suc_z is not None:
            x, y, z = suc_z
            # if open_list[x][y][z] == 0:
            thiscost = 2 * s_gap
            h = heuristic3d(suc_z, seq1, seq2, seq3, scores, option)
            g = cur_cost + thiscost
            f = g + h
            q.push(suc_z, 'z', cur_point, g, f)
            # open_list[x][y][z] = 'pass'
        
        if suc_xy is not None:
            x, y, z = suc_xy
            # if open_list[x][y][z] == 0:
            # seq1, seq2 step, seq3 gap
            __score_xy = s_match if (seq1[x - 1] == seq2[y - 1]) else s_sub
            thiscost = 2 * s_gap + __score_xy
            h = heuristic3d(suc_xy, seq1, seq2, seq3, scores, option)
            g = cur_cost + thiscost
            f = g + h
            q.push(suc_xy, 'xy', cur_point, g, f)
            # open_list[x][y][z] = 'pass'
        
        if suc_xz is not None:
            x, y, z = suc_xz
            # if open_list[x][y][z] == 0:
            __score_xz = s_match if (seq1[x - 1] == seq3[z - 1]) else s_sub
            thiscost = 2 * s_gap + __score_xz
            h = heuristic3d(suc_xz, seq1, seq2, seq3, scores, option)
            g = cur_cost + thiscost
            f = g + h
            q.push(suc_xz, 'xz', cur_point, g, f)
            # open_list[x][y][z] = 'pass'
        
        if suc_yz is not None:
            x, y, z = suc_yz
            # if open_list[x][y][z] == 0:
            __score_yz = s_match if (seq2[y - 1] == seq3[z - 1]) else s_sub
            thiscost = 2 * s_gap + __score_yz
            h = heuristic3d(suc_yz, seq1, seq2, seq3, scores, option)
            g = cur_cost + thiscost
            f = g + h
            q.push(suc_yz, 'yz', cur_point, g, f)
            # open_list[x][y][z] = 'pass'
        
        if suc_xyz is not None:
            x, y, z = suc_xyz
            # if open_list[x][y][z] == 0:
            __score_xy = s_match if (seq1[x - 1] == seq2[y - 1]) else s_sub
            __score_xz = s_match if (seq1[x - 1] == seq3[z - 1]) else s_sub
            __score_yz = s_match if (seq2[y - 1] == seq3[z - 1]) else s_sub
            thiscost = __score_xy + __score_yz + __score_xz
            h = heuristic3d(suc_xyz, seq1, seq2, seq3, scores, option)
            g = cur_cost + thiscost
            f = g + h
            q.push(suc_xyz, 'xyz', cur_point, g, f)
            # open_list[x][y][z] = 'pass'
                
    # print("%d points traversed\n" % traverse)
    
    return path, opt_res

def test():
    seq1 = 'KJXXJAJKPXKJJXJKPXKJXXJAJKPXKJJXJKPXKJXXJAJKPXKJXXJAJKHXKJXXJAJKPXKJXXJAJKHXKJXX'
    seq2 = 'VXTLKZOKMOKAPHXHMLOWZHTPPHKPKIAXPOXKSKSWJSTSGNSHIOTTLPLLMZKUJHXTPWOWHZGAHLWKKPKMPXOTMZJUOPJ'
    seq3 = 'PJJAPJJPPJJPJJAPJJPPJJPJJAPJJPPJJPJJAPJJPPJJPJJAPJJPPJJPJJAPJJPJJKJJP'
    s = time.time()
    # path, score = AStar2d(seq1, seq2, scores)
    # seq1_o, seq2_o = resolve_path2d(seq1, seq2, path)
    # print("Seq1: %s\nSeq2: %s\nScore: %d" % (seq1, seq2, score))
    # print("Seq1_aligned %s\nSeq2_aligned %s" % (seq1_o, seq2_o))
    
    path, score = AStar3d(seq1, seq2, seq3, scores)
    seq1_o, seq2_o, seq3_o = resolve_path3d(seq1, seq2, seq3, path)
    e = time.time()
    print("Seq1: %s\nSeq2: %s\nSeq3: %s\nScore: %d" % (seq1, seq2, seq3, score))
    print("Seq1_aligned %s\nSeq2_aligned %s\nSeq3_aligned %s" % (seq1_o, seq2_o, seq3_o))
    print("TIME SPENT %.6f" % (e - s))
    
if __name__ == '__main__':
    test()
    pass