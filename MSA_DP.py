"""
Algorithms implemented for Multiple Sequence Alignment (MSA)
* Dynamic Programming
    * 2-D DP and 3-D DP are implemented in this file
* A* Algorithm (in other file)
* Genetic Algorithm (in other file)
"""
import time
from util.trace_path import resolve_path2d, resolve_path3d

scores = {
        'match': 0,
        'sub': 3,
        'gap': 2,
    }


def DP2d(seq1, seq2, scores: dict):
    '''
    Using dynamic programming to solve 2d case of sequence alignment
    Finding the best alignment between **2** sequences
    
    @param:
    seq1 & seq2: first sequence and second sequence waiting for alignment
    #! NOTE: no gap shall exists in original sequence unaligned
    
    @output:
    dp (array)  : the aligning array
    path (dict) : the path taken
    #! NOTE: path is a dictionary, using (i, j) to index, outputs the previous move and coordinate
    '''
    
    m, n = len(seq1), len(seq2)
    path = {}
    # seq1: x
    # seq2: y
    s_match = scores.get('match', 0)
    s_sub = scores.get('sub', 3)
    s_gap = scores.get('gap', 2)
    # since min() method shall be used, dp is initialized to +inf
    # each grid represents a letter to be matched
    dp = [[float('inf') for _ in range(n + 1)] for _ in range(m + 1)]
    dp[0][0] = 0
    # initialize the first row and column
    for x in range(m + 1): 
        dp[x][0] = x * s_gap
        path[(x, 0)] = ('x', (x - 1, 0))
    for y in range(n + 1): 
        dp[0][y] = y * s_gap
        path[(0, y)] = ('y', (0, y - 1))
    
    # left to right, up to down
    for x in range(1, m + 1):
        for y in range(1, n + 1):
            __score_xy = s_match if (seq1[x - 1] == seq2[y - 1]) else s_sub
            val_x = dp[x - 1][y] + s_gap
            val_y = dp[x][y - 1] + s_gap
            val_xy = dp[x - 1][y - 1] + __score_xy
            
            dp[x][y] = min(val_x, val_y, val_xy)
            # save the move taken
            if dp[x][y] == val_x:
                path[(x, y)] = ('x', (x - 1, y))
            elif dp[x][y] == val_y:
                path[(x, y)] = ('y', (x, y - 1))
            elif dp[x][y] == val_xy:
                path[(x, y)] = ('xy', (x - 1, y - 1))
            else: pass
            
    return (path, dp[-1][-1])

def output2d(seq1, seq2, seq1_o, seq2_o, score):
    # output function
    print("==============================")
    print("Sequence 1: %s" % seq1)
    print("Sequence 2: %s" % seq2)
    print("After alignment")
    print("Aligned sequence 1: %s" % seq1_o)
    print("Aligned sequence 2: %s" % seq2_o)
    print("The score is %d" % score)
    print("==============================")

def DP3d(seq1, seq2, seq3, scores: dict=None):
    '''
    Using dynamic programming to solve 3d case of sequence alignment
    Finding the best alignment between **3** sequences 
    
    @param:
    seq1 & seq2 & seq3: seqs to be aligned
    #! NOTE: no gap shall exists in original sequence unaligned
    
    @output:
    path (dict) : the path taken
    score: alignment score
    #! NOTE: path is a dictionary, using (x, y, z) to index, outputs tuple (prev_move, prev_point)
    '''
    m, n, p = len(seq1), len(seq2), len(seq3)
    s_match, s_sub, s_gap = scores.get('match', 0), scores.get('sub', 3), scores.get('gap', 2)
    path = {}
    
    dp = [[[float('inf') for _ in range(p + 1)] for _ in range(n + 1)]for _ in range(m + 1)] # (m + 1) x (n + 1) x (p + 1)
    dp[0][0][0] = 0
    
    ######## init axis ########
    # x axis
    for x in range(1, m + 1):
        dp[x][0][0] = 2 * x * s_gap
        path[(x, 0, 0)] = ('x', (x - 1, 0, 0))
    # y axis
    for y in range(1, n + 1):
        dp[0][y][0] = 2 * y * s_gap
        path[(0, y, 0)] = ('y', (0, y - 1, 0))
    # z axis
    for z in range(1, p + 1):
        dp[0][0][z] = 2 * z * s_gap
        path[(0, 0, z)] = ('z', (0, 0, z - 1))
    ############################
    
    ######## init planes ########
    # xy plane
    for x in range(1, m + 1):
        for y in range(1, n + 1):
            __score_xy = s_match if (seq1[x - 1] == seq2[y - 1]) else s_sub
            # step in seq1, seq2 and seq3 are gaps
            val_x = dp[x - 1][y][0] + 2 * s_gap
            # step in seq2, seq1 and seq3 are gaps
            val_y = dp[x][y - 1][0] + 2 * s_gap
            # step in both seq1 and seq2, seq3 is gap
            val_xy = dp[x - 1][y - 1][0] + __score_xy + 2 * s_gap
            # select the minimal value
            dp[x][y][0] = min(val_x, val_y, val_xy)
            # record the path
            if dp[x][y][0] == val_x:
                path[(x, y, 0)] = ('x', (x - 1, y, 0))
            elif dp[x][y][0] == val_y:
                path[(x, y, 0)] = ('y', (x, y - 1, 0))
            elif dp[x][y][0] == val_xy:
                path[(x, y, 0)] = ('xy', (x - 1, y - 1, 0))
            else: pass # not gonna happen, code in this form for clarity
    # xz plane
    for x in range(1, m + 1):
        for z in range(1, p + 1):
            __score_xz = s_match if (seq1[x - 1] == seq3[z - 1]) else s_sub
            val_x = dp[x - 1][0][z] + 2 * s_gap
            val_z = dp[x][0][z - 1] + 2 * s_gap
            val_xz = dp[x - 1][0][z - 1] + __score_xz + 2 * s_gap
            # select the minimal value
            dp[x][0][z] = min(val_x, val_z, val_xz)
            # record the path
            if dp[x][0][z] == val_x:
                path[(x, 0, z)] = ('x', (x - 1, 0, z))
            elif dp[x][0][z] == val_z:
                path[(x, 0, z)] = ('z', (x, 0, z - 1))
            elif dp[x][0][z] == val_xz:
                path[(x, 0, z)] = ('xz', (x - 1, 0, z - 1))
            else: pass # not gonna happen, code in this form for clarity
    # yz plane
    for y in range(1, n + 1):
        for z in range(1, p + 1):
            __score_yz = s_match if (seq2[y - 1] == seq3[z - 1]) else s_sub
            val_y = dp[0][y - 1][z] + 2 * s_gap
            val_z = dp[0][y][z - 1] + 2 * s_gap
            val_yz = dp[0][y - 1][z - 1] + __score_yz + 2 * s_gap
            # select the minimal value
            dp[0][y][z] = min(val_y, val_z, val_yz)
            # record the path
            if dp[0][y][z] == val_y:
                path[(0, y, z)] = ('y', (0, y - 1, z))
            elif dp[0][y][z] == val_z:
                path[(0, y, z)] = ('z', (0, y, z - 1))
            elif dp[0][y][z] == val_yz:
                path[(0, y, z)] = ('yz', (0, y - 1, z - 1))
            else: pass # not gonna happen, code in this form for clarity
    #############################
    # end init in 1d and 2d, do 3d dp now
    for x in range(1, m + 1):
        for y in range(1, n + 1):
            for z in range(1, p + 1):
                __score_xy = s_match if (seq1[x - 1] == seq2[y - 1]) else s_sub
                __score_xz = s_match if (seq1[x - 1] == seq3[z - 1]) else s_sub
                __score_yz = s_match if (seq2[y - 1] == seq3[z - 1]) else s_sub
                # val_x: step in seq1, seq2 and seq3 gaps
                val_x = dp[x - 1][y][z] + 2 * s_gap
                val_y = dp[x][y - 1][z] + 2 * s_gap
                val_z = dp[x][y][z - 1] + 2 * s_gap
                # val_xy: step in seq1 and seq2, seq3 gaps
                val_xy = dp[x - 1][y - 1][z] + __score_xy + 2 * s_gap
                val_xz = dp[x - 1][y][z - 1] + __score_xz + 2 * s_gap
                val_yz = dp[x][y - 1][z - 1] + __score_yz + 2 * s_gap
                # val_xyz: step in both 3 axises
                val_xyz = dp[x - 1][y - 1][z - 1] + __score_xy + __score_xz + __score_yz
                dp[x][y][z] = min(val_x, val_y, val_z, val_xy, val_xz, val_yz, val_xyz)
                # record in path
                if dp[x][y][z] == val_x:
                    path[(x, y, z)] = ('x', (x - 1, y , z))
                elif dp[x][y][z] == val_y:
                    path[(x, y, z)] = ('y', (x, y - 1, z))
                elif dp[x][y][z] == val_z:
                    path[(x, y, z)] = ('z', (x, y, z - 1))
                elif dp[x][y][z] == val_xy:
                    path[(x, y, z)] = ('xy', (x - 1, y - 1, z))
                elif dp[x][y][z] == val_xz:
                    path[(x, y, z)] = ('xz', (x - 1, y, z - 1))
                elif dp[x][y][z] == val_yz:
                    path[(x, y, z)] = ('yz', (x, y - 1, z - 1))
                elif dp[x][y][z] == val_xyz:
                    path[(x, y, z)] = ('xyz', (x - 1, y - 1, z - 1))
                else: pass
    return path, dp[-1][-1][-1]

def output3d(seq1, seq2, seq3, seq1_o, seq2_o, seq3_o, score):
    print("==============================")
    print("Sequence 1: %s" % seq1)
    print("Sequence 2: %s" % seq2)
    print("Sequence 3: %s" % seq3)
    print("After alignment")
    print("Aligned sequence 1: %s" % seq1_o)
    print("Aligned sequence 2: %s" % seq2_o)
    print("Aligned sequence 3: %s" % seq3_o)
    print("The score is %d" % score)
    print("==============================")

def test():
    seq1 = 'KJXXJAJKPXKJJXJKPXKJXXJAJKPXKJJXJKPXKJXXJAJKPXKJXXJAJKHXKJXXJAJKPXKJXXJAJKHXKJXX'
    seq2 = 'VXTLKZOKMOKAPHXHMLOWZHTPPHKPKIAXPOXKSKSWJSTSGNSHIOTTLPLLMZKUJHXTPWOWHZGAHLWKKPKMPXOTMZJUOPJ'
    seq3 = 'PJJAPJJPPJJPJJAPJJPPJJPJJAPJJPPJJPJJAPJJPPJJPJJAPJJPPJJPJJAPJJPJJKJJP'
    # seq3 ='IWTJBGTJGJTWGBJTPKHAXHAGJJXJJKPJTPJHJHJHJHJHJHJHJHJHKUTJJUWXHGHHGALKLPJTPJPGVXPLBJHH'
    # start = time.time()
    # path, score = DP2d(seq1, seq2, scores)
    # seq1_o, seq2_o = resolve_path2d(seq1, seq2, path)
    # end = time.time()
    # output2d(seq1, seq2, seq1_o, seq2_o, score)
    # print("TIME SPENT: %.6f second" % (end - start))

    start = time.time()
    path, score = DP3d(seq1, seq2, seq3, scores)
    seq1_o, seq2_o, seq3_o = resolve_path3d(seq1, seq2, seq3, path)
    end = time.time()
    output3d(seq1, seq2, seq3, seq1_o, seq2_o, seq3_o, score)
    print("TIME SPENT: %.6f second" % (end - start))

if __name__ == '__main__':
    test()
    pass