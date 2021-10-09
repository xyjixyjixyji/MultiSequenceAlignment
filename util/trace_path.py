def resolve_path2d(seq1, seq2, path):
    '''
    Given path dictionary, resolve the aligned sequence pair
    (which means it supports only DP2d)
    #! NOTE: path[(i, j)] = ('prev_step', (prev_x, prev_y))
    
    @param:
    seq1 & seq2: sequences used to generate the path dictionary in method DP2d()
    #! seq1: x
    #! seq2: y
    path: The dictionary saving the path taken each step in dp
    
    @output:
    seq1_out & seq2_out: Aligned Sequences
    '''
    
    m, n = len(seq1), len(seq2)
    point = (m, n)
    moves = []
    while point != (0, 0):
        prev_step, point = path[point]
        moves.insert(0, prev_step)
    # now generate the aligned sequence based on "moves"
    seq1_out, seq2_out = "", ""
    # ptr1 and ptr2 respectively points at the position waiting to be appended
    ptr1 = ptr2 = 0
    for move in moves:
        if move == 'xy': 
            seq1_out += seq1[ptr1]
            seq2_out += seq2[ptr2]
            ptr1 += 1
            ptr2 += 1
        elif move == 'x':
            seq1_out += seq1[ptr1]
            seq2_out += '-'
            ptr1 += 1
        elif move == 'y':
            seq1_out += '-'
            seq2_out += seq2[ptr2]
            ptr2 += 1
    return (seq1_out, seq2_out)

def resolve_path3d(seq1, seq2, seq3, path):
    '''
    Given path dictionary, resolve the aligned sequence pair
    (which means it supports only DP2d)
    #! NOTE: path[(x, y, z)] = ('prev_step', (prev_x, prev_y, prev_z))
    
    @param:
    seq1 & seq2 & seq3: sequences used to generate the path dictionary in method DP3d()
    #! seq1: x  seq2:y  seq3:z
    path: The dictionary saving the path taken each step in dp
    
    @output:
    seq1_out & seq2_out & seq3_out: Aligned sequences
    '''
    m, n, p = len(seq1), len(seq2), len(seq3)
    point = (m, n, p)
    moves = []
    while point != (0, 0, 0):
        prev_step, point = path[point]
        moves.insert(0, prev_step)
    seq1_o, seq2_o, seq3_o = "", "", ""
    ptr1, ptr2, ptr3 = 0, 0, 0
    for move in moves:
        if move == 'x':
            seq1_o += seq1[ptr1]
            ptr1 += 1
            seq2_o += '-'
            seq3_o += '-'
        elif move == 'y':
            seq2_o += seq2[ptr2]
            ptr2 += 1
            seq1_o += '-'
            seq3_o += '-'
        elif move == 'z':
            seq3_o += seq3[ptr3]
            ptr3 += 1
            seq1_o += '-'
            seq2_o += '-'
        elif move == 'xy':
            seq1_o += seq1[ptr1]
            ptr1 += 1
            seq2_o += seq2[ptr2]
            ptr2 += 1
            seq3_o += '-'
        elif move == 'xz':
            seq1_o += seq1[ptr1]
            ptr1 += 1
            seq3_o += seq3[ptr3]
            ptr3 += 1
            seq2_o += '-'
        elif move == 'yz':
            seq2_o += seq2[ptr2]
            ptr2 += 1
            seq3_o += seq3[ptr3]
            ptr3 += 1
            seq1_o += '-'
        elif move == 'xyz':
            seq1_o += seq1[ptr1]
            ptr1 += 1
            seq2_o += seq2[ptr2]
            ptr2 += 1
            seq3_o += seq3[ptr3]
            ptr3 += 1
        else: pass
    return seq1_o, seq2_o, seq3_o