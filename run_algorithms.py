from util.load_data import load_base, load_query
from util.trace_path import resolve_path2d, resolve_path3d
from MSA_DP import DP2d, DP3d
from MSA_AS import AStar2d, AStar3d
from MSA_gen import genetic2d, genetic3d
from tqdm import tqdm
import time
import json

# load json files, mainly for genetic's hyperparams
json_file_path = "hparam.json"
with open(json_file_path, 'r') as f:
    params = json.load(f)

scores = params['scores']

def run_DP2d(query2, database):
    count = 1
    start = time.time()
    totaltime = len(query2) * len(database)
    
    best_seq = ""
    best_seq_a = ""
    best_query_a = ""
    
    filename = "./results/DP2d_res.txt"
    filename_opt = "./results/DP2d_opt.txt"
    # print_count = 0
    
    with open(filename, 'w') as f1:
        with open(filename_opt, 'w') as f2:
            for query in query2:
                bestscore = float('inf')
                for data in tqdm(database):
                    start_each = time.time()
                    path, score = DP2d(query, data, scores)
                    query_a, data_a = resolve_path2d(query, data, path)
                    end = time.time()
                    
                    if score <= bestscore:
                        bestscore = score
                        best_seq = data
                        best_query_a = query_a
                        best_seq_a = data_a
                
                    # print some info for convenience
                    # if print_count % print_every2d == 0:
                        # print("%d/%d time done DP2d, TOTAL TIME SPENT %.6f second\n" % (count, totaltime, (end - start)))
                        
                    count += 1
                    # print_count += 1
                    
                    # write the matching info to file
                    content = "Query: %s\nSeq: %s\nScore: %d\nTime Expired: %.6f second\n" \
                    % (query_a, data_a, score, (end - start_each))
                    
                    f1.write(content + '\n')
            
                # print("QUERY: %s\n SEQ: %s\n BEST SCORE: %d\n" % (query, best_seq, bestscore))
                
                content = "QUERY: %s\nSEQ: %s\nBEST SCORE: %d\nTOTAL TIME EXPIRED: %.6f\nALIGNED:\n%s\n%s\n" \
                % (query, best_seq, bestscore, (end - start), best_seq_a, best_query_a)
                
                f2.write(content + '\n')

def run_DP3d(query3, database):
    count = 1
    start = time.time()
    len_data = len(database)
    # C^(2)_(database)
    totaltime = len(query3) * (len_data * (len_data - 1) / 2)
    
    best_seq1 = ""
    best_seq2 = ""
    best_seq1_a = ""
    best_seq2_a = ""
    best_query_a = ""
    
    filename = "./results/DP3d_res.txt"
    filename_opt = "./results/DP3d_opt.txt"
    
    with open(filename, 'w') as f1:
        with open(filename_opt, 'w') as f2:
            for query in query3:
                best_score = float('inf')
                for i in tqdm(range(len_data)):
                    for j in range(i + 1, len_data):
                        start_each = time.time()
                        seq1 = database[i]
                        seq2 = database[j]
                        path, score = DP3d(query, seq1, seq2, scores)
                        query_a, seq1_a, seq2_a = resolve_path3d(query, seq1, seq2, path)
                        end = time.time()
                        
                        if score <= best_score:
                            best_score = score
                            best_seq1 = seq1
                            best_seq2 = seq2
                            best_query_a = query_a
                            best_seq1_a = seq1_a
                            best_seq2_a = seq2_a
                        
                        # print("%d/%d time done DP3d, TOTAL TIME SPENT %.6f second\n" % (count, totaltime, (end - start)))
                        count += 1
                        
                        # write the info to file
                        content = "Query: %s\nSeq1: %s\nSeq2: %s\nScore: %d\nTime Expired: %.6f second\n" % (
                            query_a, seq1_a, seq2_a, score, (end - start_each)
                        )
                        f1.write(content + '\n')
                content = "Query: %s\nSeq1: %s\nSeq2: %s\nScore: %d\nTotal Time Expired: %.6f second\nALIGNED:\n\%s\n%s\n%s\n" % (
                            query, best_seq1, best_seq2, best_score, (end - start),
                            best_query_a, best_seq1_a, best_seq2_a)
                f2.write(content + '\n')
                                      
def run_Astar2d(query2, database):
    count = 1
    start = time.time()
    totaltime = len(query2) * len(database)
    
    best_seq = ""
    best_seq_a = ""
    best_query_a = ""
    
    print_count = 0
    filename = "./results/AStar2d_res.txt"
    filename_opt = "./results/AStar2d_opt.txt"
    
    with open(filename, 'w') as f1:
        with open(filename_opt, 'w') as f2:
            for query in query2:
                bestscore = float('inf')
                for data in tqdm(database):
                    start_each = time.time()
                    path, score = AStar2d(query, data, scores)
                    query_a, data_a = resolve_path2d(query, data, path)
                    end = time.time()
                    
                    if score <= bestscore:
                        bestscore = score
                        best_seq = data
                        best_query_a = query_a
                        best_seq_a = data_a
                        
                    # if print_count % print_every2d == 0:
                    #     print("%d/%d time done AStar2d, TOTAL TIME SPENT %.6f second\n" % (count, totaltime, (end - start)))
                    
                    count += 1
                    print_count += 1
                    
                    # write the matching info to file
                    content = "Query: %s\nSeq: %s\nScore: %d\nTime Expired: %.6f second\n" % (query_a, data_a, score, (end - start_each))
                    f1.write(content + '\n')
                
                # print("QUERY: %s\n SEQ: %s\n BEST SCORE: %d\n" % (query, best_seq, bestscore))
                
                content = "QUERY: %s\nSEQ: %s\nBEST SCORE: %d\nTOTAL TIME EXPIRED: %.6f\nALIGNED:\n%s\n%s\n" \
                % (query, best_seq, bestscore, (end - start), best_seq_a, best_query_a)
                
                f2.write(content + '\n')

def run_Astar3d(query3, database):
    count = 1
    start = time.time()
    len_data = len(database)
    # C^(2)_(database)
    totaltime = len(query3) * (len_data * (len_data - 1) / 2)
    
    best_seq1 = ""
    best_seq2 = ""
    best_seq1_a = ""
    best_seq2_a = ""
    best_query_a = ""
    
    filename = "./results/AStar3d_res.txt"
    filename_opt = "./results/AStar3d_opt.txt"
    
    with open(filename, 'w') as f1:
        with open(filename_opt, 'w') as f2:
            for query in query3:
                best_score = float('inf')
                for i in tqdm(range(len_data)):
                    for j in range(i + 1, len_data):
                        start_each = time.time()
                        seq1 = database[i]
                        seq2 = database[j]
                        path, score = AStar3d(query, seq1, seq2, scores)
                        query_a, seq1_a, seq2_a = resolve_path3d(query, seq1, seq2, path)
                        end = time.time()
                        
                        if score <= best_score:
                            best_score = score
                            best_seq1 = seq1
                            best_seq2 = seq2
                            best_query_a = query_a
                            best_seq1_a = seq1_a
                            best_seq2_a = seq2_a
                        
                        # print("%d/%d time done DP3d, TOTAL TIME SPENT %.6f second\n" % (count, totaltime, (end - start)))
                        count += 1
                        
                        # write the info to file
                        content = "Query: %s\nSeq1: %s\nSeq2: %s\nScore: %d\nTime Expired: %.6f second\n" % (
                            query_a, seq1_a, seq2_a, score, (end - start_each)
                        )
                        f1.write(content + '\n')
                content = "Query: %s\nSeq1: %s\nSeq2: %s\nScore: %d\nTotal Time Expired: %.6f second\nALIGNED:\n\%s\n%s\n%s\n" % (
                            query, best_seq1, best_seq2, best_score, (end - start),
                            best_query_a, best_seq1_a, best_seq2_a)
                f2.write(content + '\n')
        
def run_gen2d(query2, database):
    count = 1
    start = time.time()
    totaltime = len(query2) * len(database)
    
    best_seq = ""
    best_seq_a = ""
    best_query_a = ""
    
    filename = "./results/Genetic2d_res.txt"
    filename_opt = "./results/Genetic2d_opt.txt"
    
    # parameters for genetic2d()
    kwargs = params["genetic2d"]
    
    with open(filename, 'w') as f1:
        with open(filename_opt, 'w') as f2:
            for query in query2:
                bestscore = float('inf')
                for seq in tqdm(database):
                    start_each = time.time()
                    score, query_a, seq_a = genetic2d(query, seq, kwargs, scores)
                    end = time.time()
                    
                    if score <= bestscore:
                        bestscore = score
                        best_seq = seq
                        best_query_a = query_a
                        best_seq_a = seq_a
                        
                    # print("%d/%d time done Genetic2d, TOTAL TIME SPENT %.6f second\n" % (count, totaltime, (end - start)))
                    count += 1
                    
                    # write matching info to file
                    content = "Query: %s\nSeq: %s\nScore: %d\nTime Expired: %.6f second\n" \
                    % (query_a, seq_a, score, (end - start_each))
                    
                    f1.write(content + '\n')
                # print("QUERY: %s\n SEQ: %s\n BEST SCORE: %d\n" % (query, best_seq, bestscore))
                
                content = "QUERY: %s\nSEQ: %s\nBEST SCORE: %d\nTOTAL TIME EXPIRED: %.6f\nALIGNED:\n%s\n%s\n" \
                % (query, best_seq, bestscore, (end - start), best_seq_a, best_query_a)
                
                f2.write(content + '\n')
                
def run_gen3d(query3, database):
    count = 1
    start = time.time()
    len_data = len(database)
    # C^(2)_(database)
    totaltime = len(query3) * (len_data * (len_data - 1) / 2)
    
    best_seq1 = ""
    best_seq2 = ""
    best_seq1_a = ""
    best_seq2_a = ""
    best_query_a = ""
    
    filename = "./results/Genetic3d_res.txt"
    filename_opt = "./results/Genetic3d_opt.txt"
    
    kwargs = params["genetic3d"]
    
    with open(filename, 'w') as f1:
        with open(filename_opt, 'w') as f2:
            for query in query3:
                best_score = float('inf')
                for i in tqdm(range(len_data)):
                    for j in range(i + 1, len_data):
                        start_each = time.time()
                        seq1 = database[i]
                        seq2 = database[j]
                        score, query_a, seq1_a, seq2_a = genetic3d(query, seq1, seq2, kwargs, scores)
                        end = time.time()
                        
                        if score <= best_score:
                            best_score = score
                            best_seq1 = seq1
                            best_seq2 = seq2
                            best_query_a = query_a
                            best_seq1_a = seq1_a
                            best_seq2_a = seq2_a
                        
                        # print("%d/%d time done DP3d, TOTAL TIME SPENT %.6f second\n" % (count, totaltime, (end - start)))
                        count += 1
                        
                        # write the info to file
                        content = "Query: %s\nSeq1: %s\nSeq2: %s\nScore: %d\nTime Expired: %.6f second\n" % (
                            query_a, seq1_a, seq2_a, score, (end - start_each)
                        )
                        f1.write(content + '\n')
                content = "Query: %s\nSeq1: %s\nSeq2: %s\nScore: %d\nTotal Time Expired: %.6f second\nALIGNED:\n\%s\n%s\n%s\n" % (
                            query, best_seq1, best_seq2, best_score, (end - start),
                            best_query_a, best_seq1_a, best_seq2_a)
                f2.write(content + '\n')
