from util.load_data import load_base, load_query
from util.trace_path import resolve_path2d, resolve_path3d
from MSA_DP import DP2d, DP3d
from MSA_AS import AStar2d
from MSA_gen import genetic2d
import time

scores = {
        'match': 0,
        'sub': 3,
        'gap': 2,
    }

def load_data():
    query2, query3 = load_query()
    database = load_base()
    return query2, query3, database

def run_DP2d(query2, database):
    count = 1
    start = time.time()
    totaltime = len(query2) * len(database)
    best_seq = ""
    for query in query2:
        bestscore = float('inf')
        for data in database:
            path, score = DP2d(query, data, scores)
            query_a, data_a = resolve_path2d(query, data, path)
            end = time.time()
            if score <= bestscore:
                bestscore = score
                best_seq = data
            # print("%d/%d time done DP2d, TOTAL TIME SPENT %.6f second\n" % (count, totaltime, (end - start)))
            count += 1
        print("QUERY: %s\n SEQ: %s\n BEST SCORE: %d\n" % (query, best_seq, bestscore))

def run_DP3d(query3, database):
    count = 1
    start = time.time()
    len_data = len(database)
    # C^(2)_(database)
    totaltime = len(query3) * (len_data * (len_data - 1) / 2)
    for query in query3:
        for i in range(len_data):
            for j in range(i + 1, len_data):
                seq1 = database[i]
                seq2 = database[j]
                path, score = DP3d(query, seq1, seq2, scores)
                query_a, seq1_a, seq2_a = resolve_path3d(query, seq1, seq2, path)
                end = time.time()
                print("%d/%d time done DP3d, TOTAL TIME SPENT %.6f second\n" % (count, totaltime, (end - start)))
                count += 1
                
def run_Astar2d(query2, database):
    count = 1
    start = time.time()
    totaltime = len(query2) * len(database)
    best_seq = ""
    for query in query2:
        bestscore = float('inf')
        for data in database:
            path, score = AStar2d(query, data, scores)
            query_a, data_a = resolve_path2d(query, data, path)
            end = time.time()
            if score <= bestscore:
                bestscore = score
                best_seq = data
            print("%d/%d time done AStar2d, TOTAL TIME SPENT %.6f second\n" % (count, totaltime, (end - start)))
            count += 1           

def run_gen2d(query2, database):
    count = 1
    start = time.time()
    totaltime = len(query2) * len(database)
    best_seq = ""
    for query in query2:
        bestscore = float('inf')
        for data in database:
            score, query_a, seq_a = genetic2d(query, data, scores)
            end = time.time()
            if score <= bestscore:
                bestscore = score
                best_seq = data
            print("%d/%d time done Genetic2d, TOTAL TIME SPENT %.6f second\n" % (count, totaltime, (end - start)))
            count += 1
        print("QUERY: %s\n SEQ: %s\n BEST SCORE: %d\n" % (query, best_seq, bestscore))

if __name__ == '__main__':
    query2, query3, database = load_data()
    
    # run_DP2d(query2, database)
    
    # run_DP3d(query3, database)
    
    # run_Astar2d(query2, database)
    
    # run_gen2d(query2, database)