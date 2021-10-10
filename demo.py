from run_algorithms import run_DP2d, run_DP3d, run_Astar2d, run_Astar3d, run_gen2d, run_gen3d
from util.load_data import load_data
if __name__ == '__main__':
    
    query2, query3, database = load_data()
    
    # run_DP2d(query2, database)
    
    # run_DP3d(query3, database)
    
    # run_Astar2d(query2, database)
    
    # run_Astar3d(query3, database)
    
    # run_gen2d(query2, database)
    
    run_gen3d(query3, database)