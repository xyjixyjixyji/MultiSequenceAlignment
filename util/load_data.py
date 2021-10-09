path_query = "./data/MSA_query.txt"
path_database = "./data/MSA_database.txt"

def load_query():
    query2 = []
    query3 = []
    flag = 0
    with open(path_query, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line == "2\n":   flag = 2
            elif line == "3\n": flag = 3
            else:
                if flag == 2:
                    query2.append(line[:-1])
                elif flag == 3:
                    query3.append(line[:-1])
    return query2, query3

def load_base():
    database = []
    with open(path_database, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            else:
                database.append(line[:-1])
    return database