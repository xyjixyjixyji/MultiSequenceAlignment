import json

json_file_path = "hparam.json"
with open(json_file_path, 'r') as f:
    params = json.load(f)

scores = params['scores']
print(scores)