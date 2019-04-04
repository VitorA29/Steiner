import sys, os, json
from os.path import exists

def main():
    used_time_stamp_folder = "./output/" + sys.argv[1]
    best_solution_dict = dict()
    with open(used_time_stamp_folder + "/results.json") as f:
        results = json.load(f)
        for group in results.keys():
            if group.startswith("#"):
                continue
            best_solution_dict[group] = dict()
            for instance in results[group].keys():
                if instance.startswith("#"):
                    continue
                best_seed = -1
                best_elite_index = -1
                best_cost = -1
                for seed in results[group][instance]:
                    for index, cost in enumerate(results[group][instance][seed]["elite"]):
                        if best_cost == -1 or ( cost < best_cost and cost < results[group][instance][seed]["best"] ):
                            best_seed = seed
                            best_elite_index = index + 1
                            best_cost = cost
                best_solution_dict[group][instance] = set()
                with open(used_time_stamp_folder + "/" + group + "/" + instance + "/" + best_seed + "/EliteE.txt") as elite_file:
                    index = 0
                    elite_edges_value = ""
                    while index != best_elite_index:
                        elite_edges_value = elite_file.readline()
                        index += 1
                    for edge_index in elite_edges_value.split(" "):
                        if edge_index == "\n":
                            continue
                        best_solution_dict[group][instance].add(edge_index)

    for root, dirs, files in os.walk(used_time_stamp_folder):
        if root == used_time_stamp_folder:
            continue
        for f in files:
            print(root, end=" - ")
            print(f)

    for group in best_solution_dict.keys():
        for instance in best_solution_dict[group].keys():
            print(group, end=" - ")
            print(instance, end=" - ")
            print(best_solution_dict[group][instance])

main()