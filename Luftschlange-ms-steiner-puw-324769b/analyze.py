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
                
                best_solution_dict[group][instance] = dict()
                best_solution_dict[group][instance]["edges"] = list()
                vertex_set = set()

                elite_match = False
                with open(used_time_stamp_folder + "/" + group + "/" + instance + "/" + best_seed + "/EliteD.txt") as elite_file:
                    for line in elite_file:
                        if "Index" in line:
                            if elite_match:
                                break
                            else:
                                splited_line = line.split(" ")
                                actual_index = splited_line[1][:-1]
                                elite_match = best_elite_index == int(actual_index)
                        elif elite_match and "e" in line:
                            edge_description = line.split(" ")
                            vertex_set.add(edge_description[1])
                            vertex_set.add(edge_description[2])
                            best_solution_dict[group][instance]["edges"].append("X"+edge_description[1]+"_"+edge_description[2])
                best_solution_dict[group][instance]["vertex"] = list(vertex_set)

                best_solution_dict[group][instance]["#viability"] = 0
                for root, dirs, files in os.walk(used_time_stamp_folder + "/" + group + "/" + instance):
                    if "patternE.txt" in files:
                        seed = root.split("/")[-1]
                        best_solution_dict[group][instance][seed] = dict()
                        with open(root + "/patternE.txt") as pattern_file:
                            line_index = 1
                            pattern_process_dict = dict()
                            for line in pattern_file:
                                A_diff_B = set()
                                B_diff_A = set()
                                intercession = set()
                                elite_edges = set(best_solution_dict[group][instance]["edges"].copy())
                                elite_vertex = set(best_solution_dict[group][instance]["vertex"].copy())

                                pattern_edges = line.split(";")[-1][:-1].split(" ")
                                for edge in elite_edges:
                                    if not edge in pattern_edges:
                                        A_diff_B.add(edge)
                                    else:
                                        intercession.add(edge)

                                for edge in pattern_edges:
                                    if not edge in elite_edges:
                                        B_diff_A.add(edge)
                                        pattern_vertex = edge[1:].split("_")
                                        if pattern_vertex[0] in elite_vertex and pattern_vertex[1] in elite_vertex:
                                            best_solution_dict[group][instance]["#viability"] += 1
                                pattern_process_dict[line_index] = dict()
                                pattern_process_dict[line_index]["#A-B"] = list(A_diff_B)
                                pattern_process_dict[line_index]["#B-A"] = list(B_diff_A)
                                pattern_process_dict[line_index]["#intercession"] = list(intercession)
                                line_index += 1

                            best_solution_dict[group][instance][seed] = pattern_process_dict
                        print("done " + seed)
                print("done " + instance)

    with open(used_time_stamp_folder + "/analyze.json", 'w') as f:
        json.dump(best_solution_dict, f, indent=4, sort_keys=True)

main()