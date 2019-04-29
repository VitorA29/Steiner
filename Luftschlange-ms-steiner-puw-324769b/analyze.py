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

                best_solution_dict[group][instance] = dict()
                for root, dirs, files in os.walk(used_time_stamp_folder + "/" + group + "/" + instance):
                    if len(root.split('/')) == 5:
                        continue

                    seed = root.split("/")[-1]
                    best_solution_dict[group][instance][seed] = dict()

                    best_elite_index = -1
                    best_cost = -1
                    for index, cost in enumerate(results[group][instance][seed]["elite"]):
                        if best_cost == -1 or ( cost < best_cost and cost < results[group][instance][seed]["best"] ):
                            best_elite_index = index + 1
                            best_cost = cost

                    best_solution_dict[group][instance][seed]["elite_edges"] = list()
                    elite_vertex_set = set()
                    elite_edges_set = set()
                    with open(root + "/EliteD.txt") as elite_file:
                        elite_match = False
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
                                elite_vertex_set.add(edge_description[1])
                                elite_vertex_set.add(edge_description[2])
                                edge = "X"+edge_description[1]+"_"+edge_description[2]
                                elite_edges_set.add(edge)
                                best_solution_dict[group][instance][seed]["elite_edges"].append(edge)

                    best_solution_dict[group][instance][seed]["elite_vertex"] = list(elite_vertex_set)
                    best_solution_dict[group][instance][seed]["#viability"] = 0
                    if "patternE.txt" in files:
                        with open(root + "/patternE.txt") as pattern_file:
                            line_index = 1
                            pattern_process_dict = dict()
                            pattern_union = set()
                            for line in pattern_file:
                                pattern_edges = set(line.split(";")[-1][:-1].split(" "))

                                for edge in pattern_edges:
                                    pattern_union.add(edge)
                                    if not edge in elite_edges_set:
                                        pattern_vertex = edge[1:].split("_")
                                        if pattern_vertex[0] in elite_vertex_set and pattern_vertex[1] in elite_vertex_set:
                                            best_solution_dict[group][instance][seed]["#viability"] += 1
                                pattern_process_dict[line_index] = dict()
                                pattern_process_dict[line_index]["#Elite-Pattern"] = list(elite_edges_set.difference(pattern_edges))
                                pattern_process_dict[line_index]["#Pattern-Elite"] = list(pattern_edges.difference(elite_edges_set))
                                pattern_process_dict[line_index]["#Intercession"] = list(elite_edges_set.intersection(pattern_edges))
                                line_index += 1
                            best_solution_dict[group][instance][seed]["pattern"] = pattern_process_dict
                            best_solution_dict[group][instance][seed]["pattern_union"] = list(pattern_union)
                        print("done " + seed)
                print("done " + instance)

    with open(used_time_stamp_folder + "/analyze.json", 'w') as f:
        json.dump(best_solution_dict, f, indent=4, sort_keys=True)

main()