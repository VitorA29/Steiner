import sys, os, json
from os.path import exists

def main():
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
                    best_solution_dict[group][instance][seed]["#inviability"] = 0
                    if "patternE.txt" in files:
                        with open(root + "/patternE.txt") as pattern_file:
                            line_index = 1
                            pattern_process_dict = dict()
                            pattern_union = set()
                            vertices = elite_vertex_set.copy()

                            for line in pattern_file:
                                pattern_edges = set(line.split(";")[-1][:-1].split(" "))

                                for edge in pattern_edges:
                                    pattern_union.add(edge)
                                    if not edge in elite_edges_set:
                                        pattern_vertex = edge[1:].split("_")
                                        if pattern_vertex[0] in vertices and pattern_vertex[1] in vertices:
                                            best_solution_dict[group][instance][seed]["#inviability"] += 1
                                            log.write("inviability " + "group " + group + " Instance " + instance + " edge " + edge + "\n")
                                        else:
                                            try:
                                                vertices.add(pattern_vertex[0])
                                                vertices.add(pattern_vertex[1])

                                            except:
                                                #print(pattern_vertex)
                                                log.write("group " + group + "Instance " + instance + "edge " + edge)
                                pattern_process_dict[line_index] = dict()
                                pattern_process_dict[line_index]["#Elite-Pattern"] = list(elite_edges_set.difference(pattern_edges))
                                pattern_process_dict[line_index]["#Pattern-Elite"] = list(pattern_edges.difference(elite_edges_set))
                                pattern_process_dict[line_index]["#Intercession"] = list(elite_edges_set.intersection(pattern_edges))
                                line_index += 1
                            best_solution_dict[group][instance][seed]["pattern"] = pattern_process_dict
                            best_solution_dict[group][instance][seed]["pattern_union"] = list(pattern_union)
                        print("done " + seed)
                        exportDot(best_solution_dict[group][instance][seed], group, instance, seed)
                print("done " + instance)

    with open(used_time_stamp_folder + "/analyze.json", 'w') as f:
        json.dump(best_solution_dict, f, indent=4, sort_keys=True)


def exportDot(solution, group, instance, seed):
    cores = ["magenta", "aqua", "green", "yellow", "darkblue", "lime", "red", "darkcyan", "darkred", "orange"];
    folder_name = used_time_stamp_folder + "/" + group + "/" + instance + "/" + seed
    with open(folder_name + "/solution.dot", 'w') as f:
        f.write("graph{\n")
        f.write("node [shape=circle fontsize=16]\nedge [length=100, color=gray, fontcolor=black]\n")
        for edge in solution["elite_edges"]:
            (a, b) = edge[1:].split('_')
            f.write(a + " -- " + b + ";\n")

        with open(folder_name + "/terminals.txt") as t:
            line = t.read()
            line = line[:-1].split(" ")
            solution["terminals"] = line
            for n in line:
                f.write(n + " [fontcolor=white, color=red];\n")

        pattern = solution["pattern_union"]
        for p in pattern:
            try:
                (a, b) = p[1:].split('_')
            except:
                print("Erro!" + p)
            f.write(a + " -- " + b + "[color=" + cores[0] + "];\n")

        f.write("}")

def gerar_estatistica():
    total_instances = 0
    total_solution_seeds = 0
    solutions_with_inviability = 0
    num_optimum = 0
    avg_edge_patterns = 0
    avg_solution_size_vertice = 0
    avg_steiner_vertices = 0
    avg_diff_steiner_vertices = 0
    avg_diff_edges = 0
    avg_equals_edges = 0
    avg_equals_steiner_vertices = 0

    with open(used_time_stamp_folder + "/analyze.json") as f:
        results = json.load(f)
        for group in results.keys():
            for instance in results[group].keys():
                total_solution_seeds += 1
                for seed in results[group][instance].keys():
                    result = results[group][instance][seed]
                    total_instances += 1
                    if result["#inviability"] > 0:
                        solutions_with_inviability += 1
                    try:
                        avg_edge_patterns += len(result["pattern"])
                        avg_solution_size_vertice += len(result["elite_vertex"])
                        avg_steiner_vertices += len(result["elite_vertex"]) - len(result["terminals"])
                        avg_diff_edges += len(result["#Elite-Pattern"])
                        avg_equals_edges += len(result["elite_edges"]) - avg_diff_edges

                        vertices = set(result["elite_vertex"])
                        steiner_vertices = vertices.difference(result["terminals"])
                        avg_diff_steiner_vertices += len(steiner_vertices)
                    except:
                        pass

    with open(used_time_stamp_folder + "/deviation.json") as f:
        deviation = json.load(f)
        for group in deviation.keys():
            if group.startswith("#"):
                continue
            for instance in results[group].keys():
                if deviation[group][instance] == 0.0:
                    num_optimum += 1

    with open(used_time_stamp_folder + "/statistics.txt", 'w') as f:
        f.write("Total solutions: " + str(total_solution_seeds))
        f.write("\nNumber of instances: " + str(total_instances))
        f.write("\nSolutions with inviability: " + str(solutions_with_inviability))
        f.write("\nOptimum solutions: " + str(num_optimum))
        f.write("\nAverage number of patterns: " + str(avg_edge_patterns/total_instances))
        f.write("\nAverage length of solution: " + str(avg_solution_size_vertice/total_instances))
        f.write("\nAverage number of steiner vertices: " + str(avg_steiner_vertices/total_instances))
        f.write("\nAverage difference steiner vertices: " + str(avg_diff_steiner_vertices/total_instances))
        f.write("\nAverage equals edges: " + str(avg_equals_steiner_vertices/total_instances))


used_time_stamp_folder = "./output/" + sys.argv[1]
with open(used_time_stamp_folder + "/log.log", 'w') as log:
    main()
    gerar_estatistica()