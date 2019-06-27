import os, json
import xlsxwriter as xlsxwriter


def main():
    excel = xlsxwriter.Workbook('ResultadosLP.xlsx')
    worksheet = excel.add_worksheet()
    colunas = ["Instância",
               "Mineração",
               "Interseção",
               "Interseção fracionária",
               "Está no ótimo",
               "Tamanho da RL"
               ]
    row = 0

    for i in range(len(colunas)):
        worksheet.write(row, i, colunas[i])

    lpdict = dict()
    with open("solutions_pool.json") as s, open("analyze.json") as a:
        solutions = json.load(s)
        analyze = json.load(a)
        groups = next(os.walk("lp"))[1]
        for group in groups:
            instances = next(os.walk("lp/" + group))[2]
            for instance in instances:
                instance = instance.split(".")[0]
                elite_edges = analyze[group][instance]["1"]["elite_edges"]
                solution_edges = solutions[group][instance]["solution"]
                rl_set = set()
                list_cost = []
                row += 1
                edges_in_elite = 0
                edges_in_elite_frac = 0
                edges_in_solution = 0
                edges_in_rl = 0
                vertex_solution_set = set()
                vertex_RL_set = set()

                for e in solution_edges:
                    (a, b) = e[1:].split('_')
                    vertex_solution_set.add(a)
                    vertex_solution_set.add(b)

                with open("lp/" + group + "/" + instance + ".out") as f:
                    for line in f:
                        x, y, cost = line[2:].replace("_", " ").replace(" - ", " ").split()
                        x = int(x) + 1
                        y = int(y) + 1
                        cost = float(cost)
                        edges_in_rl += 1
                        edge = "X" + str(min(x, y)) + "_" + str(max(x, y))
                        vertex_RL_set.add(x)
                        vertex_RL_set.add(y)
                        if edge not in rl_set:
                            rl_set.add(edge)
                            list_cost.append([edge, cost])
                            if edge in elite_edges:
                                edges_in_elite += 1
                                if cost < 0.999:
                                    edges_in_elite_frac += 1
                            if edge in solution_edges:
                                edges_in_solution += 1

                lpdict[instance] = list(rl_set)
                exportDot(solution_edges, list_cost, group, instance, True)
                worksheet.write(row, 0, instance)
                worksheet.write(row, 1, len(elite_edges))
                worksheet.write(row, 2, edges_in_elite)
                worksheet.write(row, 3, edges_in_elite_frac)
                worksheet.write(row, 4, edges_in_solution)
                worksheet.write(row, 5, edges_in_rl)
    excel.close()
    with open("lp/analyze.json", 'w') as f:
        json.dump(lpdict, f, indent=4, sort_keys=True)


def exportDot(list1, list2, group, instance, frac=False):
    directory = "./lp/dot/" + group
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(directory + "/" + instance + ".dot", 'w') as f:
        f.write("graph{\n")
        f.write("node [shape=circle fontsize=16]\nedge [length=100, color=gray, fontcolor=black]\n")
        for edge in list1:
            (a, b) = edge[1:].split('_')
            f.write(a + " -- " + b + ";\n")

        f.write("edge [color=magenta, style=dashed];\n")
        for p in list2:
            if frac:
                (a, b) = p[0][1:].split('_')
                f.write(a + " -- " + b + "[label=" + str(p[1]) + "];\n")
            else:
                (a, b) = p[1:].split('_')
                f.write(a + " -- " + b + ";\n")

        f.write("}")


main()
