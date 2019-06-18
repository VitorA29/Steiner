import os, json
import xlsxwriter as xlsxwriter


def main():
    excel = xlsxwriter.Workbook('ResultadosLP.xlsx')
    worksheet = excel.add_worksheet()
    colunas = ["Instância",
               "Mineração",
               "Interseção",
               "Interseção fracionária",
               "Está no ótimo"
               ]
    row = 0

    for i in range(len(colunas)):
        worksheet.write(row, i, colunas[i])

    with open("solutions_pool.json") as s, open("analyze.json") as a:
        solutions = json.load(s)
        analyze = json.load(a)
        groups = next(os.walk("lp"))[1]
        results = dict()
        for group in groups:
            instances = next(os.walk("lp/" + group))[2]
            for instance in instances:
                instance = instance.split(".")[0]
                elite_edges = analyze[group][instance]["1"]["elite_edges"]
                solution_edges = solutions[group][instance]["solution"]

                row += 1
                edges_in_elite = 0
                edges_in_elite_frac = 0
                edges_in_solution = 0

                with open("lp/" + group + "/" + instance + ".out") as f:
                    for line in f:
                        x, y, cost = line[2:].replace("_", " ").replace(" - ", " ").split()
                        x = int(x) + 1
                        y = int(y) + 1
                        cost = float(cost)
                        if "X" + str(x) + "_" + str(y) in elite_edges or "X" + str(y) + "_" + str(x) in elite_edges:
                            edges_in_elite += 1
                            if cost < 0.999:
                                edges_in_elite_frac += 1
                        if "X" + str(x) + "_" + str(y) in solution_edges or "X" + str(y) + "_" + str(x) in solution_edges:
                            edges_in_solution += 1

                worksheet.write(row, 0, instance)
                worksheet.write(row, 1, len(elite_edges))
                worksheet.write(row, 2, edges_in_elite)
                worksheet.write(row, 3, edges_in_elite_frac)
                worksheet.write(row, 4, edges_in_solution)

    excel.close()

main()