import sys, os, json


def main():
    best_solution_dict = dict()
    groups = next(os.walk(used_time_stamp_folder))[1]
    for group in groups:
        instances = next(os.walk(used_time_stamp_folder + "/" + group))[1]
        best_solution_dict[group] = dict()
        for instance in instances:
            best_sol_cost = 999999999
            folder = used_time_stamp_folder + "/" + group + "/" + instance
            seeds = next(os.walk(folder))[1]
            best_solution_dict[group][instance] = dict()
            for seed in seeds:
                folder = used_time_stamp_folder + "/" + group + "/" + instance + "/" + seed
                for root, dirs, files in os.walk(folder):
                    if "bestSol.txt" in files:
                        with open(folder + "/bestSol.txt", "r") as f:
                            line = f.readline()
                            size = int(line.split(' ')[1])
                            line = f.readline()
                            best_sol_cost = int(line.split(' ')[1])
                            sol = list()
                            for i in range(size):
                                line = f.readline()
                                line = line.split(' ')
                                sol.append("X" + line[1] + "_" + line[2])
                        best_solution_dict[group][instance]["solution"] = sol
                        best_solution_dict[group][instance]["cost"] = best_sol_cost
                        best_solution_dict[group][instance]["optimun"] = True

                    if "timeoutSol.txt" in files:
                        with open(folder + "/timeoutSol.txt", "r") as f:
                            line = f.readline()
                            size = int(line.split(' ')[1])
                            line = f.readline()
                            cost = int(line.split(' ')[1])
                            if cost < best_sol_cost:
                                best_sol_cost = cost
                                sol = list()
                                for i in range(size):
                                    line = f.readline()
                                    line = line.split(' ')
                                    sol.append("X" + line[1] + "_" + line[2])
                                best_solution_dict[group][instance]["solution"] = sol
                                best_solution_dict[group][instance]["cost"] = best_sol_cost
                                best_solution_dict[group][instance]["optimun"] = False

            with open(used_time_stamp_folder + "/solutions_pool.json", 'w') as f:
                    json.dump(best_solution_dict, f, indent=4, sort_keys=True)


used_time_stamp_folder = "./output/" + sys.argv[1]
with open(used_time_stamp_folder + "/log.log", 'w') as log:
    main()
