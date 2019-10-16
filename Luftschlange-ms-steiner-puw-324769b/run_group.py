import os
from os.path import isfile, join
import datetime, sys, time, json, zipfile
from pprint import pprint
from collections import OrderedDict

def main():
    #Getting timestamp
    dts = datetime.datetime.utcnow()
    timestamp = round(time.mktime(dts.timetuple()) + dts.microsecond/1e6)
    os.system("mkdir output/" + str(timestamp))

    #Getting start time and iterations
    details_report = dict()
    time_format = "%d:%H:%M:%S"
    start_time = time.time()
    details_report["iterations"] = int(sys.argv[2])
    details_report["max_seed"] = int(sys.argv[3])

    # Some parameters initialization
    max_seed = sys.argv[3]
    data_mining = sys.argv[4]
    early_stop = False
    if not sys.argv[5] == 0:
        early_stop = True

    # prepare model folder
    if data_mining == '1':
        os.system("mkdir output/" + str(timestamp) + "/model")

    #Opening instances folder
    report_dict = dict()
    deviation_dict = dict()
    elite_dict = dict()
    for folder in os.listdir("instances/" + sys.argv[1]):
        start_folder_time = time.time()

        path = sys.argv[1] + "/" + folder

        #Preparing output folder
        output_folder = str(timestamp) + "/" + folder
        os.system("mkdir output/" + output_folder)

        #Getting instances files
        onlyfiles = [f for f in os.listdir("instances/" + path) if isfile(join("instances/" + path, f))]

        # Getting the best value information
        opt = dict()
        with open("instances/" + path + "/opt.json") as fopt:
                opt = json.load(fopt)

        #Executing instances
        report = dict()
        for instance in onlyfiles:
            instance_name = instance.split(".")[0]
            if ".json" in instance:
                continue

            bash_cmd = "./bin/run_instance " + path + "/" + instance + " -seed " + max_seed + " -maxit " + sys.argv[2] + " -mine " + data_mining + " -folder " + output_folder

            # Check if should early stop
            if early_stop:
                bash_cmd += " -bestbound " + str(opt[instance_name])

            os.system(bash_cmd)
            print("output/" + output_folder + "/" + instance_name + "/" + instance_name + ".json")
            with open("output/" + output_folder + "/" + instance_name + "/" + instance_name + ".json") as f:
                data = json.load(f)
                report[instance_name] = data
            os.system("rm output/" + output_folder + "/" + instance_name + "/" + instance_name + ".json")

            if data_mining == '1':
                # prepare the model file name
                model_file_name = "output/" + str(timestamp) + "/model/" + instance_name + ".eu"

                # find the best cost among all seeds
                best_cost = -1
                for seed in report[instance_name].keys():
                    for cost in report[instance_name][str(seed)]["elite"]:
                        if best_cost == -1 or ( cost < best_cost and cost < report[instance_name][str(seed)]["best"] ):
                            best_cost = cost

                with open(model_file_name, 'a') as f:
                    f.write(str(best_cost) + "\n")

                keep_file = False

                #Building pattern instance file
                for pattern_file_marker in ["E", "V"]:
                    pattern_dict = dict()

                    # creating the array of all patterns
                    pattern_array = list()

                    for seed in range (1, int(max_seed) + 1):
                        try:
                            with open("output/" + str(timestamp) + "/" + folder + "/" + instance_name + "/" + str(seed) + "/pattern" + pattern_file_marker + ".txt") as pattern_file:
                                line = pattern_file.readline()
                                while line:
                                    line_split = line.split(";")
                                    if(len(line_split) == 3):
                                        # adding the line to the patterns array(removing the line break)
                                        pattern_array.append(line)
                                        
                                        # preparing new dict entry
                                        pattern_execution_dict = dict()
                                        pattern_execution_dict["seed"] = seed
                                        pattern_execution_helper = line_split[2].split(" ")
                                        pattern_execution_helper[-1] = pattern_execution_helper[-1][:-1]
                                        if pattern_file_marker == "V":
                                            pattern_execution_dict["pattern"] = [int(vertex) for vertex in pattern_execution_helper]
                                        else:
                                            pattern_execution_dict["pattern"] = pattern_execution_helper

                                        if( str(line_split[1]) in pattern_dict ):
                                            if( str(line_split[0]) in pattern_dict[str(line_split[1])] ):
                                                #adding new entry to the array
                                                pattern_dict[str(line_split[1])][str(line_split[0])].append(pattern_execution_dict)
                                            else:
                                                #adding new entry to the created array
                                                pattern_dict[str(line_split[1])][str(line_split[0])] = []
                                                pattern_dict[str(line_split[1])][str(line_split[0])].append(pattern_execution_dict)

                                        else:
                                            #adding new entry for first time
                                            pattern_dict[str(line_split[1])] = dict()
                                            pattern_dict[str(line_split[1])][str(line_split[0])] = []
                                            pattern_dict[str(line_split[1])][str(line_split[0])].append(pattern_execution_dict)
                                    line = pattern_file.readline()
                        except IOError:
                            pass
                    #saving the instance pattern dict
                    if len(pattern_dict) > 0:
                        pattern_dict_json_name = "output/" + str(timestamp) + "/" + folder + "/" + instance_name + "/pattern" + pattern_file_marker + ".json"
                        with open(pattern_dict_json_name, 'w') as f:
                            json.dump(pattern_dict, f, indent=4, sort_keys=True)
                        #Correcting the keys order
                        with open(pattern_dict_json_name, 'r') as f:
                            pattern_dict = json.load(f)
                        reversed_pattern_dict = OrderedDict()
                        for key, value in sorted(pattern_dict.items(), reverse=True):
                            reversed_pattern_dict[key] = {key_: value_ for key_, value_ in sorted(value.items(), reverse=True)}
                        with open(pattern_dict_json_name, 'w') as f:
                            json.dump(reversed_pattern_dict, f, indent=4)

                    # create a single pattern file
                    if (len(pattern_array) > 0):
                        keep_file = True
                    with open(model_file_name, 'a') as f:
                        f.write(str(len(pattern_array)) + "\n")
                        for entry in pattern_array:
                            f.write(entry)
                
                # delete the file if there was no pattern mining
                if not keep_file:
                    os.remove(model_file_name)


        #Validating group
        deviation_data = dict()
        elite_data = dict()
        for key in report:
            if not len(report[key]) > 0:
                break
            average = sum(seed_execution["best"] for seed_execution in report[key].values())/len(report[key])
            deviation_data[key] = ((average-opt[key])/opt[key])
            elite_data[key] = [len(seed_execution["elite"]) for seed_execution in report[key].values()]

        #Saving execution reports
        elapsed_folder_time = time.time() - start_folder_time
        report["#elapsed"] = time.strftime(time_format, time.localtime(elapsed_folder_time))
        deviation_dict[folder] = deviation_data
        elite_dict[folder] = elite_data
        report_dict[folder] = report
            
        #Marking iteration time
        elapsed_time = time.time() - start_time
        details_report["elapsed"] = time.strftime(time_format, time.localtime(elapsed_time))

        #Saving execution
        deviation_dict["#details"] = details_report
        elite_dict["#details"] = details_report
        report_dict["#details"] = details_report
        with open("output/" + str(timestamp) + "/deviation.json", 'w') as f:
            json.dump(deviation_dict, f, indent=4, sort_keys=True)
        with open("output/" + str(timestamp) + "/elite.json", 'w') as f:
            json.dump(elite_dict, f, indent=4, sort_keys=True)
        with open("output/" + str(timestamp) + "/results.json", 'w') as f:
            json.dump(report_dict, f, indent=4, sort_keys=True)

    # zipping the pattern files
    if data_mining == '1':
        model_path_name = "output/" + str(timestamp) + "/model"
        zip_file = zipfile.ZipFile(model_path_name +'.zip', 'w')
        for root, directories, files in os.walk(model_path_name):
            for filename in files:
                # get all files in the model folder.
                write_path = root.split(model_path_name)[1][1:] + "/" + filename
                zip_file.write(join(root, filename), write_path)
        os.system("rm -r output/" + str(timestamp) + "/model")

if __name__ == "__main__":
    main()
