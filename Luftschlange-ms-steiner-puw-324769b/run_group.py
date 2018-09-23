import os
from os.path import isfile, join
import sys
import datetime, time
import json
from pprint import pprint

def main():
    #Getting timestamp
    dts = datetime.datetime.utcnow()
    timestamp = round(time.mktime(dts.timetuple()) + dts.microsecond/1e6)
    os.system("mkdir output/" + str(timestamp))

    #Getting start time and iterations
    details_report = dict()
    time_format = '%Y-%m-%d %H:%M:%S'
    details_report["start"] = datetime.datetime.now().strftime(time_format)
    details_report["iterations"] = int(sys.argv[2])

    #Opening instances folder
    report_dict = dict()
    deviation_dict = dict()
    elite_dict = dict()
    for folder in os.listdir("instances/" + sys.argv[1]):
        path = sys.argv[1] + "/" + folder

        #Preparing output folder
        output_folder = str(timestamp) + "/" + folder
        os.system("mkdir output/" + output_folder)

        #Getting instances files
        onlyfiles = [f for f in os.listdir("instances/" + path) if isfile(join("instances/" + path, f))]

        #Executing instances
        report = dict()
        for instance in onlyfiles:
            instance_name = instance.split(".")[0]
            if ".json" in instance:
                continue
            os.system("./bin/run_instance " + path + "/" + instance + " -seed 1 -maxit " + sys.argv[2] + " -mine 0 -folder " + output_folder)
            with open("output/" + output_folder + "/" + instance_name + "/" + instance_name + ".json") as f:
                data = json.load(f)
                report[instance_name] = data
            os.system("rm output/" + output_folder + "/" + instance_name + "/" + instance_name + ".json")

        #Validating group
        deviation_data = dict()
        elite_data = dict()
        with open("instances/" + path + "/opt.json") as fopt:
                opt = json.load(fopt)
                for key in report:
                    average = sum(seed_execution["best"] for seed_execution in report[key].values())/len(report[key])
                    deviation_data[key] = ((average-opt[key])/opt[key])
                    elite_data[key] = [len(seed_execution["elite"]) for seed_execution in report[key].values()]

        #Saving execution reports
        deviation_dict[folder] = deviation_data
        elite_dict[folder] = elite_data
        report_dict[folder] = report
            
    #Ending execution
    details_report["stop"] = datetime.datetime.now().strftime(time_format)
    deviation_dict["#details"] = details_report
    elite_dict["#details"] = details_report
    report_dict["#details"] = details_report
    with open("output/" + str(timestamp) + "/deviation.json", 'w') as f:
        json.dump(deviation_dict, f, indent=4, sort_keys=True)
    with open("output/" + str(timestamp) + "/elite.json", 'w') as f:
        json.dump(elite_dict, f, indent=4, sort_keys=True)
    with open("output/" + str(timestamp) + "/results.json", 'w') as f:
        json.dump(report_dict, f, indent=4, sort_keys=True)

if __name__ == "__main__":
    main()