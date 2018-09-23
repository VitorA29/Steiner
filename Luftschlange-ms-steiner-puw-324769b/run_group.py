import os
from os.path import isfile, join
import sys
import datetime, time
import json
from pprint import pprint

def main():
    dts = datetime.datetime.utcnow()
    timestamp = round(time.mktime(dts.timetuple()) + dts.microsecond/1e6)
    os.system("mkdir output/" + str(timestamp))

    final_dictionary = dict()
    for folder in os.listdir("instances/luidi"):
        path = "luidi/" + folder

        output_folder = str(timestamp) + "/" + folder
        os.system("mkdir output/" + output_folder)

        onlyfiles = [f for f in os.listdir("instances/" + path) if isfile(join("instances/" + path, f))]

        execution_dictionary = dict()

        for instance in onlyfiles:
            instance_name = instance.split(".")[0]
            if ".json" in instance:
                continue
            os.system("./bin/run_instance " + path + "/" + instance + " -seed 10 -maxit " + sys.argv[1] + " -mine 0 -folder " + output_folder)
            with open("output/" + output_folder + "/" + instance_name + "/" + instance_name + ".json") as f:
                data = json.load(f)
                execution_dictionary[instance_name] = data
            os.system("rm output/" + output_folder + "/" + instance_name + "/" + instance_name + ".json")

        with open("output/" + output_folder + "/execution.json", 'w') as f:
            json.dump(execution_dictionary, f, indent=4, sort_keys=True)

        data = dict()
        with open("output/" + output_folder + "/execution.json") as f:
            report = json.load(f)
            with open("instances/" + path + "/opt.json") as fopt:
                opt = json.load(fopt)
                for key in report:
                    counter = 0
                    average = 0
                    instance = report[key]
                    for seed in instance:
                        average += instance[seed]["best"]
                        counter += 1
                    average /= counter
                    data[key] = ((average-opt[key])/opt[key])

        with open("output/" + output_folder + "/results.json", 'w') as f:
            json.dump(data, f, indent=4, sort_keys=True)
            final_dictionary[folder] = data
        os.system("rm output/" + output_folder + "/results.json")
    with open("output/" + str(timestamp) + "/results.json", 'w') as f:
        json.dump(final_dictionary, f, indent=4, sort_keys=True)

if __name__ == "__main__":
    main()