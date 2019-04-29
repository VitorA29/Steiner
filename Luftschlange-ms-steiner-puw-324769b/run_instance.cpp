#include <sstream>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using namespace std;

static void ShowUsage () {
    fprintf (stderr, "Usage: run <stp_file> [-seed] [-mine] [-maxit] [-folder]\n");
    exit (-1);
}

const vector<string> explode(const string& s, const char& c)
{
	string buff{""};
	vector<string> v;
	
	for(auto n:s)
	{
		if(n != c) buff+=n; else
		if(n == c && buff != "") { v.push_back(buff); buff = ""; }
	}
	if(buff != "") v.push_back(buff);
	
	return v;
}

int main (int argc, char **argv) {
    if(argc < 2) ShowUsage();

    int max_seed = 10;
    int iterations = 50;
    bool MINING = false;
    bool DEFAULTFOLDER = true;
    stringstream ss;

    for (int i=2; i<argc; i+=2) {
        if (i == argc-1) ShowUsage();
        if (strcmp(argv[i], "-seed")==0) {
            max_seed = atoi(argv[i+1]);
            continue;
        }
        if (strcmp(argv[i], "-mine")==0) {
            if (atoi(argv[i+1])!=0) {
                MINING = true;
            }
            continue;
        }
        if (strcmp(argv[i], "-maxit")==0) {
            iterations = atoi(argv[i+1]);
            continue;
        }
        if (strcmp(argv[i], "-folder")==0) {
            DEFAULTFOLDER = false;
            ss << argv[i+1];
            continue;
        }
    }

	ostringstream buffer;
    string instance_path = argv[1];

    //getting instance name
    vector<string> v{explode(instance_path, '/')};
    string instance = v[v.size()-1];
    vector<string> n{explode(instance, '.')};
    string instance_name = n[0];

    //starting output folder
    if(DEFAULTFOLDER){
        time_t time_marker = time(0);
        ss << time_marker;
        buffer.str("");
        buffer << "mkdir output/" << ss.str().c_str();
        system(buffer.str().c_str());
    }
    ss << "/";
    ss << instance_name;
    buffer.str("");
    buffer << "mkdir output/" << ss.str().c_str();
    system(buffer.str().c_str());
    string output_folder = ss.str();
    string output_folder_path = "output/";
    output_folder_path += output_folder;

    //starting execution
    printf("Executing %s instance\n", instance_name.c_str());
    string mining = "";
    if(MINING){
        mining = "-mine 1 ";
    }
    int elite_cap = 10;
    string name_json = output_folder_path;
    name_json += "/";
    name_json += instance_name;
    name_json += ".json";
    FILE *json = fopen(name_json.c_str(), "w");
    fprintf(json, "{\n");

    //preparing the instance, making the preprocess
    buffer.str("");
    buffer << "./bin/steiner instances/"<< instance_path << " -preprocess " << output_folder_path << "/" << instance_name << "_preproceded.stp";
    system(buffer.str().c_str());

    for(int seed=1;seed<=max_seed;seed++){
        buffer.str("");
        buffer << "./bin/steiner instances/"<< instance_path << " -seed " << seed << " -msit " << iterations << " -elite "<< elite_cap << " " << mining <<"-folder " << output_folder;
        // buffer << "./bin/steiner " << output_folder_path << "/" << instance_name << "_preproceded.stp " << " -seed " << seed << " -msit " << iterations << " -elite "<< elite_cap << " " << mining <<"-folder " << output_folder;
        system(buffer.str().c_str());

        //preparing file name best
        stringstream ss;
        ss << output_folder_path << "/" << seed << "/best.bin";
        string file_name = ss.str();
        //opening best
        FILE *fp = fopen(file_name.c_str(), "r");
        float best_cost;
        fscanf(fp, "%f", &best_cost);
        fclose(fp);
        fprintf(json, "\t\"%d\":{\n\t\t\"best\":%.0f,\n\t\t\"elite\":[\n", seed, best_cost);
        //deleting best
        buffer.str("");
        buffer << "rm " << file_name;
        system(buffer.str().c_str());

        //preparing file name elite
        ss.str("");
        ss << output_folder_path << "/" << seed << "/elite.bin";
        file_name = ss.str();
        //opening elite
        fp = fopen(file_name.c_str(), "r");
        vector<double> elite_cost;
        int read_elite_cost = fscanf(fp, "%f", &best_cost);
        while(read_elite_cost==1){
            elite_cost.push_back(best_cost);
            read_elite_cost = fscanf(fp, "%f", &best_cost);
        }
        fclose(fp);
        //deleting elite
        buffer.str("");
        buffer << "rm " << file_name;
        system(buffer.str().c_str());

        //adding elite to json
        for(int i = 0;i<elite_cost.size();i++){
            fprintf(json, "\t\t\t%.0f", elite_cost[i]);
            if(i<elite_cost.size()-1)
                fprintf(json, ",");
            fprintf(json, "\n");
        }
        fprintf(json, "\t\t\]\n\t}");
        if(seed<max_seed)
            fprintf(json, ",");
        fprintf(json, "\n");
        fflush(json);
    }
    fprintf(json, "}");
    fclose(json);
	return 0;
}