#include <sstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int main () {
	ostringstream buffer;
    int max_seed = 10;
    printf("Executing small instance\n");
    for(int seed=1;seed<=max_seed;seed++){
        buffer.str("");
        buffer << "./bin/steiner instances/geo_original/G207.stp -seed " << seed << " -msit 50 -elite 10 -mine 8";
        system(buffer.str().c_str());
    }
    printf("Executing medium instance\n");
    for(int seed=1;seed<=max_seed;seed++){
        buffer.str("");
        buffer << "./bin/steiner instances/geo_original/G102.stp -seed " << seed << " -msit 50 -elite 10 -mine 9";
        system(buffer.str().c_str());
    }
    printf("Executing big instance\n");
    for(int seed=1;seed<=max_seed;seed++){
        buffer.str("");
        buffer << "./bin/steiner instances/geo_original/G303.stp -seed " << seed << " -msit 50 -elite 10 -mine 9";
        system(buffer.str().c_str());
    }
	return 0;
}