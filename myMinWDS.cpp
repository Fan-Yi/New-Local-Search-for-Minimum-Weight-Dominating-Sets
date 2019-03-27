#include <iostream>
#include "localSearch.h"

using namespace std;

int main(int argc, char* argv[])
{
	char solver_name[1024];
	sscanf(argv[0], "%s", solver_name);

	char filename[1024];
	sscanf(argv[1], "%s", filename);

	int seed;
	sscanf(argv[2], "%d", &seed);
	srand(seed);

	LL time_limit;
	sscanf(argv[3], "%lld", &time_limit);

	StochasticLocalSearch mySLS(solver_name, filename, seed, time_limit);

	mySLS.cover_sls();

	mySLS.show_results();

	return 0;
}
