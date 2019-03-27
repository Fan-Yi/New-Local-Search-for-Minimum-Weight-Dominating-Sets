all: myMinWDS

myMinWDS: myMinWDS.cpp graph.h localSearch.h constants.h myBijection.h hugeInt.h config.h weightBuckets.h operandSets.h scenarioHash.h
	g++ -std=gnu++0x -g -O3 -static myMinWDS.cpp -o myMinWDS

clean: rm -f *~
