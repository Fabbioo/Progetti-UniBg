#include <iostream>
#include "global.hpp"
#include "queue.hpp"
#include "simulator.hpp"

using namespace std;

int main(int argc, char *argv[]) {

	simulator *eval;

	printf("**********************************************************\n\n");
	printf("               G/G/1 QUEUE SIMULATION PROGRAM\n\n");
	printf("**********************************************************\n\n");

	eval = new queue(argc, argv);

	eval->init();
	eval->run();
	eval->results();

	delete eval;
	eval = NULL;

	return 0;
}
