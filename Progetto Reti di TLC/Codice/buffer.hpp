#ifndef BUFFER_H
#define BUFFER_H

#include "packet.hpp"

class buffer {
	packet *head;
	packet *last;
public:
	buffer();
	~buffer() {}
	void insert(packet *pack);
	packet* get();
	packet* full() {
		return head;
	}
	double tot_delay;
	double tot_packs;
	int status;
	int K = 0;
	int G = 0;
	double P;
	int current;
	int tot_scarti_casuali;
	int tot_scarti_trabocco;
	double tot_arrivals;
	double tot_in_queue;
	void cancella(int);
};

#endif
