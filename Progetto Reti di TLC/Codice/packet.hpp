#ifndef PACKET_H
#define PACKET_H

#include "global.hpp"

class packet {
	double gen_time;
public:
	packet(double Gen_time);
	~packet() {}
	packet *next;
	double get_time() {
		return gen_time;
	}
};

inline packet::packet(double Gen_time) {
	gen_time = Gen_time;
	next = NULL;
}

#endif
