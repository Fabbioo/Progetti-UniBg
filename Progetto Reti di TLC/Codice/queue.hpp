#ifndef _QUEUE_H
#define _QUEUE_H

#include "buffer.hpp"
#include "calendar.hpp"
#include "event.hpp"
#include "packet.hpp"
#include "simulator.hpp"
#include "stat.hpp"

class queue: public simulator {
	virtual void input(void);
	buffer *buf;
	int traffic_model;
	double load;
	int service_model;
	double packets;
	double tot_delay;
	Sstat *delay;
	virtual void clear_counters(void);
	virtual void clear_stats(void);
	virtual void update_stats(void);
	virtual void print_trace(int Run);
	virtual void results(void);
	virtual int isconfsatisf(double perc);
	Sstat *Sstat_Pcasuali;
	Sstat *Sstat_Ptrabocco;
	double prob_scarti_casuali;
	double prob_scarti_trabocco;
public:
	queue(int argc, char *argv[]);
	virtual ~queue(void);
	virtual void init(void);
	virtual void run(void);
};

#endif
