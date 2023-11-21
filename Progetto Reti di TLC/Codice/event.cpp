#include "event.hpp"
#include "buffer.hpp"
#include "calendar.hpp"
#include "rand.hpp"
#include <iostream>

using namespace std;

extern calendar *cal;
extern double inter;
extern double duration;

void arrival::body() {
	buf->tot_arrivals++;
	event *ev;
	double esito;
	GEN_EXP(SEED, inter, esito);
	ev = new arrival(time + esito, buf);
	cal->put(ev);
	packet *pack = new packet(time);
	if (buf->current == buf->K) {
		delete pack;
		buf->tot_scarti_trabocco++;
	} else if (buf->full() || buf->status) {
		buf->insert(pack);
		buf->tot_in_queue++;
	} else {
		buf->tot_in_queue++;
		buf->tot_packs += 1.0;
		delete pack;
		GEN_RAY(SEED, duration, esito);
		ev = new service(time + esito, buf);
		cal->put(ev);
		buf->status = 1;
	}
}

void service::body() {
	if (buf->current > buf->G) {
		int p = rand() % 10 + 1;
		if (p <= (buf->P * 10)) {
			buf->tot_scarti_casuali++;
			int X = rand() % (buf->current);
			if (X == 0)
				buf->get();
			else
				buf->cancella(X);
		}
	}
	packet *pack;
	pack = buf->get();
	event *ev;
	double esito;
	GEN_RAY(SEED, duration, esito);
	if (pack != NULL) {
		ev = new service(time + esito, buf);
		cal->put(ev);
		buf->tot_delay += time - pack->get_time();
		buf->tot_packs += 1.0;
		delete pack;
	} else
		buf->status = 0;
}
