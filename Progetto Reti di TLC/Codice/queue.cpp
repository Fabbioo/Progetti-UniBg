#include "queue.hpp"
#include <stdio.h>
#include "buffer.hpp"
#include "calendar.hpp"
#include "easyio.hpp"
#include "event.hpp"
#include "global.hpp"
#include "rand.hpp"

using namespace std;

calendar *cal;
double inter;
double duration;
double Trslen;
double Runlen;
int NRUNmin;
int NRUNmax;

queue::queue(int argc, char *argv[]) :
		simulator(argc, argv) {
	cal = new calendar();
	buf = new buffer();
	delay = new Sstat();
	Sstat_Pcasuali = new Sstat();
	Sstat_Ptrabocco = new Sstat();
}

queue::~queue() {
	delete delay;
	delete cal;
	delete buf;
}

void queue::input() {
	printf("MODEL PARAMETERS:\n");
	printf("\nArrivals model: ");
	printf("Poisson>\n");

	traffic_model = read_int((char*) "", 1, 1, 1);
	load = read_double((char*) "Traffic load (Erlang)", 0.9, 0.01, 0.999);

	printf("\nService model: ");
	printf("Rayleigh>\n");

	service_model = read_int((char*) "", 1, 1, 1);
	duration = read_double((char*) "\nAverage service duration (s)", 0.4, 0.01,
			100);
	inter = duration / load;
	printf("\nSIMULATION PARAMETERS:\n\n");

	Trslen = read_double((char*) "Simulation transient len (s)", 1000, 0.01,
			10000);
	Runlen = read_double((char*) "\nSimulation RUN len (s)", 1000, 0.01, 10000);
	NRUNmin = read_int((char*) "\nSimulation number of RUNs", 50, 2, 100);

	buf->K = read_int((char*) "\nBuffer Capacity K", 6, 1, 10);
	buf->G = read_int((char*) "\nBuffer Threshold G", 4, 1, buf->K - 1);
	buf->P = read_double((char*) "\nProbabilita' di scarto casuale P", 0.8, 0.1,
			1);
}

void queue::init() {
	input();
	event *Ev;
	Ev = new arrival(0.0, buf);
	cal->put(Ev);
	buf->status = 0;
}

void queue::run() {
	extern double Trslen;
	extern double Runlen;
	extern int NRUNmin;
	double clock = 0.0;
	event *ev;
	while (clock < Trslen) {
		ev = cal->get();
		ev->body();
		clock = ev->time;
		delete (ev);
	}
	clear_stats();
	clear_counters();
	int current_run_number = 1;
	while (current_run_number <= NRUNmin) {
		while (clock < (current_run_number * Runlen + Trslen)) {
			ev = cal->get();
			ev->body();
			clock = ev->time;
			delete (ev);
		}
		update_stats();
		print_trace(current_run_number);
		clear_counters();
		current_run_number++;
	}
}

void queue::results() {
	extern double Trslen;
	extern double Runlen;
	extern int NRUNmin;
	fprintf(fpout, "\n*********************************************\n");
	fprintf(fpout, "           SIMULATION RESULTS                \n");
	fprintf(fpout, "*********************************************\n\n");
	fprintf(fpout, "Input parameters:\n");
	fprintf(fpout, "Transient length (s)         %5.3f\n", Trslen);
	fprintf(fpout, "Run length (s)               %5.3f\n", Runlen);
	fprintf(fpout, "Number of runs               %5d\n", NRUNmin);
	fprintf(fpout, "Traffic load                 %5.3f\n", load);
	fprintf(fpout, "Average service duration     %5.3f\n", duration);
	fprintf(fpout, "Results:\n");
	fprintf(fpout, "Average Delay                %2.6f   +/- %.2e  p:%3.2f\n",
			delay->mean(), delay->confidence(.95), delay->confpercerr(.95));
	fprintf(fpout, "D  %2.6f   %2.6f   %.2e %2.6f\n", load, delay->mean(),
			delay->confidence(.95), duration * (load) / (1 - load));
	fprintf(fpout, "Probabilita' media di scarto casuale: %f %% +/- %.2e\n",
			Sstat_Pcasuali->mean(), Sstat_Pcasuali->confidence(.95));
	fprintf(fpout,
			"Probabilita' media di scarto per trabocco: %f %% +/- %.2e\n",
			Sstat_Ptrabocco->mean(), Sstat_Ptrabocco->confidence(.95));
}

void queue::print_trace(int n) {
	fprintf(fptrc, "\n*********************************************\n");
	fprintf(fptrc, "                 TRACE RUN %d                \n", n);
	fprintf(fptrc, "*********************************************\n\n");
	fprintf(fptrc, "Average Delay                %2.6f   +/- %.2e  p:%3.2f\n",
			delay->mean(), delay->confidence(.95), delay->confpercerr(.95));
	fprintf(fptrc, "Tot Arrivals: 				%f\n", buf->tot_arrivals);
	fprintf(fptrc, "Tot in Queue: 				%f\n", buf->tot_in_queue);
	fprintf(fptrc, "Prob Scarti Casuali: 		%f %%\n", prob_scarti_casuali);
	fprintf(fptrc, "Prob Scarti Trabocco: 		%f %%\n", prob_scarti_trabocco);
	fflush(fptrc);
}

void queue::clear_counters() {
	buf->tot_delay = 0.0;
	buf->tot_packs = 0.0;
	buf->current = 0;
	buf->tot_arrivals = 0;
	buf->tot_in_queue = 0;
	buf->tot_scarti_casuali = 0;
	buf->tot_scarti_trabocco = 0;
}

void queue::clear_stats() {
	delay->reset();
	Sstat_Pcasuali->reset();
	Sstat_Ptrabocco->reset();
}
void queue::update_stats() {
	*delay += buf->tot_delay / buf->tot_packs;

	// Calcolo delle probabilita' per la RUN attuale
	prob_scarti_casuali = (buf->tot_scarti_casuali / buf->tot_in_queue) * 100;
	prob_scarti_trabocco = (buf->tot_scarti_trabocco / buf->tot_arrivals) * 100;

	*Sstat_Pcasuali += prob_scarti_casuali;
	*Sstat_Ptrabocco += prob_scarti_trabocco;
}

int queue::isconfsatisf(double perc) {
	return delay->isconfsatisfied(10, .95);
}
