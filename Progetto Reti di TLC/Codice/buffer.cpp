#include "buffer.hpp"

buffer::buffer() {
	head = NULL;
	last = NULL;
	status = 0;
	tot_delay = 0.0;
	tot_packs = 0.0;
	tot_scarti_casuali = 0;
	tot_scarti_trabocco = 0;
	current = 0;
	tot_arrivals = 0;
	tot_in_queue = 0;
	P = 0;
}

void buffer::insert(packet *pack) {
	current++;
	if (head == NULL) {
		head = pack;
		last = pack;
		last->next = head;
	} else {
		last->next = pack;
		last = pack;
		last->next = head;
	}
}

packet* buffer::get() {
	packet *pack;
	if (head == NULL)
		return NULL;
	current--;
	if (last == head) {
		pack = head;
		last = NULL;
		head = NULL;
	} else {
		pack = head;
		head = head->next;
		last->next = head;
	}
	return pack;
}

void buffer::cancella(int X) {
	packet *e = head;
	for (int i = 1; i < X; i++)
		if (e->next != NULL)
			e = e->next;
	e->next = e->next->next;
	if (e->next == head)
		last = e;
	current--;
}
