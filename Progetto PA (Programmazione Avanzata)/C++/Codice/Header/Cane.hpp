#ifndef HEADER_CANE_HPP_
#define HEADER_CANE_HPP_

#include "Animale.hpp"

class Cane: public Animale {
	int oreGioco;
public:
	Cane(char*, int, char, bool, int);
	~Cane();
	void stampa();
	void descrizione(char*);
};

#endif /* HEADER_CANE_HPP_ */
