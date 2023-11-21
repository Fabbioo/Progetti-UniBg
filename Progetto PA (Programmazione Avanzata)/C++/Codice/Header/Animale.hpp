#ifndef HEADER_ANIMALE_HPP_
#define HEADER_ANIMALE_HPP_

#include "EssereVivente.hpp"

class Animale: public EssereVivente {
	bool pedigree;
public:
	Animale(char*, int, char, bool);
	~Animale();
	void stampa();
	void descrizione(char*);
};

#endif /* HEADER_ANIMALE_HPP_ */
