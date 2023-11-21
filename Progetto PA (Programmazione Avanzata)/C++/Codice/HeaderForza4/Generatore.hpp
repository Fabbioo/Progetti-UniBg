#ifndef GENERATORE_HPP_
#define GENERATORE_HPP_

#include "Gettone.hpp"

class Generatore {
public:
	Generatore();
	virtual ~Generatore();
	int gettoneDaColore(char*, Gettone*&);
	int gettoneDaGiocatore(int, Gettone*&);
	int coloreDaGiocatore(int, char*);
};

#endif /* GENERATORE_HPP_ */
