#ifndef HEADER_GATTO_HPP_
#define HEADER_GATTO_HPP_

#include "Animale.hpp"

#include <memory>

class Gatto: public Animale {
	string dedica;
public:
	Gatto(char*, int, char, bool, string);
	~Gatto();
	void stampa();
	void descrizione(char*);
	void giocaCon(unique_ptr<Gatto>);
};

#endif /* HEADER_GATTO_HPP_ */
