#ifndef HEADER_UMANO_HPP_
#define HEADER_UMANO_HPP_

#include "EssereVivente.hpp"
#include "Cane.hpp"

#include <vector>

class Umano: public EssereVivente {
	char *cognome;
	string codiceFiscale;
public:
	Umano(char*, char*, int, char, string);
	~Umano();
	void stampa();
	void descrizione(char*);
	void iniziaGara(vector<Cane*>);
};

#endif /* HEADER_UMANO_HPP_ */
