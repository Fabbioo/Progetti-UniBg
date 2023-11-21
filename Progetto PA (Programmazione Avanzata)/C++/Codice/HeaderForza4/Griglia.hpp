#ifndef GRIGLIA_HPP_
#define GRIGLIA_HPP_

#include "Gettone.hpp"

class Griglia {
	Gettone *griglia[RIGHE][COLONNE];
	int piene;
	int cadutaGettone(int);
	int cercaSequenza(int, int, char*);
public:
	Griglia();
	virtual ~Griglia();
	void stampa();
	int esploraGriglia(char*);
	int inserisciGettone(char*, int);
	int piena();
};

#endif /* GRIGLIA_HPP_ */
