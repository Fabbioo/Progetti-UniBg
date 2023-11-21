#ifndef GETTONE_HPP_
#define GETTONE_HPP_

#define RIGHE 8
#define COLONNE 10

#include <iostream>

using namespace std;

class Gettone {
public:
	Gettone();
	virtual ~Gettone();
	virtual void colore(char*);
	virtual void stampa();
};

#endif /* GETTONE_HPP_ */
