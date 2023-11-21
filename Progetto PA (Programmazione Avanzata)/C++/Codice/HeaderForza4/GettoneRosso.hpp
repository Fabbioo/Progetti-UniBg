#ifndef GETTONEROSSO_HPP_
#define GETTONEROSSO_HPP_

#include "Gettone.hpp"

class GettoneRosso: public Gettone {
public:
	GettoneRosso();
	virtual ~GettoneRosso();
	void colore(char*);
	void stampa();
};

#endif /* GETTONEROSSO_HPP_ */
