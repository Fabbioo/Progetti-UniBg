#ifndef GETTONEVERDE_HPP_
#define GETTONEVERDE_HPP_

#include "Gettone.hpp"

class GettoneVerde: public Gettone {
public:
	GettoneVerde();
	virtual ~GettoneVerde();
	void colore(char*);
	void stampa();
};

#endif /* GETTONEVERDE_HPP_ */
