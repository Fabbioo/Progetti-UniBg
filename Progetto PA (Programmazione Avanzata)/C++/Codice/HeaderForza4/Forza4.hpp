#ifndef FORZA4_HPP_
#define FORZA4_HPP_

#include "Griglia.hpp"

class Forza4 {
	Griglia griglia;
	int letturaColonna(int);
public:
	Forza4();
	~Forza4();
	int gioca();
};

#endif /* FORZA4_HPP_ */
