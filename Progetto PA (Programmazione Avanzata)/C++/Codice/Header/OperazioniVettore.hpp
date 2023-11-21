#ifndef OPERAZIONIVETTORE_HPP_
#define OPERAZIONIVETTORE_HPP_

#include "EssereVivente.hpp"

#include <vector>

class Print {
	vector<EssereVivente*> vettore;
public:
	Print(vector<EssereVivente*>);
	~Print();
	void printing(vector<EssereVivente*>);
};

class Sort {
	vector<EssereVivente*> vettore;
public:
	Sort(vector<EssereVivente*>);
	~Sort();
	void sorting(vector<EssereVivente*>&);
};

class VettoreFacade {
	vector<EssereVivente*> vettore;
	Sort *sort;
	Print *print;
public:
	VettoreFacade(vector<EssereVivente*>);
	~VettoreFacade();
	void startSortPrint();
};

#endif /* OPERAZIONIVETTORE_HPP_ */
