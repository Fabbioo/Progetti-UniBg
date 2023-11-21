#include "../Header/EssereVivente.hpp"
#include "../Header/OperazioniVettore.hpp"

#include <algorithm>

Print::Print(vector<EssereVivente*> vettore) {
	this->vettore = vettore;
}

Print::~Print() {}

void Print::printing(vector<EssereVivente*> vettore) {
	for (int i = 0; i < vettore.size(); i++)
		cout << "\tNome: " << vettore.at(i)->getNome() << " -> età: " << vettore.at(i)->getAnni() << endl;
}

Sort::Sort(vector<EssereVivente*> vettore) {
	this->vettore = vettore;
}

Sort::~Sort() {}

void Sort::sorting(vector<EssereVivente*> &vettore) {
	sort(vettore.begin(), vettore.end(), EssereVivente::compareByAnni);
}

VettoreFacade::VettoreFacade(vector<EssereVivente*> vettore) {
	this->vettore = vettore;
	this->sort = new Sort(vettore);
	this->print = new Print(vettore);
}

VettoreFacade::~VettoreFacade() {
	delete sort;
	delete print;
}

void VettoreFacade::startSortPrint() {
	sort->sorting(vettore);
	print->printing(vettore);
}
