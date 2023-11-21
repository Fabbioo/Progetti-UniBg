#include "../Header/Animale.hpp"

#define valBool(x) ((x) ? "True" : "False")

Animale::Animale(char *nome, int anni, char sesso, bool pedigree) :
		EssereVivente(nome, anni, sesso) {
	this->pedigree = pedigree;
}

Animale::~Animale() {}

void Animale::stampa() {
	EssereVivente::stampa();
	cout << "\tPedigree: " << valBool(pedigree) << endl;
}

void Animale::descrizione(char *nome) {
	cout << "\t" << nome << " è un oggetto della classe Animale" << endl;
}
