#include "../Header/Cane.hpp"

Cane::Cane(char *nome, int anni, char sesso, bool pedigree, int oreGioco) :
		Animale(nome, anni, sesso, pedigree) {
	this->oreGioco = oreGioco;
}

Cane::~Cane() {}

void Cane::stampa() {
	Animale::stampa();
	cout << "\toreGioco: " << oreGioco << endl;
}

void Cane::descrizione(char *nome) {
	cout << "\t" << nome << " è un oggetto della classe Cane" << endl;
}
