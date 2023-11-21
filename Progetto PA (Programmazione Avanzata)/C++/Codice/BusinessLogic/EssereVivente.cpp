#include "../Header/EssereVivente.hpp"

EssereVivente::EssereVivente(char *nome, int anni, char sesso) {
	this->nome = new char[strlen(nome) + 1];
	strcpy(this->nome, nome);
	this->anni = anni;
	this->sesso = sesso;
}

EssereVivente::~EssereVivente() {
	delete[] nome;
}

char* EssereVivente::getNome() {
	return nome;
}

int EssereVivente::getAnni() {
	return anni;
}

void EssereVivente::stampa() {
	cout << "\tNome: " << nome << endl;
	cout << "\tEtà: " << anni << " anni" << endl;
	cout << "\tSesso: " << sesso << endl;
}

void EssereVivente::descrizione(char *nome) {
	cout << "\t" << nome << " è un oggetto della classe EssereVivente" << endl;
}

bool EssereVivente::compareByAnni(EssereVivente *a, EssereVivente *b) {
	return a->getAnni() < b->getAnni();
}
