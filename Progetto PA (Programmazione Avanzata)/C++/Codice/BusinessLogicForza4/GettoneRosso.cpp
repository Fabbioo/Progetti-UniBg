#include "../HeaderForza4/GettoneRosso.hpp"

GettoneRosso::GettoneRosso() {
}

GettoneRosso::~GettoneRosso() {
}

void GettoneRosso::colore(char *s) {
	strcpy(s, "Rosso");
}

void GettoneRosso::stampa() {
	cout << "R";
}
