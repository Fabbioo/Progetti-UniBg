#include "../HeaderForza4/GettoneVerde.hpp"

GettoneVerde::GettoneVerde() {
}

GettoneVerde::~GettoneVerde() {
}

void GettoneVerde::colore(char *s) {
	strcpy(s, "Verde");
}

void GettoneVerde::stampa() {
	cout << "V";
}
