#include "../HeaderForza4/Gettone.hpp"
#include "../HeaderForza4/Generatore.hpp"
#include "../HeaderForza4/GettoneRosso.hpp"
#include "../HeaderForza4/GettoneVerde.hpp"

#include <string.h>

Generatore::Generatore() {
}

Generatore::~Generatore() {
}

int Generatore::gettoneDaColore(char *colore, Gettone *&gettone) {
	Gettone *p;
	if (strcmp(colore, "Rosso") == 0) {
		p = new GettoneRosso();
		gettone = p;
		return 0;
	}
	if (strcmp(colore, "Verde") == 0) {
		p = new GettoneVerde();
		gettone = p;
		return 0;
	}
	cout << "\tColore sconosciuto" << endl;
	return 1;
}

int Generatore::gettoneDaGiocatore(int giocatore, Gettone *&gettone) {
	Gettone *p;
	if (giocatore == 1) {
		p = new GettoneRosso();
		gettone = p;
		return 0;
	}
	if (giocatore == 2) {
		p = new GettoneVerde();
		gettone = p;
		return 0;
	}
	cout << "\tColore sconosciuto" << endl;
	return 1;
}

int Generatore::coloreDaGiocatore(int giocatore, char *colore) {
	if (giocatore == 1) {
		strcpy(colore, "Rosso");
		return 0;
	}
	if (giocatore == 2) {
		strcpy(colore, "Verde");
		return 0;
	}
	cout << "\tColore sconosciuto" << endl;
	return 1;
}
