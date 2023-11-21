#include "../HeaderForza4/Forza4.hpp"
#include "../HeaderForza4/Generatore.hpp"

Forza4::Forza4() {
}

Forza4::~Forza4() {
}

int Forza4::letturaColonna(int giocatore) {
	int colonnaScelta;
	Generatore g;
	char colore[10];
	if (g.coloreDaGiocatore(giocatore, colore))
		return 1;
	do {
		cout << "\tGiocatore " << colore << ". Scegli colonna: ";
		cin >> colonnaScelta;
		if (colonnaScelta <= 0 || colonnaScelta > COLONNE) {
			cout << "\tHai scelto una colonna sbagliata. Riprova" << endl;
			griglia.stampa();
		}
	} while (colonnaScelta <= 0 || colonnaScelta > COLONNE);
	return colonnaScelta;
}

int Forza4::gioca() {
	int colonnaScelta, giocatore = 1, fine = 0, risp;
	Generatore g;
	char colore[10];
	do {
		griglia.stampa();
		colonnaScelta = letturaColonna(giocatore);
		if (g.coloreDaGiocatore(giocatore, colore))
			return 1;
		risp = griglia.inserisciGettone(colore, colonnaScelta);
		if (risp == 1)
			continue;
		if (risp == -1)
			return 1;
		if (griglia.esploraGriglia(colore)) {
			griglia.stampa();
			cout << "\tHa vinto il giocatore " << colore << "!!!" << endl;
			fine = 1;
		}
		if (griglia.piena() && !fine) {
			cout << "\tStallo: la griglia e' piena." << endl;
			fine = 1;
		}
		if (fine)
			continue;
		if (giocatore == 1)
			giocatore = 2;
		else
			giocatore = 1;
	} while (!fine);
	return 0;
}
