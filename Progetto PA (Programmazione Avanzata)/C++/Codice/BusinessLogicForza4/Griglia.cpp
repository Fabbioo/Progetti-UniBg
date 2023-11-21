#include "../HeaderForza4/Gettone.hpp"
#include "../HeaderForza4/Griglia.hpp"
#include "../HeaderForza4/Generatore.hpp"

Griglia::Griglia() {
	for (int i = 0; i < RIGHE; i++)
		for (int j = 0; j < COLONNE; j++)
			griglia[i][j] = NULL;
	piene = 0;
}

Griglia::~Griglia() {
	for (int i = 0; i < RIGHE; i++)
		for (int j = 0; j < COLONNE; j++)
			if (griglia[i][j])
				delete griglia[i][j];
}

int Griglia::piena() {
	if (piene == RIGHE * COLONNE)
		return 1;
	return 0;
}

void Griglia::stampa() {
	cout << "        \t";
	for (int i = 0; i < COLONNE; i++)
		cout << "--";
	cout << endl;
	for (int i = RIGHE - 1; i >= 0; i--) {
		cout << "\tRiga: " << i + 1 << "| ";
		for (int j = 0; j < COLONNE; j++) {
			if (griglia[i][j] != NULL) {
				griglia[i][j]->stampa();
				cout << " ";
			} else
				cout << "  ";
		}
		cout << "|" << endl;
	}
	cout << "\tColonne: ";
	for (int j = 0; j < COLONNE; j++)
		cout << j + 1 << " ";
	cout << endl << endl;
}

int Griglia::cadutaGettone(int colonna) {
	int posizione = -1;
	for (int i = RIGHE - 1; i >= 0; i--)
		if (griglia[i][colonna - 1] == 0)
			posizione = i;
	return posizione;
}

int Griglia::esploraGriglia(char *colore) {
	for (int i = RIGHE - 1; i >= 0; i--)
		for (int j = 0; j < COLONNE; j++)
			if (cercaSequenza(i, j, colore) != 0)
				return 1;
	return 0;
}

int Griglia::cercaSequenza(int riga, int colonna, char *colore) {
	int r, c, drighe, dcolonne, conta;
	char coloreGettone[10];
	for (int i = 1; i <= 4; i++) {
		switch (i) {
		case 1:
			drighe = 0;
			dcolonne = 1;
			break;
		case 2:
			drighe = -1;
			dcolonne = 1;
			break;
		case 3:
			drighe = -1;
			dcolonne = 0;
			break;
		case 4:
			drighe = -1;
			dcolonne = -1;
			break;
		}
		r = riga;
		c = colonna;
		conta = 0;
		for (int i = 0; i < 4; i++) {
			if (r < 0 || r >= RIGHE)
				break;
			if (c < 0 || c >= COLONNE)
				break;
			if (griglia[r][c] != NULL) {
				griglia[r][c]->colore(coloreGettone);
				if (strcmp(coloreGettone, colore) == 0)
					conta++;
			}
			r += drighe;
			c += dcolonne;
		}
		if (conta == 4)
			return 1;
	}
	return 0;
}

int Griglia::inserisciGettone(char *colore, int colonna) {
	Gettone *p;
	int posizione;
	Generatore g;
	posizione = cadutaGettone(colonna);
	if (posizione < 0) {
		cout << "\tLa colonna scelta e' piena. Ritenta" << endl;
		return 1;
	}
	if (g.gettoneDaColore(colore, p))
		return -1;
	griglia[posizione][colonna - 1] = p;
	piene++;
	return 0;
}
