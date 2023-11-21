#include "../Header/Umano.hpp"

#include <map>
#include <utility>

Umano::Umano(char *nome, char *cognome, int anni, char sesso, string codiceFiscale) :
		EssereVivente(nome, anni, sesso) {
	this->cognome = new char[strlen(cognome) + 1];
	strcpy(this->cognome, cognome);
	this->codiceFiscale = codiceFiscale;
}

Umano::~Umano() {
	delete[] cognome;
}

void Umano::stampa() {
	EssereVivente::stampa();
	cout << "\tCognome: " << cognome << endl;
	cout << "\tCodice fiscale: " << codiceFiscale << endl;
}

void Umano::descrizione(char *nome) {
	cout << "\t" << nome << " è un oggetto della classe Umano" << endl;
}

void Umano::iniziaGara(vector<Cane*> caniPartecipantiGara) {
	srand(time(0));
	map<int, string> classificaPunteggi;
	for (int i = 0; i < caniPartecipantiGara.size(); i++)
		classificaPunteggi.insert(pair<int, string>(rand() % 101, caniPartecipantiGara.at(i)->getNome()));
	vector<pair<int, string> > arrayAppoggio;
	for(auto item: classificaPunteggi)
		arrayAppoggio.push_back(item);
	cout << "\tClassifica generale:" << endl;
	for(int i = 0; i < arrayAppoggio.size(); i++)
		cout << "\t" << (i + 1) << "^ candidato: " << arrayAppoggio.at(i).second << "\tPunteggio: " << arrayAppoggio.at(i).first << endl;
	cout << endl << "\tPodio:" << endl;
	for(int i = 0; i < 3; i++)
		cout << "\t" << (i + 1) << "^ classificato: " << arrayAppoggio.at(4 - i).second << "\tPunteggio: " << arrayAppoggio.at(4 - i).first << endl;
}
