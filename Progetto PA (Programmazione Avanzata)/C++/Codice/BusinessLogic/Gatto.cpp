#include "../Header/Gatto.hpp"

Gatto::Gatto(char *nome, int anni, char sesso, bool pedigree, string dedica) :
		Animale(nome, anni, sesso, pedigree) {
	this->dedica = dedica;
}

Gatto::~Gatto() {}

void Gatto::stampa() {
	Animale::stampa();
	cout << "\tDedica: " << dedica << endl;
}

void Gatto::descrizione(char *nome) {
	cout << "\t" << nome << " è un oggetto della classe Gatto" << endl;
}

void Gatto::giocaCon(unique_ptr<Gatto> gatto) {
	srand(time(0));
	int c = (char)(rand() % 256);
	string s = ((c > 0) ? "" : " non");
	cout << "\t" << this->getNome() << " e " << gatto->getNome() << s << " si stanno divertendo!" << endl;
}
