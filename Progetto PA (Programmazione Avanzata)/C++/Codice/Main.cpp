#include "Header/EssereVivente.hpp"
#include "Header/Umano.hpp"
#include "Header/Animale.hpp"
#include "Header/Cane.hpp"
#include "Header/Gatto.hpp"
#include "Header/OperazioniVettore.hpp"
#include "HeaderForza4/Forza4.hpp"

int main() {

	vector<EssereVivente*> vettoreUmani;

	Umano *u1 = new Umano((char*) "Fabio", (char*) "Mandalari", 25, 'M', (char*) "Codice1Fiscale2");
	vettoreUmani.push_back(u1);
	Umano *u2 = new Umano((char*) "Marcelo", (char*) "Rossi", 5, 'M', (char*) "Codice2Fiscale3");
	vettoreUmani.push_back(u2);
	Umano *u3 = new Umano((char*) "Antonia", (char*) "Verdi", 20, 'F', (char*) "Codice3Fiscale4");
	vettoreUmani.push_back(u3);

	vector<Cane*> caniPartecipantiGara;

	Cane c1 = Cane((char*) "Maya", 6, 'F', false, 8);
	caniPartecipantiGara.push_back(&c1);
	Cane c2 = Cane((char*) "Merlin", 12, 'M', true, 4);
	caniPartecipantiGara.push_back(&c2);
	Cane c3 = Cane((char*) "Elvis", 3, 'M', true, 6);
	caniPartecipantiGara.push_back(&c3);
	Cane c4 = Cane((char*) "Loki", 9, 'M', false, 10);
	caniPartecipantiGara.push_back(&c4);
	Cane c5 = Cane((char*) "Kira", 5, 'F', false, 7);
	caniPartecipantiGara.push_back(&c5);

	unique_ptr<Gatto> g1(new Gatto((char*) "Pasqualino", 2, 'M', true, (char*) "Sei un dormiglione"));
	unique_ptr<Gatto> g2(new Gatto((char*) "Ares", 10, 'M', true, (char*) "Distruggi tutti i divani"));

	cout << endl << "----------------------------------------------------------"
		 << endl << endl;
	u1->descrizione(u1->getNome());
	c1.descrizione(c1.getNome());
	g1->descrizione(g1->getNome());

	cout << endl << "----------------------------------------------------------"
		 << endl << endl;
	u2->stampa();
	cout << endl;
	c2.stampa();
	cout << endl;
	g2->stampa();

	cout << endl << "----------------------------------------------------------"
		 << endl << endl;
	cout << "\tElenco di persone in ordine crescente di età:" << endl << endl;
	VettoreFacade vetFac1 = VettoreFacade(vettoreUmani);
	vetFac1.startSortPrint();

	cout << endl << "----------------------------------------------------------"
		 << endl << endl;
	Umano giudice = Umano((char*) "Marco", (char*) "Montemagno", 38, 'M', (char*) "Codice4Fiscale5");
	giudice.iniziaGara(caniPartecipantiGara);

	cout << endl << "----------------------------------------------------------"
		 << endl << endl;
	g1->giocaCon(move(g2));

	cout << endl << "----------------------------------------------------------"
		 << endl << endl;

	cout << "\tEstrazione dei partecipanti per una partita a Forza 4" << endl;

	char *giocatore1;
	int rand1 = rand() % 101;
	if (rand1 >= 0 && rand1 < 33)
		giocatore1 = u1->getNome();
	else if (rand1 >= 33 && rand1 < 66)
		giocatore1 = u2->getNome();
	else
		giocatore1 = u3->getNome();

	char *giocatore2;
	int rand2;
	int uguali = false;
	do {
		rand2 = rand() % 101;
		if (rand2 >= 0 && rand2 < 33)
			giocatore2 = u1->getNome();
		else if (rand2 >= 33 && rand2 < 66)
			giocatore2 = u2->getNome();
		else
			giocatore2 = u3->getNome();
		uguali = giocatore1 == giocatore2 ? true : false;
	} while (uguali);

	cout << endl << "\tGiocatori estratti: " << giocatore1 << " e " << giocatore2;
	cout << endl << endl << "\t" << giocatore1 << " sarà il giocatore Rosso" << endl;
	cout << "\t" <<giocatore2 << " sarà il giocatore Verde" << endl << endl;

	Forza4 gioco;
	gioco.gioca();

	delete u1;
	u1 = NULL;
	delete u2;
	u2 = NULL;
	delete u3;
	u3 = NULL;

	return 0;
}
