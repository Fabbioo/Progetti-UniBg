#ifndef HEADER_ESSEREVIVENTE_HPP_
#define HEADER_ESSEREVIVENTE_HPP_

#include <iostream>

using namespace std;

class EssereVivente {
	char *nome;
	int anni;
	char sesso;
public:
	EssereVivente(char*, int, char);
	virtual ~EssereVivente();
	char* getNome();
	int getAnni();
	void stampa();
	virtual void descrizione(char*);
	static bool compareByAnni(EssereVivente*, EssereVivente*);
};

#endif /* HEADER_ESSEREVIVENTE_HPP_ */
