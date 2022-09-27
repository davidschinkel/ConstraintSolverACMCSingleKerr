#include "main.h"
#define pI parInit
#define pG parGoal


int main(int argc, char **argv)
{
	int Newton_iter, i_Seq;
	parameters parInit, par, parGoal;
	ftype *X0, *X, h, sig_min, sig_max, j_BH;
	int ns_dat, nt_dat;
	FILE *stream;
	char name_config[50];

	cout <<" ********************************************************************" << endl;
	cout <<" *                        ACMC_Data   27/07/12                      *" << endl;
	cout <<" ********************************************************************" << endl;

	//lesen der Anfangsdaten und der Config
	ScanFiles(&parInit, &X0, &parGoal);
	par = parGoal;

	//ändern der Auflösung der Anfangsdaten
	X = new ftype[par.ntotal];
	GetNewResolution(parInit, par, X0, X);
	delete[] X0;

	//lesen der config der Koeffizientenfunktionen
	sprintf(name_config, "../%s/config", par.which_data);
	read_config(name_config, &ns_dat, &nt_dat, &sig_min, &sig_max, &j_BH);

	//erzeigen von Anfangsdaten und beenden des Programms
	if (par.Create_ID == 1)
	{
		par.sig1 = 2. / (1. + sqrt(1. - sqr(j_BH)));
		par.sig0 = 0.5 * par.sig1;
		CreateInitialData(par);
		cout.precision(4);
		cout << fixed;
		cout << "Create InitialData with"<< endl;
		cout << "\t sig0    = "<< par.sig0 <<"\tsig1    = "<<par.sig1<<" for data in "<<par.which_data<<", j_BH =" << j_BH<<endl;
		cout << "\t sig_min = "<<sig_min<< "\t   sig_max = "<<sig_max << endl;
		exit(1);
	}

	//Newton-Raphson-Löser
	for (i_Seq = 0; i_Seq <= par.n_Seq; i_Seq++)
	{
		//Berechnung der Schrittweite h
		if (par.n_Seq > 0)
			h = 1. * i_Seq / par.n_Seq;
		else
			h = 0.;

		cout << " -----------------------------------------------------" << endl;
		cout << " Calculation of Sequence element No. \t"<<i_Seq<< endl;
		cout << " -----------------------------------------------------" << endl;

		//Berechnung der Zielparameter
		par.sig0 = pI.sig0 + h * (pG.sig0 - pI.sig0);
		par.sig1 = pI.sig1 + h * (pG.sig1 - pI.sig1);
		par.j_BH = j_BH;

		//Lesen der Koeffizientenfunktionen
		Get_Coefficients(&par);

		//Newton-Verfahren
		cout << scientific;
		Newton_iter = newton(par, X);
		cout << fixed;

		//Speichern des Newton-Schrittes
		cout << "\n Writing File 'SequenceElement_"<<par.which_data;
		if (i_Seq < 10)
			cout<<"00"<<i_Seq << endl;
		if (i_Seq >= 10 && i_Seq < 100)
			cout<<"0"<<i_Seq << endl;
		if (i_Seq >= 100)
			cout <<i_Seq << endl;
		cout << endl;
		PrintToFile(par, X, i_Seq);

		//Analyse der Daten (Bondi-Masse und Drehimpuls)
		BondiMass(par, X);

		if (i_Seq < par.n_Seq)
			Free_Coefficients(&par);
	}

	delete[] X;
	Free_Arrays_n_And_n_v(&parGoal);
	Free_Arrays_n_And_n_v(&parInit);
	Free_Coefficients(&par);

	return 1;
}
#undef pI
#undef pG

