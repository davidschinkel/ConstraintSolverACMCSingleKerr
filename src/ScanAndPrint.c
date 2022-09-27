#include "main.h"

// Scanning and printing data files.

void ScanInitialData(parameters *par, ftype **X, char file[])
{
	FILE *stream;
	char string[NMAX], number[NMAX];
	int j_aux, ns, nt, n, n_v, idom, ipot, i, j;
	ftype d_aux;	// d_aux: a ftype auxiliary variable

	stream = fopen(file, "r");

	for (idom = 0; idom < NDOM; idom++)
	{
		j_aux = fscanf(stream, " %d %s ", &ns, string);
		(*par).ns[idom] = ns;
	}
	j_aux = fscanf(stream, " %d %s ", &nt, string);
	(*par).nt = nt;
	(*par).nst[0] = 0;
	for (idom = 0; idom < NDOM; idom++)
	{
		ns = (*par).ns[idom];
		(*par).nst[idom + 1] = (*par).nst[idom] + ns * nt;
	}
	Get_Arrays_n_And_n_v(par);

	j_aux = fscanf(stream, " %s %s ",  number, string); 	d_aux=cast<ftype>(number);	(*par).sig0 = d_aux;
	j_aux = fscanf(stream, " %s %s ",  number, string);     d_aux=cast<ftype>(number);	(*par).sig1 = d_aux;

	*X = new ftype[(*par).ntotal];

	for (idom = 0; idom < NDOM; idom++)
	{ // Reading the Potential values
		ns = (*par).ns[idom];
		nt = (*par).nt;
		for (ipot = 0; ipot < NPOT; ipot++)
		{
			for (j = 0; j < nt; j++)
			{
				for (i = 0; i < ns; i++)
				{
					n_v = Index(*par, idom, ipot, i, j);
					n = (*par).n_of_n_v[n_v];
					if (n >= 0)
					{
						j_aux = fscanf(stream, " %s %s ", number, string);
						d_aux=cast<ftype>(number);
						(*X)[n] = d_aux;
					}
				}
			}
		}
	}
	fclose(stream);
}

void ScanConfig(parameters *parInit, parameters *parGoal)
{
	FILE *stream_config;
	char string[NMAX], number[NMAX], which_data[Ndat], ref_data[NMAX];
	int i_aux, j_aux, idom, ns, nt;  // i_aux: an int auxiliary variable
	ftype d_aux;                    // d_aux: a ftype auxiliary variable

	stream_config = fopen("Config", "r");
	j_aux = fscanf(stream_config, " %d  %s ",  &i_aux, string);		(*parInit).Newton_itmin            = i_aux;
	j_aux = fscanf(stream_config, " %d  %s ",  &i_aux, string);		(*parInit).Newton_itmax            = i_aux;
	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*parInit).Newton_tol = d_aux;
	j_aux = fscanf(stream_config, " %d  %s ",  &i_aux, string);		(*parInit).Newton_verb             = i_aux;

	j_aux = fscanf(stream_config, " %d  %s ",  &i_aux, string);		(*parInit).bicgstab_itmax          = i_aux;
	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*parInit).bicgstab_decr = d_aux;
	j_aux = fscanf(stream_config, " %d  %s ",  &i_aux, string);		(*parInit).bicgstab_verb           = i_aux;

	j_aux = fscanf(stream_config, " %s %s ",  number, string);		d_aux=cast<ftype>(number);	(*parInit).PreCond_Improve = d_aux;
	j_aux = fscanf(stream_config, " %d  %s ",  &i_aux, string);		(*parInit).PreCond_itmax           = i_aux;
	j_aux = fscanf(stream_config, " %d  %s ",  &i_aux, string);		(*parInit).PreCond_FD_Order        = i_aux;

	j_aux = fscanf(stream_config,  "%s %s", which_data, string);
	sprintf((*parInit).which_data, "%s", which_data);
	sprintf((*parGoal).which_data, "%s", which_data);

	*parGoal = *parInit;
	for (idom = 0; idom < NDOM; idom++)
	{
		j_aux = fscanf(stream_config, " %d %s ", &ns, string);
		(*parGoal).ns[idom] = ns;
	}
	j_aux = fscanf(stream_config, " %d %s ", &nt, string);
	(*parGoal).nt = nt;

	(*parGoal).nst[0] = 0;
	for (idom = 0; idom < NDOM; idom++)
	{
		ns = (*parGoal).ns[idom];
		(*parGoal).nst[idom + 1] = (*parGoal).nst[idom] + ns * nt;
	}
	Get_Arrays_n_And_n_v(parGoal);

	j_aux = fscanf(stream_config, " %s %s ", number, string);     d_aux=cast<ftype>(number);	(*parGoal).sig0 = d_aux;
	j_aux = fscanf(stream_config, " %s %s ", number, string);     d_aux=cast<ftype>(number);	(*parGoal).sig1 = d_aux;

	j_aux = fscanf(stream_config, " %d  %s ", &i_aux, string);
	(*parInit).n_Seq = (*parGoal).n_Seq = i_aux;
	j_aux = fscanf(stream_config, " %d  %s ", &i_aux, string);
	(*parInit).Create_ID = (*parGoal).Create_ID = i_aux;

	fclose(stream_config);
}

void ScanFiles(parameters *parInit, ftype **X, parameters *parGoal)
{
	char file[NMAX];

	sprintf(file, "InitialData");
	ScanInitialData(parInit, X, file);
	ScanConfig(parInit, parGoal);
}

void PrintToFile(parameters par, ftype *X, int i_Seq)
{
	char name[NMAX];
	//FILE *stream;
	ofstream fp;
	int idom, ipot, i, j, ns, nt, n, n_v;

	if (i_Seq <  10)              sprintf(name, "./SequenceElements/SequenceElement_%s.00%1d", par.which_data, i_Seq);
	if (i_Seq >= 10 && i_Seq<100) sprintf(name, "./SequenceElements/SequenceElement_%s.0%2d",  par.which_data, i_Seq);
	if (i_Seq >= 100)             sprintf(name, "./SequenceElements/SequenceElement_%s.%3d",   par.which_data, i_Seq);

	fp.open(name);

	for (idom = 0; idom < NDOM; idom++)
	{
		fp << par.ns[idom] << "\t ns_domain_" << idom << endl;
	}
	fp << endl << par.nt << "\t nt" << endl;
	fp.precision(PREC_OUTPUT);
	fp << par.sig0 << "\t sig0" << endl;
	fp << par.sig1 << "\t sig1" << endl;
	fp << endl;

	for (idom = 0; idom < NDOM; idom++)
	{       // Writing the Potential values
		ns = par.ns[idom];
		nt = par.nt;
		for (ipot = 0; ipot < NPOT; ipot++)
		{
			for (j = 0; j < nt; j++)
			{
				for (i = 0; i < ns; i++)
				{
					n_v = Index(par, idom, ipot, i, j);
					n = par.n_of_n_v[n_v];
					if (n >= 0)
						fp << X[n] << "\t (idom,ipot,i,j)=("<<idom<<","<<ipot<<","<<i<<","<<j<<")"<<endl;
				}
			}
		}
	}
	fp.close();
}

void CreateInitialData(parameters par)
{
	int idom, ipot, i, j, ns, nt, Ns, Nt, n, n_v;
	ofstream fp;

	fp.open("Created_InitialData");
	fp.precision(PREC_OUTPUT);

	for (idom = 0; idom < NDOM; idom++)
	{
		fp << par.ns[idom] << "\t ns_domain_" << idom << endl;
	}
	fp << par.nt <<"\t nt" << endl;
	fp<< par.sig0 << "\t sig0" << endl;
	fp << par.sig1 << "\t sig1" << endl;
	fp << endl;

	for (idom = 0; idom < NDOM; idom++)
	{       // Writing the Potential values
		ns = par.ns[idom];
		Ns = ns - 1;
		nt = par.nt;
		Nt = nt - 1;
		for (ipot = 0; ipot < NPOT; ipot++)
		{
			for (j = 0; j < nt; j++)
			{
				for (i = 0; i < ns; i++)
				{

					ftype X_n;

					switch (ipot)
					{
					case 0:
						X_n = 2.;
						break;  // phi   = 2 for Kerr black hole
					case 1:
						X_n = 0.;
						break;  // V^sig = 0 for Kerr black hole
					case 2:
						X_n = 0.;
						break;  // V^mu  = 0 for Kerr black hole
					case 3:
						X_n = 0.;
						break;  // V^phi = 0 for Kerr black hole
					}
					n_v = Index(par, idom, ipot, i, j);
					n = par.n_of_n_v[n_v];
					if (n >= 0)
					{
						fp<< X_n << "\t" <<"(idom,ipot,i,j)=("<<idom<<","<<ipot<<","<<i<<","<<j<<")"<<endl;
					}
				}
			}
		}
	}
	fp.close();
}
