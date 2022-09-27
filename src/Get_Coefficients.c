#include "main.h"

#define hmax 3

#define n1_4 3
#define n2_4 3
#define n3_4 2
#define n4_4 3

#define n1_3 3
#define n2_3 3
#define n3_3 3

#define n1_2 3
#define n2_2 3

#define n1_1 3


void read_config(char name[], int *ns_dat, int *nt_dat, ftype *sig_min,
		ftype *sig_max, ftype *j_BH)
{
	char string[NMAX], number[NMAX];
	int j_aux;
	FILE *stream;

 	stream = fopen(name, "r");
	j_aux = fscanf(stream, " %s\n",     string);
	j_aux = fscanf(stream, " %d %s\n",  ns_dat,  string);
	j_aux = fscanf(stream, " %d %s\n",  nt_dat,  string);
	j_aux = fscanf(stream, " %s %s\n", number, string);
	*sig_min=cast<ftype>(number);
	j_aux = fscanf(stream, " %s %s\n", number, string);
	*sig_max=cast<ftype>(number);
	j_aux = fscanf(stream, " %s %s\n", number, string);
	*j_BH=cast<ftype>(number);

	fclose(stream);
}

int scan_coefficient(char name[], ftype **a, int ns_dat, int nt_dat)
{
	char string[NMAX], number[NMAX];
	int i, j, j_aux;
	ftype d_aux;
	FILE *stream;

	stream = fopen(name, "r");

	if (stream == NULL)
		return 0;
	else
	{
		j_aux = fscanf(stream, " %s ", string);

		for (i = 0; i < ns_dat; i++)
		{
			for (j = 0; j < nt_dat; j++)
			{
				j_aux = fscanf(stream, " %s %s ", number, string);
				d_aux=cast<ftype>(number);
				a[i][j] = d_aux;
			}
		}
		fclose(stream);
		chebft_Zeros_2D(a, ns_dat, nt_dat);

		return 1;
	}
}

int indx_4(parameters par, int i1, int i2, int i3, int i4, int idom, int i,
		int j)
{
	int ns0 = par.ns[0], nt = par.nt;
	return i1 - 1
			+ n1_4
					* (i2 - 1
							+ n2_4
									* (i3 - 1
											+ n3_4
													* (i4 - 1
															+ n4_4
																	* (j
																			+ nt
																					* (i
																							+ ns0
																									* idom)))));
}

void Get_Coefficient_4(parameters *par, ftype *par_a, char prefix[])
{
	int i1, i2, i3, i4, i, j, idom, ns, Ns, nt, Nt, ns_dat, nt_dat, isignum,
			file_exists, na_max = n1_4 * n2_4 * n3_4 * n4_4 * (*par).nt
					* ((*par).ns[0] + (*par).ns[1]);
	ftype sig_min, sig_max, j_BH, sig0 = (*par).sig0, sig1 = (*par).sig1, **a;
	char name[50], name_config[50];

	sprintf(name_config, "../%s/config", (*par).which_data);
	read_config(name_config, &ns_dat, &nt_dat, &sig_min, &sig_max, &j_BH);

	a = dmatrix(0, ns_dat - 1, 0, nt_dat - 1);
	fill0_dvector(par_a, 0, na_max);

	for (i1 = 1; i1 <= n1_4; i1++)
	{
		for (i2 = 1; i2 <= n2_4; i2++)
		{
			for (i3 = 1; i3 <= n3_4; i3++)
			{
				for (i4 = 1; i4 <= n4_4; i4++)
				{

					if ((i1 + i2 + i3 + i4) % 2 == 0)
						isignum = 1;
					else
						isignum = -1;

					sprintf(name, "%s%1d%1d%1d%1d.dat", prefix, i1, i2, i3, i4);
					file_exists = scan_coefficient(name, a, ns_dat, nt_dat);

					if (file_exists == 1)
					{

						for (idom = 0; idom < NDOM; idom++)
						{
							ns = (*par).ns[idom];
							Ns = ns - 1;
							nt = (*par).nt;
							Nt = nt - 1;
							for (j = 0; j < nt; j++)
							{
								for (i = 0; i < ns; i++)
								{
									int indx = indx_4(*par, i1, i2, i3, i4,
											idom, i, j);
									ftype s = sqr(sin(Pih * i / Ns)), t =
											-cos(Pi * j / Nt), sig, mu = t;

									if (idom == 0)
										sig = sig0 * s;
									else
										sig = sig1 + s * (sig0 - sig1);

									if (mu > 0)
										par_a[indx] = chebevxy(0., sig_max, 0.,
												1., a, ns_dat, nt_dat, sig, mu);
									else
										par_a[indx] = chebevxy(0., sig_max, 0.,
												1., a, ns_dat, nt_dat, sig, -mu)
												* isignum;
								}
							}
						}
					}
				}
			}
		}
	}
	free_dmatrix(a, 0, ns_dat - 1, 0, nt_dat - 1);
}

int indx_3(parameters par, int i1, int i2, int i3, int idom, int i, int j)
{
	int ns0 = par.ns[0], nt = par.nt;
	return i1 - 1
			+ n1_3
					* (i2 - 1
							+ n2_3
									* (i3 - 1
											+ n3_3 * (j + nt * (i + ns0 * idom))));
}

void Get_Coefficient_3(parameters *par, ftype *par_b, char prefix[])
{
	int i1, i2, i3, i, j, idom, ns, Ns, nt, Nt, ns_dat, nt_dat, isignum,
			file_exists, nb_max = n1_3 * n2_3 * n3_3 * (*par).nt
					* ((*par).ns[0] + (*par).ns[1]);
	ftype sig_min, sig_max, j_BH, sig0 = (*par).sig0, sig1 = (*par).sig1, **b;
	char name[50], name_config[50];

	sprintf(name_config, "../%s/config", (*par).which_data);
	read_config(name_config, &ns_dat, &nt_dat, &sig_min, &sig_max, &j_BH);

	b = dmatrix(0, ns_dat - 1, 0, nt_dat - 1);
	fill0_dvector(par_b, 0, nb_max);

	for (i1 = 1; i1 <= n1_3; i1++)
	{
		for (i2 = 1; i2 <= n2_3; i2++)
		{
			for (i3 = 1; i3 <= n3_3; i3++)
			{

				if ((i1 + i2 + i3) % 2 == 0)
					isignum = -1;
				else
					isignum = 1;

				sprintf(name, "%s%1d%1d%1d.dat", prefix, i1, i2, i3);
				file_exists = scan_coefficient(name, b, ns_dat, nt_dat);

				if (file_exists == 1)
				{

					for (idom = 0; idom < NDOM; idom++)
					{
						ns = (*par).ns[idom];
						Ns = ns - 1;
						nt = (*par).nt;
						Nt = nt - 1;
						for (j = 0; j < nt; j++)
						{
							for (i = 0; i < ns; i++)
							{
								int indx = indx_3(*par, i1, i2, i3, idom, i, j);
								ftype s = sqr(sin(Pih * i / Ns)), t = -cos(
										Pi * j / Nt), sig, mu = t;

								if (idom == 0)
									sig = sig0 * s;
								else
									sig = sig1 + s * (sig0 - sig1);

								if (mu > 0)
									par_b[indx] = chebevxy(0., sig_max, 0., 1.,
											b, ns_dat, nt_dat, sig, mu);
								else
									par_b[indx] = chebevxy(0., sig_max, 0., 1.,
											b, ns_dat, nt_dat, sig, -mu)
											* isignum;
							}
						}
					}
				}
			}
		}
	}
	free_dmatrix(b, 0, ns_dat - 1, 0, nt_dat - 1);
}

int indx_2(parameters par, int i1, int i2, int idom, int i, int j)
{
	int ns0 = par.ns[0], nt = par.nt;
	return i1 - 1 + n1_2 * (i2 - 1 + n2_2 * (j + nt * (i + ns0 * idom)));
}

void Get_Coefficient_2(parameters *par, ftype *par_c, char prefix[])
{
	int i1, i2, i, j, idom, ns, Ns, nt, Nt, ns_dat, nt_dat, isignum,
			file_exists, nc_max = n1_2 * n2_2 * (*par).nt
					* ((*par).ns[0] + (*par).ns[1]);
	ftype sig_min, sig_max, j_BH, sig0 = (*par).sig0, sig1 = (*par).sig1, **c;
	char name[50], name_config[50];

	sprintf(name_config, "../%s/config", (*par).which_data);
	read_config(name_config, &ns_dat, &nt_dat, &sig_min, &sig_max, &j_BH);

	c = dmatrix(0, ns_dat - 1, 0, nt_dat - 1);
	fill0_dvector(par_c, 0, nc_max);

	for (i1 = 1; i1 <= n1_2; i1++)
	{
		for (i2 = 1; i2 <= n2_2; i2++)
		{

			if ((i1 + i2) % 2 == 0)
				isignum = 1;
			else
				isignum = -1;

			sprintf(name, "%s%1d%1d.dat", prefix, i1, i2);
			file_exists = scan_coefficient(name, c, ns_dat, nt_dat);

			if (file_exists == 1)
			{

				for (idom = 0; idom < NDOM; idom++)
				{
					ns = (*par).ns[idom];
					Ns = ns - 1;
					nt = (*par).nt;
					Nt = nt - 1;
					for (j = 0; j < nt; j++)
					{
						for (i = 0; i < ns; i++)
						{
							int indx = indx_2(*par, i1, i2, idom, i, j);
							ftype s = sqr(sin(Pih * i / Ns)), t = -cos(
									Pi * j / Nt), sig, mu = t;

							if (idom == 0)
								sig = sig0 * s;
							else
								sig = sig1 + s * (sig0 - sig1);

							if (mu > 0)
								par_c[indx] = chebevxy(0., sig_max, 0., 1., c,
										ns_dat, nt_dat, sig, mu);
							else
								par_c[indx] = chebevxy(0., sig_max, 0., 1., c,
										ns_dat, nt_dat, sig, -mu) * isignum;
						}
					}
				}
			}
		}
	}
	free_dmatrix(c, 0, ns_dat - 1, 0, nt_dat - 1);
}

int indx_1(parameters par, int i1, int idom, int i, int j)
{
	int ns0 = par.ns[0], nt = par.nt;
	return i1 - 1 + n1_1 * (j + nt * (i + ns0 * idom));
}

void Get_Coefficient_1(parameters *par, ftype *par_d, char prefix[])
{
	int i1, i, j, idom, ns, Ns, nt, Nt, ns_dat, nt_dat, isignum, file_exists,
			nd_max = n1_1 * (*par).nt * ((*par).ns[0] + (*par).ns[1]);
	ftype sig_min, sig_max, j_BH, sig0 = (*par).sig0, sig1 = (*par).sig1, **d;
	char name[50], name_config[50];

	sprintf(name_config, "../%s/config", (*par).which_data);
	read_config(name_config, &ns_dat, &nt_dat, &sig_min, &sig_max, &j_BH);

	d = dmatrix(0, ns_dat - 1, 0, nt_dat - 1);
	fill0_dvector(par_d, 0, nd_max);

	for (i1 = 1; i1 <= n1_1; i1++)
	{

		if (i1 % 2 == 0)
			isignum = -1;
		else
			isignum = +1;

		sprintf(name, "%s%1d.dat", prefix, i1);
		file_exists = scan_coefficient(name, d, ns_dat, nt_dat);

		if (file_exists == 1)
		{

			for (idom = 0; idom < NDOM; idom++)
			{
				ns = (*par).ns[idom];
				Ns = ns - 1;
				nt = (*par).nt;
				Nt = nt - 1;
				for (j = 0; j < nt; j++)
				{
					for (i = 0; i < ns; i++)
					{
						int indx = indx_1(*par, i1, idom, i, j);
						ftype s = sqr(sin(Pih * i / Ns)), t = -cos(
								Pi * j / Nt), sig, mu = t;

						if (idom == 0)
							sig = sig0 * s;
						else
							sig = sig1 + s * (sig0 - sig1);

						if (mu > 0)
							par_d[indx] = chebevxy(0., sig_max, 0., 1., d,
									ns_dat, nt_dat, sig, mu);
						else
							par_d[indx] = chebevxy(0., sig_max, 0., 1., d,
									ns_dat, nt_dat, sig, -mu) * isignum;
					}
				}
			}
		}
	}
	free_dmatrix(d, 0, ns_dat - 1, 0, nt_dat - 1);
}

int indx_H(parameters par, int i1, int idom, int i, int j)
{
	int ns0 = par.ns[0], nt = par.nt;
	return i1 - 1 + Hmax * (j + nt * (i + ns0 * idom));
}

void Get_Coefficient_H(parameters *par, ftype *par_H, char prefix[])
// void Get_Coefficient_H(parameters *par, ftype *(*par).H, char prefix[])
{
	int i1, i, j, idom, ns, Ns, nt, Nt, ns_dat, nt_dat, isignum, file_exists,
			nH_max = Hmax * (*par).nt * ((*par).ns[0] + (*par).ns[1]);
	ftype sig_min, sig_max, j_BH, sig0 = (*par).sig0, sig1 = (*par).sig1, **H;
	char name[50], name_config[50];

	sprintf(name_config, "../%s/config", (*par).which_data);
	read_config(name_config, &ns_dat, &nt_dat, &sig_min, &sig_max, &j_BH);

	H = dmatrix(0, ns_dat - 1, 0, nt_dat - 1);
	fill0_dvector(par_H, 0, nH_max);

	for (i1 = 1; i1 <= Hmax; i1++)
	{

		if (i1 == 3 || i1 == 5)
			isignum = -1;
		else
			isignum = +1;

		sprintf(name, "%s%1d.dat", prefix, i1);
		file_exists = scan_coefficient(name, H, ns_dat, nt_dat);

		if (file_exists == 1)
		{

			for (idom = 0; idom < NDOM; idom++)
			{
				ns = (*par).ns[idom];
				Ns = ns - 1;
				nt = (*par).nt;
				Nt = nt - 1;
				for (j = 0; j < nt; j++)
				{
					for (i = 0; i < ns; i++)
					{
						int indx = indx_H(*par, i1, idom, i, j);
						ftype s = sqr(sin(Pih * i / Ns)), t = -cos(
								Pi * j / Nt), sig, mu = t;

						if (idom == 0)
							sig = sig0 * s;
						else
							sig = sig1 + s * (sig0 - sig1);

						if (mu > 0)
							par_H[indx] = chebevxy(0., sig_max, 0., 1., H,
									ns_dat, nt_dat, sig, mu);
						else
							par_H[indx] = chebevxy(0., sig_max, 0., 1., H,
									ns_dat, nt_dat, sig, -mu) * isignum;
					}
				}
			}
		}
	}
	free_dmatrix(H, 0, ns_dat - 1, 0, nt_dat - 1);
}

int indx_Horizon(parameters par, int i1, int j)
{
	return i1 - 1 + hmax * j;
}

void Get_Coefficient_Horizon(parameters *par, ftype *par_h, char prefix[],
		int odd)
{
	int i1, j, nt, Nt, ns_dat, nt_dat, isignum, file_exists, nh_max = hmax
			* (*par).nt;
	ftype sig_min, sig_max, j_BH, sig1 = (*par).sig1, **h;
	char name[50], name_config[50];

	sprintf(name_config, "../%s/config", (*par).which_data);
	read_config(name_config, &ns_dat, &nt_dat, &sig_min, &sig_max, &j_BH);

	h = dmatrix(0, ns_dat - 1, 0, nt_dat - 1);
	fill0_dvector(par_h, 0, nh_max);

	for (i1 = 1; i1 <= hmax; i1++)
	{
		if (i1 == 2)
			isignum = 1 - 2 * odd;
		else
			isignum = 1;

		sprintf(name, "%s%1d.dat", prefix, i1);
		file_exists = scan_coefficient(name, h, ns_dat, nt_dat);

		if (file_exists == 1)
		{

			nt = (*par).nt;
			Nt = nt - 1;
			for (j = 0; j < nt; j++)
			{
				int indx = indx_Horizon(*par, i1, j);
				ftype t = -cos(Pi * j / Nt), sig = sig1, mu = t;

				if (mu > 0)
					par_h[indx] = chebevxy(sig_min, sig_max, 0., 1., h, ns_dat,
							nt_dat, sig, mu);
				else
					par_h[indx] = chebevxy(sig_min, sig_max, 0., 1., h, ns_dat,
							nt_dat, sig, -mu) * isignum;
			}
		}
	}
	free_dmatrix(h, 0, ns_dat - 1, 0, nt_dat - 1);
}

void Get_Coefficients(parameters *par)
{
	char prefix[NMAX], which_data[NMAX];
	int na_max = n1_4 * n2_4 * n3_4 * n4_4 * (*par).nt	* ((*par).ns[0] + (*par).ns[1]),
		nb_max = n1_3 * n2_3 * n3_3	* (*par).nt * ((*par).ns[0] + (*par).ns[1]),
		nc_max = n1_2 * n2_2 * (*par).nt * ((*par).ns[0] + (*par).ns[1]),
		nd_max = n1_1 * (*par).nt * ((*par).ns[0] + (*par).ns[1]),
		nH_max = Hmax * (*par).nt * ((*par).ns[0] + (*par).ns[1]),
		nh_max = hmax * (*par).nt;

	cout << " Get Coefficients ... please wait ... " << flush;
	sprintf(which_data, "%s", (*par).which_data);

	(*par).a = new ftype[na_max];
	sprintf(prefix, "../%s/a", which_data);
	Get_Coefficient_4(par, (*par).a, prefix);

	(*par).kdd = new ftype[na_max];
	sprintf(prefix, "../%s/kdd", which_data);
	Get_Coefficient_4(par, (*par).kdd, prefix);

	(*par).kuu = new ftype[na_max];
	sprintf(prefix, "../%s/kuu", which_data);
	Get_Coefficient_4(par, (*par).kuu, prefix);

	(*par).b = new ftype[nb_max];
	sprintf(prefix, "../%s/b", which_data);
	Get_Coefficient_3(par, (*par).b, prefix);

	(*par).ldd = new ftype[nb_max];
	sprintf(prefix, "../%s/ldd", which_data);
	Get_Coefficient_3(par, (*par).ldd, prefix);

	(*par).luu = new ftype[nb_max];
	sprintf(prefix, "../%s/luu", which_data);
	Get_Coefficient_3(par, (*par).luu, prefix);

	(*par).c = new ftype[nc_max];
	sprintf(prefix, "../%s/c", which_data);
	Get_Coefficient_2(par, (*par).c, prefix);

	(*par).mdd = new ftype[nc_max];
	sprintf(prefix, "../%s/mdd", which_data);
	Get_Coefficient_2(par, (*par).mdd, prefix);

	(*par).muu = new ftype[nc_max];
	sprintf(prefix, "../%s/muu", which_data);
	Get_Coefficient_2(par, (*par).muu, prefix);

	(*par).d = new ftype[nd_max];
	sprintf(prefix, "../%s/d", which_data);
	Get_Coefficient_1(par, (*par).d, prefix);

	(*par).e = new ftype[nd_max];
	sprintf(prefix, "../%s/e", which_data);
	Get_Coefficient_1(par, (*par).e, prefix);

	(*par).H = new ftype[nH_max];
	sprintf(prefix, "../%s/H", which_data);
	Get_Coefficient_H(par, (*par).H, prefix);

	(*par).h = new ftype[nh_max];
	sprintf(prefix, "../%s/hh", which_data);
	Get_Coefficient_Horizon(par, (*par).h, prefix, 0);

	(*par).s = new ftype[nh_max];
	sprintf(prefix, "../%s/s", which_data);
	Get_Coefficient_Horizon(par, (*par).s, prefix, 1);

	cout << "Done.\n" << endl;
}

void Free_Coefficients(parameters *par)
{
	int na_max = n1_4 * n2_4 * n3_4 * n4_4 * (*par).nt
			* ((*par).ns[0] + (*par).ns[1]), nb_max = n1_3 * n2_3 * n3_3
			* (*par).nt * ((*par).ns[0] + (*par).ns[1]), nc_max = n1_2 * n2_2
			* (*par).nt * ((*par).ns[0] + (*par).ns[1]), nd_max = n1_1
			* (*par).nt * ((*par).ns[0] + (*par).ns[1]), nH_max = Hmax
			* (*par).nt * ((*par).ns[0] + (*par).ns[1]), nh_max = hmax
			* (*par).nt, ns_max = hmax * (*par).nt;

	delete[] (*par).a;
	delete[] (*par).kdd;
	delete[] (*par).kuu;
	delete[] (*par).b;
	delete[] (*par).ldd;
	delete[] (*par).luu;
	delete[] (*par).c;
	delete[] (*par).mdd;
	delete[] (*par).muu;

	delete[] (*par).d;
	delete[] (*par).e;
	delete[] (*par).H;
	delete[] (*par).h;
	delete[] (*par).s;
}

#undef hmax

#undef n1_4
#undef n2_4
#undef n3_4
#undef n4_4

#undef n1_3
#undef n2_3
#undef n3_3

#undef n1_2
#undef n2_2

#undef n1_1

