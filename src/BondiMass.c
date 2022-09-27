#include "main.h"

void LVmm(parameters par, derivs_2D v, ftype *LVmm_Scri)
{
	int j, ii, jj, a, nt = par.nt;
	ftype LVdd_ij, fac = 1 / (par.sig0);
	int i = 0, idom = 0; //Scri

	for (j = 0; j < nt; j++)
	{
		LVmm_Scri[j] = 0;
		for (ii = 1; ii <= 3; ii++)
		{
			for (jj = 1; jj <= 3; jj++)
			{
				LVdd_ij = 0;

				for (a = 1; a <= 3; a++)
				{
					int indx = Index(par, idom, a, i, j);

					LVdd_ij += par.kdd[indx_4(par, ii, jj, 1, a, idom, i, j)]
							* v.d1[indx] * fac
							+ par.kdd[indx_4(par, ii, jj, 2, a, idom, i, j)]
									* v.d2[indx]
							+ par.ldd[indx_3(par, ii, jj, a, idom, i, j)]
									* v.d0[indx];
				}
				LVmm_Scri[j] += 0.5 * LVdd_ij * qcinv_Scri(par, ii, jj, j); //regularisierung ist in quinv eingearbeitet
			}
		}
	}
}

ftype qcinv_Scri(parameters par, int ii, int jj, int j)
{ //conform
//output of quinv - regularisierung beachten!
	ftype result;

	if (ii != jj)
		return 0;
	else if (ii == 1)
		return 0;
	if (ii == 2)
		return 4; //q^22=(1-mu^2)Q^22 -> Q^22=4
	if (ii == 3)
		return 4; //q^33=1/(1-mu^2)Q^33 -> Q^33=4

	cout << "FEHLER in qcinv!" << endl;

	return -1;
}

ftype gammac110_Scri(ftype j)
{
	return 2000 / (4260 + 34 * j + 25 * sqr(j));
}

ftype sqrtdetqc_Scri()
{
	return 0.25;
}

void getphi3(ftype **phi, ftype **phi3, int ns, int nt, ftype ds)
{
	int j, i;
	ftype p[NMAX], dp[NMAX], d2p[NMAX], d3p[NMAX];

	for (j = 0; j < nt; j++) // Calculation of Derivatives w.r.t. s-Direction
	{
		for (i = 0; i < ns; i++)
		{	// (Chebyshev_Extremes)
			p[i] = phi[i][j];
		}
		chebft_Extremes(p, ns, 0);
		chder(0., 1., p, dp, ns);
		chder(0., 1., dp, d2p, ns);
		chder(0., 1., d2p, d3p, ns);
		chebft_Extremes(d3p, ns, 1);
		for (i = 0; i < ns; i++)
		{
			phi3[i][j] = d3p[i] / (6 * ds * ds * ds);
		}
	}
}

ftype sqrtdetqc(parameters par, ftype sig, ftype mu, ftype B, ftype Hd2)
{	//conform
	ftype detq2, a = par.j_BH;

	detq2 = ((4 + 16 * B * pow(Hd2, 2) * (-1 + pow(mu, 2))
			+ pow(a, 2) * pow(mu, 2) * pow(sig, 2))
			* (16 + pow(a, 4) * pow(mu, 2) * pow(sig, 4)
					+ 4 * pow(a, 2) * pow(sig, 2)
							* (1 - pow(mu, 2) * (-1 + sig) + sig))
			+ 16 * (-1 + pow(mu, 2)) * pow(sig, 2)
					* (4 - 4 * sig + pow(a, 2) * pow(sig, 2))
					* (4 + pow(a, 2) * pow(mu, 2) * pow(sig, 2))
					* pow(A2(sig, mu, a), 2))
			/ (256. * (4 + pow(a, 2) * pow(mu, 2) * pow(sig, 2)));
	return sqrt(detq2);
}

void GetAdd3j(parameters par, derivs_2D v, ftype **Add31, ftype **Add32,
		ftype **Add33)
{
	int i, j, ii, jj, a, ns = par.ns[1], nt = par.nt, idom = 1;
	ftype fac1 = 1 / (par.sig0 - par.sig1), aux, test;

	for (i = 0; i < ns; i++)
	{	//iteration over gridpoints
		for (j = 0; j < nt; j++)
		{	//iteration over gridpoints
			ii = 3;
			jj = 1; //tensor element ^(ii,jj)
			aux = par.mdd[indx_2(par, ii, jj, idom, i, j)];
			for (a = 1; a <= 3; a++)
			{
				int indx = Index(par, idom, a, i, j);
				aux += par.kdd[indx_4(par, ii, jj, 1, a, idom, i, j)]
						* v.d1[indx] * fac1
						+ par.kdd[indx_4(par, ii, jj, 2, a, idom, i, j)]
								* v.d2[indx]
						+ par.ldd[indx_3(par, ii, jj, a, idom, i, j)]
								* v.d0[indx];
			}
			Add31[i][j] = aux;

			ii = 3;
			jj = 2; //tensor element ^(ii,jj)
			aux = par.mdd[indx_2(par, ii, jj, idom, i, j)];
			for (a = 1; a <= 3; a++)
			{
				int indx = Index(par, idom, a, i, j);
				aux += par.kdd[indx_4(par, ii, jj, 1, a, idom, i, j)]
						* v.d1[indx] * fac1
						+ par.kdd[indx_4(par, ii, jj, 2, a, idom, i, j)]
								* v.d2[indx]
						+ par.ldd[indx_3(par, ii, jj, a, idom, i, j)]
								* v.d0[indx];
			}
			Add32[i][j] = aux;

			ii = 3;
			jj = 3; //tensor element ^(ii,jj)
			aux = par.mdd[indx_2(par, ii, jj, idom, i, j)];
			for (a = 1; a <= 3; a++)
			{
				int indx = Index(par, idom, a, i, j);
				aux += par.kdd[indx_4(par, ii, jj, 1, a, idom, i, j)]
						* v.d1[indx] * fac1
						+ par.kdd[indx_4(par, ii, jj, 2, a, idom, i, j)]
								* v.d2[indx]
						+ par.ldd[indx_3(par, ii, jj, a, idom, i, j)]
								* v.d0[indx];
			}
			Add33[i][j] = aux;
		}
	}
}

ftype gammacuuij(parameters par, ftype sig, ftype mu, ftype i, ftype j)
{
	ftype a = par.j_BH;

	if (1 == i && 1 == j)
	{
		return (4
				* (16 + pow(a, 4) * pow(mu, 2) * pow(sig, 4)
						+ 4 * pow(a, 2) * pow(sig, 2)
								* (1 - pow(mu, 2) * (-1 + sig) + sig)
						+ 16 * (-1 + pow(mu, 2)) * pow(sig, 2)
								* (4 - 4 * sig + pow(a, 2) * pow(sig, 2))
								* pow(A2(sig, mu, a), 2)))
				/ ((4 + pow(a, 2) * pow(mu, 2) * pow(sig, 2))
						* (16 - pow(a, 2) + pow(a, 2) * pow(mu, 2) + 16 * sig
								- 4 * pow(a, 2) * sig
								- 4 * pow(a, 2) * pow(sig, 2)
								+ 4
										* (4 + (-8 + pow(a, 2)) * pow(sig, 2)
												+ 2 * pow(a, 2) * pow(sig, 3))
										* A1(sig, mu, a)
								- 4 * pow(sig, 2)
										* (4 - 4 * sig
												+ pow(a, 2) * pow(sig, 2))
										* pow(A1(sig, mu, a), 2)
								+ 16 * (-1 + pow(mu, 2))
										* pow(A2(sig, mu, a), 2)));
	}
	else if ((1 == i && 2 == j) || (2 == i && 1 == j))
	{
		return (-32 * (-1 + pow(mu, 2))
				* (-4 - (-8 + pow(a, 2)) * pow(sig, 2)
						- 2 * pow(a, 2) * pow(sig, 3)
						+ 2 * pow(sig, 2)
								* (4 - 4 * sig + pow(a, 2) * pow(sig, 2))
								* A1(sig, mu, a)) * A2(sig, mu, a))
				/ ((4 + pow(a, 2) * pow(mu, 2) * pow(sig, 2))
						* (16 - pow(a, 2) + pow(a, 2) * pow(mu, 2) + 16 * sig
								- 4 * pow(a, 2) * sig
								- 4 * pow(a, 2) * pow(sig, 2)
								+ 4
										* (4 + (-8 + pow(a, 2)) * pow(sig, 2)
												+ 2 * pow(a, 2) * pow(sig, 3))
										* A1(sig, mu, a)
								- 4 * pow(sig, 2)
										* (4 - 4 * sig
												+ pow(a, 2) * pow(sig, 2))
										* pow(A1(sig, mu, a), 2)
								+ 16 * (-1 + pow(mu, 2))
										* pow(A2(sig, mu, a), 2)));
	}
	else if (1 == i && 3 == j)
	{
		return (-8 * a
				* (4 + 8 * sig + 8 * pow(sig, 2)
						+ pow(a, 2) * pow(mu, 2) * pow(sig, 2)
						- 8 * pow(sig, 3) * A1(sig, mu, a)
						+ 16 * (-1 + pow(mu, 2)) * pow(sig, 2)
								* pow(A2(sig, mu, a), 2)))
				/ ((4 + pow(a, 2) * pow(mu, 2) * pow(sig, 2))
						* (16 - pow(a, 2) + pow(a, 2) * pow(mu, 2) + 16 * sig
								- 4 * pow(a, 2) * sig
								- 4 * pow(a, 2) * pow(sig, 2)
								+ 4
										* (4 + (-8 + pow(a, 2)) * pow(sig, 2)
												+ 2 * pow(a, 2) * pow(sig, 3))
										* A1(sig, mu, a)
								- 4 * pow(sig, 2)
										* (4 - 4 * sig
												+ pow(a, 2) * pow(sig, 2))
										* pow(A1(sig, mu, a), 2)
								+ 16 * (-1 + pow(mu, 2))
										* pow(A2(sig, mu, a), 2)));
	}
	else if (2 == i && 2 == j)
	{
		return (-16 * (-1 + pow(mu, 2))
				* (16 * (1 + sig)
						+ pow(a, 2) * (pow(mu, 2) - pow(1 + 2 * sig, 2))
						+ 4
								* (4 + (-8 + pow(a, 2)) * pow(sig, 2)
										+ 2 * pow(a, 2) * pow(sig, 3))
								* A1(sig, mu, a)
						- 4 * pow(sig, 2)
								* (4 - 4 * sig + pow(a, 2) * pow(sig, 2))
								* pow(A1(sig, mu, a), 2)))
				/ ((4 + pow(a, 2) * pow(mu, 2) * pow(sig, 2))
						* (16 - pow(a, 2) + pow(a, 2) * pow(mu, 2) + 16 * sig
								- 4 * pow(a, 2) * sig
								- 4 * pow(a, 2) * pow(sig, 2)
								+ 4
										* (4 + (-8 + pow(a, 2)) * pow(sig, 2)
												+ 2 * pow(a, 2) * pow(sig, 3))
										* A1(sig, mu, a)
								- 4 * pow(sig, 2)
										* (4 - 4 * sig
												+ pow(a, 2) * pow(sig, 2))
										* pow(A1(sig, mu, a), 2)
								+ 16 * (-1 + pow(mu, 2))
										* pow(A2(sig, mu, a), 2)));
	}
	else if (2 == i && 3 == j)
	{
		return (64 * a * (-1 + pow(mu, 2))
				* (-1 - 2 * sig + 2 * pow(sig, 2) * A1(sig, mu, a))
				* A2(sig, mu, a))
				/ ((4 + pow(a, 2) * pow(mu, 2) * pow(sig, 2))
						* (16 - pow(a, 2) + pow(a, 2) * pow(mu, 2) + 16 * sig
								- 4 * pow(a, 2) * sig
								- 4 * pow(a, 2) * pow(sig, 2)
								+ 4
										* (4 + (-8 + pow(a, 2)) * pow(sig, 2)
												+ 2 * pow(a, 2) * pow(sig, 3))
										* A1(sig, mu, a)
								- 4 * pow(sig, 2)
										* (4 - 4 * sig
												+ pow(a, 2) * pow(sig, 2))
										* pow(A1(sig, mu, a), 2)
								+ 16 * (-1 + pow(mu, 2))
										* pow(A2(sig, mu, a), 2)));
	}
	else
	{
		cout <<"Fehler in gammacuuij!"<<endl;
		exit(1);
	}
}

ftype gammacdd3i(ftype sig, ftype mu, ftype j, int i)
{
	if (i == 1)
	{
		return -(j * (-1 + pow(mu, 2))
				* (4 + sig * (8 + (8 + pow(j, 2) * pow(mu, 2)) * sig)
						- 8 * pow(sig, 3) * A1(sig, mu, j)))
				/ (8. * (4 + pow(j, 2) * pow(mu, 2) * pow(sig, 2)));
	}
	if (i == 2)
	{
		return (j * (-1 + pow(mu, 2)) * pow(sig, 3) * A2(sig, mu, j))
				/ (4 + pow(j, 2) * pow(mu, 2) * pow(sig, 2));
	}
	if (i == 3)
	{
		return -((-1 + pow(mu, 2))
				* (16
						+ pow(j, 2) * pow(sig, 2)
								* (4 * (1 + sig)
										+ pow(mu, 2)
												* (4
														+ sig
																* (-4
																		+ pow(
																				j,
																				2)
																				* sig)))))
				/ (16. * (4 + pow(j, 2) * pow(mu, 2) * pow(sig, 2)));
	}
	cout << "Fehler in gammacdd3i!" << endl;
}

ftype A1(ftype sig, ftype mu, ftype j)
{ //dA/ds
	return (3
			* (-2 * sig
					* (3515625 * pow(j, 8) * pow(6 - 5 * pow(mu, 2), 4) * sig
							- 19125000 * pow(j, 7)
									* pow(-6 + 5 * pow(mu, 2), 3) * sig
							+ 576000 * j * (133329750 + 268101149 * sig)
							+ 160000 * (7082911500 + 28585055449 * sig)
							- 15000 * pow(j, 6) * pow(6 - 5 * pow(mu, 2), 2)
									* (-93750
											+ (-461101 + 798750 * pow(mu, 2))
													* sig)
							+ 40800 * pow(j, 5) * (-6 + 5 * pow(mu, 2))
									* (-93750
											+ (-938617 + 1198125 * pow(mu, 2))
													* sig)
							- 1920 * pow(j, 3)
									* (-846431250 + 3355017843 * sig
											+ 625 * pow(mu, 2)
													* (2688750 + 7111457 * sig))
							- 3200 * pow(j, 2)
									* (20808990000 + 76399756193 * sig
											+ 11250 * pow(mu, 2)
													* (1744250 + 10373597 * sig))
							+ 16 * pow(j, 4)
									* (-725724937500 + 3371453243689 * sig
											+ 2343750 * pow(mu, 4)
													* (-3750 + 418321 * sig)
											- 22500 * pow(mu, 2)
													* (-23881250
															+ 69264921 * sig)))
					* (1171875 * pow(j, 8) * pow(-6 + 5 * pow(mu, 2), 3)
							* (-6 + 5 * pow(mu, 2) - 20 * sig)
							- 6375000 * pow(j, 7)
									* pow(6 - 5 * pow(mu, 2), 2)
									* (-6 + 5 * pow(mu, 2) - 15 * sig)
							+ 64000 * j * (212878097 + 1042922485 * sig)
							+ 160000 * (1950390283 + 4472788780 * sig)
							- 640 * pow(j, 3)
									* (-2379695207 - 7867637100 * sig
											+ 625 * pow(mu, 2)
													* (4022557 + 15168270 * sig))
							- 3200 * pow(j, 2)
									* (-13540744969 - 9981770940 * sig
											+ 1250 * pow(mu, 2)
													* (12522241 + 19999355 * sig))
							+ 800 * pow(j, 5)
									* (102751434 + 31949695 * sig
											+ 9375 * pow(mu, 4)
													* (8313 + 17075 * sig)
											- 5 * pow(mu, 2)
													* (33173239 + 50814750 * sig))
							+ 5000 * pow(j, 6)
									* (28749636 - 70250940 * sig
											+ 5 * pow(mu, 2)
													* (-12109212
															+ pow(mu, 2)
																	* (9765505
																			- 9341250
																					* sig)
															+ 25367990 * sig
															+ 18750
																	* pow(mu,
																			4)
																	* (-163
																			+ 25
																					* sig)))
							+ 16 * pow(j, 4)
									* (207420014563 - 648467123500 * sig
											+ 1250 * pow(mu, 2)
													* (-280529426
															+ 170372400 * sig
															+ 625 * pow(mu, 2)
																	* (236621
																			+ 160410
																					* sig))))
					+ 500
							* (1431292409600 + 133494078080 * j
									+ 191250 * pow(j, 7)
											* pow(6 - 5 * pow(mu, 2), 2)
									- 46875 * pow(j, 8)
											* pow(-6 + 5 * pow(mu, 2), 3)
									- 192 * pow(j, 3)
											* (-52450914
													+ 63201125 * pow(mu, 2))
									- 64 * pow(j, 2)
											* (-998177094
													+ 2499919375 * pow(mu, 2))
									+ 8 * pow(j, 5)
											* (6389939 - 50814750 * pow(mu, 2)
													+ 32015625 * pow(mu, 4))
									+ 16 * pow(j, 4)
											* (-1296934247
													+ 425931000 * pow(mu, 2)
													+ 250640625 * pow(mu, 4))
									+ 100 * pow(j, 6)
											* (-7025094 + 12683995 * pow(mu, 2)
													- 4670625 * pow(mu, 4)
													+ 234375 * pow(mu, 6)))
							* sig
							* (3515625 * pow(j, 8)
									* pow(6 - 5 * pow(mu, 2), 4)
									* pow(sig, 2)
									- 19125000 * pow(j, 7)
											* pow(-6 + 5 * pow(mu, 2), 3)
											* pow(sig, 2)
									+ 192000 * j
											* (123471425 + 799978500 * sig
													+ 804303447 * pow(sig, 2))
									+ 160000
											* (5178027300
													+ sig
															* (14165823000
																	+ 28585055449
																			* sig))
									- 15000 * pow(j, 6)
											* pow(6 - 5 * pow(mu, 2), 2)
											* (-112500
													- sig
															* (187500
																	+ 461101
																			* sig)
													+ 3750 * pow(mu, 2)
															* (25
																	+ 213
																			* pow(
																					sig,
																					2)))
									+ 40800 * pow(j, 5)
											* (-6 + 5 * pow(mu, 2))
											* (-168750
													- sig
															* (187500
																	+ 938617
																			* sig)
													+ 5625 * pow(mu, 2)
															* (25
																	+ 213
																			* pow(
																					sig,
																					2)))
									- 1920 * pow(j, 3)
											* (-818705975
													+ 3 * sig
															* (-564287500
																	+ 1118339281
																			* sig)
													+ 625 * pow(mu, 2)
															* (1374450
																	+ sig
																			* (5377500
																					+ 7111457
																							* sig)))
									- 3200 * pow(j, 2)
											* (-18653709450
													+ sig
															* (41617980000
																	+ 76399756193
																			* sig)
													+ 3750 * pow(mu, 2)
															* (7263025
																	+ 3 * sig
																			* (3488500
																					+ 10373597
																							* sig)))
									+ 16 * pow(j, 4)
											* (3371453243689 * pow(sig, 2)
													- 75000
															* (-4021053
																	+ 19352665
																			* sig)
													+ 2343750 * pow(mu, 4)
															* (80850
																	+ sig
																			* (-7500
																					+ 418321
																							* sig))
													- 7500 * pow(mu, 2)
															* (48202525
																	+ 3 * sig
																			* (-47762500
																					+ 69264921
																							* sig))))
					+ (1171875 * pow(j, 8) * pow(-6 + 5 * pow(mu, 2), 3)
							* (-6 + 5 * pow(mu, 2) - 20 * sig)
							- 6375000 * pow(j, 7)
									* pow(6 - 5 * pow(mu, 2), 2)
									* (-6 + 5 * pow(mu, 2) - 15 * sig)
							+ 64000 * j * (212878097 + 1042922485 * sig)
							+ 160000 * (1950390283 + 4472788780 * sig)
							- 640 * pow(j, 3)
									* (-2379695207 - 7867637100 * sig
											+ 625 * pow(mu, 2)
													* (4022557 + 15168270 * sig))
							- 3200 * pow(j, 2)
									* (-13540744969 - 9981770940 * sig
											+ 1250 * pow(mu, 2)
													* (12522241 + 19999355 * sig))
							+ 800 * pow(j, 5)
									* (102751434 + 31949695 * sig
											+ 9375 * pow(mu, 4)
													* (8313 + 17075 * sig)
											- 5 * pow(mu, 2)
													* (33173239 + 50814750 * sig))
							+ 5000 * pow(j, 6)
									* (28749636 - 70250940 * sig
											+ 5 * pow(mu, 2)
													* (-12109212
															+ pow(mu, 2)
																	* (9765505
																			- 9341250
																					* sig)
															+ 25367990 * sig
															+ 18750
																	* pow(mu,
																			4)
																	* (-163
																			+ 25
																					* sig)))
							+ 16 * pow(j, 4)
									* (207420014563 - 648467123500 * sig
											+ 1250 * pow(mu, 2)
													* (-280529426
															+ 170372400 * sig
															+ 625 * pow(mu, 2)
																	* (236621
																			+ 160410
																					* sig))))
							* (3515625 * pow(j, 8)
									* pow(6 - 5 * pow(mu, 2), 4)
									* pow(sig, 2)
									- 19125000 * pow(j, 7)
											* pow(-6 + 5 * pow(mu, 2), 3)
											* pow(sig, 2)
									+ 192000 * j
											* (123471425 + 799978500 * sig
													+ 804303447 * pow(sig, 2))
									+ 160000
											* (5178027300
													+ sig
															* (14165823000
																	+ 28585055449
																			* sig))
									- 15000 * pow(j, 6)
											* pow(6 - 5 * pow(mu, 2), 2)
											* (-112500
													- sig
															* (187500
																	+ 461101
																			* sig)
													+ 3750 * pow(mu, 2)
															* (25
																	+ 213
																			* pow(
																					sig,
																					2)))
									+ 40800 * pow(j, 5)
											* (-6 + 5 * pow(mu, 2))
											* (-168750
													- sig
															* (187500
																	+ 938617
																			* sig)
													+ 5625 * pow(mu, 2)
															* (25
																	+ 213
																			* pow(
																					sig,
																					2)))
									- 1920 * pow(j, 3)
											* (-818705975
													+ 3 * sig
															* (-564287500
																	+ 1118339281
																			* sig)
													+ 625 * pow(mu, 2)
															* (1374450
																	+ sig
																			* (5377500
																					+ 7111457
																							* sig)))
									- 3200 * pow(j, 2)
											* (-18653709450
													+ sig
															* (41617980000
																	+ 76399756193
																			* sig)
													+ 3750 * pow(mu, 2)
															* (7263025
																	+ 3 * sig
																			* (3488500
																					+ 10373597
																							* sig)))
									+ 16 * pow(j, 4)
											* (3371453243689 * pow(sig, 2)
													- 75000
															* (-4021053
																	+ 19352665
																			* sig)
													+ 2343750 * pow(mu, 4)
															* (80850
																	+ sig
																			* (-7500
																					+ 418321
																							* sig))
													- 7500 * pow(mu, 2)
															* (48202525
																	+ 3 * sig
																			* (-47762500
																					+ 69264921
																							* sig))))))
			/ pow(
					3515625 * pow(j, 8) * pow(6 - 5 * pow(mu, 2), 4)
							* pow(sig, 2)
							- 19125000 * pow(j, 7)
									* pow(-6 + 5 * pow(mu, 2), 3)
									* pow(sig, 2)
							+ 192000 * j
									* (123471425 + 799978500 * sig
											+ 804303447 * pow(sig, 2))
							+ 160000
									* (5178027300
											+ sig
													* (14165823000
															+ 28585055449 * sig))
							- 15000 * pow(j, 6) * pow(6 - 5 * pow(mu, 2), 2)
									* (-112500 - sig * (187500 + 461101 * sig)
											+ 3750 * pow(mu, 2)
													* (25 + 213 * pow(sig, 2)))
							+ 40800 * pow(j, 5) * (-6 + 5 * pow(mu, 2))
									* (-168750 - sig * (187500 + 938617 * sig)
											+ 5625 * pow(mu, 2)
													* (25 + 213 * pow(sig, 2)))
							- 1920 * pow(j, 3)
									* (-818705975
											+ 3 * sig
													* (-564287500
															+ 1118339281 * sig)
											+ 625 * pow(mu, 2)
													* (1374450
															+ sig
																	* (5377500
																			+ 7111457
																					* sig)))
							- 3200 * pow(j, 2)
									* (-18653709450
											+ sig
													* (41617980000
															+ 76399756193 * sig)
											+ 3750 * pow(mu, 2)
													* (7263025
															+ 3 * sig
																	* (3488500
																			+ 10373597
																					* sig)))
							+ 16 * pow(j, 4)
									* (3371453243689 * pow(sig, 2)
											- 75000
													* (-4021053 + 19352665 * sig)
											+ 2343750 * pow(mu, 4)
													* (80850
															+ sig
																	* (-7500
																			+ 418321
																					* sig))
											- 7500 * pow(mu, 2)
													* (48202525
															+ 3 * sig
																	* (-47762500
																			+ 69264921
																					* sig))),
					2);
}

ftype A2(ftype sig, ftype mu, ftype j)
{ //dA/dmu
	return (3000 * pow(j, 2) * mu * sig
			* (-3
					* (-191250 * pow(j, 5) * pow(6 - 5 * pow(mu, 2), 2)
							* pow(sig, 2)
							+ 46875 * pow(j, 6) * pow(-6 + 5 * pow(mu, 2), 3)
									* pow(sig, 2)
							- 800 * j
									* (1374450 + 5377500 * sig
											+ 7111457 * pow(sig, 2))
							- 8000
									* (7263025 + 10465500 * sig
											+ 31120791 * pow(sig, 2))
							- 100 * pow(j, 4) * (-6 + 5 * pow(mu, 2))
									* (-168750 - 187500 * sig
											- 940351 * pow(sig, 2)
											+ 5625 * pow(mu, 2)
													* (25 + 213 * pow(sig, 2)))
							+ 136 * pow(j, 3)
									* (-337500 - 187500 * sig
											- 2376367 * pow(sig, 2)
											+ 11250 * pow(mu, 2)
													* (25 + 213 * pow(sig, 2)))
							+ 80 * pow(j, 2)
									* (-48202525 + 143287500 * sig
											- 207794763 * pow(sig, 2)
											+ 625 * pow(mu, 2)
													* (80850 - 7500 * sig
															+ 418321
																	* pow(sig,
																			2))))
					* (1171875 * pow(j, 8) * pow(-6 + 5 * pow(mu, 2), 3)
							* (-6 + 5 * pow(mu, 2) - 20 * sig)
							- 6375000 * pow(j, 7)
									* pow(6 - 5 * pow(mu, 2), 2)
									* (-6 + 5 * pow(mu, 2) - 15 * sig)
							+ 64000 * j * (212878097 + 1042922485 * sig)
							+ 160000 * (1950390283 + 4472788780 * sig)
							- 640 * pow(j, 3)
									* (-2379695207 - 7867637100 * sig
											+ 625 * pow(mu, 2)
													* (4022557 + 15168270 * sig))
							- 3200 * pow(j, 2)
									* (-13540744969 - 9981770940 * sig
											+ 1250 * pow(mu, 2)
													* (12522241 + 19999355 * sig))
							+ 800 * pow(j, 5)
									* (102751434 + 31949695 * sig
											+ 9375 * pow(mu, 4)
													* (8313 + 17075 * sig)
											- 5 * pow(mu, 2)
													* (33173239 + 50814750 * sig))
							+ 5000 * pow(j, 6)
									* (28749636 - 70250940 * sig
											+ 5 * pow(mu, 2)
													* (-12109212
															+ pow(mu, 2)
																	* (9765505
																			- 9341250
																					* sig)
															+ 25367990 * sig
															+ 18750
																	* pow(mu,
																			4)
																	* (-163
																			+ 25
																					* sig)))
							+ 16 * pow(j, 4)
									* (207420014563 - 648467123500 * sig
											+ 1250 * pow(mu, 2)
													* (-280529426
															+ 170372400 * sig
															+ 625 * pow(mu, 2)
																	* (236621
																			+ 160410
																					* sig))))
					+ (46875 * pow(j, 6) * pow(6 - 5 * pow(mu, 2), 2)
							* (-6 + 5 * pow(mu, 2) - 15 * sig)
							- 191250 * pow(j, 5) * (-6 + 5 * pow(mu, 2))
									* (-6 + 5 * pow(mu, 2) - 10 * sig)
							- 800 * j * (4022557 + 15168270 * sig)
							- 8000 * (12522241 + 19999355 * sig)
							+ 100 * pow(j, 4)
									* (-6054606
											+ pow(mu, 2)
													* (9765505 - 9341250 * sig)
											+ 12683995 * sig
											+ 28125 * pow(mu, 4)
													* (-163 + 25 * sig))
							+ 8 * pow(j, 3)
									* (-33173239 - 50814750 * sig
											+ 3750 * pow(mu, 2)
													* (8313 + 17075 * sig))
							+ 80 * pow(j, 2)
									* (-140264713 + 85186200 * sig
											+ 625 * pow(mu, 2)
													* (236621 + 160410 * sig)))
							* (3515625 * pow(j, 8)
									* pow(6 - 5 * pow(mu, 2), 4)
									* pow(sig, 2)
									- 19125000 * pow(j, 7)
											* pow(-6 + 5 * pow(mu, 2), 3)
											* pow(sig, 2)
									+ 192000 * j
											* (123471425 + 799978500 * sig
													+ 804303447 * pow(sig, 2))
									+ 160000
											* (5178027300
													+ sig
															* (14165823000
																	+ 28585055449
																			* sig))
									- 15000 * pow(j, 6)
											* pow(6 - 5 * pow(mu, 2), 2)
											* (-112500
													- sig
															* (187500
																	+ 461101
																			* sig)
													+ 3750 * pow(mu, 2)
															* (25
																	+ 213
																			* pow(
																					sig,
																					2)))
									+ 40800 * pow(j, 5)
											* (-6 + 5 * pow(mu, 2))
											* (-168750
													- sig
															* (187500
																	+ 938617
																			* sig)
													+ 5625 * pow(mu, 2)
															* (25
																	+ 213
																			* pow(
																					sig,
																					2)))
									- 1920 * pow(j, 3)
											* (-818705975
													+ 3 * sig
															* (-564287500
																	+ 1118339281
																			* sig)
													+ 625 * pow(mu, 2)
															* (1374450
																	+ sig
																			* (5377500
																					+ 7111457
																							* sig)))
									- 3200 * pow(j, 2)
											* (-18653709450
													+ sig
															* (41617980000
																	+ 76399756193
																			* sig)
													+ 3750 * pow(mu, 2)
															* (7263025
																	+ 3 * sig
																			* (3488500
																					+ 10373597
																							* sig)))
									+ 16 * pow(j, 4)
											* (3371453243689 * pow(sig, 2)
													- 75000
															* (-4021053
																	+ 19352665
																			* sig)
													+ 2343750 * pow(mu, 4)
															* (80850
																	+ sig
																			* (-7500
																					+ 418321
																							* sig))
													- 7500 * pow(mu, 2)
															* (48202525
																	+ 3 * sig
																			* (-47762500
																					+ 69264921
																							* sig))))))
			/ pow(
					3515625 * pow(j, 8) * pow(6 - 5 * pow(mu, 2), 4)
							* pow(sig, 2)
							- 19125000 * pow(j, 7)
									* pow(-6 + 5 * pow(mu, 2), 3)
									* pow(sig, 2)
							+ 192000 * j
									* (123471425 + 799978500 * sig
											+ 804303447 * pow(sig, 2))
							+ 160000
									* (5178027300
											+ sig
													* (14165823000
															+ 28585055449 * sig))
							- 15000 * pow(j, 6) * pow(6 - 5 * pow(mu, 2), 2)
									* (-112500 - sig * (187500 + 461101 * sig)
											+ 3750 * pow(mu, 2)
													* (25 + 213 * pow(sig, 2)))
							+ 40800 * pow(j, 5) * (-6 + 5 * pow(mu, 2))
									* (-168750 - sig * (187500 + 938617 * sig)
											+ 5625 * pow(mu, 2)
													* (25 + 213 * pow(sig, 2)))
							- 1920 * pow(j, 3)
									* (-818705975
											+ 3 * sig
													* (-564287500
															+ 1118339281 * sig)
											+ 625 * pow(mu, 2)
													* (1374450
															+ sig
																	* (5377500
																			+ 7111457
																					* sig)))
							- 3200 * pow(j, 2)
									* (-18653709450
											+ sig
													* (41617980000
															+ 76399756193 * sig)
											+ 3750 * pow(mu, 2)
													* (7263025
															+ 3 * sig
																	* (3488500
																			+ 10373597
																					* sig)))
							+ 16 * pow(j, 4)
									* (3371453243689 * pow(sig, 2)
											- 75000
													* (-4021053 + 19352665 * sig)
											+ 2343750 * pow(mu, 4)
													* (80850
															+ sig
																	* (-7500
																			+ 418321
																					* sig))
											- 7500 * pow(mu, 2)
													* (48202525
															+ 3 * sig
																	* (-47762500
																			+ 69264921
																					* sig))),
					2);
}

void BondiMass(parameters par, ftype *X)
{

	ftype *F, *LVmm_Scri, *H, *Hd2, H0, **phi, **phi3, *dI, *Int, dMB, MB, J,
			J_IB, J_NB, dJ, j_B, Area, eps, Omega, Delta_sigma,
			Penrose_inequality;
	derivs_2D v;
	int ns[2], nt = par.nt, i, j, iaux;
	ns[0] = par.ns[0];
	ns[1] = par.ns[1];
	ftype ds[2], fac[2], mu, phi_j, Add31_j, Add32_j, Add33_j, aux, B, epsm1;
	ds[0] = par.sig0;
	fac[0] = 1 / ds[0];
	ds[1] = par.sig0 - par.sig1;
	fac[1] = 1 / ds[1];
	ofstream fp;
	FILE *fp_test;

	allocate_derivs_2D(&v, par.num_v);
	fill0_derivs_2D(v, par.num_v);
	F = dvector(0, par.ntotal - 1);
	fill0_dvector(F, 0, par.ntotal);
	F_of_X(par, X, v, F);

	H = dvector(0, nt);
	Hd2 = dvector(0, nt);
	phi = dmatrix(0, ns[0], 0, nt);
	phi3 = dmatrix(0, ns[0], 0, nt);
	LVmm_Scri = dvector(0, nt);

	//get the apparent Horizon H(mu)
	//GetApparentHorizon(par,H,Hd2);
	for (j = 0; j < nt; j++)
	{
		H[j] = par.sig1;
		Hd2[j] = 0;
	}
	//calculate Bondi-Mass
	dI = dvector(0, nt);
	Int = dvector(0, nt);

	LVmm(par, v, LVmm_Scri);

	for (i = 0; i < ns[0]; i++)
	{
		for (j = 0; j < nt; j++)
		{
			phi[i][j] = v.d0[Index(par, 0, 0, i, j)];
		}
	}

	getphi3(phi, phi3, ns[0], nt, ds[0]);

	for (j = 0; j < nt; j++)
	{
		dI[j] = sqrtdetqc_Scri()
				* (0.125 * sqrt(gammac110_Scri(par.j_BH)) * LVmm_Scri[j]
						+ 8 * gammac110_Scri(par.j_BH) * phi3[0][j]);
	}
	chebft_Extremes(dI, nt, 0);
	chint(-1, 1, dI, Int, nt);
	dMB = 0.5 * chebev(-1, 1, Int, nt, 1);
	MB = 1. + dMB;
	free_dvector(Int, 0, nt);
	free_dvector(dI, 0, nt);

	//calculate the Area
	Int = dvector(0, nt);
	dI = dvector(0, nt);
	for (j = 0; j < nt; j++)
	{
		mu = -cos(Pi * j / ((ftype) nt - 1));
		//phi_j=chebevxy(par.sig0, par.sig1,-1, 1,phi, ns,nt,H[j],mu);
		phi_j = v.d0[Index(par, 1, 0, 0, j)];
		Omega = H[j] / sqr(phi_j);
		aux = gammacuuij(par, H[j], mu, 1, 1)
				- 2 * Hd2[j] * gammacuuij(par, H[j], mu, 1, 2)
				+ sqr(Hd2[j]) * gammacuuij(par, H[j], mu, 2, 2);
		B = 1 / sqrt(aux);
		dI[j] = sqrtdetqc(par, H[j], mu, B, Hd2[j]) / sqr(Omega);
	}
	chebft_Extremes(dI, nt, 0);
	chint(-1, 1, dI, Int, nt);
	Area = 2 * Pi * chebev(-1, 1, Int, nt, 1);
	free_dvector(Int, 0, nt);
	free_dvector(dI, 0, nt);

	//calculate the angular momentum at the inner boundary

	J = CalculateAngularMomentum(par, v, H, Hd2);

	//calculate j_B

	j_B = J / sqr(MB);

	//calculate deviation

	Delta_sigma = par.sig1 - 2 / (1 + sqrt(1 - sqr(par.j_BH)));

	//calculate eps
	eps = 0.125 * Area / (sqr(MB) * Pi * (1 + sqrt(1 - sqr(J / (sqr(MB))))));
	epsm1 = 1 - eps;

	//calculate penrose inequality

	Penrose_inequality = Area / (16 * Pi * sqr(MB));

	ftype JoA, AoMB;

	JoA = 8 * Pi * J / (Area);

	AoMB = Area / sqr(MB);

	cout << "dMB="<<dMB<<endl;
	fp_test = fopen("./BondiMass/MB","a");
	size_t size;
	//check whether the file is empty or not
	iaux = fseek(fp_test, 0, SEEK_END); //go to the end of the file
	size = ftell(fp_test); //measure the filesize
	fclose(fp_test);
	fp.precision(PREC_OUTPUT);
	fp.open("./BondiMass/MB",ios_base::app);
	if (size == 0)
	{
		fp<<"#1:j_ID \t 2:sig1 3:Delta_sigma \t 4:MB \t 5:J \t 6:Area \t 7:j_B \t 8:eps \t 9: Penrose_inequality(A/(16*Pi*MB^2)) \t 10:8*pi*J/A \t 11: A/MB^2 \t12:1-eps "<<endl;
	}
	fp<<par.j_BH << "\t" << par.sig1 << "\t" << Delta_sigma << "\t" << MB << "\t" << "\t" << J << "\t" << Area << "\t" << j_B << "\t" << eps << "\t" <<	Penrose_inequality << "\t" << JoA << "\t" << AoMB << "\t" << epsm1 <<endl;
	fp.close();

	free_dvector(H, 0, nt);
	free_dvector(Hd2, 0, nt);
	free_dmatrix(phi, 0, ns[0], 0, nt);
	free_dmatrix(phi3, 0, ns[0], 0, nt);
	free_dvector(LVmm_Scri, 0, nt);
	free_derivs_2D(&v, par.num_v);
	free_dvector(F, 0, par.ntotal - 1);
}

ftype CalculateAngularMomentum(parameters par, derivs_2D v, ftype *H,
		ftype *Hd2)
{
	int nt = par.nt, ns = par.ns[1], j, i;
	ftype *Int, *dI, **phi, J, phi_j, Add31_j, Add32_j, Add33_j, Omega, aux, mu,
			B, **Add31, **Add32, **Add33;

	Int = dvector(0, nt);
	dI = dvector(0, nt);
	phi = dmatrix(0, ns, 0, nt);
	Add31 = dmatrix(0, ns, 0, nt);
	Add32 = dmatrix(0, ns, 0, nt);
	Add33 = dmatrix(0, ns, 0, nt);

	GetAdd3j(par, v, Add31, Add32, Add33);

	chebft_Extremes_2D(Add31, ns, nt);
	chebft_Extremes_2D(Add32, ns, nt);
	chebft_Extremes_2D(Add33, ns, nt);

	for (j = 0; j < nt; j++)
	{
		mu = -cos(Pi * j / ((ftype) nt - 1));
		//phi_j=chebevxy(par.sig0, par.sig1,-1, 1,phi, ns,nt,H[j],mu);
		phi_j = v.d0[Index(par, 1, 0, 0, j)];

		Add31_j = chebevxy(par.sig1, par.sig0, -1, 1, Add31, ns, nt, H[j], mu);
		Add32_j = chebevxy(par.sig1, par.sig0, -1, 1, Add32, ns, nt, H[j], mu);
		Add33_j = chebevxy(par.sig1, par.sig0, -1, 1, Add33, ns, nt, H[j], mu)
				* (1 - sqr(mu));
		Omega = H[j] / sqr(phi_j);
		aux = gammacuuij(par, H[j], mu, 1, 1)
				- 2 * Hd2[j] * gammacuuij(par, H[j], mu, 1, 2)
				+ sqr(Hd2[j]) * gammacuuij(par, H[j], mu, 2, 2);
		B = 1 / sqrt(aux);
		dI[j] = sqrtdetqc(par, H[j], mu, B, Hd2[j]) * B
				* (Add31_j
						* (gammacuuij(par, H[j], mu, 1, 1)
								- gammacuuij(par, H[j], mu, 2, 1) * Hd2[j])
						+ Add32_j
								* (gammacuuij(par, H[j], mu, 1, 2)
										- gammacuuij(par, H[j], mu, 2, 2)
												* Hd2[j])
						+ Add33_j
								* (gammacuuij(par, H[j], mu, 1, 3)
										- gammacuuij(par, H[j], mu, 2, 3)
												* Hd2[j]));
	}
	chebft_Extremes(dI, nt, 0);
	chint(-1, 1, dI, Int, nt);
	J = 0.25 * chebev(-1, 1, Int, nt, 1);
	cout <<"J="<<J<<endl;

	free_dvector(Int, 0, nt);
	free_dvector(dI, 0, nt);
	free_dmatrix(phi, 0, ns, 0, nt);
	free_dmatrix(Add31, 0, ns, 0, nt);
	free_dmatrix(Add32, 0, ns, 0, nt);
	free_dmatrix(Add33, 0, ns, 0, nt);

	return J;
}
