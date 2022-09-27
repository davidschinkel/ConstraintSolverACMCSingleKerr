#include "main.h"

#define iH4 1

void func_AijAij(parameters par, int idom, int i, int j, derivs_2D v,
		ftype fac1, ftype *AijAij)
{
	int ii, jj, a;

	*AijAij = 0.;
	for (ii = 1; ii <= 3; ii++)
	{
		for (jj = 1; jj <= 3; jj++)
		{
			ftype Add_ij = par.mdd[indx_2(par, ii, jj, idom, i, j)], Auu_ij =
					par.muu[indx_2(par, ii, jj, idom, i, j)];

			for (a = 1; a <= 3; a++)
			{
				int indx = Index(par, idom, a, i, j);

				Add_ij += par.kdd[indx_4(par, ii, jj, 1, a, idom, i, j)]
						* v.d1[indx] * fac1
						+ par.kdd[indx_4(par, ii, jj, 2, a, idom, i, j)]
								* v.d2[indx]
						+ par.ldd[indx_3(par, ii, jj, a, idom, i, j)]
								* v.d0[indx];

				Auu_ij += par.kuu[indx_4(par, ii, jj, 1, a, idom, i, j)]
						* v.d1[indx] * fac1
						+ par.kuu[indx_4(par, ii, jj, 2, a, idom, i, j)]
								* v.d2[indx]
						+ par.luu[indx_3(par, ii, jj, a, idom, i, j)]
								* v.d0[indx];
			}
			*AijAij += Add_ij * Auu_ij;
		}
	}
}

void func_DAijAij(parameters par, int idom, int i, int j, derivs_2D v,
		derivs_2D Dv, ftype fac1, ftype *AijAij, ftype *DAijAij)
{
	int ii, jj, a;

	*AijAij = *DAijAij = 0.;
	for (ii = 1; ii <= 3; ii++)
	{
		for (jj = 1; jj <= 3; jj++)
		{
			ftype DAdd_ij = 0., DAuu_ij = 0., Add_ij = par.mdd[indx_2(par, ii,
					jj, idom, i, j)], Auu_ij = par.muu[indx_2(par, ii, jj, idom,
					i, j)];

			for (a = 1; a <= 3; a++)
			{
				int indx = Index(par, idom, a, i, j);

				Add_ij += par.kdd[indx_4(par, ii, jj, 1, a, idom, i, j)]
						* v.d1[indx] * fac1
						+ par.kdd[indx_4(par, ii, jj, 2, a, idom, i, j)]
								* v.d2[indx]
						+ par.ldd[indx_3(par, ii, jj, a, idom, i, j)]
								* v.d0[indx];

				Auu_ij += par.kuu[indx_4(par, ii, jj, 1, a, idom, i, j)]
						* v.d1[indx] * fac1
						+ par.kuu[indx_4(par, ii, jj, 2, a, idom, i, j)]
								* v.d2[indx]
						+ par.luu[indx_3(par, ii, jj, a, idom, i, j)]
								* v.d0[indx];

				DAdd_ij += par.kdd[indx_4(par, ii, jj, 1, a, idom, i, j)]
						* Dv.d1[indx] * fac1
						+ par.kdd[indx_4(par, ii, jj, 2, a, idom, i, j)]
								* Dv.d2[indx]
						+ par.ldd[indx_3(par, ii, jj, a, idom, i, j)]
								* Dv.d0[indx];

				DAuu_ij += par.kuu[indx_4(par, ii, jj, 1, a, idom, i, j)]
						* Dv.d1[indx] * fac1
						+ par.kuu[indx_4(par, ii, jj, 2, a, idom, i, j)]
								* Dv.d2[indx]
						+ par.luu[indx_3(par, ii, jj, a, idom, i, j)]
								* Dv.d0[indx];
			}
			*AijAij += Add_ij * Auu_ij;
			*DAijAij += DAdd_ij * Auu_ij + Add_ij * DAuu_ij;
		}
	}
}

void func_h3(parameters par, int idom, int i, int j, derivs_2D v, ftype fac1,
		ftype sin2th, ftype *h3)
{
	int ii, jj, a;

	*h3 = 0.;
	for (ii = 1; ii <= 3; ii++)
	{
		for (jj = 1; jj <= 3; jj++)
		{
			ftype Add_ij = par.mdd[indx_2(par, ii, jj, idom, i, j)], sisj;

			for (a = 1; a <= 3; a++)
			{
				int indx = Index(par, idom, a, i, j);

				Add_ij += par.kdd[indx_4(par, ii, jj, 1, a, idom, i, j)]
						* v.d1[indx] * fac1
						+ par.kdd[indx_4(par, ii, jj, 2, a, idom, i, j)]
								* v.d2[indx]
						+ par.ldd[indx_3(par, ii, jj, a, idom, i, j)]
								* v.d0[indx];
			}
			sisj = par.s[indx_Horizon(par, ii, j)]
					* par.s[indx_Horizon(par, jj, j)];
			if (ii == 2 || jj == 2 || (ii == 3 && jj == 3))
				*h3 += Add_ij * sisj * sin2th;
			else
				*h3 += Add_ij * sisj;
		}
	}
}

void func_Dh3(parameters par, int idom, int i, int j, derivs_2D v, derivs_2D Dv,
		ftype fac1, ftype sin2th, ftype *h3, ftype *Dh3)
{
	int ii, jj, a;

	*h3 = *Dh3 = 0.;
	for (ii = 1; ii <= 3; ii++)
	{
		for (jj = 1; jj <= 3; jj++)
		{
			ftype Add_ij = par.mdd[indx_2(par, ii, jj, idom, i, j)], DAdd_ij =
					0., sisj;

			for (a = 1; a <= 3; a++)
			{
				int indx = Index(par, idom, a, i, j);

				Add_ij += par.kdd[indx_4(par, ii, jj, 1, a, idom, i, j)]
						* v.d1[indx] * fac1
						+ par.kdd[indx_4(par, ii, jj, 2, a, idom, i, j)]
								* v.d2[indx]
						+ par.ldd[indx_3(par, ii, jj, a, idom, i, j)]
								* v.d0[indx];
				DAdd_ij += par.kdd[indx_4(par, ii, jj, 1, a, idom, i, j)]
						* Dv.d1[indx] * fac1
						+ par.kdd[indx_4(par, ii, jj, 2, a, idom, i, j)]
								* Dv.d2[indx]
						+ par.ldd[indx_3(par, ii, jj, a, idom, i, j)]
								* Dv.d0[indx];
			}
			sisj = par.s[indx_Horizon(par, ii, j)]
					* par.s[indx_Horizon(par, jj, j)];
			if (ii == 2 || jj == 2 || (ii == 3 && jj == 3))
			{
				*h3 += Add_ij * sisj * sin2th;
				*Dh3 += DAdd_ij * sisj * sin2th;
			}
			else
			{
				*h3 += Add_ij * sisj;
				*Dh3 += DAdd_ij * sisj;
			}
		}
	}
}

void NonLinFieldEqns(parameters par, int idom, int i, int j, derivs_2D v,
		ftype *F)
{	// Nonlinear field equations

	int ns = par.ns[idom], Ns = ns - 1, nt = par.nt, Nt = nt - 1, ipot, i1, i4,
			indx0 = Index(par, idom, 0, i, j);
	ftype sig0 = par.sig0, sig1 = par.sig1, sig, fac1, fac2, H[Hmax + 1], s =
			sqr(sin(Pih * i / Ns));

	if (idom == 0)
	{
		sig = sig0 * s;
		fac1 = 1. / sig0;
	}            // fac1 = 1/D[sig,s]
	else
	{
		sig = sig1 + s * (sig0 - sig1);
		fac1 = 1. / (sig0 - sig1);
	}
	fac2 = sqr(fac1);

	ftype phi = v.d0[indx0], phid1 = v.d1[indx0] * fac1, phid11 = v.d11[indx0]
			* fac2, phid2 = v.d2[indx0], phid22 = v.d22[indx0], phid12 =
			v.d12[indx0] * fac1, phi2 = sqr(phi), phi5 = sqr(phi2) * phi, phi6 =
			phi5 * phi, phim7 = 1. / (phi5 * phi2);

	if (i > 0 && i < Ns)
	{ // Field equations
		ftype sig6 = sqr(sig * sig * sig), AijAij, H8;

		for (i1 = 1; i1 <= Hmax; i1++)
			H[i1] = par.H[indx_H(par, i1, idom, i, j)];

		func_AijAij(par, idom, i, j, v, fac1, &AijAij);
		H8 = sig6 * AijAij;
// 			H8 = 0.;

// 			F[indx0] = H[1]*phid11 + H[2]*phid22 + H[3]*phid12 + H[4]*phid1 + H[5]*phid2 + H[6]*phi + H[7]*phi5 + H8*phim7;
		F[indx0] = H[1] * phid11 + H[2] * phid22 + H[3] * phid12
				+ H[4] * phid1 * iH4 + H[5] * phid2 + H[6] * phi + H[7] * phi5
				+ H8 * phim7;

		for (i1 = 1; i1 <= 3; i1++)
		{
			F[Index(par, idom, i1, i, j)] = par.d[indx_1(par, i1, idom, i, j)]
					+ par.e[indx_1(par, i1, idom, i, j)] * phi6;
			for (i4 = 1; i4 <= 3; i4++)
			{
				F[Index(par, idom, i1, i, j)] += par.a[indx_4(par, i1, 1, 1, i4,
						idom, i, j)] * v.d11[Index(par, idom, i4, i, j)] * fac2
						+ (par.a[indx_4(par, i1, 1, 2, i4, idom, i, j)]
								+ par.a[indx_4(par, i1, 2, 1, i4, idom, i, j)])
								* v.d12[Index(par, idom, i4, i, j)] * fac1
						+ par.a[indx_4(par, i1, 2, 2, i4, idom, i, j)]
								* v.d22[Index(par, idom, i4, i, j)]

						+ par.b[indx_3(par, i1, 1, i4, idom, i, j)]
								* v.d1[Index(par, idom, i4, i, j)] * fac1
						+ par.b[indx_3(par, i1, 2, i4, idom, i, j)]
								* v.d2[Index(par, idom, i4, i, j)]

						+ par.c[indx_2(par, i1, i4, idom, i, j)]
								* v.d0[Index(par, idom, i4, i, j)];
			}
		}
	}
	else
	{
		if (i == 0)
		{ // Dirichlet-Boundary values for the momentum constraints

			for (i1 = 1; i1 <= 3; i1++)
				F[Index(par, idom, i1, i, j)] =
						v.d0[Index(par, idom, i1, i, j)];

			if (idom == 1)
			{ // Apparent Horizon Boundary condition
				ftype phi3 = phi2 * phi, phim3 = 1. / phi3, h[3], h3, sig2 =
						sqr(sig), mu = -cos(Pi * j / Nt), sin2th = 1.
						- sqr(mu);

				for (i1 = 1; i1 <= 2; i1++)
					h[i1] = par.h[indx_Horizon(par, i1, j)];

				func_h3(par, idom, i, j, v, fac1, sin2th, &h3);
				//h3 *= 0.25*sig2;
				h3 *= -0.25 * sig2;

				F[indx0] = h[1] * phi + h[2] * phi3 + h3 * phim3
						+ par.s[indx_Horizon(par, 1, j)] * phid1
						+ par.s[indx_Horizon(par, 2, j)] * sin2th * phid2;

			}
			else
			{ // Require Hamiltonian Constraint (Yamabe type field equation) at Scri+ (sig = 0)

				for (i1 = 1; i1 <= Hmax; i1++)
					H[i1] = par.H[indx_H(par, i1, idom, i, j)];

//  						F[indx0] = H[1]*phid11 + H[2]*phid22 + H[3]*phid12 + H[4]*phid1 + H[5]*phid2 + H[6]*phi + H[7]*phi5;
				F[indx0] = H[1] * phid11 + H[2] * phid22 + H[3] * phid12
						+ H[4] * phid1 * iH4 + H[5] * phid2 + H[6] * phi
						+ H[7] * phi5;
			}
		}
		else
		{  // transition conditions
			if (idom == 0)
			{
				ftype fac1 = 1. / sig0;           // fac1 = 1/D[sig,s]
				for (ipot = 0; ipot < NPOT; ipot++)
				{
					int indx = Index(par, idom, ipot, i, j);
					F[indx] += v.d1[indx] * fac1;
				}
			}
			else
			{
				int NS_0 = par.ns[0] - 1;       // NS_0 = Ns[idom=0]
				ftype fac1 = 1. / (sig0 - sig1); // fac1 = 1/D[sig,s]
				for (ipot = 0; ipot < NPOT; ipot++)
				{
					int indx0 = Index(par, 0, ipot, NS_0, j), indx1 = Index(par,
							1, ipot, i, j);
					F[indx0] -= v.d1[indx1] * fac1;
				}
			}
		}
	}
}

void LinFieldEqns(parameters par, int idom, int i, int j, derivs_2D v,
		derivs_2D Dv, ftype *JDX)
{	// Linear version of 'Nonlinear field equations'

	int ns = par.ns[idom], Ns = ns - 1, nt = par.nt, Nt = nt - 1, ipot, i1, i4,
			indx0 = Index(par, idom, 0, i, j);
	ftype sig0 = par.sig0, sig1 = par.sig1, sig, fac1, fac2, H[9], s = sqr(
			sin(Pih * i / Ns));

	if (idom == 0)
	{
		sig = sig0 * s;
		fac1 = 1. / sig0;
	}            // fac1 = 1/D[sig,s]
	else
	{
		sig = sig1 + s * (sig0 - sig1);
		fac1 = 1. / (sig0 - sig1);
	}
	fac2 = sqr(fac1);

	ftype phi = v.d0[indx0], phi2 = sqr(phi), phi4 = sqr(phi2), phi5 = phi4
			* phi, phim4 = 1. / phi4, phim7 = 1. / (phi5 * phi2), phim8 = sqr(
			phim4),

	Dphi = Dv.d0[indx0], Dphid1 = Dv.d1[indx0] * fac1, Dphid11 = Dv.d11[indx0]
			* fac2, Dphid2 = Dv.d2[indx0], Dphid22 = Dv.d22[indx0], Dphid12 =
			Dv.d12[indx0] * fac1,

	Dphi5 = 5. * phi4 * Dphi, Dphi6 = 6. * phi5 * Dphi, Dphim7 = -7. * phim8
			* Dphi;

	if (i > 0 && i < Ns)
	{ // Field equations
		ftype sig6 = sqr(sig * sig * sig), AijAij, DAijAij, H8, DH8;

		for (i1 = 1; i1 <= Hmax; i1++)
			H[i1] = par.H[indx_H(par, i1, idom, i, j)];

		func_DAijAij(par, idom, i, j, v, Dv, fac1, &AijAij, &DAijAij);
		H8 = sig6 * AijAij;
		DH8 = sig6 * DAijAij;
// 			H8 = DH8 = 0.;

		/*			JDX[indx0] = H[1]*Dphid11 + H[2]*Dphid22 + H[3]*Dphid12 + H[4]*Dphid1
		 + H[5]*Dphid2 + H[6]*Dphi + H[7]*Dphi5 + DH8*phim7 + H8*Dphim7;*/
		JDX[indx0] = H[1] * Dphid11 + H[2] * Dphid22 + H[3] * Dphid12
				+ H[4] * Dphid1 * iH4 + H[5] * Dphid2 + H[6] * Dphi
				+ H[7] * Dphi5 + DH8 * phim7 + H8 * Dphim7;

		for (i1 = 1; i1 <= 3; i1++)
		{
			JDX[Index(par, idom, i1, i, j)] = par.e[indx_1(par, i1, idom, i, j)]
					* Dphi6;
			for (i4 = 1; i4 <= 3; i4++)
			{
				JDX[Index(par, idom, i1, i, j)] += par.a[indx_4(par, i1, 1, 1,
						i4, idom, i, j)] * Dv.d11[Index(par, idom, i4, i, j)]
						* fac2
						+ (par.a[indx_4(par, i1, 1, 2, i4, idom, i, j)]
								+ par.a[indx_4(par, i1, 2, 1, i4, idom, i, j)])
								* Dv.d12[Index(par, idom, i4, i, j)] * fac1
						+ par.a[indx_4(par, i1, 2, 2, i4, idom, i, j)]
								* Dv.d22[Index(par, idom, i4, i, j)]

						+ par.b[indx_3(par, i1, 1, i4, idom, i, j)]
								* Dv.d1[Index(par, idom, i4, i, j)] * fac1
						+ par.b[indx_3(par, i1, 2, i4, idom, i, j)]
								* Dv.d2[Index(par, idom, i4, i, j)]

						+ par.c[indx_2(par, i1, i4, idom, i, j)]
								* Dv.d0[Index(par, idom, i4, i, j)];
			}
		}
	}
	else
	{
		if (i == 0)
		{ // Dirichlet-Boundary values for the momentum constraints

			for (i1 = 1; i1 <= 3; i1++)
				JDX[Index(par, idom, i1, i, j)] = Dv.d0[Index(par, idom, i1, i,
						j)];

			if (idom == 1)
			{ // Apparent Horizon Boundary condition
				ftype Dphi3 = 3. * phi2 * Dphi, phim3 = 1. / (phi2 * phi),
						Dphim3 = -3. * phim4 * Dphi, h[3], h3, Dh3, sig2 = sqr(
								sig), mu = -cos(Pi * j / Nt), sin2th = 1.
								- sqr(mu);

				for (i1 = 1; i1 <= 2; i1++)
					h[i1] = par.h[indx_Horizon(par, i1, j)];

				func_Dh3(par, idom, i, j, v, Dv, fac1, sin2th, &h3, &Dh3);
				//h3  *= 0.25*sig2;
				//Dh3 *= 0.25*sig2;

				h3 *= -0.25 * sig2;
				Dh3 *= -0.25 * sig2;

				JDX[indx0] = h[1] * Dphi + h[2] * Dphi3 + Dh3 * phim3
						+ h3 * Dphim3 + par.s[indx_Horizon(par, 1, j)] * Dphid1
						+ par.s[indx_Horizon(par, 2, j)] * sin2th * Dphid2;
			}
			else
			{ // Require Hamiltonian Constraint (Yamabe type field equation) at Scri+ (sig = 0)

				for (i1 = 1; i1 <= Hmax; i1++)
					H[i1] = par.H[indx_H(par, i1, idom, i, j)];

// 						JDX[indx0] = H[1]*Dphid11 + H[2]*Dphid22 + H[3]*Dphid12 + H[4]*Dphid1 + H[5]*Dphid2 + H[6]*Dphi + H[7]*Dphi5;
				JDX[indx0] = H[1] * Dphid11 + H[2] * Dphid22 + H[3] * Dphid12
						+ H[4] * Dphid1 * iH4 + H[5] * Dphid2 + H[6] * Dphi
						+ H[7] * Dphi5;
			}
		}
		else
		{  // transition conditions
			if (idom == 0)
			{
				ftype fac1 = 1. / sig0;           // fac1 = 1/D[sig,s]
				for (ipot = 0; ipot < NPOT; ipot++)
				{
					int indx = Index(par, idom, ipot, i, j);
					JDX[indx] += Dv.d1[indx] * fac1;
				}
			}
			else
			{
				int NS_0 = par.ns[0] - 1;          // NS_0 = Ns[idom=0]
				ftype fac1 = 1. / (sig0 - sig1);    // fac1 = 1/D[sig,s]
				for (ipot = 0; ipot < NPOT; ipot++)
				{
					int indx0 = Index(par, 0, ipot, NS_0, j), indx1 = Index(par,
							1, ipot, i, j);
					JDX[indx0] -= Dv.d1[indx1] * fac1;
				}
			}
		}
	}
}

#undef iH4

