#include "main.h"


void Derivatives_st(parameters par, derivs_2D v)
{
	int idom, ipot, ns, nt, i, j, indx[NMAX];
	ftype p[NMAX], dp[NMAX], d2p[NMAX], q[NMAX], dq[NMAX];

	for (idom = 0; idom < NDOM; idom++)
	{

		ns = par.ns[idom];
		nt = par.nt;

		for (ipot = 0; ipot < NPOT; ipot++)
		{

			for (j = 0; j < nt; j++)
			{		// Calculation of Derivatives w.r.t. s-Direction
				for (i = 0; i < ns; i++)
				{	// (Chebyshev_Extremes)
					indx[i] = Index(par, idom, ipot, i, j);
					p[i] = v.d0[indx[i]];
				}
				chebft_Extremes(p, ns, 0);
				chder(0., 1., p, dp, ns);
				chder(0., 1., dp, d2p, ns);
				chebft_Extremes(dp, ns, 1);
				chebft_Extremes(d2p, ns, 1);
				for (i = 0; i < ns; i++)
				{
					v.d1[indx[i]] = dp[i];
					v.d11[indx[i]] = d2p[i];
				}
			}

			for (i = 0; i < ns; i++)
			{		// Calculation of Derivatives w.r.t. t-Direction
				for (j = 0; j < nt; j++)
				{	// (Chebyshev_Extremes)
					indx[j] = Index(par, idom, ipot, i, j);
					p[j] = v.d0[indx[j]];
					q[j] = v.d1[indx[j]];
				}
				chebft_Extremes(p, nt, 0);
				chebft_Extremes(q, nt, 0);
				chder(-1., 1., p, dp, nt);
				chder(-1., 1., dp, d2p, nt);
				chder(-1., 1., q, dq, nt);
				chebft_Extremes(dp, nt, 1);
				chebft_Extremes(d2p, nt, 1);
				chebft_Extremes(dq, nt, 1);
				for (j = 0; j < nt; j++)
				{
					v.d2[indx[j]] = dp[j];
					v.d22[indx[j]] = d2p[j];
					v.d12[indx[j]] = dq[j];
				}
			}

		}
	}

}

void Get_v_From_X(parameters par, ftype *X, derivs_2D v)
{
	int n_2D = par.n_2D, n, n_v, idom, ipot, i, j;

	for (n = 0; n < n_2D; n++)
	{
		n_v = par.n_v_of_n[n];
		Get_Indices_From_n_v(par, n_v, &idom, &ipot, &i, &j);
		v.d0[n_v] = X[n];
		if (idom == 0 && i == par.ns[0] - 1)
			v.d0[Index(par, 1, ipot, par.ns[1] - 1, j)] = X[n];
	}
	Derivatives_st(par, v);
}

void F_of_X(parameters par, ftype *X, derivs_2D v, ftype *F)
{
	int n_2D = par.n_2D, num_v = par.num_v, idom, i, j, ns, Ns, nt, Nt, n;
	ftype *FF;

	FF = new ftype[num_v];
	fill0_dvector(FF, 0, num_v);

	Get_v_From_X(par, X, v);

	for (idom = 0; idom < NDOM; idom++)
	{
		ns = par.ns[idom];
		Ns = ns - 1;
		nt = par.nt;
		Nt = nt - 1;
		for (j = 0; j < nt; j++)
		{
			for (i = 0; i < ns; i++)
			{
				NonLinFieldEqns(par, idom, i, j, v, FF);
			}
		}
	}

// reordering
	for (n = 0; n < n_2D; n++)
		F[n] = FF[par.n_v_of_n[n]];

	delete[] FF;
}

void J_times_DX(parameters par, ftype *X, ftype *DX, derivs_2D v, ftype *JDX)
{
	int n_2D = par.n_2D, num_v = par.num_v, idom, i, j, ns, Ns, nt, Nt, n;
	ftype *JJDX;
	derivs_2D Dv;

	allocate_derivs_2D(&Dv, num_v);
	JJDX = new ftype[num_v];
	fill0_dvector(JJDX, 0, num_v);

	Get_v_From_X(par, DX, Dv);

	for (idom = 0; idom < NDOM; idom++)
	{
		ns = par.ns[idom];
		Ns = ns - 1;
		nt = par.nt;
		Nt = nt - 1;
		for (j = 0; j < nt; j++)
		{
			for (i = 0; i < ns; i++)
			{
				LinFieldEqns(par, idom, i, j, v, Dv, JJDX);
			}
		}
	}

// reordering
	for (n = 0; n < n_2D; n++)
		JDX[n] = JJDX[par.n_v_of_n[n]];

	free_derivs_2D(&Dv, num_v);
	delete[] JJDX;
}
