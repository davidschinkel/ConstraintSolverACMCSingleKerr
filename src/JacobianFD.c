#include "main.h"


// 2nd and 4th order finite differencing


void GetFD_Dv_a(parameters par, derivs_2D Dv, int idom, int ipot, int i, int j)
{	// Derivatives of (Dv) with respect to (a):
	int i00 = Index(par, idom, ipot, i, j), ns = par.ns[idom], Ns = ns - 1,
			ip1 = Index(par, idom, ipot, i + 1, j), ip2 = Index(par, idom, ipot,
					i + 2, j), im1 = Index(par, idom, ipot, i - 1, j), im2 =
					Index(par, idom, ipot, i - 2, j);
	ftype ha = Pih / Ns, ga = 1. / ha, ga2 = ga * ga;// ha: Stepsize with respect to (a)

	if (i * (Ns - i) != 0)
	{ // here: Dv.d1  = Dv_a,  Dv.d11 = Dv_(aa)
	  // Dv.d1 = Dv_a, Dv.d11 = Dv_(aa):
		if (par.PreCond_FD_Order == 2)
		{
			Dv.d1[i00] = 0.5 * ga * (Dv.d0[ip1] - Dv.d0[im1]);
			Dv.d11[i00] = ga2 * (Dv.d0[ip1] + Dv.d0[im1] - 2. * Dv.d0[i00]);
		}
		if (par.PreCond_FD_Order == 4)
		{
			Dv.d1[i00] = ga
					* (2. * (Dv.d0[ip1] - Dv.d0[im1]) / 3.
							- (Dv.d0[ip2] - Dv.d0[im2]) / 12.);
			Dv.d11[i00] = ga2
					* (-2.5 * Dv.d0[i00] + 4. * (Dv.d0[ip1] + Dv.d0[im1]) / 3.
							- (Dv.d0[ip2] + Dv.d0[im2]) / 12.);
		}
	}
	else
	{ // here: Dv.d1  = Dv_(aa),  Dv.d11 = Dv_(aaaa)
		ftype ga4 = ga2 * ga2;
		if (par.PreCond_FD_Order == 2)
		{
			Dv.d1[i00] = 2. * ga2 * (Dv.d0[ip1] - Dv.d0[i00]);
			Dv.d11[i00] = 2. * ga4
					* (3. * Dv.d0[i00] - 4. * Dv.d0[ip1] + Dv.d0[ip2]);
		}
		if (par.PreCond_FD_Order == 4)
		{
			Dv.d1[i00] = ga2
					* (-2.5 * Dv.d0[i00] + 4. * (Dv.d0[ip1] + Dv.d0[im1]) / 3.
							- (Dv.d0[ip2] + Dv.d0[im2]) / 12.);
			Dv.d11[i00] = ga4
					* (6. * Dv.d0[i00] - 4. * (Dv.d0[ip1] + Dv.d0[im1])
							+ (Dv.d0[ip2] + Dv.d0[im2]));
		}
	}
}

void GetFD_Dv_b(parameters par, derivs_2D Dv, int idom, int ipot, int i, int j)
{	// Derivatives of (Dv) with respect to (b):
	int i00 = Index(par, idom, ipot, i, j), nt = par.nt, Nt = nt - 1, ip1 =
			Index(par, idom, ipot, i, j + 1), ip2 = Index(par, idom, ipot, i,
			j + 2), im1 = Index(par, idom, ipot, i, j - 1), im2 = Index(par,
			idom, ipot, i, j - 2);
	ftype hb = Pih / Nt, gb = 1. / hb, gb2 = gb * gb;// hb: Stepsize with respect to (b)

	if (j * (Nt - j) != 0)
	{ // here: Dv.d2  = Dv_b,  Dv.d22 = Dv_(bb)
	  // Dv.d2      = Dv_b, Dv.d22 = Dv_(bb):
		if (par.PreCond_FD_Order == 2)
		{
			Dv.d2[i00] = 0.5 * gb * (Dv.d0[ip1] - Dv.d0[im1]);
			Dv.d22[i00] = gb2 * (Dv.d0[ip1] + Dv.d0[im1] - 2. * Dv.d0[i00]);
		}
		if (par.PreCond_FD_Order == 4)
		{
			Dv.d2[i00] = gb
					* (2. * (Dv.d0[ip1] - Dv.d0[im1]) / 3.
							- (Dv.d0[ip2] - Dv.d0[im2]) / 12.);
			Dv.d22[i00] = gb2
					* (-2.5 * Dv.d0[i00] + 4. * (Dv.d0[ip1] + Dv.d0[im1]) / 3.
							- (Dv.d0[ip2] + Dv.d0[im2]) / 12.);
		}
	}
	else
	{ // here: Dv.d2  = Dv_(bb),  Dv.d22 = Dv_(bbbb)
		ftype gb4 = gb2 * gb2;
		// Dv.d2  = Dv_bb, Dv.d22 = Dv_(bbbb):
		if (par.PreCond_FD_Order == 2)
		{
			Dv.d2[i00] = 2. * gb2 * (Dv.d0[ip1] - Dv.d0[i00]);
			Dv.d22[i00] = 2. * gb4
					* (3. * Dv.d0[i00] - 4. * Dv.d0[ip1] + Dv.d0[ip2]);
		}
		if (par.PreCond_FD_Order == 4)
		{
			Dv.d2[i00] = gb2
					* (-2.5 * Dv.d0[i00] + 4. * (Dv.d0[ip1] + Dv.d0[im1]) / 3.
							- (Dv.d0[ip2] + Dv.d0[im2]) / 12.);
			Dv.d22[i00] = gb4
					* (6. * Dv.d0[i00] - 4. * (Dv.d0[ip1] + Dv.d0[im1])
							+ (Dv.d0[ip2] + Dv.d0[im2]));
		}
	}
}

void GetFD_Dv_ab_FD2(parameters par, derivs_2D Dv, int idom, int ipot, int i,
		int j)
{	// Derivatives of (Dv) with respect to (ab), Finite Difference Order 2:
	int ns = par.ns[idom], Ns = ns - 1, nt = par.nt, Nt = nt - 1, icc = Index(
			par, idom, ipot, i, j), ipc = Index(par, idom, ipot, i + 1, j),
			icp = Index(par, idom, ipot, i, j + 1), imc = Index(par, idom, ipot,
					i - 1, j), icm = Index(par, idom, ipot, i, j - 1), ipp =
					Index(par, idom, ipot, i + 1, j + 1), ipm = Index(par, idom,
					ipot, i + 1, j - 1), imp = Index(par, idom, ipot, i - 1,
					j + 1), imm = Index(par, idom, ipot, i - 1, j - 1);
	ftype // ha: Stepsize with respect to (a), hb: Stepsize with respect to (b)
	ha = Pih / Ns, ga = 1. / ha, ga2 = ga * ga, hb = Pih / Nt, gb = 1. / hb,
			gb2 = gb * gb;

	if (i * (Ns - i) * j * (Nt - j) != 0) // here: Dv.d12 = Dv_(ab)
		Dv.d12[icc] = 0.25 * ga * gb
				* (Dv.d0[ipp] + Dv.d0[imm] - Dv.d0[ipm] - Dv.d0[imp]);
	else
	{
		if (i * (Ns - i) != 0)  // here: j*(Nt-j) = 0: Dv.d12 = Dv_(abb)
			Dv.d12[icc] = ga * gb2
					* ((Dv.d0[ipp] - Dv.d0[ipc]) - (Dv.d0[imp] - Dv.d0[imc]));
		else
		{
			if (j * (Nt - j) != 0)  // here: i*(Nt-i) = 0: Dv.d12 = Dv_(aab)
				Dv.d12[icc] =
						ga2 * gb
								* ((Dv.d0[ipp] - Dv.d0[icp])
										- (Dv.d0[ipm] - Dv.d0[icm]));
			else
				// here: i*(Ns-i)*j*(Nt-j) == 0: Dv.d12 = Dv_(aabb)
				Dv.d12[icc] = 4. * ga2 * gb2
						* (Dv.d0[ipp] + Dv.d0[icc] - Dv.d0[ipc] - Dv.d0[icp]);
		}
	}
}

void GetFD_Dv_ab_FD4(parameters par, derivs_2D Dv, int idom, int ipot, int i,
		int j)
{	// Derivatives of (Dv) with respect to (ab), Finite Difference Order 4:
	int ns = par.ns[idom], Ns = ns - 1, nt = par.nt, Nt = nt - 1,

	im2m2 = Index(par, idom, ipot, i - 2, j - 2), im1m2 = Index(par, idom, ipot,
			i - 1, j - 2), i00m2 = Index(par, idom, ipot, i, j - 2), ip1m2 =
			Index(par, idom, ipot, i + 1, j - 2), ip2m2 = Index(par, idom, ipot,
			i + 2, j - 2),

	im2m1 = Index(par, idom, ipot, i - 2, j - 1), im1m1 = Index(par, idom, ipot,
			i - 1, j - 1), i00m1 = Index(par, idom, ipot, i, j - 1), ip1m1 =
			Index(par, idom, ipot, i + 1, j - 1), ip2m1 = Index(par, idom, ipot,
			i + 2, j - 1),

	im200 = Index(par, idom, ipot, i - 2, j), im100 = Index(par, idom, ipot,
			i - 1, j), i0000 = Index(par, idom, ipot, i, j), ip100 = Index(par,
			idom, ipot, i + 1, j), ip200 = Index(par, idom, ipot, i + 2, j),

	im2p1 = Index(par, idom, ipot, i - 2, j + 1), im1p1 = Index(par, idom, ipot,
			i - 1, j + 1), i00p1 = Index(par, idom, ipot, i, j + 1), ip1p1 =
			Index(par, idom, ipot, i + 1, j + 1), ip2p1 = Index(par, idom, ipot,
			i + 2, j + 1),

	im2p2 = Index(par, idom, ipot, i - 2, j + 2), im1p2 = Index(par, idom, ipot,
			i - 1, j + 2), i00p2 = Index(par, idom, ipot, i, j + 2), ip1p2 =
			Index(par, idom, ipot, i + 1, j + 2), ip2p2 = Index(par, idom, ipot,
			i + 2, j + 2);
	ftype // ha: Stepsize with respect to (a), hb: Stepsize with respect to (b)
	ha = Pih / Ns, ga = 1. / ha, ga2 = ga * ga, hb = Pih / Nt, gb = 1. / hb,
			gb2 = gb * gb;

	if (i * (Ns - i) != 0)
	{
		ftype // a(jj) = [Dv_a at Position (i,jj) with jj \in {j-2, j-1, j, j+1, j+2}]
		am2 = ga
				* (2. * (Dv.d0[ip1m2] - Dv.d0[im1m2]) / 3.
						- (Dv.d0[ip2m2] - Dv.d0[im2m2]) / 12.), am1 = ga
				* (2. * (Dv.d0[ip1m1] - Dv.d0[im1m1]) / 3.
						- (Dv.d0[ip2m1] - Dv.d0[im2m1]) / 12.), a00 = ga
				* (2. * (Dv.d0[ip100] - Dv.d0[im100]) / 3.
						- (Dv.d0[ip200] - Dv.d0[im200]) / 12.), ap1 = ga
				* (2. * (Dv.d0[ip1p1] - Dv.d0[im1p1]) / 3.
						- (Dv.d0[ip2p1] - Dv.d0[im2p1]) / 12.), ap2 = ga
				* (2. * (Dv.d0[ip1p2] - Dv.d0[im1p2]) / 3.
						- (Dv.d0[ip2p2] - Dv.d0[im2p2]) / 12.);
		if (j * (Nt - j) != 0) // here: i*(Ns-i)*j*(Nt-j) != 0; Dv.d12 = Dv_(ab)
			Dv.d12[i0000] = gb * (2. * (ap1 - am1) / 3. - (ap2 - am2) / 12.);
		else
			// here: i*(Ns-i) != 0, j*(Nt-j) = 0: Dv.d12 = Dv_(abb)
			Dv.d12[i0000] = gb2
					* (-2.5 * a00 + 4. * (ap1 + am1) / 3. - (ap2 + am2) / 12.);
	}
	else
	{
		ftype // a(jj) = [Dv_(aa) at Position (i,jj) with jj \in {j-2, j-1, j, j+1, j+2}]
		am2 = ga2
				* (-2.5 * Dv.d0[i00m2] + 4. * (Dv.d0[ip1m2] + Dv.d0[im1m2]) / 3.
						- (Dv.d0[ip2m2] + Dv.d0[im2m2]) / 12.), am1 = ga2
				* (-2.5 * Dv.d0[i00m1] + 4. * (Dv.d0[ip1m1] + Dv.d0[im1m1]) / 3.
						- (Dv.d0[ip2m1] + Dv.d0[im2m1]) / 12.), a00 = ga2
				* (-2.5 * Dv.d0[i0000] + 4. * (Dv.d0[ip100] + Dv.d0[im100]) / 3.
						- (Dv.d0[ip200] + Dv.d0[im200]) / 12.), ap1 = ga2
				* (-2.5 * Dv.d0[i00p1] + 4. * (Dv.d0[ip1p1] + Dv.d0[im1p1]) / 3.
						- (Dv.d0[ip2p1] + Dv.d0[im2p1]) / 12.), ap2 = ga2
				* (-2.5 * Dv.d0[i00p2] + 4. * (Dv.d0[ip1p2] + Dv.d0[im1p2]) / 3.
						- (Dv.d0[ip2p2] + Dv.d0[im2p2]) / 12.);
		if (j * (Nt - j) != 0) // here: i*(Ns-i) = 0, j*(Nt-j) != 0; Dv.d12 = Dv_(aab)
			Dv.d12[i0000] = gb * (2. * (ap1 - am1) / 3. - (ap2 - am2) / 12.);
		else
			// here: i*(Ns-i) = 0, j*(Nt-j) = 0: Dv.d12 = Dv_(aabb)
			Dv.d12[i0000] = gb2
					* (-2.5 * a00 + 4. * (ap1 + am1) / 3. - (ap2 + am2) / 12.);
	}
}

void GetFD_Dv_st(parameters par, derivs_2D Dv, int idom, int ipot, int i, int j)
{	// Derivatives of (Dv) with respect to (s,t):
	int ns = par.ns[idom], Ns = ns - 1, nt = par.nt, Nt = nt - 1, indx = Index(
			par, idom, ipot, i, j);
	ftype a = Pih * i / Ns, sa = sin(a), ca = cos(a), s2a = 2. * sa * ca, s =
			sa * sa, c2a = 1. - 2. * s, b = Pih * j / Nt, sb = sin(b), cb =
			cos(b), s2b = 2. * sb * cb, t = sb * sb, c2b = 1. - 2. * t;

	GetFD_Dv_a(par, Dv, idom, ipot, i, j);
	GetFD_Dv_b(par, Dv, idom, ipot, i, j);
	if (par.PreCond_FD_Order == 2)
		GetFD_Dv_ab_FD2(par, Dv, idom, ipot, i, j);
	if (par.PreCond_FD_Order == 4)
		GetFD_Dv_ab_FD4(par, Dv, idom, ipot, i, j);

	if (i * (Ns - i) != 0)
	{
		Dv.d1[indx] = Dv.d1[indx] / s2a; // Dv_s = Dv_a/s2a
		// Dv_(ss) = (Dv_(aa) - 2*c2a*Dv_s)/(s2a^2)
		Dv.d11[indx] = (Dv.d11[indx] - 2. * c2a * Dv.d1[indx])
				/ (4. * s * (1. - s));
	}
	else
	{
		Dv.d1[indx] = 0.5 * c2a * Dv.d1[indx]; // lim_{s2a->0}[Dv_s] = c2a*Dv_(aa)/2
		// lim_{s2a->0}[Dv_(ss)] = [Dv_(aaaa)+8*c2a*Dv_s]/12
		Dv.d11[indx] = (Dv.d11[indx] + 8. * c2a * Dv.d1[indx]) / 12.;
	}

	if (j * (Nt - j) != 0)
	{
		Dv.d2[indx] = Dv.d2[indx] / s2b; // Dv_t = Dv_b/s2b
		// Dv_(tt) = (Dv_(bb) - 2*c2b*Dv_t)/(s2b^2)
		Dv.d22[indx] = (Dv.d22[indx] - 2. * c2b * Dv.d2[indx])
				/ (4. * t * (1. - t));
	}
	else
	{
		Dv.d2[indx] = 0.5 * c2b * Dv.d2[indx]; // lim_{s2b->0}[Dv_t] = c2b*Dv_(bb)/2
		// lim_{s2b->0}(Dv_tt) = (Dv_(bbbb)+8*c2b*Dv_t)/12
		Dv.d22[indx] = (Dv.d22[indx] + 8. * c2b * Dv.d2[indx]) / 12.;
	}

	if (i * (Ns - i) * j * (Nt - j) != 0)
		Dv.d12[indx] = Dv.d12[indx] / (s2a * s2b);// Dv_(st) = Dv_(ab)/(s2a*s2b)
	else
	{
		if (i * (Ns - i) != 0)  // lim_{s2b->0}(Dv_(st)) = c2b*Dv_(abb)/(2*s2a)
			Dv.d12[indx] = 0.5 * c2b * Dv.d12[indx] / s2a;
		else
		{
			if (j * (Nt - j) != 0) // lim_{s2a->0}(Dv_(st)) = c2a*Dv_(aab)/(2*s2b)
				Dv.d12[indx] = 0.5 * c2a * Dv.d12[indx] / s2b;
			else
				// lim_{s2a->0, s2b->0}(Dv_(st)) = c2a*c2b*Dv_(aabb)/4
				Dv.d12[indx] = 0.25 * c2a * c2b * Dv.d12[indx];
		}
	}
}

void GetFD_Dv(parameters par, int n, ftype *DX, derivs_2D v, derivs_2D Dv)
{  // Calculating the finite difference strcuture Dv from a unit vector DX
// with the non-vanishing entry n_v[n]
	int n_v = par.n_v_of_n[n], idom, ipot, i, j, ns, Ns, nt, Nt;

	Get_Indices_From_n_v(par, n_v, &idom, &ipot, &i, &j);
	ns = par.ns[idom];
	nt = par.nt;
	Ns = ns - 1;
	Nt = nt - 1;

	int ii, jj, i0 = maximum2(0, i - 2) - i, j0 = maximum2(0, j - 2) - j, i1 =
			minimum2(i + 2, Ns) - i, j1 = minimum2(j + 2, Nt) - j;

	Dv.d0[n_v] = 1.;
	if (idom == 0 && i == Ns)
		Dv.d0[Index(par, 1, ipot, par.ns[1] - 1, j)] = 1.;

	for (ii = i0; ii <= i1; ii++)
	{
		for (jj = j0; jj <= j1; jj++)
		{
			GetFD_Dv_st(par, Dv, idom, ipot, i + ii, j + jj);
			if (idom == 0 && i == Ns)
				GetFD_Dv_st(par, Dv, 1, ipot, par.ns[1] - 1 + ii, j + jj);
		}
	}
}

void fill0_Dv_of_indx(int indx, derivs_2D Dv)
{
	Dv.d0[indx] = Dv.d1[indx] = Dv.d2[indx] = Dv.d11[indx] = Dv.d12[indx] =
			Dv.d22[indx] = 0.;
}

void fill0_Dv(parameters par, int n, derivs_2D Dv)
{
	int n_v = par.n_v_of_n[n], idom, ipot, i, j, ns, Ns, nt, Nt;
	Get_Indices_From_n_v(par, n_v, &idom, &ipot, &i, &j);
	ns = par.ns[idom];
	nt = par.nt;
	Ns = ns - 1;
	Nt = nt - 1;

	int ii, jj, i0 = maximum2(0, i - 2) - i, j0 = maximum2(0, j - 2) - j, i1 =
			minimum2(i + 2, Ns) - i, j1 = minimum2(j + 2, Nt) - j;

	for (ii = i0; ii <= i1; ii++)
	{
		for (jj = j0; jj <= j1; jj++)
		{
			int indx = Index(par, idom, ipot, i + ii, j + jj);
			fill0_Dv_of_indx(indx, Dv);
			if (idom == 0 && i == Ns)
			{
				indx = Index(par, 1, ipot, par.ns[1] - 1 + ii, j + jj);
				fill0_Dv_of_indx(indx, Dv);
			}
		}
	}
}

void Get_boundary_i1(parameters par, int jdom, int *adom, int *ai0, int *ai1,
		int *aj0, int *aj1, int *mdom)
{
	adom[*mdom] = jdom;
	ai0[*mdom] = par.ns[jdom] - 3;
	ai1[*mdom] = par.ns[jdom] - 1;
	aj0[*mdom] = aj0[0];
	aj1[*mdom] = aj1[0];
	*mdom += 1;
}

void Get_gridpoints(parameters par, int n, int *mdom, int *adom, int *ai0,
		int *ai1, int *aj0, int *aj1)
{ // Getting the relevant gridpoints at which for a given unit vector DX the evaluation of the
// field equations and boundary conditions yields a non-vanishing entry for J*DX
// At the exit:
//      mdom: number of domains within which relevant gridpoints are located
//      for (jdom = 0; jdom < mdom; jdom ++) :
//          adom[jdom]: contains the number of the domain (jdom), ranging from 0 to 1
//           ai0[jdom]: contains the minimal s-index of gridpoints in the domain adom[jdom]
//           ai1[jdom]: contains the maximal s-index of gridpoints in the domain adom[jdom]
//           aj0[jdom]: contains the minimal t-index of gridpoints in the domain adom[jdom]
//           aj1[jdom]: contains the maximal t-index of gridpoints in the domain adom[jdom]

	int j, ns, Ns, nt, Nt, idom, ipot, i, n_v = par.n_v_of_n[n];

	Get_Indices_From_n_v(par, n_v, &idom, &ipot, &i, &j);
	ns = par.ns[idom];
	nt = par.nt;
	Ns = ns - 1;
	Nt = nt - 1;

	*mdom = 1;
	adom[0] = idom;

	ai0[0] = maximum2(0, i - 2);
	ai1[0] = minimum2(i + 2, Ns);
	if (idom == 0 && i <= 2)
	{
		aj0[0] = 0;
		aj1[0] = par.nt - 1;
	}
	else
	{
		aj0[0] = maximum2(0, j - 2);
		aj1[0] = minimum2(j + 2, Nt);
		if (idom == 0 && ai1[0] == Ns)
			Get_boundary_i1(par, 1, adom, ai0, ai1, aj0, aj1, mdom);
		if (idom == 1 && ai1[0] == Ns)
			Get_boundary_i1(par, 0, adom, ai0, ai1, aj0, aj1, mdom);
	}
}

void Get_JFD_Matrix(parameters par, ftype *X, derivs_2D v, JFD_Components *JFD)
{	// Calculating JFD
	int ntotal = par.ntotal, num_v = par.num_v, ipot, idom, i, i0, i1, j, j0,
			j1, ns, Ns, nt, Nt, row, row_v, column, mcol, mdom, jdom,
			adom[NDOM], ai0[NDOM], aj0[NDOM], ai1[NDOM], aj1[NDOM];
	ftype *DX, *JJDX;
	derivs_2D Dv;

	allocate_derivs_2D(&Dv, num_v);
	fill0_derivs_2D(Dv, num_v);

	DX = new ftype[ntotal];
	fill0_dvector(DX, 0, ntotal);
	JJDX = new ftype[num_v];
	fill0_dvector(JJDX, 0, num_v);

	(*JFD).m1 = (*JFD).m2 = 0;

	for (column = 0; column < ntotal; column++)
	{
		DX[column] = 1.;	// setting the unit vector DX
		fill0_Dv(par, column, Dv);
		GetFD_Dv(par, column, DX, v, Dv);
		Get_gridpoints(par, column, &mdom, adom, ai0, ai1, aj0, aj1);
		for (jdom = 0; jdom < mdom; jdom++)
		{
			idom = adom[jdom];
			i0 = ai0[jdom];
			i1 = ai1[jdom];
			j0 = aj0[jdom];
			j1 = aj1[jdom];
			ns = par.ns[idom];
			Ns = ns - 1;
			nt = par.nt;
			Nt = nt - 1;
			for (j = j0; j <= j1; j++)
			{
				for (i = i0; i <= i1; i++)
				{
					LinFieldEqns(par, idom, i, j, v, Dv, JJDX);
				}
			}
		}
		for (jdom = 0; jdom < mdom; jdom++)
		{
			idom = adom[jdom];
			i0 = ai0[jdom];
			i1 = ai1[jdom];
			j0 = aj0[jdom];
			j1 = aj1[jdom];
			for (ipot = 0; ipot < NPOT; ipot++)
			{
				for (j = j0; j <= j1; j++)
				{
					for (i = i0; i <= i1; i++)
					{
						row_v = Index(par, idom, ipot, i, j);
						row = par.n_of_n_v[row_v];
						if (row >= 0)
						{
							if (fabs(JJDX[row_v]) > TINY)
							{
								mcol = (*JFD).ncols_J[row];
//									cout<<"row = "<<row<<" \t row_max="<<n_2D-1<<"\t mcol="<<mcol<<endl;
								(*JFD).cols_J[row][mcol] = column;
								(*JFD).J[row][mcol] = JJDX[row_v];
								(*JFD).ncols_J[row] += 1;
								if (column > row) // determining the numbers m1 and m2
									(*JFD).m2 = maximum2((*JFD).m2,
											column - row);
								else
									(*JFD).m1 = maximum2((*JFD).m1,
											row - column);
							}
						}
						JJDX[row_v] = 0.;
					}
				}
			}
		}
		fill0_Dv(par, column, Dv);
		DX[column] = 0.;
	}
	free_derivs_2D(&Dv, num_v);
	delete[] JJDX;
	delete[] DX;
}

void Get_JFD_Components(parameters par, JFD_Components *JFD)
{ // Calculates the remaining components of JFD, see pages 14-18
// Note: For the band-matrix operations we use the 'Numerical Recipes in C'
// which provide the corresponding routines for matrices and vectors with indices
// starting from 1.
	int n_2D = par.n_2D, m1 = (*JFD).m1, m2 = (*JFD).m2, j, k, row_K, col_K,
			col_J;
	ftype d = -1., *aux_F, *aux_X;

	aux_X = new ftype[n_2D];
	aux_F = new ftype[n_2D];

	for (j = 0; j < n_2D; j++)
	{
		row_K = j + 1; // see Note
		for (k = 0; k < (*JFD).ncols_J[j]; k++)
		{
			col_J = (*JFD).cols_J[j][k];
			col_K = col_J + 1; // see Note
			if (row_K - col_K <= m1 && col_K - row_K <= m2) // writing K
				(*JFD).K[row_K][m1 + 1 + col_K - row_K] = (*JFD).J[j][k];
		}
	}

	// Calculating the LU-decomposition of K; introducing Kl and iK,
	// see 'Numerical Recipes in C' pages 51-54
	bandec((*JFD).K, n_2D, m1, m2, (*JFD).Kl, (*JFD).iK, &d);

	delete[] aux_X;
	delete[] aux_F;
}


