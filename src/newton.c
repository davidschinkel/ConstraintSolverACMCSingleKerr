#include "main.h"
#define STENCILSIZE 32


void resid(parameters par, ftype *DX, ftype *rhs, JFD_Components JFD,
		ftype *res)
{ // Calculates res = rhs - JFD*DX
	int j, k, col, n_2D = par.n_2D;

	for (j = 0; j < n_2D; j++)
	{
		res[j] = 0.;
		for (k = 0; k < JFD.ncols_J[j]; k++)
		{
			col = JFD.cols_J[j][k];
			res[j] -= JFD.J[j][k] * DX[col];
		}
		res[j] += rhs[j];
	}
}

void Solve_JFD(parameters par, ftype *F, ftype *X, JFD_Components JFD)
{ // Solves the submatrix-problem J*X = F
	int j, n_2D = par.n_2D, m1 = JFD.m1, m2 = JFD.m2;
	ftype *y;

	y = new ftype[n_2D + 1];

	// Calculating the vector y of eqn. 2.7.20 with A=K, B=F
	for (j = 1; j <= n_2D; j++)
		y[j] = F[j - 1];
	banbks(JFD.K, n_2D, m1, m2, JFD.Kl, JFD.iK, y);
	for (j = 1; j <= n_2D; j++)
		X[j - 1] = y[j];
	delete[] y;
}

int PreCond(parameters par, ftype *F, ftype *DX, ftype *norm,
		JFD_Components JFD)
{ // Performs and improves the preconditioning step up to <PreCond_itmax> times
// using the methods of iterative improvements, 'Numerical Recipes in C', pages 55-58
	int iter = 0, n, ntotal = par.ntotal;
	ftype norm0 = norm2(F, ntotal), *res, *XX;

	res = new ftype[ntotal];
	XX = new ftype[ntotal];

	fill0_dvector(DX, 0, ntotal);
	copy_dvector(res, F, 0, ntotal);
	do
	{
		iter += 1;
		Solve_JFD(par, res, XX, JFD);
		for (n = 0; n < ntotal; n++)
			DX[n] += XX[n];
		resid(par, DX, F, JFD, res);
		*norm = norm2(res, ntotal);
	} while (iter < par.PreCond_itmax && *norm > par.PreCond_Improve * norm0);

	delete[] res;
	delete[] XX;

	return iter;
}

int bicgstab_Iterations(parameters par, ftype *X, derivs_2D v, ftype *F,
		ftype *DX, ftype *normres, JFD_Components JFD)
{	// Iterations for the Biconjugate Stabilized Gradient Method
	int ntotal = par.ntotal, iter, n, N_Pre;
	ftype norm, alpha = 0., beta = 0., rho = 0., rho1 = 1., rhotol = 1e-50,
			omega = 0., omegatol = 1e-50, *p, *ph, *q, *r, *rt, *s, *sh, *t;

	p = new ftype[ntotal];
	ph = new ftype[ntotal];
	q = new ftype[ntotal];
	r = new ftype[ntotal];
	s = new ftype[ntotal];
	sh = new ftype[ntotal];
	rt = new ftype[ntotal];
	t = new ftype[ntotal];

	// 	compute initial residual rt = r = p = F - J*DX
	J_times_DX(par, X, DX, v, r);
	for (n = 0; n < ntotal; n++)
		rt[n] = r[n] = p[n] = F[n] - r[n];

	*normres = norm2(r, ntotal);
	if (*normres <= par.bicgstab_tol)
		return 0;

	for (iter = 0; iter < par.bicgstab_itmax; iter++)
	{
		rho = scalarproduct(rt, r, ntotal);
		if (fabs(rho) < rhotol)
			break;
		if (iter > 0)
		{ // compute direction vector p
			beta = (rho / rho1) * (alpha / omega);
			for (n = 0; n < ntotal; n++)
				p[n] = r[n] + beta * (p[n] - omega * q[n]);
		}
		if (par.bicgstab_verb == 1)
			 cout <<" ("<< norm1(p, ntotal) << ",";
		// compute direction adjusting vector ph and scalar alpha :
		N_Pre = PreCond(par, p, ph, &norm, JFD);
		if (par.bicgstab_verb == 1)
			 cout<<norm<<") N_Pre="<<N_Pre<<" | ";
		J_times_DX(par, X, ph, v, q);	// q = J*ph
		alpha = rho / scalarproduct(rt, q, ntotal);
		for (n = 0; n < ntotal; n++)
			s[n] = r[n] - alpha * q[n];
		// early check of tolerance:
		*normres = norm2(s, ntotal);
		if (*normres <= par.bicgstab_tol)
		{
			for (n = 0; n < ntotal; n++)
				DX[n] += alpha * ph[n];
			if (par.bicgstab_verb == 1)
			{
				cout<<iter+1<<" "<<*normres<<flush;
			}

			break;
		}
		// compute stabilizer vector sh and scalar omega:
		if (par.bicgstab_verb == 1)
			cout <<" ("<< norm1(s, ntotal)<<",";
		N_Pre = PreCond(par, s, sh, &norm, JFD);
		if (par.bicgstab_verb == 1)
			cout<<norm<<") N_Pre="<<N_Pre<<" | ";
		J_times_DX(par, X, sh, v, t);	// t = J*sh
		omega = scalarproduct(t, s, ntotal) / scalarproduct(t, t, ntotal);
		// compute new solution approximation:
		for (n = 0; n < ntotal; n++)
		{
			DX[n] += alpha * ph[n] + omega * sh[n];
			r[n] = s[n] - omega * t[n];
		}
		// are we done? 
		*normres = norm2(r, ntotal);
		if (par.bicgstab_verb == 1)
		{
			cout<<iter+1<<" "<<*normres<<flush;
		}
		if (*normres <= par.bicgstab_tol)
			break;
		rho1 = rho;
		if (fabs(omega) < omegatol)
			break;
	}
	if (par.bicgstab_verb == 1)
		cout<<endl;

	delete[] p;
	delete[] ph;
	delete[] q;
	delete[] r;
	delete[] s;
	delete[] sh;
	delete[] rt;
	delete[] t;

	/* iteration failed */
	if (iter > par.bicgstab_itmax)
		return -1;

	/* breakdown */
	if (fabs(rho) < rhotol)
		return -10;
	if (fabs(omega) < omegatol)
		return -11;

	/* success! */
	return iter + 1;
}

int bicgstab(parameters par, ftype *X, derivs_2D v, ftype *F, ftype *DX,
		ftype *normres)
{// Biconjugate Stabilized Gradient Method, see pages 1, 2, http://www.netlib.org/templates/
	int n_2D = par.n_2D, maxcol = NPOT * STENCILSIZE, iter;
	JFD_Components JFD;

	// allocating the components of JFD
	JFD.J = dmatrix(0, n_2D - 1, 0, maxcol - 1);
	fill0_dmatrix(JFD.J, 0, n_2D, 0, maxcol);
	JFD.cols_J = imatrix(0, n_2D - 1, 0, maxcol - 1);
	fill0_imatrix(JFD.cols_J, 0, n_2D, 0, maxcol);
	JFD.ncols_J = ivector(0, n_2D - 1);
	fill0_ivector(JFD.ncols_J, 0, n_2D);

	Get_JFD_Matrix(par, X, v, &JFD);

	JFD.K = dmatrix(1, n_2D, 1, JFD.m1 + JFD.m2 + 1);
	fill0_dmatrix(JFD.K, 1, n_2D + 1, 1, JFD.m1 + JFD.m2 + 2);
	JFD.Kl = dmatrix(1, n_2D, 1, JFD.m1);
	fill0_dmatrix(JFD.Kl, 1, n_2D + 1, 1, JFD.m1 + 1);
	JFD.iK = ivector(1, n_2D);
	fill0_ivector(JFD.iK, 1, n_2D + 1);

	Get_JFD_Components(par, &JFD);

	if (par.bicgstab_verb == 1)
		cout<<" bicgstab:  itmax = "<<par.bicgstab_itmax<<", tol = "<<par.bicgstab_tol<<endl;

	iter = bicgstab_Iterations(par, X, v, F, DX, normres, JFD);

	// free the components of JFD
	free_dmatrix(JFD.J, 0, n_2D - 1, 0, maxcol - 1);
	free_imatrix(JFD.cols_J, 0, n_2D - 1, 0, maxcol - 1);
	free_ivector(JFD.ncols_J, 0, n_2D - 1);

	free_dmatrix(JFD.K, 1, n_2D, 1, JFD.m1 + JFD.m2 + 1);
	free_dmatrix(JFD.Kl, 1, n_2D, 1, JFD.m1);
	free_ivector(JFD.iK, 1, n_2D);

	return iter;
}

int newton(parameters par, ftype *X)
{	// Newton Raphson Method, see pages 1, 2
	int ntotal = par.ntotal, num_v = par.num_v, iter = 0, j, bicgstab_iter;
	ftype *F, *DX, norm, normres;
	derivs_2D v;

	F = new ftype[ntotal];
	DX = new ftype[ntotal];

	allocate_derivs_2D(&v, num_v);

	F_of_X(par, X, v, F);
	norm = norm2(F, ntotal);

	if (par.Newton_verb == 1)
	{
		cout<<" Newton Raphson Method: Initial Residual: \t |F| = "<<norm<<endl;
		cout<<" --------------------------------------------------------"<<endl;
	}
	while (iter < par.Newton_itmin || (norm > par.Newton_tol && iter < par.Newton_itmax))
	{
		iter += 1;
		fill0_dvector(DX, 0, ntotal);

		par.bicgstab_tol = par.bicgstab_decr * norm;
		bicgstab_iter = bicgstab(par, X, v, F, DX, &normres);

		for (j = 0; j < ntotal; j++)
		{
			X[j] -= DX[j];
		}

		F_of_X(par, X, v, F);
		norm = norm2(F, ntotal);
		if (par.Newton_verb == 1)
			if (boost::math::isinf(norm) || boost::math::isnan(norm) || norm > 1.0e+10)
			{
				cout<<"\n No Convergence of the Newton Raphson Method. Now exiting to system."<<endl;
 				exit(1);
			}
		cout<<" Newton: iter ="<<iter<<" \t |F| = "<<norm<<endl;
	}
	if (norm > par.Newton_tol)
	{
		cout<<" Newton Raphson Method failed to converge to prescribed tolerance. Now move on to the next sequence element."<<endl;
	}
	delete[] F;
	delete[] DX;
	free_derivs_2D(&v, num_v);

	return iter;
}
