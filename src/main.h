#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include "utilities.h"
#include "cast.h"

#define NDOM 		2
#define NPOT 		4
#define Hmax 		7
#define Ndat        10

typedef struct PARAMETERS
{
	ftype sig0, sig1, j_BH, *a, *kdd, *kuu, *b, *ldd, *luu, *c, *mdd, *muu, *d,
			*e, *h, *s, *H, Newton_tol, bicgstab_decr, bicgstab_tol,
			PreCond_Improve, xmin, xmax, rhomax;
	int ntotal, num_v, n_2D, ns_dat, nt_dat, ns[NDOM], nt, nst[NDOM + 1],
			*n_of_n_v, *n_v_of_n, n_Seq, Create_ID,
			Newton_itmin, Newton_itmax, Newton_verb, bicgstab_itmax,
			bicgstab_verb, PreCond_itmax, PreCond_FD_Order, test_wr_Jac;
	char which_data[Ndat];
} parameters;

typedef struct COORDINATE
{
	ftype d0, d1, d2, d11, d12, d22, d122;
} coordinate;

typedef struct DERIVS_2D
{
	ftype *d0, *d1, *d2, *d11, *d12, *d22;
} derivs_2D;

typedef struct JFD_COMPONENTS
{
	ftype **J, **K, **Kl;
	int *ncols_J, **cols_J, *iK, m1, m2;
} JFD_Components;

// Subroutines in "ACMC_Data.c"
int main(int argc, char **argv);

// Subroutines in "FuncAndJacobian.c"
void Derivatives_st(parameters par, derivs_2D v);
void Get_v_From_X(parameters par, ftype *X, derivs_2D v);
void F_of_X(parameters par, ftype *X, derivs_2D v, ftype *F);
void J_times_DX(parameters par, ftype *X, ftype *DX, derivs_2D v, ftype *JDX);

// Subroutines in "IndexRoutines.c"
int Index_i(parameters par, int idom, int i);
int Index_j(parameters par, int idom, int j);
int Index(parameters par, int idom, int ipot, int i, int j);
void Get_Arrays_n_And_n_v_domain_0(parameters *par, int *n);
void Get_Arrays_n_And_n_v_domain_1(parameters *par, int *n);
void Get_Arrays_n_And_n_v(parameters *par);
void Get_Indices_From_n_v(parameters par, int n_v, int *idom, int *ipot, int *i,
		int *j);
void Free_Arrays_n_And_n_v(parameters *par);

// Subroutines in "FieldAndBoundEqns.c"
void func_AijAij(parameters par, int idom, int i, int j, derivs_2D v,
		ftype fac1, ftype *AijAij);
void func_DAijAij(parameters par, int idom, int i, int j, derivs_2D v,
		derivs_2D Dv, ftype fac1, ftype *AijAij, ftype *DAijAij);
void func_h3(parameters par, int idom, int i, int j, derivs_2D v, ftype fac1,
		ftype sin2th, ftype *h3);
void func_Dh3(parameters par, int idom, int i, int j, derivs_2D v, derivs_2D Dv,
		ftype fac1, ftype sin2th, ftype *h3, ftype *Dh3);
void NonLinFieldEqns(parameters par, int idom, int i, int j, derivs_2D v,
		ftype *F);
void LinFieldEqns(parameters par, int idom, int i, int j, derivs_2D v,
		derivs_2D Dv, ftype *JDX);

// Subroutines in "ScanAndPrint.c"
void ScanInitialData(parameters *par, ftype **X, char file[]);
void ScanConfig(parameters *parInit, parameters *parGoal);
void ScanFiles(parameters *parInit, ftype **X, parameters *parGoal);
void PrintToFile(parameters par, ftype *X, int i_Seq);
void CreateInitialData(parameters par);
void TestData(parameters par, ftype *X);

// Subroutines in "GetNewResolution.c"
void v1_To_v2(parameters par1, derivs_2D v1, parameters par2, derivs_2D v2);
void GetNewResolution(parameters par0, parameters par, ftype *X0, ftype *X);

// Subroutines in "allocate_and_free_derivs.c"
void allocate_derivs_2D(derivs_2D *v, int n);
void fill0_derivs_2D(derivs_2D v, int n);
void free_derivs_2D(derivs_2D *v, int n);

// Subroutines in "newton.c"
void resid(parameters par, ftype *DX, ftype *rhs, JFD_Components JFD,
		ftype *res);
void Solve_JFD(parameters par, ftype *F, ftype *X, JFD_Components JFD);
void Solve_JFD(parameters par, ftype *F, ftype *DX, JFD_Components JFD);
int PreCond(parameters par, ftype *F, ftype *DX, ftype *norm,
		JFD_Components JFD);
int bicgstab_Iterations(parameters par, ftype *X, derivs_2D v, ftype *F,
		ftype *DX, ftype *normres, JFD_Components JFD);
int bicgstab(parameters par, ftype *X, derivs_2D v, ftype *F, ftype *DX,
		ftype *normres);
int newton(parameters par, ftype *X);

// Subroutines in "JacobianFD.c"
void GetFD_Dv_a(parameters par, derivs_2D Dv, int idom, int ipot, int i, int j);
void GetFD_Dv_b(parameters par, derivs_2D Dv, int idom, int ipot, int i, int j);
void GetFD_Dv_ab_FD2(parameters par, derivs_2D Dv, int idom, int ipot, int i,
		int j);
void GetFD_Dv_ab_FD4(parameters par, derivs_2D Dv, int idom, int ipot, int i,
		int j);
void GetFD_Dv_st(parameters par, derivs_2D Dv, int idom, int ipot, int i,
		int j);
void GetFD_Dv(parameters par, int n, ftype *DX, derivs_2D v, derivs_2D Dv);
void fill0_Dv_of_indx(int indx, derivs_2D Dv);
void fill0_Dv(parameters par, int n, derivs_2D Dv);
void Get_boundary_i1(parameters par, int jdom, int *adom, int *ai0, int *ai1,
		int *aj0, int *aj1, int *mdom);
void Get_gridpoints(parameters par, int n, int *mdom, int *adom, int *ai0,
		int *ai1, int *aj0, int *aj1);
void Get_JFD_Matrix(parameters par, ftype *X, derivs_2D v, JFD_Components *JFD);
void Get_JFD_Components(parameters par, JFD_Components *JFD);

// Subroutines in "Get_Coefficients.c"
void read_config(char name[], int *ns_dat, int *nt_dat, ftype *sig_min,
		ftype *sig_max, ftype *j_BH);
int scan_coefficient(char name[], ftype **a, int ns_dat, int nt_dat);
int indx_4(parameters par, int i1, int i2, int i3, int i4, int idom, int i,
		int j);
void Get_Coefficient_4(parameters *par, ftype *par_a, char prefix[]);
int indx_3(parameters par, int i1, int i2, int i3, int idom, int i, int j);
void Get_Coefficient_3(parameters *par, ftype *par_b, char prefix[]);
int indx_2(parameters par, int i1, int i2, int idom, int i, int j);
void Get_Coefficient_2(parameters *par, ftype *par_c, char prefix[]);
int indx_1(parameters par, int i1, int idom, int i, int j);
void Get_Coefficient_1(parameters *par, ftype *par_d, char prefix[]);
int indx_H(parameters par, int i1, int idom, int i, int j);
void Get_Coefficient_H(parameters *par, ftype *par_H, char prefix[]);
int indx_Horizon(parameters par, int i1, int j);
void Get_Coefficient_Horizon(parameters *par, ftype *par_h, char prefix[],
		int odd);
void Get_Coefficients(parameters *par);
void Free_Coefficients(parameters *par);

//Subroutines in "BondiMass.c"
void LVmm(parameters par, derivs_2D v, ftype *LVmm_Scri);
ftype qcinv_Scri(parameters par, int ii, int jj, int j);
ftype gammac110_Scri(ftype j);
ftype sqrtdetqc_Scri();
void getphi3(ftype **phi, ftype **phi3, int ns, int nt, ftype ds);
ftype sqrtdetqc(parameters par, ftype sig, ftype mu, ftype B, ftype Hd2);
void GetAdd3j(parameters par, derivs_2D v, ftype **Add31, ftype **Add32,
		ftype **Add33);
ftype gammacuuij(parameters par, ftype sig, ftype mu, ftype i, ftype j);
ftype gammacdd3i(ftype sig, ftype mu, ftype j, int i);
ftype A1(ftype sig, ftype mu, ftype j);
ftype A2(ftype sig, ftype mu, ftype j);
void BondiMass(parameters par, ftype *X);
ftype CalculateAngularMomentum(parameters par, derivs_2D v, ftype *H,
		ftype *Hd2);

//Subroutines in "functions.c"
ftype cast_char(const char *number);
