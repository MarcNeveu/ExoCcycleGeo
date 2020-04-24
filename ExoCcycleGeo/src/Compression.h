/*
 * Compression.h
 *
 *  Created on: Jun 2, 2014, modified for ExoCcycleGeo Jun 13, 2019
 *      Author: Marc Neveu (marc.f.neveu@nasa.gov)
 *
 *      Mass-radius relationship code developed by Lorenzo, Desch and Shim, modified to include water and ammonia.
 *
 *      This code does not use the bulk modulus, but solves the equation of state (EOS) directly. It uses
 *      an EOS rho = rho0 + c*P^n, as described in Seager et al. 2007.
 *      This provides a smooth transition from the experimental regime to the theoretical regime, and avoids
 *      unphysical results arising from extrapolating Birch-Murnaghan outside its range of validity.  Parameters
 *      taken from Seager et al. (2007).
 *      This code allows for three layers: inner & outer core, & mantle.
 *
 *      The code is in SI throughout.
 *
 *      Reference: Lorenzo A., Desch S.J. & Shim S.-H. (2014) On the Lower Radius Limit of Exoplanets
 *      45th LPSC, abstract 1636
 */

#include "ExoCcycleGeo.h"

#ifndef COMPRESSION_H_
#define COMPRESSION_H_

//-------------------------------------------------------------------
// SUBROUTINE DECLARATIONS
//-------------------------------------------------------------------

int compression(int NR, double m_p, double m_c, double Tsurf, int dbincore, int dboutcore, int dbmantle, double **r,
		double **P, double **rho, double **g, int *cmbindex, char path[1024]);

int planmat(int ncomp, int **dbindex, int **eos, double **rho0, double **c, double **nn, double **Ks0, double **Ksp,
		double **V0, double **Tref, double **a0, double **a1, double **b0, double **b1, double **b2, char path[1024]);

int planmat_index(int mat, int ncomp, int *dbindex);

/*--------------------------------------------------------------------
 *
 * Subroutine Compression
 *
 * Calculates pressure ad density radial profiles inside a planetary
 * body by allowing the materials to compress. Handles up to 3
 * compositional layers (e.g. core, mantle, crust).
 *
 *--------------------------------------------------------------------*/

int compression(int NR, double m_p, double m_c, double Tsurf, int dbincore, int dboutcore, int dbmantle, double **r,
		double **P, double **rho, double **g, int *cmbindex, char path[1024]) {

	//-------------------------------------------------------------------
	// Declarations and initializations
	//-------------------------------------------------------------------

	// Planetary parameters
	double Mincore = 0.0; // Inner core mass
	double Moutcore = 0.0; // Outer core mass
	double Rp = 0.0; // Planetary radius
	double Rincore = 0.0; // Inner core radius
	double Routcore = 0.0; // Outer core radius
	double Mmantle = 0.0; // Mantle mass

	// Grid variables
	int ir = 0;
	int nic = 0; // N inner core
	int noc = 0; // N outer core
	double Vol = 0.0; // Volume
	double Pavg = 0.0; // Average pressure

	double *M = (double*) malloc((NR+1)*sizeof(double)); // Mass
	if (M == NULL) printf("Compression: Not enough memory to create M[NR]\n");

	double *dM = (double*) malloc((NR+1)*sizeof(double)); // dMass
	if (dM == NULL) printf("Compression: Not enough memory to create dM[NR]\n");

	double *rhonew = (double*) malloc((NR+1)*sizeof(double)); // New density
	if (rhonew == NULL) printf("Compression: Not enough memory to create rhonew[NR]\n");

	// Compositional variables
	int ncomp = 50;
	int iincore = 0;
	int ioutcore = 0;
	int imantle = 0;

	int *icomp = (int*) malloc((NR+1)*sizeof(int));
	if (icomp == NULL) printf("Compression: Not enough memory to create icomp[NR]\n");

	int *dbindex = (int*) malloc(ncomp*sizeof(int)); // Planetary material database index
	if (dbindex == NULL) printf("Compression: Not enough memory to create dbindex[ncomp]\n");

	int *eos = (int*) malloc(ncomp*sizeof(int)); // Switch between types of equation of state. 1: 3rd order Birch-Murnaghan EOS, 2: rho0+cP^n
	if (eos == NULL) printf("Compression: Not enough memory to create eos[ncomp]\n");

	double *rho0 = (double*) malloc(ncomp*sizeof(double)); // Constant density term of rho = rho0 + c*P^nn
	if (rho0 == NULL) printf("Compression: Not enough memory to create rho0[ncomp]\n");

	double *c = (double*) malloc(ncomp*sizeof(double)); // c of rho = rho0 + c*P^nn
	if (c == NULL) printf("Compression: Not enough memory to create c[ncomp]\n");

	double *nn = (double*) malloc(ncomp*sizeof(double)); // Exponent of rho = rho0 + c*P^nn
	if (nn == NULL) printf("Compression: Not enough memory to create nn[ncomp]\n");

	double *Ks0 = (double*) malloc(ncomp*sizeof(double)); // Constant bulk modulus term
	if (Ks0 == NULL) printf("Compression: Not enough memory to create Ks0[ncomp]\n");

	double *Ksp = (double*) malloc(ncomp*sizeof(double)); // Derivative w.r.t. P of the bulk modulus term
	if (Ks0 == NULL) printf("Compression: Not enough memory to create Ks0[ncomp]\n");

	double *V0 = (double*) malloc(ncomp*sizeof(double)); // V0 term (m3/kg) in eqs. (6) & (7) of Choukroun & Grasset (2010)
	if (V0 == NULL) printf("Compression: Not enough memory to create V0[ncomp]\n");

	double *Tref = (double*) malloc(ncomp*sizeof(double)); // Tref term (K) in eqs. (6) & (7) of Choukroun & Grasset (2010)
	if (Tref == NULL) printf("Compression: Not enough memory to create Tref[ncomp]\n");

	double *a0 = (double*) malloc(ncomp*sizeof(double)); // a0 term (no dim) in eqs. (6) & (7) of Choukroun & Grasset (2010)
	if (a0 == NULL) printf("Compression: Not enough memory to create a0[ncomp]\n");

	double *a1 = (double*) malloc(ncomp*sizeof(double)); // a1 term (K-1) in eqs. (6) & (7) of Choukroun & Grasset (2010)
	if (a1 == NULL) printf("Compression: Not enough memory to create a1[ncomp]\n");

	double *b0 = (double*) malloc(ncomp*sizeof(double)); // b0 term (no dim) in eqs. (6) & (7) of Choukroun & Grasset (2010)
	if (b0 == NULL) printf("Compression: Not enough memory to create b0[ncomp]\n");

	double *b1 = (double*) malloc(ncomp*sizeof(double)); // b1 term (no dim) in eqs. (6) & (7) of Choukroun & Grasset (2010)
	if (b1 == NULL) printf("Compression: Not enough memory to create b1[ncomp]\n");

	double *b2 = (double*) malloc(ncomp*sizeof(double)); // b2 term (Pa-1) in eqs. (6) & (7) of Choukroun & Grasset (2010)
	if (b2 == NULL) printf("Compression: Not enough memory to create b2[ncomp]\n");

	// Iteration control
	int iter = 0;
	int itermax = 100;
	int j = 0;
	int jmax = 20;
	double delta = 0.0; // Difference in density change
	double mix = 1.0; // Memory of old iteration
	double dy = 0.0; // Birch-Murnaghan equation (right term - left term) to solve by Newton-Raphson
	double dydx = 0.0; // Derivative of Birch-Murnaghan equation
	double x = 0.0; // rho/rho0 or eta
	double z = 0.0;

	for (ir=0;ir<=NR;ir++) {
		(*r)[ir] = 0.0;
		M[ir] = 0.0;
		dM[ir] = 0.0;
		(*g)[ir] = 0.0;
		(*rho)[ir] = 0.0;
		rhonew[ir] = 0.0;
		(*P)[ir] = 0.0;
	}
	for (ir=0;ir<ncomp;ir++) {
		icomp[ir] = 0;
		dbindex[ir] = 0;
		eos[ir] = 0;
		rho0[ir] = 0.0;
		c[ir] = 0.0;
		nn[ir] = 0.0;
		Ks0[ir] = 0.0;
		Ksp[ir] = 0.0;
		V0[ir] = 0.0;
		Tref[ir] = 0.0;
		a0[ir] = 0.0;
		a1[ir] = 0.0;
		b0[ir] = 0.0;
		b1[ir] = 0.0;
		b2[ir] = 0.0;
	}

	//-------------------------------------------------------------------
	// Setup
	//-------------------------------------------------------------------

	// Load planetary materials database
	planmat(ncomp, &dbindex, &eos, &rho0, &c, &nn, &Ks0, &Ksp, &V0, &Tref, &a0, &a1, &b0, &b1, &b2, path);

	iincore = planmat_index(dbincore, ncomp, dbindex);
	ioutcore = planmat_index(dboutcore, ncomp, dbindex);
	imantle = planmat_index(dbmantle, ncomp, dbindex);

	// Initialize grids
	Moutcore = m_c;
	Mmantle = m_p-Mincore-Moutcore; // Here, assume no inner core

	Vol = Mincore/rho0[iincore];
	Rincore = pow(0.75*Vol/PI_greek,1.0/3.0);
	Vol = Moutcore/rho0[ioutcore];
	Routcore = pow(0.75*Vol/PI_greek + Rincore*Rincore*Rincore,1.0/3.0);
	Vol = Mmantle/rho0[imantle];
	Rp = pow(0.75*Vol/PI_greek + Routcore*Routcore*Routcore,1.0/3.0);

	nic = floor(Rincore/Rp*(double)NR);
	noc = floor(Routcore/Rp*(double)NR);

	for (ir=0;ir<=nic;ir++) {
		icomp[ir] = iincore;
		(*r)[ir] = Rincore*(double)ir/(double)nic;
		(*rho)[ir] = rho0[icomp[ir]];
	}
	if (nic == 0) {
		ir = 0;
		(*r)[ir] = 0.0;
		(*rho)[ir] = rho0[icomp[ir]];
	}

	for (ir=nic+1;ir<=noc;ir++) {
		icomp[ir] = ioutcore;
		(*r)[ir] = Rincore+(Routcore-Rincore)*(double)(ir-nic)/(double)(noc-nic);
		(*rho)[ir] = rho0[icomp[ir]];
	}
	if (noc == 0) {
		ir = 0;
		(*r)[ir] = 0.0;
		(*rho)[ir] = rho0[icomp[ir]];
	}

	for (ir=noc+1;ir<=NR;ir++) {
		icomp[ir] = imantle;
		(*r)[ir] = Routcore+(Rp-Routcore)*(double)(ir-noc)/(double)(NR-noc);
		(*rho)[ir] = rho0[icomp[ir]];
	}

	for (ir=1;ir<=NR;ir++) {
		Vol = (PI_greek/0.75)*(pow((*r)[ir],3)-pow((*r)[ir-1],3));
		dM[ir] = (*rho)[ir]*Vol;
		M[ir] = M[ir-1]+dM[ir];
	}

	//-------------------------------------------------------------------
	// Begin iterations
	//-------------------------------------------------------------------

	while (iter<itermax) {
		iter++;

		// Calculate gravity at zone boundaries
		for (ir=1;ir<=NR;ir++) (*g)[ir] = G*M[ir]/(*r)[ir]/(*r)[ir];

		// Calculate pressure at zone boundaries, for a given density structure, by integrating hydrostatic equilibrium
		for (ir=NR-1;ir>=0;ir--) (*P)[ir] = (*P)[ir+1] + 0.5*((*g)[ir+1]+(*g)[ir])*((*r)[ir+1]-(*r)[ir])*(*rho)[ir+1];

		// Calculate what compression is needed to produce this pressure, using rho=rho0+c*P^nn (Seager et al. 2007)
		for (ir=1;ir<=NR;ir++) {
			Pavg = 0.5*((*P)[ir]+(*P)[ir-1]);

			switch (eos[icomp[ir]]) {
				case 1:
					z = 0.75*(Ksp[icomp[ir]]-4.0);
					x = 1.0;
					for (j=0;j<jmax;j++) { // Newton-Raphson
						dy = 1.5*Ks0[icomp[ir]]*x*x*x*x*x*(z*x*x*x*x + (1.0-2.0*z)*x*x + (z-1.0)) - Pavg;
						dydx = 1.5*Ks0[icomp[ir]]*x*x*x*x*(9.0*z*x*x*x*x + 7.0*(1.0-2.0*z)*x*x + 5.0*(z-1));
						x = x - dy/dydx;
					}
					rhonew[ir] = rho0[icomp[ir]]*x*x*x;
					break;
				case 2:
					rhonew[ir] = rho0[icomp[ir]] + c[icomp[ir]]*pow(Pavg,nn[icomp[ir]]);
					break;
				case 3: // Note: switch Tsurf to local T[ir] if available.
					rhonew[ir] = V0[icomp[ir]] * (1.0 + a0[icomp[ir]]*tanh(a1[icomp[ir]]*(Tsurf-Tref[icomp[ir]]))) *
									(b0[icomp[ir]] + b1[icomp[ir]]*(1.0-tanh(b2[icomp[ir]]*(*P)[ir])));
					rhonew[ir] = 1.0/rhonew[ir];
					break;
				default:
					printf("Compression: Error: specify EOS type in database");
					exit(0);
			}
		}

		// Update the density using a mix of old density and predicted new density
		delta = fabs(rhonew[1]/(*rho)[1]-1.0);

		for (ir=1;ir<=NR;ir++) (*rho)[ir] = rhonew[ir]*mix + (*rho)[ir]*(1.0-mix);

		// Update radial grid so that each shell, with known mass dM, now has updated density rho
		(*r)[0] = 0.0;
		for (ir=1;ir<=NR;ir++) {
			Vol = dM[ir]/(*rho)[ir];
			(*r)[ir] = pow(0.75*dM[ir]/(*rho)[ir]/PI_greek + pow((*r)[ir-1],3),1.0/3.0);
		}

		// If the central density change from one iteration to the next does not exceed a threshold, then we're done
		if (delta < 1.0e-9) break;
	}
	*cmbindex = noc; // Index of core-mantle boundary

	//-------------------------------------------------------------------
	// Free mallocs
	//-------------------------------------------------------------------

	free(M);
	free(dM);
	free(rhonew);
	free(icomp);
	free(dbindex);
	free(eos);
	free(rho0);
	free(c);
	free(nn);
	free(Ks0);
	free(Ksp);
	free(V0);
	free(Tref);
	free(a0);
	free(a1);
	free(b0);
	free(b1);
	free(b2);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine planmat
 *
 * Reads database file of mass-radius relationships
 *
 *--------------------------------------------------------------------*/

int planmat(int ncomp, int **dbindex, int **eos, double **rho0, double **c, double **nn, double **Ks0, double **Ksp,
		double **V0, double **Tref, double **a0, double **a1, double **b0, double **b1, double **b2, char path[1024]) {

	FILE *fid;
	int i = 0;
	char str[1024];

	char *planmatdb = (char*)malloc(1024);       // Don't forget to free!
	planmatdb[0] = '\0';
	if (cmdline == 1) strncat(planmatdb,path,strlen(path)-20);
	else strncat(planmatdb,path,strlen(path)-18);
	strcat(planmatdb,"Data/Compression_planmat.txt");

	fid = fopen (planmatdb,"r");
	if (fid == NULL) printf("ExoCcycleGeo/Compression: Missing Compression_planmat.txt file.\n");
	else {
		printf("Reading planmat database file...\n");
		if (fgets(str, 1024, fid) != NULL) puts(str);
		if (fgets(str, 1024, fid) != NULL) puts(str);
		if (fgets(str, 1024, fid) != NULL) puts(str);
		for (i=0;i<ncomp;i++) {
			if (fgets(str, 1024, fid) != NULL) puts(str);
			int scan = fscanf(fid, "%d %d %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &(*dbindex)[i], &(*eos)[i],
						&(*rho0)[i], &(*c)[i], &(*nn)[i], &(*Ks0)[i], &(*Ksp)[i], &(*V0)[i], &(*Tref)[i], &(*a0)[i],
						&(*a1)[i], &(*b0)[i], &(*b1)[i], &(*b2)[i]);
			if (scan != 14) {                                                         // If scanning error
				printf("%d components found in planmat database\n",i);
				break;
			}
			else {
				printf("%d \t %d \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg \t %lg", (*dbindex)[i], (*eos)[i],
						(*rho0)[i], (*c)[i], (*nn)[i], (*Ks0)[i], (*Ksp)[i], (*V0)[i], (*Tref)[i], (*a0)[i],
						(*a1)[i], (*b0)[i], (*b1)[i], (*b2)[i]);
			}
			if (fgets(str, 1024, fid) != NULL) puts(str);
		}
	}

	fclose(fid);
	free(planmatdb);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine planmat_index
 *
 * Returns rank in database corresponding to input index
 *
 *--------------------------------------------------------------------*/

int planmat_index(int mat, int ncomp, int *dbindex) {
	int dbentry = 0;
	int i = 0;

	for (i=0;i<ncomp;i++) if (dbindex[i] == mat) dbentry = i;

	return dbentry;
}

#endif /* COMPRESSION_H_ */
