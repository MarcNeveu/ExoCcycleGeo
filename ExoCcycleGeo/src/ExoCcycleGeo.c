/*
 ============================================================================
 Name        : ExoCcycleGeo.c
 Author      : Marc Neveu
 Computes net C fluxes (in mol C m-2 s-1) at the surface-atmosphere interface
 of a terrestrial planet due to geophysical and geochemical processes.
 ============================================================================

 *  Copyright (C) 2017-2024 Marc Neveu (marc.f.neveu@nasa.gov)
 *
 *  This program is free software: you can redistribute it and/or modify it under the terms of the
 *  GNU General Public License as published by the Free Software Foundation, either version 3 of
 *  the License, or (at your option) any later version. This program is distributed in the hope
 *  that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 *  more details. You should have received a copy of the GNU General Public License along with this
 *  program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "ExoCcycleGeo.h"
#include "CHNOSZcmds.h"
#include "Compression.h"
#include "Geochem.h"
#include "Geodynamics.h"
#include <signal.h>

//-------------------------------------------------------------------
// MAIN PROGRAM
//-------------------------------------------------------------------

int main(int argc, char *argv[]) {

	int pid = getpid();

	//-------------------------------------------------------------------
	// Variable declarations and initializations
	//-------------------------------------------------------------------

	int i = 0;
	int recover = 0;
//	int t = 0;
//	int t_reset = 0;

	int n_inputs = 25;
	double *input = (double*) malloc(n_inputs*sizeof(double));
	if (input == NULL) printf("IcyDwarf: Not enough memory to create input[%d]\n", n_inputs);
	for (i=0;i<n_inputs;i++) input[i] = 0.0;

// USER-SPECIFIED INPUTS

	// User-specified grid parameters
	double dtime = 0.0;                // Time step (s)
	double dtime0 = 0.0;			   // Input time step (s)
	double dtime_old = 0.0;
	double tstart = 0.0;               // Time at which to start carbon cycling calculations (s)
	double tend = 0.0;                 // Time at which to end carbon cycling calculations (s)

	// User-specified planet interior parameters
	double m_p = 0.0;                  // Planet mass (kg)
	double Fe_FeMg = 0.0;			   // Fe/(Fe+Mg) ratio (mol mol-1)
	double Mg_Si = 0.0;                // Bulk Mg/Si ratio (mol mol-1)
	double FeO_Fe = 0.0;               // Bulk FeO/Fe ratio (mol mol-1)
	double C_OH = 1000.0;        	   // ppm H/Si in the mantle (Korenaga & Karato 2008, citing Hirth & Kohlstedt 1996)
    double magmaCmassfrac = 0.0;       // Mass fraction of C in magmas. Default 0.004 = 0.4±0.25% H2O and CO2 in MORB and OIB parent magmas (Jones et al. 2018; HArtley et al. 2014; Hekinian et al. 2000; Gerlach & Graeber 1985; Anderson 1995)
	int radionuclides = 0;             // 0 = Custom (LK07 lo), 1 = High (TS02), 2 = Intermediate (R91), 3 = Low (LK07), default = Intermediate (McDS95)
    double radioMult = 0.0;            // Multiplier of radionuclide abundances, 0.5-2 (Unterborn et al. 2020, p. 5-19)
    double fCH4 = 0.0;                 // Mole fraction of C outgassed as CH4, relative to CH4+CO2

    // User-specified planet surface parameters
	double Mocean = 0.0;               // Mass of ocean (kg, default: Earth=1.4e21)
	double L = 0.0;                    // Areal fraction of planet surface covered by land
	double Lnow = 0.0;                 // Present-day areal fraction covered by land
	double Tsurf0 = 0.0;               // Initial surface temperature (K)
	double Psurf = 0.0;                // Surface pressure (bar)
	double runoff = 0.0;               // Global mean river runoff rate (m s-1), default runoff_Earth = 0.67e-3 m day-1 (consistent, with 65% of evaporation, with Broeker & Peng 1982 p. 247 who have 0.7 m of rainfall per year)
	double runoff0 = 0.0;			   // Input global mean runoff rate (m s-1)
	double tResLand = 0.0;             // Residence time of water on continents (s)
	double tResLandNow = 0.0;		   // Present-day residence time of water on continents (s)
	double WplateRdg = 50000.0;        // Width of tectonic ridges and rifts (m)

    // User-specified planet atmosphere parameters
	double *xgas = (double*) malloc(nAtmSpecies*sizeof(double));
	if (xgas == NULL) printf("ExoCcycleGeo: Not enough memory to create xgas[nAtmSpecies]\n"); // Mixing ratios by volume (or by mol since all gases are pretty much ideal and have the same molar volume) of atmospheric gases
    for (i=0;i<nAtmSpecies;i++) xgas[i] = 0.0;

// OTHER QUANTITIES COMPUTED BY THE CODE

	int ir = 0;
	int j = 0;
	int NR = 10000;                    // Number of radial grid zones used in determining planetary structure at setup. 10000 seems to achieve numerical convergence.
	int iLith = 0;                     // Index of grid zone that corresponds to base of lithosphere (geotherm inflexion to mantle adiabat)
	int imin = 0;                      // Index of innermost grid zone of melting in MELTS output
	int imax = 0;                      // Index of outermost grid zone of melting in MELTS output
	int iCMB = 0;                      // Index of grid zone that corresponds to core-mantle boundary

//	int nslopeAvg = 3;                 // Target number of data points used for averaging melt fraction slope in the mantle
//	int islope = 0;                    // Actual numbers of data points used for averaging melt fraction slope in the mantle

//	double slope = 0.0;                // Average melt fraction slope in the mantle (bar-1)
	double meltfrac = 0.0;             // Melt fraction extrapolated at higher pressures than MELTS can handle

	double realtime = 0.0;             // Real time since birth of planetary system (s)

	// Planet parameters
	double m_c = 0.0;                  // Core mass (kg), default ≈0.3*m_p for Earth, Kite et al. (2009) use 0.325*m_p
	double fPx = 0.0;       		   // Fraction of upper mantle that is pyroxene (mol/mol-1, the rest is assumed to be olivine)
	int layerCode1 = 0;                // Inner layer material code used to compute a compressed planet structure (see Compression_planmat.txt database)
	int layerCode2 = 0;                // Middle layer material code (see Compression_planmat.txt database)
	int layerCode3 = 0;                // Outer layer material code (see Compression_planmat.txt database)
	double r_p = 0.0;                  // Planet radius (m)
	double r_c = 0.0;				   // Planet core radius (m)
	double rho_p = 0.0;                // Planet bulk density (kg m-3)
	double Asurf = 0.0;                // Planet surface area (m2)
	double Tsurf = 0.0;                // Surface temperature (K)

	double *r = (double*) malloc((NR+1)*sizeof(double)); // Radius (m)
	if (r == NULL) printf("Compression: Not enough memory to create r[NR+1]\n");
	for (ir=0;ir<NR+1;ir++) r[ir] = 0.0;

	double *rho = (double*) malloc((NR+1)*sizeof(double)); // Density (kg m-3)
	if (rho == NULL) printf("Compression: Not enough memory to create rho[NR]\n");
	for (ir=0;ir<NR+1;ir++) rho[ir] = 0.0;

	double *g = (double*) malloc((NR+1)*sizeof(double)); // Gravitational acceleration (m s-2)
	if (g == NULL) printf("Compression: Not enough memory to create g[NR]\n");
	for (ir=0;ir<NR+1;ir++) g[ir] = 0.0;

	double *P = (double*) malloc((NR+1)*sizeof(double)); // Pressure (Pa)
	if (P == NULL) printf("Compression: Not enough memory to create P[NR]\n");
	for (ir=0;ir<NR+1;ir++) P[ir] = 0.0;

	double *T = (double*) malloc((NR+1)*sizeof(double)); // Temperature in lithosphere and upper mantle convective boundary layer (K)
	if (T == NULL) printf("Compression: Not enough memory to create T[NR]\n");
	for (ir=0;ir<NR+1;ir++) T[ir] = 0.0;

	double *Meltfrac = (double*) malloc((NR+1)*sizeof(double*));
	if (Meltfrac == NULL) printf("ExoCcycleGeo: Not enough memory to create Meltfrac_geoth[%d]\n", NR); // Melt fraction by mass
	for (ir=0;ir<NR+1;ir++) Meltfrac[ir] = 0.0;

	// Reservoirs
//	double RCplate = 0.0;              // Plate/crust C reservoir (mol)
	double RCmantle = 0.0;             // Mantle C reservoir (mol)
	double RCmantle_old = 0.0;
	double RCatm = 0.0;                // Atmospheric C reservoir (mol)
	double RCocean = 0.0;              // Ocean C reservoir (mol)
	double RCatmoc = RCatm + RCocean;  // Combined atmospheric and ocean C reservoir (mol)
	double RCatmoc_old = RCatmoc;
	double RCorg = 0.0;			   	   // Organic carbon (biomass)
	double hChazeFallout = 0.0;	       // Cumulative thickness of haze fallout deposit (m)

	// Fluxes
	double FCoutgas = 0.0;             // C flux from outgassing (subaerial+submarine) (mol s-1)
	double FCcontW = 0.0;              // C flux from continental weathering (mol s-1)
	double FCseafsubd = 0.0;           // C flux from seafloor weathering and subduction (mol s-1)
	double netFC = 0.0;                // Net surface C flux from all geological processes (mol s-1)
	double farc = 0.0;                 // Fraction of subducted C that makes it back out through arc volcanism (as opposed to into the mantle)

	// Geochem parameters
	double pH = 0.0;                   // pH of the surface ocean (default 8.22)
	double pHout = 0.0;				   // pH output (calculated from ocean-atmosphere equilibrium)
	double rainpH = 0.0;               // pH of rainwater
	double pe = 0.0;                   // pe (-log activity e-) corresponding to logfO2
	double logfO2 = 0.0;               // log O2 fugacity
	double logKO2H2O = 0.0;            // log K for reaction 4 H+ + 4 e- + O2 = 2 H2O, from CHNOSZ: subcrt(c("H+","e-","O2","H2O"),c(-4,-4,-1,2),c("aq","aq","g","liq"),T=25,P=1)
	double Tfreeze = 273.15;           // Temperature at which water freezes at surface pressure (K)

	// Atmosphere parameters
	double nAir = 0.0;                 // Number of mol in atmosphere (mol)
	double S = 0.0;                    // Stellar flux at planet (W m-2)
	double S0 = 0.0;                   // Initial stellar flux at planet (W m-2)
	double albedo = 0.0;               // Average surface albedo (dimensionless)
	double Teff = 0.0;                 // Effective temperature (K)
	double DeltaTghe = 0.0;            // Temperature increase over effective temperature due to greenhouse effect (K)
	double psi = 0.0;                  // log10(pCO2) (bar)
	double Vatm1 = 0.0;                // Initial atmospheric volume (m3), determined at second time step to avoid skew of starting values
	double Tsurf1 = 0.0;               // Surface temperature used to determine Vatm1 (K)
	double Vatm = 0.0;                 // Atmospheric volume (m3)

	// Quantities to be computed by thermal/geodynamic model
	int staglid = 1;                   // 1 if stagnant-lid, 0 if mobile-lid
    double Tmantle = 0.0;              // Temperature of the upper mantle (K), initially set by accretion
	double H = 0.0;                    // Specific radiogenic heating rate (J s-1 kg-1)
	double x40K = 0.0;                 // Abundance of the radionuclide 40 K by mass
	double x232Th = 0.0;               // Abundance of the radionuclide 232 Th by mass
	double x235U = 0.0;                // Abundance of the radionuclide 235 U by mass
	double x238U = 0.0;                // Abundance of the radionuclide 238 U by mass
	double kappa = 0.0;                // Mantle thermal diffusivity (m2 s-1)
	double nu = 0.0;                   // Mantle kinematic viscosity (m2 s-1)
	double zCrust = 0.0;               // Crustal thickness (m), i.e. depth of layer crystallized from a melt
	double Ra = 0.0;                   // Rayleigh number for mantle convection (no dim)
	double Nu = 1.0;                   // Nusselt number (no dim)
	double Tref = 0.0;                 // Temperature at outer boundary of convective zone (surface or base of stagnant lid)
	double meltmass = 0.0;             // Total mass of outgassing mantle melt (kg)
	double rhomelt = 0.0;              // Density of the melt (kg m-3)
    double dMmelt_dt = 0.0;            // Rate of mantle mass melting (kg s-1)
	double dzCrust_dt = 0.0;           // Rate of crustal thickening due to mantle melting (m s-1)
	double zLith = 0.0;                // Depth of lithosphere (both thermal: geotherm inflexion to adiabat and mechanical: brittle-ductile transition) (m)
	double tConv = 0.0;                // Convection timescale (s), 200 Myr for deep Earth mantle today
	double vConv = 0.0;                // Convective velocity (m s-1)
	double p = 0.0;                    // Frank-Kamenetskii parameter (Solomatov 1995; denoted theta in Foley et al. 2020 chapter) (no dim)

	// Viscosity law constants (can't be #define'd because they are used as exponents)
//	double DryOlDiff[nParamFlowLaw];   // Dry diffusion creep flow law (Korenaga & Karato 2008, applicable for Fe/(Mg+Fe) = 0.1, modified for Fe content by Zhao et al. (2009))
//	DryOlDiff[0] = 261.0e3;            // Activation energy (J mol-1), default 261.0±28e3 (KK08)
//	DryOlDiff[1] = 6.0e-6;             // Activation volume (m3 mol-1), default 6±5e-6 (KK08)
//	DryOlDiff[2] = 5.25;               // Exponent of pre-exponential factor, default 5.25±0.03 (KK08)
//	DryOlDiff[3] = 1.0;                // Stress exponent (KK08)
//	DryOlDiff[4] = 2.98;               // Grain size exponent, default 2.98±0.02 (KK08)
//	DryOlDiff[5] = 0.0;                // Water content exponent, default 0.0 (this is a dry law) (KK08)
//	DryOlDiff[6] = -45.0e3;            // Activation energy dependence on Fe/(Fe+Mg) (Zhao et al. 2009 by analogy with gbs-accommodated dislocation creep, as did Ruedas et al. 2013 eq. 13)
//	DryOlDiff[7] = 1.5;                // Fe/(Fe+Mg) exponent (Zhao et al. 2009 by analogy with gbs-accommodated dislocation creep, as did Ruedas et al. 2013)

	double WetOlDiff[nParamFlowLaw];   // Wet diffusion creep flow law (Korenaga & Karato 2008, applicable for Fe/(Mg+Fe) = 0.1, modified for Fe content by Zhao et al. (2018))
	WetOlDiff[0] = 387.0e3;            // Activation energy (J mol-1), default 387±53e3
	WetOlDiff[1] = 25.0e-6;            // Activation volume (m3 mol-1), default 25±4e-6
	WetOlDiff[2] = 4.32;               // Exponent of pre-exponential factor, default 4.32±0.38
	WetOlDiff[3] = 1.0;                // Stress exponent
	WetOlDiff[4] = 2.56;               // Grain size exponent, default 2.56±0.24
	WetOlDiff[5] = 1.93;               // Water content exponent, default 1.93±0.07
	WetOlDiff[6] = -120.0e3;           // Activation energy dependence on Fe/(Fe+Mg) (Zhao et al. 2009 by analogy with dislocation creep)
	WetOlDiff[7] = 0.5;                // Fe/(Fe+Mg) exponent (Zhao et al. 2009  by analogy with dislocation creep)

//	double DryOlDisl[nParamFlowLaw];   // Dry dislocation creep flow law (Korenaga & Karato 2008, applicable for Fe/(Mg+Fe) = 0.1, modified for Fe content by Zhao et al. (2009))
//	DryOlDisl[0] = 610.0e3;            // Activation energy (J mol-1), default 610±30e3
//	DryOlDisl[1] = 13.0e-6;            // Activation volume (m3 mol-1), default 13±8e-6
//	DryOlDisl[2] = 6.09;               // Exponent of pre-exponential factor, default 6.09±0.11
//	DryOlDisl[3] = 4.94;               // Stress exponent, default 4.94±0.05
//	DryOlDisl[4] = 0.0;                // Grain size exponent
//	DryOlDisl[5] = 0.0;                // Water content exponent, default 0.0 (this is a dry law)
//	DryOlDisl[6] = -45.0e3;            // Activation energy dependence on Fe/(Fe+Mg) (Zhao et al. 2009, really for grain boundary sliding-accommodated dislocation)
//	DryOlDisl[7] = 1.5;                // Fe/(Fe+Mg) exponent (Zhao et al. 2009)

	double WetOlDisl[nParamFlowLaw];   // Wet dislocation creep flow law (Korenaga & Karato 2008, applicable for Fe/(Mg+Fe) = 0.1, modified for Fe content by Zhao et al. (2018))
	WetOlDisl[0] = 523.0e3;            // Activation energy (J mol-1), default 523±100e3
	WetOlDisl[1] = 4.0e-6;             // Activation volume (m3 mol-1), default 4±3e-6
	WetOlDisl[2] = 0.6;                // Exponent of pre-exponential factor, default 0.6±0.5
	WetOlDisl[3] = 3.60;               // Stress exponent, default 3.60±0.24
	WetOlDisl[4] = 0.0;                // Grain size exponent
	WetOlDisl[5] = 1.95;               // Water content exponent, default 1.95±0.05
	WetOlDisl[6] = -120.0e3;           // Activation energy dependence on Fe/(Fe+Mg) (Zhao et al. 2009)
	WetOlDisl[7] = 0.5;                // Fe/(Fe+Mg) exponent (Zhao et al. 2009)

	double flowLawDiff[nParamFlowLaw];
	for (i=0;i<nParamFlowLaw;i++) flowLawDiff[i] = 0.0;

	double flowLawDisl[nParamFlowLaw];
	for (i=0;i<nParamFlowLaw;i++) flowLawDisl[i] = 0.0;

	const double grainSize = 1.0e4;	   // Grain size for computation of viscosity and creep laws (µm)

//	printf("fPx \t Viscosity (Pa s)\n");
//	for (i=0;i<100;i++) {
//		printf("%g \t %g\n", (double)i/100.0, viscosity(1360.0+273.15, 0.3e9, WetOlDisl, grainSize, C_OH, Fe_FeMg, (double)i/100.0, 1.0/(2.0e-6)));
//	}
//	exit(0);
//
//	int tK = 1000;
//	printf("T(K) \t Viscosity(Pa s)\n");
//	for (tK=500;tK<2500;tK=tK+10) {
//		printf("%g \t %g\n", (double)tK, viscosity(tK, 3.0e9, WetOlDiff, grainSize, C_OH, Fe_FeMg, fPx, 0.1*Gyr2sec));
//	}
//	exit(0);

	// Quantities to be computed by continental weathering model
    int iResTime = 0;
	int kinsteps = 1000;                  // Number of time steps of PHREEQC kinetic simulation

	double **xriver = (double**) malloc(kinsteps*sizeof(double*)); // Molalities of aqueous species and moles (remaining or consumed since last step) of minerals at different times (mol (kg H2O)-1 or mol, xriver[][1]: time in s)
	if (xriver == NULL) printf("ExoCcycleGeo: Not enough memory to create xriver[kinsteps]\n");
    for (i=0;i<kinsteps;i++) {
    	xriver[i] = (double*) malloc(nvarKin*sizeof(double));
    	if (xriver[i] == NULL) printf("ExoCcycleGeo: Not enough memory to create xriver[kinsteps][nvarKin]\n");
    }
	for (i=0;i<kinsteps;i++) {
		for (j=0;j<nvarKin;j++) xriver[i][j] = 0.0;
	}

	double kintime = 0.0;              // Total time of PHREEQC kinetic simulation (s)
	double WRcontW = 0.0;              // Water:rock mass ratio for continental weathering reactions (kg/kg)

	double fracEvap = 0.65;            // Fraction of rain/river water that evaporates before reaching the ocean (Berner et al. 1983)
	double massH2Oriver = 0.0;         // Mass of water after residence time elapsed in PHREEQC calculation
	double Mriver = 0.0;               // Mass of river water delivered to the ocean after a time step

//	double xriver_Mg_evap = 0.0;       // River abundance of Mg accounting for evaporation (mol/kg)
//	double xriver_Ca_evap = 0.0;       // River abundance of Ca accounting for evaporation (mol/kg)
//	double xriver_Si_evap = 0.0;       // River abundance of Si accounting for evaporation (mol/kg)
//	double xriver_Na_evap = 0.0;       // River abundance of Na accounting for evaporation (mol/kg)
//	double xriver_Fe_evap = 0.0;       // River abundance of Fe accounting for evaporation (mol/kg)

	double Mg_carb_consumed = 0.0;     // Amount of crustal magnesium carbonate dissolved (mol)
	double Ca_carb_consumed = 0.0;     // Amount of crustal calcium carbonate dissolved (mol)
	double Ca_sulf_consumed = 0.0;     // Amount of crustal calcium sulfate dissolved (mol)
	double Fe_sulf_consumed = 0.0;     // Amount of crustal iron sulfide dissolved (mol)

	double FC_Mg = 0.0;                // Contribution to weathering flux (HCO3- trapping capacity) from Mg ions (mol/s)
	double FC_Mg_carb = 0.0;           // Contribution to weathering flux (HCO3- trapping capacity) from Mg carbonate dissolution (mol/s)
	double FC_Mg_sil = 0.0;            // Contribution to weathering flux (HCO3- trapping capacity) from Mg silicate dissolution (mol/s)

	double FC_Ca = 0.0;                // Contribution to weathering flux (HCO3- trapping capacity) from Ca ions (mol/s)
	double FC_Ca_carb = 0.0;           // Contribution to weathering flux (HCO3- trapping capacity) from Ca carbonate dissolution (mol/s)
	double FC_Ca_sulf = 0.0;           // Contribution to weathering flux (HCO3- trapping capacity) from Ca sulfate dissolution (mol/s)
	double FC_Ca_sil = 0.0;            // Contribution to weathering flux (HCO3- trapping capacity) from Ca silicate dissolution (mol/s)

	double FC_Fe = 0.0;                // Contribution to weathering flux (HCO3- trapping capacity) from Fe ions (mol/s)
	double FC_Fe_sulf = 0.0;           // Contribution to weathering flux (HCO3- trapping capacity) from Fe sulfide dissolution (mol/s)
	double FC_Fe_sil = 0.0;            // Contribution to weathering flux (HCO3- trapping capacity) from Fe silicate dissolution (mol/s)


	// Quantities to be computed by seafloor weathering model
	double LplateRdg = 0.0;            // Length of tectonic plate boundaries, today 6.5*Earth circumference (m)
	double zCrack = 6000.0;            // Depth of fracturing below seafloor (m), default 6000 m (Vance et al. 2007)
	double Pseaf = 0.0;                // Avg. seafloor pressure calculated assuming avg. ocean depth
	double WRseafW = 1.0;              // Water to rock ratio of seafloor weathering, here multiplied by rock mass
	double mix = 0.0;                  // Mocean/Mriver, optionally scaled by sub-timestep
	double deltaCreac = 0.0;           // Net C leached/precipitated per kg of rock (mol kg-1) over one sub-timestep
	double xaq0 = 0.0;                 // Memory of xaq[0] (mol/kg) before seafloor weathering calculation to avoid removing oxidized C twice (during PHREEQC calculation and from netFC)
	double xaq1 = 0.0;                 // Memory of xaq[1] (mol/kg) before seafloor weathering calculation to avoid removing reduced C twice (during PHREEQC calculation and from netFC)
	double volSeafCrust = 0.0;         // Volume of seafloor reacting with ocean water (m3)

	double Tadiab[NR]; // Mantle adiabat temperature profile (K) TODO dynamically allocate memory
	for (i=0;i<NR;i++) Tadiab[i] = 0.0;

	double *xaq = (double*) malloc(nAqSpecies*sizeof(double));
	if (xaq == NULL) printf("ExoCcycleGeo: Not enough memory to create xaq[nAqSpecies]\n"); // Molalities of aqueous species (mol (kg H2O)-1)
    for (i=0;i<nAqSpecies;i++) xaq[i] = 0.0;

	double **sys_tbl = (double**) malloc(NR*sizeof(double*));
	if (sys_tbl == NULL) printf("ExoCcycleGeo: Not enough memory to create sys_tbl[%d]\n", NR); // Storage of System_main_tbl.txt alphaMELTS output
	for (i=0;i<NR;i++) {
		sys_tbl[i] = (double*) malloc(18*sizeof(double));
		if (sys_tbl[i] == NULL) printf("ExoCcycleGeo: Not enough memory to create sys_tbl[%d][18]\n", NR); // 18 columns in System_main_tbl.txt
	}
	for (i=0;i<NR;i++) {
		for (j=0;j<18;j++) sys_tbl[i][j] = 0.0;
	}

	FILE *fout;
	char *title = (char*)malloc(1024*sizeof(char)); title[0] = '\0';

	//-------------------------------------------------------------------
	// Startup
	//-------------------------------------------------------------------

	// Initialize the R environment. We do it here, in the main loop, because this can be done only once.
	setenv("R_HOME","/Library/Frameworks/R.framework/Resources",1); // Specify R home directory
	Rf_initEmbeddedR(argc, argv);                                   // Launch R
	CHNOSZ_init(1);                                                 // Launch CHNOSZ

	// Get current directory. Works for Mac only! To switch between platforms, see:
	// http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe
	char path[1024];
	unsigned int size = sizeof(path);
	path[0] = '\0';
	if (_NSGetExecutablePath(path, &size) == 0) printf("\n");
	else printf("ExoCcycleGeo: Couldn't retrieve executable directory. Buffer too small; need size %u\n", size);

	// Set up path of executable in alphaMELTS files
	if (alphaMELTS_init(path) == 1) {
		printf("ExoCcycleGeo: Couldn't successfully set up path of executable in alphaMELTS files. Exiting.\n");
		exit(0);
	}
	else printf("Successfully set up path of executable in alphaMELTS files.\n");

	// Read input file
	input = exoCinput(input, path);

	i = 0;
	// Grid inputs
	dtime0 = input[i]*Myr2sec; i++;      // Timestep
	tstart = input[i]*Gyr2sec; i++;      // Simulation start time
	tend = input[i]*Gyr2sec; i++;        // Simulation end time
	// Interior inputs
	m_p = input[i]*mEarth; i++;
	Fe_FeMg = input[i]; i++;             // Fe/Mg, to be converted into Fe/(Fe+Mg) ratio
	Mg_Si = input[i]; i++;               // Mg/Si, to be converted into pyroxene/(pyroxene + olivine)
	FeO_Fe = input[i]; i++;              // Core mass, TODO convert from FeO/Fe
	C_OH = input[i]; i++;                // Mantle H2O (ppm by mass)
	magmaCmassfrac = input[i]; i++;      // Mass fraction of C in magmas; range 0.03-1 wt.%, higher than mass fraction of C in upper mantle, which is 350 ppm; range 115-670 ppm (Aiuppa et al. 2021)
	radionuclides = (int) input[i]; i++; // 0 = Custom (LK07 lo), 1 = High (TS02), 2 = Intermediate (R91), 3 = Low (LK07), default = Intermediate (McDS95)
	radioMult = input[i]; i++;           // Radionuclide content multiplier
	staglid = (int) input[i]; i++;       // 0 = plate tectonics, 1 = stagnant lid
    fCH4 = input[i]; i++;				 // Mole fraction of C outgassed as CH4, relative to CH4+CO2
    // Surface inputs
	Mocean = input[i]; i++;              // Mass of ocean (kg, default: Earth=1.4e21)
	Lnow = input[i]; i++;                // Areal fraction of planet surface covered by land
	Tsurf0 = input[i]; i++;              // Initial surface temperature (K)
	Psurf = input[i]; i++;               // Surface pressure (bar)
	runoff0 = input[i]/86400.0; i++;     // River runoff (m s-1)
	tResLandNow = input[i]*Yr2sec; i++;  // Residence time of rain- and river water on continents (s)
	// Atmospheric inputs
	xgas[0] = input[i]; i++;             // CO2 mixing ratio
    xgas[1] = input[i]; i++;             // CH4 mixing ratio
    xgas[2] = input[i]; i++;             // O2 mixing ratio
    xgas[3] = input[i]; i++;             // N2 mixing ratio
    xgas[4] = input[i]; i++;             // H2O mixing ratio

	printf("\n");
	printf("ExoCcycleGeo v26.2\n");
	printf("Copyright (C) 2017-2026 Marc Neveu (marc.f.neveu@nasa.gov)\n\n"       );
	printf("This program is free software: you can redistribute it and/or\n"      );
	printf("modify it under the terms of the GNU General Public License as\n"     );
	printf("published by the Free Software Foundation, either version 3 of the\n" );
	printf("License, or (at your option) any later version. This program is\n"    );
	printf("distributed in the hope that it will be useful, but without any\n"    );
	printf("warranty. See the GNU General Public License for more details:\n"     );
	printf("<http://www.gnu.org/licenses/>.\n\n"                                  );
	if (cmdline == 1) printf("Command line mode\n");

	printf("|--------------------------------------------------------------|\n");
	printf("| Grid |||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
	printf("|-----------------------------------------------|--------------|\n");
	printf("| Simulation time step (Myr, default 1, adapts) | %g \n", dtime0/Myr2sec);
	printf("| Simulation start time (Gyr, default 0.6)      | %g \n", tstart/Gyr2sec);
	printf("| Simulation end   time (Gyr, default 5)        | %g \n", tend/Gyr2sec);
	printf("|-----------------------------------------------|--------------|\n");
	printf("| Planet Interior ||||||||||||||||||||||||||||||||||||||||||||||\n");
	printf("|-----------------------------------------------|--------------|\n");
	printf("| Planet mass (Earth masses, min 0.5, max 2)    | %g \n", m_p/mEarth);
	printf("|        Fe/Mg (0.6 to 1.5, Earth 1.0)          | %g \n", Fe_FeMg);
	printf("|        Mg/Si (0.5 to 2, Earth 0.95)           | %g \n", Mg_Si);
	printf("|        FeO/(FeO+Fe) (0 to 0.2, Earth 0.08)    | %g \n", FeO_Fe);
	printf("| Mantle H2O (ppm by mass, Earth 300-1000)      | %g \n", C_OH);
	printf("| Magma C (ppm by mass, Earth 300-10000)        | %g \n", magmaCmassfrac);
	printf("| Radionuclides (1 hi, 2 int, 3 lo, def. int)   | %d \n", radionuclides);
	printf("| Radionuclide content multiplier (0.5-2)       | %g \n", radioMult);
	printf("| Tectonic mode (0 plate tecton, 1 stagnant lid)| %d \n", staglid);
	printf("| Mole fraction of C outgassed as CH4/(CH4+CO2) | %g \n", fCH4);
	printf("|-----------------------------------------------|--------------|\n");
	printf("| Surface ||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
	printf("|-----------------------------------------------|--------------|\n");
	printf("| Mass of surface ocean (kg, modrn Earth 1.4e21)| %g \n", Mocean);
	printf("| Areal land fraction (default 0.29)            | %g \n", Lnow);
	printf("| Initial temperature (K, default 288)          | %g \n", Tsurf0);
	printf("| Initial pressure (bar)                        | %g \n", Psurf);
	printf("| Runoff rate (m/day, default 0.67e-3)          | %g \n", runoff0*86400.0);
	printf("| Water residence time, continents (yr, def. 10)| %g \n", tResLandNow/Yr2sec);
	printf("|-----------------------------------------------|--------------|\n");
	printf("| Initial Atmosphere |||||||||||||||||||||||||||||||||||||||||||\n");
	printf("|-----------------------------------------------|--------------|\n");
	printf("| CO2 (default 280e-6)                          | %g \n", xgas[0]);
	printf("| CH4                                           | %g \n", xgas[1]);
	printf("| O2  (default 0.2)                             | %g \n", xgas[2]);
	printf("| N2  (default 0.79)                            | %g \n", xgas[3]);
	printf("| H2O                                           | %g \n", xgas[4]);
	printf("|--------------------------------------------------------------|\n");

	printf("\n");
	
	// Core mass fraction, based on Fe/Mg and FeO/Fe
	m_c = m_p * (0.25 + (0.45-0.25)/(1.5-0.6)*(Fe_FeMg-0.6)); // With Fe/FeO = Earth, m_c/m_p = 0.25 for Fe/Mg = 0.6, 0.45 for Fe/Mg = 1.5 (Unterborn and Panero 2019), 0.325 canonical (Fe/Mg = 1). Linear relationship gives m_c/m_p = 0.338 for Fe/Mg = 1 and 0.325 for Fe/Mg = 0.94
	
	/* Earth has 6-8% FeO, 33% CMF
       Mars has 18% FeO, 15-20% CMF. CMF actually closer to 25% and so FeO likely lower (Stahler et al. 2021 Science).
       Mercury has no FeO, 70% CMF.
       Little is known about Venus’ core. But CMF (Fe/FeO) can’t be decided by present-day oxidation state which is much more oxidized. So, keep the two independent. 
       (Unterborn et al. 2020 p. 5-31) */
    if (FeO_Fe < 0.0 || FeO_Fe > 0.15) {
		printf ("ExoCcycleGeo: FeO fraction = %g is not between 0 and 0.15. Adjust accordingly and restart\n", fPx);
		exit(0);
	}
    m_c = (m_c / 0.338) * (0.7 + (0.33-0.7)/0.08*FeO_Fe); // = (m_c / m_c_Earth) * (FeO-based modifier). Linear modifier pegged at Earth and Mercury, this gives Mars 15% CMF for 12% FeO
	printf("Core mass fraction = %.3g\n", m_c/m_p);
	
	// Set planmat core structure layers
	layerCode1 = 101; // Default 101, alpha (bcc) Fe
	layerCode2 = 107; // Default 107, liquid Fe + 10% S
	
	// Set mantle layerCode3 based on composition. 203 (MgO) provides the best structural match to modern Earth (because crust has extracted Si?) but is not an option below.
	
	// Adjust Fe/Mg of the mantle = FeO/Mg, assuming all pure Fe is in the core and all FeO in the mantle. Fe and silicates stay well-mixed if Fe is oxidized (Foley et al. 2020; Elkins-Tanton & Seager 2008)
	Fe_FeMg = FeO_Fe*Fe_FeMg;
	printf("Mantle Fe/Mg %g (Primitive Mantle is 0.12);\t Mg/Si %g (PM is 1.2, higher due to segregation of Si in the core)\n", Fe_FeMg, Mg_Si); // Values lower than primitive mantle

	// Pyroxene fraction in upper mantle
	fPx = 2.0 - (Fe_FeMg + 1.0)*Mg_Si;     // Conversion from Mg/Si to fPx = (Fe,Mg)SiO3 / ((Fe,Mg)SiO3 + (Fe,Mg)2SiO4). 
	                                       // Since (Fe+Mg)/Si = (Fe/Mg) (Mg/Si) + (Mg/Si) = 1 for fPx = 1 (all pyroxene), 2 for fPx = 0 (all olivine), then fPx = 2 - (Fe+Mg)/Si = 2 - ((Fe/Mg) + 1) (Mg/Si)
	if (fPx < 0.0 || fPx > 1.0) {
		printf ("ExoCcycleGeo: Pyroxene fraction = %g is not between 0 and 1. Adjust Fe/Mg and Mg/Si accordingly and restart\n", fPx);
		exit(0);
	}
	
	if (Mg_Si*(1.0+Fe_FeMg) < 0.5) { // Stishovite SiO2 dominates, unreachable in practice since this corresponds to fPx > 1
		layerCode3 = 204; // Stishovite SiO2, high pressure		
		printf ("Mantle material code in Data/Compression_planmat: %d, Stishovite SiO2, high pressure\n", layerCode3);
	}
	else {
		if (fPx > 0.5) { // Pyroxene dominates
			if (Fe_FeMg < 0.25) {
				layerCode3 = 201; // MgSiO3 perovskite, Earth-like case, reachable if FeO/Fe = 0.15, Fe/Mg = 0.6-1.5 and Mg/Si = 0.9-0.95
				printf ("Mantle material code in Data/Compression_planmat: %d, MgSiO3 perovskite\n", layerCode3);
			}
			else {
				layerCode3 = 202; // (Mg,Fe)SiO3, reachable if FeO/Fe = 0.15, Fe/Mg = 7.5 and Mg/Si = 0.5-0.6
				printf ("Mantle material code in Data/Compression_planmat: %d, (Mg,Fe)SiO3\n", layerCode3);
			}
		}
		else { // Olivine dominates
			if (Fe_FeMg < 0.25) {
					layerCode3 = 205; // Olivine, forsterite Mg2SiO4 0.1-760 MPa, reachable if FeO/Fe = 0.15, Fe/Mg = 0.6-1.6 and Mg/Si = 1.5-1.8
					printf ("Mantle material code in Data/Compression_planmat: %d, Olivine, forsterite Mg2SiO4 0.1-760 MPa\n", layerCode3);
			}
			else if (Fe_FeMg < 0.75) {
					layerCode3 = 206; // Olivine, 50% Fo/50% Fa MgFeSiO4, 0.1-760 MPa, reachable if FeO/Fe = 0.15, Fe/Mg = 1.7-2.5 and Mg/Si = 1.4-1.5
					printf ("Mantle material code in Data/Compression_planmat: %d, Olivine, 50%% Fo/50%% Fa MgFeSiO4, 0.1-760 MPa\n", layerCode3);
			}
			else {
				layerCode3 = 207; // Olivine, fayalite Fe2SiO4 0.1-760 MPa, reachable if FeO/Fe = 0.15, Fe/Mg = 8.5 and Mg/Si = 0.7
				printf ("Mantle material code in Data/Compression_planmat: %d, Olivine, fayalite Fe2SiO4 0.1-760 MP\n", layerCode3);
			}
		}
	}
	
	Fe_FeMg = Fe_FeMg/(Fe_FeMg+1.0);       // Conversion from Fe/Mg to Fe/(Fe+Mg) = Fe/Mg / (Fe/Mg + Mg/Mg)
	
	// Conversion from ppm of magma C mass fraction
	magmaCmassfrac = magmaCmassfrac/1.0e6; // Conversion from ppm
	
	printf("fPx = %g; Fe/(Fe+Mg) = %g\n", fPx, Fe_FeMg);
	
	char infile[2048];
	char *meltsfile = (char*) malloc(2048*sizeof(char));
	
	infile[0] = '\0';
	meltsfile[0] = '\0';
	
	if (cmdline == 1) strncat(infile ,path, strlen(path)-20);
	else strncat(infile, path, strlen(path)-18);
	strcat(infile,"alphaMELTS-1.9/ExoC/ExoCcycleGeo_primMantleMS95.melts");
	
	if (cmdline == 1) strncat(meltsfile ,path, strlen(path)-20);
	else strncat(meltsfile, path, strlen(path)-18);
	strcat(meltsfile,"alphaMELTS-1.9/ExoC/ExoCcycleGeo.melts");

	WriteMELTSinput(infile, Fe_FeMg, Mg_Si, &meltsfile);

	free (meltsfile);

	int timesteps = (int) ((tend-tstart)/dtime0);
	if (timesteps <= 0) timesteps = 1;
//	double FCoutgas_mem[timesteps][2];        // Memorized outgassing flux in case MELTS calculations fluctuate too much (mol s-1) TODO this only works if timestep doesn't change, but it will
//	for (i=0;i<timesteps;i++) FCoutgas_mem[i][0] = 0.0; FCoutgas_mem[i][1] = 0.0;

	printf("Computing geo C fluxes through time...\n");

	// Get pressure and density profiles with depth, accounting for compression and self-gravity
	compression(NR, m_p, m_c, Tsurf, layerCode1, layerCode2, layerCode3, &r, &P, &rho, &g, &iCMB, path);

	r_p = r[NR];
	r_c = r[iCMB];
	kappa = k/(rho[(int)(NR+iCMB)/2]*Cp);
	rho_p = m_p / (4.0/3.0*PI_greek*pow(r_p,3));

	printf("\nPlanet radius = %g km, Core radius = %g km, Mid-mantle density = %g kg m-3, avg density = %g kg m-3\n", r_p/km2m, r_c/km2m, rho[(NR+iCMB)/2], m_p/(4.0/3.0*PI_greek*pow(r_p,3.0)));
	Tsurf = Tsurf0;
	runoff = runoff0;
	Asurf = 4.0*PI_greek*r_p*r_p;
	nAir = Psurf*bar2Pa*Asurf/g[NR]/molmass_atm(xgas);

	zLith = (r_p-r_c)/Nu;              // Initial value, thickness of entire mantle since Nu is initiated to 1
	for (i=0;i<NR;i++) {
		if (r[i] <= r_p - zLith && r[i+1] > r_p-zLith) iLith = i;
	}
	tConv = 10.0*Myr2sec;              // Initial value
	vConv = (r_p-r_c)/tConv;

	// Calculate Tmantle from accretion (Canup et al. 2021 Pluto chapter; volume averaging leads to 3/5 term)
	Tmantle = Tsurf0 + 0.6*4.0*h_frac*chi*G*PI_greek*rho_p*r_p*r_p/(3.0*Cp);

	// Calculate temperatures (not computed in the core)
	for (i=0   ;i<iCMB;i++) T[i] = 0.0; // Temperatures not computed in the core
	for (i=iCMB;i<NR  ;i++) {
//		if (r[i] < r_p - zLith) T[i] = Tmantle - alpha*g[(int)(0.5*(NR+iCMB))]*Tmantle/Cp * (r[i] - (r_p - zLith));
		if (r[i] < r_p - zLith) T[i] = Tmantle*exp( alpha*g[i]*(r_p - r[i] - zLith)/Cp );
		else T[i] = Tsurf + (Tmantle-Tsurf)*(r_p-r[i])/zLith; // Linear profile (Foley et al. 2020)
	}

	// Radionuclides
	switch (radionuclides) {
	case 0: // Custom: LK07 (case 3), low endmember
		x40K   =  18.0e-9;
		x232Th =  51.9e-9;
		x235U  =  0.10e-9;
		x238U  =  14.2e-9;
		break;
	case 1: // Abundances from Turcotte & Schubert (2002) as reported in Table 1 of Kite et al. (2009). Highest heating.
		x40K   =  36.9e-9;
		x232Th = 124.0e-9;
		x235U  =  0.22e-9;
		x238U  =  30.8e-9;
		break;
	case 2: // Abundances derived from Ringwood (1991) as reported in Table 1 of Kite et al. (2009). Intermediate heating.
		x40K   =  30.7e-9;
		x232Th =  84.1e-9;
		x235U  =  0.15e-9;
		x238U  =  21.0e-9;
		break;
	case 3: // Abundances from Table 2 of Lyubetskaya & Korenaga 2007. Lowest heating.
		x40K   =  22.8e-9; // = 0.012/100.0 (isotopic abundance) *190.0±40e3*1.0e-9
		x232Th =  62.6e-9; // = 0.9998 (isotopic abundance) *62.6±10.7e-9
		x235U  =  0.12e-9; // = 0.720/100.0 (isotopic abundance) *17.3±3.0e-9
		x238U  =  17.2e-9; // = 0.99274 (isotopic abundance) *17.3±3.0*1.0e-9
		break;
	default: // Abundances of McDonough & Sun (1995) as reported in Table 2 of Lyubetskaya & Korenaga 2007. Intermediate heating similar to R91.
		x40K   =  28.8e-9; // = 0.012/100.0 (isotopic abundance) *240.0±48e3*1.0e-9
		x232Th =  79.5e-9; // = 0.9998 (isotopic abundance) *79.5±11.9e-9
		x235U  =  0.15e-9; // = 0.720/100.0 (isotopic abundance) *20.3±4.1e-9
		x238U  = 20.15e-9; // = 0.99274 (isotopic abundance) *20.3±4.1*1.0e-9
		break;
	}

	// Rheology
//	switch (rheology) {
//	case 0: // Dry olivine rheology of Korenaga & Karato (2008)
//		for (i=0;i<6;i++) {
//			flowLawDiff[i] = DryOlDiff[i];
//			flowLawDisl[i] = DryOlDisl[i];
//		}
//		break;
//	case 1: // Wet olivine rheology of Korenaga & Karato (2008)
		for (i=0;i<6;i++) {
			flowLawDiff[i] = WetOlDiff[i];
			flowLawDisl[i] = WetOlDisl[i];
		}
//		break;
//	default:
//		printf("ExoCcycleGeo: Choice of rheology is not recognized. Exiting.\n");
//		exit(0);
//	}

	// Redox (from most oxidized to most reduced)
	int redox = 1; // TODO deprecated, remove
	switch(redox) {
	case 1: // Present-day Earth surface
		printf("Redox set to present-day Earth surface\n");
		logfO2 = log(0.2)/log(10.0);
		break;
	case 2: // Use CHNOSZ to get log fO2 for hematite-magnetite (HM) buffer at given T and P.
		printf("Redox set to HM buffer\n");
		logfO2 = -6.0*CHNOSZ_logK("hematite", "cr", Tsurf - Kelvin, Psurf, "SUPCRT92")
			     +4.0*CHNOSZ_logK("magnetite", "cr", Tsurf - Kelvin, Psurf, "SUPCRT92")
			     +1.0*CHNOSZ_logK("O2", "g", Tsurf - Kelvin, Psurf, "SUPCRT92");
		break;
	case 3: // Use CHNOSZ to get log fO2 for fayalite-magnetite-quartz (FMQ) buffer at given T and P.
		printf("Redox set to FMQ buffer\n");
		logfO2 = -3.0*CHNOSZ_logK("quartz", "cr", Tsurf - Kelvin, Psurf, "SUPCRT92")
			     -2.0*CHNOSZ_logK("magnetite", "cr", Tsurf - Kelvin, Psurf, "SUPCRT92")
		         +3.0*CHNOSZ_logK("fayalite", "cr", Tsurf - Kelvin, Psurf, "SUPCRT92")
			     +1.0*CHNOSZ_logK("O2", "g", Tsurf - Kelvin, Psurf, "SUPCRT92");
		break;
	case 4: // Use CHNOSZ to get log fO2 for iron-wustite (IW) buffer at given T and P.
		printf("Redox set to IW buffer\n");
		logfO2 = +2.0*CHNOSZ_logK("Fe", "cr", Tsurf - Kelvin, Psurf, "SUPCRT92")
			     -2.0*CHNOSZ_logK("FeO", "cr", Tsurf - Kelvin, Psurf, "SUPCRT92")
			     +1.0*CHNOSZ_logK("O2", "g", Tsurf - Kelvin, Psurf, "SUPCRT92");
		break;
	default:
		printf("ExoCcycleGeo: Redox switch incorrectly specified, should be 1 (current Earth surface), 2 (hematite-magnetite), "
				"3 (fayalite-magnetite-quartz), or 4 (iron-wustite). Exiting.\n");
		exit(0);
	}
	logKO2H2O = -4.0*CHNOSZ_logK("H+", "aq", Tsurf - Kelvin, Psurf, "SUPCRT92")
				-4.0*CHNOSZ_logK("e-", "aq", Tsurf - Kelvin, Psurf, "SUPCRT92")
				-1.0*CHNOSZ_logK("O2", "g", Tsurf - Kelvin, Psurf, "SUPCRT92")
				+2.0*CHNOSZ_logK("H2O", "liq", Tsurf - Kelvin, Psurf, "SUPCRT92");

	printf("log f(O2) = %g\n", logfO2);
	pe = -pH + 0.25*(logfO2+logKO2H2O);

	//-------------------------------------------------------------------
	// Initialize parameters not specified by user
	//-------------------------------------------------------------------

	// Radiative transfer
	S0 = 1368.0; // W m-2
	albedo = 0.3;

    // Ocean
	pH = 8.22;

	//-------------------------------------------------------------------
	// Initialize reservoirs
	//-------------------------------------------------------------------

	if (!recover) { // New simulation start
		RCatm = (xgas[0]+xgas[1])*nAir;

		// Equilibrate ocean and atmosphere at input pressure and atmospheric C/N ratio
		if (Psurf < 0.01) {
			printf("ExoCcycleGeo: Pressure = %g bar too close to the triple point of H2O, oceans not stable at the surface. Exiting.\n", Psurf);
			exit(0);
		}

		printf("Equilibrating ocean and atmosphere at input pressure and atmospheric C/N ratio...\n");
		for (i=0;i<nAtmSpecies;i++) {
			if (xgas[i] > 0.0 && xaq[i] == 0.0) xaq[i] = 1.0; // xaq must be >0 otherwise PHREEQC ignores it, set to 1 mol/kgw (initial guess).
		}

		AqueousChem(path, "io/OceanStart.txt", Tsurf, &Psurf, &Vatm, &nAir, &pH, &pe, &Mocean, &xgas, &xaq, &xriver, 0, 1, 0.0, 1, nvarEq, 0.0, 0.0, &deltaCreac, staglid, dtime);

		RCocean = (xaq[0]+xaq[1])*Mocean;
		RCatmoc = RCatm + RCocean;
		RCmantle = 1.0e4*(RCatmoc); // Dasgupta and Hirschmann (2010): 0.8 to 12.5 e23 g = 0.07 to 1.04 e23 mol
		RCmantle = (m_p-m_c)*magmaCmassfrac/10.0/0.012; // Dasgupta and Hirschmann (2010): parent mantle of magma is depleted in C relative to magma, has 20-1800 ppm. Aiuppa et al. (2021) have 117-669 ppm. Choosing 200 ppm by default (0.002/10)

		// Print first line of outputs
		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"Outputs/ReservoirsFluxes.txt");
		fout = fopen(title,"w");
		if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
		else {
			fprintf(fout, "'Reservoirs in mol, fluxes in mol s-1'\n");
			fprintf(fout, "'Time (Gyr)' \t 'Tmantle (K)' \t RCmantle \t RCatm \t RCocean \t RCatm+RCoc \t RCatmoc \t FCoutgas \t FCcontw \t FCseafsubd \t 'Net C flux' \t "
					"'Total N g+aq' \t RCorg\n");
			fprintf(fout, "0 \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n",
					Tmantle, RCmantle, RCatm, RCocean, RCatm+RCocean, RCatmoc, FCoutgas, FCcontW, FCseafsubd, netFC, xaq[3]*Mocean + xgas[3]*2.0*nAir, RCorg);
		}
		fclose (fout);

		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"Outputs/CompoAtmosph.txt");
		fout = fopen(title,"w");
		if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
		else {
			fprintf(fout, "'Atmospheric species in mixing ratio (by mol, i.e. by volume for ideal gases), total C and N in mol'\n");
			fprintf(fout, "'Time (Gyr)' \t CO2(g) \t CH4(g) \t O2(g) \t N2(g) \t H2O(g) \t 'P_surface (bar)' \t 'T_surface (K)' \t 'DeltaT GHE (K)' \t 'Mean rainfall (m/yr) \t 'nAir (mol)'\n");
			fprintf(fout, "0 \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n",
					xgas[0], xgas[1], xgas[2], xgas[3], xgas[4], Psurf, Tsurf, DeltaTghe, runoff/(1-fracEvap)*Yr2sec, nAir);
		}
		fclose (fout);

		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"Outputs/Geotherm.txt");
		fout = fopen(title,"w");
		if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
		fclose (fout);

		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"Outputs/Outgassing.txt");
		fout = fopen(title,"w");
		if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
		else fprintf(fout, "Time (Gyr) \t Tmantle (K) \t Ra \t Visc (Pa s) \t Heat flux (mW m-2) \t Lithospheric thickness (km) \t Frank-Kamenetskii parameter \t Outgassing flux (mol C s-1) \t Convective velocity (m yr-1) \t "
				"Convective timescale (s) \t Crustal thickness (m) \t Crust generation timescale (s) \t Stagnant lid?\n");
		fclose (fout);

		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"Outputs/ContWeatherFluxes.txt");
		fout = fopen(title,"w");
		if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
		else {
			fprintf(fout, "'Continental weathering fluxes in mol s-1 = bicarbonate trapping capacity due to dissolved cations'; river abundances in mol kg-1\n");
			fprintf(fout, "'Time (Gyr)' \t 'Land areal fraction' \t 'Residence time of water on land' \t 'Total Mg' \t 'Mg from silicates' \t 'Mg from carbonates' \t 'Total Ca' \t 'Ca from silicates' \t 'Ca from carbonates' \t 'Ca from sulfates'"
					" \t 'Total Fe' \t 'Fe from silicates' \t 'Fe from sulfides' \t 'Total C flux' \t 'River pH' \t HCO3-(aq) \t CO3-2(aq) \t CO2(aq) \t Na+(aq) \t Mg+2(aq) \t SiO2(aq) \t Ca+2(aq) \t Fe+2(aq) \t 'Haze fallout height (m)'\n");
			fprintf(fout, "0 \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n",
					L, tResLand, FC_Mg, FC_Mg_sil, FC_Mg_carb, FC_Ca, FC_Ca_sil, FC_Ca_carb, FC_Ca_sulf, FC_Fe, FC_Fe_sil, FC_Fe_sulf, FC_Mg+FC_Ca+FC_Fe,
					xriver[iResTime][3], xriver[iResTime][17], xriver[iResTime][18], xriver[iResTime][19], xriver[iResTime][20], xriver[iResTime][21], xriver[iResTime][22], xriver[iResTime][23], xriver[iResTime][24], hChazeFallout);

		}
		fclose (fout);

		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"Outputs/CompoOcean.txt");
		fout = fopen(title,"w");
		if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
		else {
			fprintf(fout, "'Ocean concentrations in mol/(kg H2O)'\n");
			fprintf(fout, "'Time (Gyr)' \t 'Mocean (kg)' \t 'Mriver (kg)' \t 'Seafloor weathering volume (m3)' \t 'Ocean pH' \t 'Seawater-seafloor equil pH' \t 'Ocean log f(O2) at Tsurf0' \t 'Rain pH' \t 'Ox C(aq)' \t 'Red C(aq)' \t Mg(aq) \t "
					"Ca(aq) \t Fe(aq) \t Si(aq) \t Na(aq) \t S(aq) \t Cl(aq)\n");
			fprintf(fout, "0 \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n",
					Mocean, Mriver, volSeafCrust, pHout, pH, 4.0*(pe+pH)-logKO2H2O, rainpH, xaq[0], xaq[1], xaq[6], xaq[7], xaq[8], xaq[9], xaq[10], xaq[11], xaq[12]);
		}
		fclose (fout);
	}
	else { // Recover where simulation left off
		printf("Recovering previous output from PHREEQC-3.1.2/io/OceanDiss.txtExec.txt...\n");
		int line_no = 0;
		int line_length = 300;
		char line[line_length];
		char buffer[30]; buffer[0] = '\0';
		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"PHREEQC-3.1.2/io/OceanDiss.txtExec.txt");
		fout = fopen(title,"r");
		if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
		else {
			while (fgets(line, line_length, fout)) {
				line_no++;
				if (line[1] == 'p' && line[2] == 'H') {
					for (i=0;i<10;i++) buffer[i] = line[i+7];
					pH = strtod((const char*)buffer, NULL); buffer[0] = '\0'; printf("pH \t %g\n", pH);
				}
				if (line[1] == 't' && line[2] == 'e' && line[3] == 'm' && line[4] == 'p') {
					for (i=0;i<10;i++) buffer[i] = line[i+9];
					Tsurf = strtod((const char*)buffer, NULL)+Kelvin; buffer[0] = '\0'; printf("Tsurf(K) \t %g\n", Tsurf);
				}
				if (line[1] == 'p' && line[2] == 'r' && line[3] == 'e' && line[4] == 's') {
					for (i=0;i<10;i++) buffer[i] = line[i+11];
					Psurf = strtod((const char*)buffer, NULL)*1.01325; buffer[0] = '\0'; printf("Psurf(bar) \t %g\n", Psurf);
				}
				if (line[1] == 'p' && line[2] == 'e') {
					for (i=0;i<10;i++) buffer[i] = line[i+7];
					pe = strtod((const char*)buffer, NULL); buffer[0] = '\0'; printf("pe \t %g\n", pe);
				}
				if (line[1] == 'C' && line[2] == 'a') {
					for (i=0;i<10;i++) buffer[i] = line[i+6];
					xaq[7] = strtod((const char*)buffer, NULL); buffer[0] = '\0'; printf("Ca \t %g\n", xaq[7]);
				}
				if (line[1] == 'M' && line[2] == 'g') {
					for (i=0;i<10;i++) buffer[i] = line[i+6];
					xaq[6] = strtod((const char*)buffer, NULL); buffer[0] = '\0'; printf("Mg \t %g\n", xaq[6]);
				}
				if (line[1] == 'N' && line[2] == 'a') {
					for (i=0;i<10;i++) buffer[i] = line[i+6];
					xaq[10] = strtod((const char*)buffer, NULL); buffer[0] = '\0'; printf("Na \t %g\n", xaq[10]);
				}
				if (line[1] == 'F' && line[2] == 'e') {
					for (i=0;i<10;i++) buffer[i] = line[i+6];
					xaq[8] = strtod((const char*)buffer, NULL); buffer[0] = '\0'; printf("Fe \t %g\n", xaq[8]);
				}
				if (line[1] == 'S' && line[2] == 'i') {
					for (i=0;i<10;i++) buffer[i] = line[i+6];
					xaq[9] = strtod((const char*)buffer, NULL); buffer[0] = '\0'; printf("Si \t %g\n", xaq[9]);
				}
				if (line[1] == 'C' && line[2] == 'l') {
					for (i=0;i<10;i++) buffer[i] = line[i+6];
					xaq[12] = strtod((const char*)buffer, NULL); buffer[0] = '\0'; printf("Cl \t %g\n", xaq[12]);
				}
				if (line[1] == 'S' && line[2] == '(' && line[3] == '6' && line[4] == ')') {
					for (i=0;i<10;i++) buffer[i] = line[i+8];
					xaq[11] = strtod((const char*)buffer, NULL); buffer[0] = '\0'; printf("S(6) \t %g\n", xaq[11]);
				}
				if (line[1] == 'C' && line[2] == '(' && line[3] == '4' && line[4] == ')') {
					for (i=0;i<10;i++) buffer[i] = line[i+8];
					xaq[0] = strtod((const char*)buffer, NULL); buffer[0] = '\0'; printf("C(4) \t %g\n", xaq[0]);
				}
				if (line[1] == 'N' && line[2] == 't' && line[3] == 'g' && line[4] == ' ') {
					for (i=0;i<10;i++) buffer[i] = line[i+7];
					xaq[3] = strtod((const char*)buffer, NULL); buffer[0] = '\0'; printf("Ntg \t %g\n", xaq[3]);
				}
				if (line[1] == 'N' && line[2] == 't' && line[3] == 'g' && line[4] == '(') {
					for (i=0;i<10;i++) buffer[i] = line[i+10];
					xgas[3] = strtod((const char*)buffer, NULL)*1.01325/Psurf; buffer[0] = '\0'; printf("Ntg(g) \t %g\n", xgas[3]);
				}
				if (line[1] == 'O' && line[2] == '2' && line[3] == '(' && line[4] == 'g') {
					for (i=0;i<10;i++) buffer[i] = line[i+9];
					xgas[2] = strtod((const char*)buffer, NULL)*1.01325/Psurf; buffer[0] = '\0'; printf("O2(g) \t %g\n", xgas[2]);
				}
				if (line[1] == 'C' && line[2] == 'O' && line[3] == '2' && line[4] == '(') {
					for (i=0;i<10;i++) buffer[i] = line[i+10];
					xgas[0] = strtod((const char*)buffer, NULL)*1.01325/Psurf; buffer[0] = '\0'; printf("CO2(g) \t %g\n", xgas[0]);
				}
				if (line[1] == 'C' && line[2] == 'H' && line[3] == '4' && line[4] == '(') {
					for (i=0;i<10;i++) buffer[i] = line[i+10];
					xgas[1] = strtod((const char*)buffer, NULL)*1.01325/Psurf; buffer[0] = '\0'; printf("CH4(g) \t %g\n", xgas[1]);
				}
				if (line[1] == 'H' && line[2] == '2' && line[3] == 'O' && line[4] == '(') {
					for (i=0;i<10;i++) buffer[i] = line[i+10];
					xgas[4] = strtod((const char*)buffer, NULL)*1.01325/Psurf; buffer[0] = '\0'; printf("H2O(g) \t %g\n", xgas[4]);
				}
				nAir = Psurf*bar2Pa*Asurf/g[NR]/molmass_atm(xgas);
			}
		}
		fclose (fout);
		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"Outputs/ReservoirsFluxes.txt");
		fout = fopen(title,"r");
		if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
		else {
			int scan = 0;
			i = 0;
			while (fgets(line, line_length, fout)) { // Divert realtime to tstart and Tmantle, all fluxes, and total N to dtime since dtime will be reset momentarily
				fgets(line, 1, fout);
				scan = fscanf(fout, "%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
						&tstart, &dtime, &RCmantle, &RCatm, &RCocean, &RCatmoc, &RCatmoc, &dtime, &dtime, &dtime, &dtime, &dtime, &RCorg);
//				if (scan != 1) printf("tail(): Error scanning %s output file at entry i = %d\n", title, i);
				i++;
			}
			printf("RCmantle \t %g\n", RCmantle);
			printf("Tmantle \t %g\n", Tmantle);
			printf("netFC \t %g\n", netFC);
			printf("RCorg \t %g\n", RCorg);
		}
		fclose (fout);
	}

	//-------------------------------------------------------------------
	// Start time loop
	//-------------------------------------------------------------------

//	double TotalN0 = xgas[3]*2.0*nAir + xaq[3]*Mocean;

	printf("Starting time loop...\n");
	while (realtime < tend) {

		dtime = fmin(dtime0, 0.1*fmin(RCatmoc, RCmantle)/fabs(netFC));
		if (realtime >= tstart && realtime < tstart + 1.0*Myr2sec) dtime = fmin(dtime, 0.2*Myr2sec); // Start slow

		if (realtime < 0.5*Gyr2sec) dtime = 0.1*Myr2sec;                           // Sufficient to achieve numerical convergence for geodynamics alone
		realtime += dtime;
//		realtime = 4.55*Gyr2sec; // Present day

		printf("\nTime: %g Gyr\n", realtime/Gyr2sec);

		//-------------------------------------------------------------------
		// Update surface temperature
		//-------------------------------------------------------------------

		// Parameterization of Caldeira & Kasting (1992) equation (4), valid between psi=-8 to -2 and T=0 to 100C
		S = S0/(1.0-0.38*(realtime/Gyr2sec/4.55-1.0));
//		S = S0;
		Teff = pow((1.0-albedo)*S/(4.0*sigStefBoltz),0.25);

		psi = log(xgas[0]*Psurf)/log(10.0);
		if (psi <= -2) DeltaTghe = 815.17 + 4.895e7/Tsurf/Tsurf - 3.9787e5/Tsurf - 6.7084/psi/psi + 73.221/psi - 30882.0/Tsurf/psi;
		else           DeltaTghe = 815.17 + 4.895e7/Tsurf/Tsurf - 3.9787e5/Tsurf - 6.7084/4.0     - 73.221/2.0 + 30882.0/Tsurf/2.0
								 + (psi+2.0)*10.0 + (Tsurf - Tfreeze)*0.05; // Coarse extrapolation beyond pCO2=0.01 bar

		// Fancier parameterization including the effects of CH4 (GHE and haze), needs work
//		if (xgas[1] < 0.1*xgas[0]) { // GHE due to CO2 and CH4 compound
//			psi = log((xgas[0]+25.0*xgas[1])*Psurf)/log(10.0); // Equivalent log pCO2 = log (pCO2 + 25*pCH4)
//
//			if (psi <= -2) DeltaTghe = 815.17 + 4.895e7/Tsurf/Tsurf - 3.9787e5/Tsurf - 6.7084/psi/psi + 73.221/psi - 30882.0/Tsurf/psi;
//			else           DeltaTghe = 815.17 + 4.895e7/Tsurf/Tsurf - 3.9787e5/Tsurf - 6.7084/4.0     - 73.221/2.0 + 30882.0/Tsurf/2.0
//									 + (psi+2.0)*10.0 + (Tsurf - Tfreeze)*0.05;
//		}
//		else DeltaTghe = 0.0; // Pale orange dot, anti-GHE haze

		Tsurf = Teff + DeltaTghe;

		if (Tsurf <= Tfreeze+0.01) printf("ExoCcycleGeo: Surface temperature = %g K <= 0.01 C. DeltaTghe=%g K.\n", Tsurf, DeltaTghe);

		//-------------------------------------------------------------------
		// Update atmosphere
		//-------------------------------------------------------------------

//		if (itime < 2) { // First 2 time steps: re-compute Vatm to avoid it being too influenced by ballpark initial conditions
			Vatm1 = nAir*R_G*Tsurf/(Psurf*bar2Pa);
			Tsurf1 = Tsurf;
//		}
		Vatm = Vatm1*Tsurf/Tsurf1; // After 2nd time step, keep Vatm constant

		// Redox of outgassing (see e.g. Gaillard and Scaillet 2014 EPSL or Scaillet and Gaillard 2011 Nature comment)

		if (xgas[0]*nAir + 1.0*dtime*netFC < 0.0) {
			netFC = -xgas[0]*nAir/(1.0*dtime)*(1.0-1.0e-8);
			RCmantle = RCmantle_old - dtime_old*netFC; // Re-adjust reservoirs accordingly
			RCatmoc  = RCatmoc_old  + dtime_old*netFC;
		}

		// TODO update xgas[0] and xgas[1] with Photochem calculation provided Foutgas as input; in that case update below only with netFC-Foutgas
		xgas[0] = (xgas[0]*nAir + (1.0-fCH4)*dtime*netFC)/(nAir + dtime*netFC); // (1-fCH4) fraction of added gas is CO2. Let equilibration with ocean speciate accurately
		xgas[1] = (xgas[1]*nAir +      fCH4 *dtime*netFC)/(nAir + dtime*netFC); // Dilute CH4 (or concentrate if netFC < 0)

		for (i=2;i<nAtmSpecies;i++) xgas[i] = xgas[i]*nAir/(nAir + dtime*netFC); // Dilute other gases accordingly

		nAir = nAir + dtime*netFC; // Update nAir
		Psurf = nAir*R_G*Tsurf/Vatm/bar2Pa; // Update Psurf assuming ideal gas law. Used to define Psurf = nAir*molmass_atm(xgas)*gsurf/Asurf/bar2Pa, but that introduces a drift
		P[NR] = Psurf*bar2Pa;

		if (Psurf < 0.01) {
			printf("ExoCcycleGeo: Pressure = %g bar too close to the triple point of H2O, oceans not stable at the surface. Exiting.\n", Psurf);
			exit(0);
		}

		//-------------------------------------------------------------------
		// Equilibrate ocean and atmosphere
		//-------------------------------------------------------------------

//		if (Tsurf > Tfreeze+0.01) {
//			printf("Equilibrating ocean and atmosphere... ");
//
//			AqueousChem(path, "io/OceanDiss.txt", Tsurf, &Psurf, &Vatm, &nAir, &pH, &pe, &Mocean, &xgas, &xaq, &xriver, 0, 0, 0.0, 1, nvarEq, 0.0, 0.0, &deltaCreac, staglid, dtime);
//
//			if (Psurf < 0.01) {
//				printf("ExoCcycleGeo: Pressure = %g bar too close to the triple point of H2O, oceans not stable at the surface. Exiting.\n", Psurf);
//				exit(0);
//			}
//
//			// If CH4 mixing ratio > 1.0 CO2, organic haze (Haqq-Misra et al. 2008, https://doi.org/10.1089/ast.2007.0197) (really Photochem should be telling ExoCcycleGeo this)
//			if (xgas[1] > 1.0*xgas[0]) {
//				hChazeFallout += 1.0e-4*xgas[1]*nAir/128.0*1320.7*1.0e-6/Asurf; // 1320.7 cm3/mol is the molar volume of KerogenC128 (Helgeson et al. 2009, https://doi.org/10.1016/j.gca.2008.03.004)
//				// Thickness of aerosol deposit, output, in what reservoir to put it?, decrease riverine flux accordingly
//				xgas[1] -= 1.0e-4*xgas[1];
//
//				AqueousChem(path, "io/OceanDiss.txt", Tsurf, &Psurf, &Vatm, &nAir, &pH, &pe, &Mocean, &xgas, &xaq, &xriver, 0, 0, 0.0, 1, nvarEq, 0.0, 0.0, &deltaCreac, staglid, dtime);
//
//				if (Psurf < 0.01) {
//					printf("ExoCcycleGeo: Pressure = %g bar too close to the triple point of H2O, oceans not stable at the surface. Exiting.\n", Psurf);
//					exit(0);
//				}
//			}
//
//			RCatm = (xgas[0]+xgas[1])*nAir;
//			RCocean = (xaq[0]+xaq[1])*Mocean;
//			pHout = pH;
//
//			cleanup(path); // Remove PHREEQC selected output file
//		}

		//-------------------------------------------------------------------
		// Calculate surface C flux from outgassing
		//-------------------------------------------------------------------

		// 1. Determine the base of the lithosphere, temperature profile, and vigor of convection as a function of tectonic mode
		// In stagnant lid mode, melting at base of lithosphere.
		// In plate tectonics mode, melting at base of lithosphere away from mid-ocean ridges (MOR)
		//                        + melting at surface at MOR.

//		for (i=0;i<NR;i++) Tadiab[i] = Tmantle - alpha*g[(int)(0.5*(NR+iCMB))]*Tmantle/Cp * (r[i] - (r_p - zLith));
		for (i=0;i<NR;i++) Tadiab[i] = Tmantle*exp( alpha*g[i]*(r_p - r[i] - zLith)/Cp );
//		double Tp = Tadiab[NR-1]; // Mantle potential temperature
//		printf("Tmantle=%g T_CMB=%g Tp=%g\n", Tmantle, Tadiab[iCMB], Tp);

		// Solomatov (1995) equation (24) stagnant lid scaling, see also Foley et al. (2020) chapter just before equation (4.5)
		p = (flowLawDiff[0] + flowLawDiff[6]*(Fe_FeMg - Fe_FeMg_ref))/(R_G*Tmantle*Tmantle)*(Tmantle-Tsurf); // Taking diffusion creep as the dominant mechanism at low stress

		//  Upper mantle viscosity
		nu = combVisc(Tmantle, P[iLith], flowLawDiff, flowLawDisl, grainSize, C_OH, Fe_FeMg, fPx, tConv)/rho[iLith];
//		nu = 1.0e19/3500.0;
//		nu = 1.0e16*exp((2.0e5 + P[(int)((iLith+iCMB)/2.0)]*1.1e-6)/(R_G*Tmantle))/rho[(int)((iLith+iCMB)/2)]; // Cízková et al. (2012)
//		nu = 1.0e16*exp((2.0e5 + P[iCMB]*1.1e-6)/(R_G*(Tmantle - alpha*g[iCMB]*Tmantle/Cp * (r[iCMB] - 0.5*(r_p+r_c)))))/rho[iCMB]; // Cízková et al. (2012) for CMB depth

		// Vigor of convection, assuming whole-mantle convection (Kite et al. 2009)
		if (staglid) {
			// Delta_Tth that drives stagnant lid convection
			for (i=iLith;i<NR;i++) {
				if ((Tmantle-T[i])/(Tmantle-Tsurf) > 1.0/p) {
					Tref = T[i];
					break;
				}
			}
		}
		else Tref = Tsurf;

		Ra = alpha*g[iLith]*(Tmantle-Tref)*pow(r_p-r_c,3.0)/(kappa*nu);
		if (Ra < Ra_c) {
			printf("No convection\n");
			Ra = Ra_c;
		}
		Nu = pow(Ra/Ra_c, beta);

		zLith = (r_p-r_c)/Nu; // zLith = d/Nu per equation 4.3 of Foley et al. (2020). If Nu = 1, Tmantle = T_CMB (initialization).
		for (i=iCMB;i<NR;i++) {
			if (r[i] < r_p-zLith && r[i+1] > r_p-zLith) iLith = i;
		}

		// Convective velocity
//		vConv = 2.0*(Nu-1.0) * k / (rho[(int)((iLith+iCMB)/2)]*Cp*(r[iLith]-r_c)) * (Tmantle-Tsurf) / (Tmantle-Tref); // Kite et al. (2009) equation (24)
		vConv = 0.0 * 0.271*kappa/(r_p-r_c)*pow(Ra,2.0/3.0) // Eq. (6.369) of Turcotte & Schubert (2002), p. 511, fluid heated from below, which is true in part for Earth's mantle, see Romanowicz et al. 2002; Lay et al. 2008
		      + 1.0 * 0.354*kappa/(r_p-r_c)*sqrt(Ra); // Eq. (6.379) of Turcotte & Schubert (2002), p. 514, fluid heated from within
		tConv = (r_p-r_c)/vConv; // Update convection timescale

		// ------------------------------------
		// 2. Melting & outgassing model (Kite et al. 2009)

		// 2a. Calculate average temperatures: conductive profile in lid, adiabat below
		for (i=0   ;i<iCMB;i++) T[i] = 0.0; // Temperatures not computed in the core
		for (i=iCMB;i<NR  ;i++) {
			if (r[i] < r_p - zLith) T[i] = Tadiab[i];
			else T[i] = Tsurf + (Tmantle-Tsurf)*(r_p-r[i])/zLith; // Linear profile (Foley et al. 2020)
		}

		// 2b. Update alphaMELTS input files
		// Output P-T profile in lithosphere+boundary layer to be fed into alphaMELTS through PTexoC.txt
		// MELTS crashes below solidus for default composition, but that's enough to capture all depths above zLith at which melting occurs.

		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"alphaMELTS-1.9/ExoC/PTexoC.txt");

		fout = fopen (title,"w");
		if (fout == NULL) printf("ExoCcycleGeo: Missing PT file path: %s\n", title);

		if (staglid) { // Stagnant lid: partial melting of passively upwelling mantle (e.g., O'Rourke and Korenaga 2015, equations 17-20)
			ir = 0;
			for (i=iCMB;i<NR;i++) {
				if (P[i] > PminMELTS && P[i] < PmaxMELTS) { // Nominally, 1 to 4 GPa
					if (T[i]-Kelvin > TminMELTS) {
						fprintf(fout, "%g %g\n", P[i]/bar2Pa, T[i]-Kelvin); // Don't input below 750 C to avoid MELTS crashing.
						ir++;
					}
					else break;
				}
			}
		}
		else { // Plate tectonics: only calculate melting near the surface at ridges or rifts
			ir = 1;
			fprintf(fout, "%g %g\n", 1000.0*g[NR-1]*(r_p - pow( pow(r_p,3.0) - Mocean/1000.0*3.0/(4.0*PI_greek) ,1.0/3.0))/bar2Pa, Tadiab[NR-1]-Kelvin);
		}

		fclose(fout);

		// 2c. Run alphaMELTS to compute melt fraction along P-T profile in lithosphere and mantle
		// Reset sys_tbl
		for (i=0;i<NR;i++) {
			for (j=0;j<18;j++) sys_tbl[i][j] = 0.0;
		}

		// Call alphaMELTS
		if (realtime > tstart && realtime <= tend) {
			printf("Outgassing... ");
			alphaMELTS(path, 0, ir, "ExoC/ExoC_env.txt", &sys_tbl);
		}

		imin = 0;
		imax = 0;
		j = 0;
		if (staglid) {
			for (i=iCMB+1;i<NR;i++) {
				Meltfrac[i] = 0.0;
				if (P[i] > PminMELTS && P[i] < PmaxMELTS) { // Nominally, 1 to 4 GPa
					if (sys_tbl[j][4] > 0.0 && (sys_tbl[j][13] < 1.5*sys_tbl[j][14] || sys_tbl[j][4] == 1.0)) { // Don't store if volume melt fraction negative or if density of melt > n*density of solid
						if (T[i]-Kelvin > TminMELTS) { // There was no input below 750 C to avoid MELTS crashing.
							if (fabs(1.0-P[i]/(sys_tbl[j][0]*bar2Pa)) > 1.0e-3) printf("ExoCcycleGeo: pressures from grid (%g bar) and MELTS (%g bar) misaligned at grid index %d\n", P[i]/bar2Pa, sys_tbl[j][0], i);
							if (fabs(1.0-T[i]/ sys_tbl[j][1]        ) > 1.0e-3) printf("ExoCcycleGeo: temperatures from grid (%g K) and MELTS (%g K) misaligned at grid index %d\n", T[i], sys_tbl[j][1], i);
							Meltfrac[i] = sys_tbl[j][3];
							rhomelt = sys_tbl[j][13]*1000.0; // g cm-3 to kg m-3, taken at depth of highest melt fraction for calculation of thickness of new crust generated
							if (imin == 0) imin = i;
						}
						imax = i;
					}
					j++;
				}
			}
		}
		else {
			if (sys_tbl[0][3] > 0.0 && (sys_tbl[0][13] < 1.5*sys_tbl[0][14] || sys_tbl[0][3] == 1.0)) { // Don't store if volume melt fraction negative or if density of melt > n*density of solid
				if (Tadiab[NR-1]-Kelvin > TminMELTS) {
					if (fabs(1.0-1000.0*g[NR-1]*(r_p - pow( pow(r_p,3.0) - Mocean/1000.0*3.0/(4.0*PI_greek) ,1.0/3.0))/(sys_tbl[0][0]*bar2Pa)) > 1.0e-3) printf("ExoCcycleGeo: pressures from grid (%g bar) and MELTS (%g bar) differ\n", P[i]/bar2Pa, sys_tbl[0][0]);
					if (fabs(1.0-Tadiab[NR-1]/ sys_tbl[0][1]) > 1.0e-3) printf("ExoCcycleGeo: temperatures from grid (%g K) and MELTS (%g K) differ\n", Tadiab[NR-1], sys_tbl[0][1]);
					Meltfrac[0] = sys_tbl[0][3];
					rhomelt = sys_tbl[0][13]*1000.0; // g cm-3 to kg m-3, taken at depth of highest melt fraction, for calculation of thickness of new crust generated
				}
				else printf("Mantle potential temperature too low for melting\n");
			}
			else {
				printf("Error: melt fraction negative or density of melt  > 1.5*density of solid\n");
			}
		}

		meltmass = 0.0;
		ir = 0;
//		slope = 0.0;
//		islope = 0;
//		if (imax > imin) {
//			if (Meltfrac[imin] > 0.0) {
////				printf("ExoCcycleGeo: alphaMELTS could only calculate melting down to depth %g km (melt fraction %.2g > 0.1), extrapolating melting curve to 0 linearly with depth\n", (r_p-r[imin])/km2m, Meltfrac[imin]);
//				// Extrapolate starting from last 5 indices, provided melt fraction of all is < 1 (otherwise, use less indices)
//				if (imax-imin > 5) { // Only do this if there are > 5 grid zones of melt, otherwise chances are the slope will be skewed.
//					for (i=imin;i<imin+nslopeAvg;i++) {
//						if (P[i+1]-P[i] == 0) printf("ExoCcycleGeo: Can't extrapolate rock melting curve with pressure because denominator is zero\n");
//						else {
//							if (Meltfrac[i] < Meltfrac[i+1]) { // Ensure we're below depth of max melt
//								slope = slope + (Meltfrac[i+1]-Meltfrac[i])/(P[i+1]-P[i]);
//								islope++;
//							}
//						}
//					}
//					slope = slope/(double) islope;
//
//					if (islope) {
//						for (i=imin-1;i>iCMB;i--) {
//							meltfrac = Meltfrac[i+1] + slope*(P[i]-P[i+1]);
//							if (meltfrac > 0.0) Meltfrac[i] = meltfrac;
//							else break;
//							ir++;
//						}
//					}
//				}
//			}
//		}

//		// 2d. Alternatively, use analytical formulation of McKenzie & Bickle (1988) at the expense of compositional versatility?
//		// McKenzie & Bickle (1988): determine Tsolidus at P, then Tliquidus at P, then T'(T, Tsol, Tliq), then X(T')
//		//2c-1 Geotherm
//		double Tsolidus = 0.0;  // Solidus temperature (Kelvin)
//		double Tliquidus = 0.0; // Liquidus temperature (Kelvin)
//		double Tprime = 0.0;    // Temperature used in calculation of Xmelt (dimensionless)
//		double Xmelt[100];      // Melt fraction (dimensionless)
//		for (i=0;i<100;i++) Xmelt[i] = 0.0;
//
//		for (i=0;i<100;i++) {
//			// Determine Tsolidus
//
//			// Initialize bounds for Tsolidus
//			T_INF = Tsurf;
//			T_SUP = Tmantle;
//
//			// Ensure that f(T_INF)<0 and f(T_SUP)>0
//			f_inf = Psolidus(T_INF) - Pgeotherm[i];
//			f_sup = Psolidus(T_SUP) - Pgeotherm[i];
//
//			if (f_inf*f_sup > 0.0) {
//				printf("ExoCcycleGeo: Solidus not between Tsurf = %g K and Tmantle = %g K, cannot determine Tsolidus and melt fraction!\n", Tsurf, Tmantle);
//			}
//			else {
//				if (f_inf > 0.0) { // Swap INF and SUP if f_inf > 0 and f_sup < 0
//					T_TEMP = T_INF;
//					T_INF = T_SUP;
//					T_SUP = T_TEMP;
//				}
//				Tsolidus = 0.5*(T_INF+T_SUP); // Initialize the guess for the root,
//				dTold = fabs(T_INF-T_SUP); // Initialize the "stepsize before last"
//				dT = dTold;                // Initialize the last stepsize
//
//				f_x = Psolidus(Tgeotherm[i]) - Pgeotherm[i];
//				f_prime_x = dPsolidusdT(Tgeotherm[i]);
//
//				// Loop over allowed iterations to find root of f
//				n_iter = 0;
//				while (fabs(f_x) > NewtRaphThresh) {
//
//					// Bisect if Newton is out of range, or if not decreasing fast enough
//					if ((((Tsolidus-T_SUP)*f_prime_x-f_x)*((Tsolidus-T_INF)*f_prime_x-f_x) > 0.0) || (fabs(2.0*f_x) > fabs(dTold*f_prime_x))) {
//						dTold = dT;
//						dT = 0.5*(T_SUP-T_INF);
//						Tsolidus = T_INF + dT;
//					}
//					else { // Do Newton-Raphson
//						dTold = dT;
//						dT = f_x/f_prime_x;
//						Tsolidus = Tsolidus - dT;
//					}
//					// Calculate updated f and f'
//					f_x = Psolidus(Tsolidus) - Pgeotherm[i];
//					f_prime_x = dPsolidusdT(Tsolidus);
//
//					if (f_x < 0.0) T_INF = Tsolidus; // Maintain the bracket on the root
//					else T_SUP = Tsolidus;
//
//					n_iter++;
//					if (n_iter>=n_iter_max) {
//						printf("ExoCcycleGeo: could not find the brittle-ductile transition after %d iterations\n",n_iter_max);
//						break;
//					}
//				}
//			}
//
//			if (Tgeotherm[i] < Tsolidus) Xmelt[i] = 0.0;
//			else {
//				// Determine Tliquidus
//				Tliquidus = Kelvin + 1736.2 + 4.343*Pgeotherm[i]/1.0e9 + 180.0 * atan(Pgeotherm[i]/1.0e9/2.2169);
//				if (Tgeotherm[i] > Tliquidus) {
//					Xmelt[i] = 1.0;
//				}
//				else {
//					Tprime = (Tgeotherm[i] - (Tsolidus+Tliquidus)/2.0) / (Tliquidus-Tsolidus); // K or C doesn't matter as the numerator and the denominator are both temp differences
//					Xmelt[i] = 0.5 + Tprime + (Tprime*Tprime - 0.25)*(0.4256+2.988*Tprime); // Again, K or C doesn't matter, Tprime is dimensionless
//				}
//			}
//			printf("%g \t %g \t %g\n", Pgeotherm[i]/bar2Pa, Xmelt[i], Tgeotherm[i]);
//		}
//
//		//2c-2 Mantle adiabat
//		double T[NR]; // Temperature (K)
//		for (ir=0;ir<NR;ir++) T[ir] = 0.0;
//
//		for (ir=NR-1;ir>=0;ir--) {
//			if (ir < iLith && P[ir] < 100000.0*bar2Pa) {
//				// Adiabat
//				T[ir] = T[iLith] + gsurf*alpha*Tmantle/Cp*(r_p-zLith-r[ir]);
//
//				// Determine Tsolidus
//
//				// Initialize bounds for Tsolidus
//				T_INF = Tsurf;
//				T_SUP = 10000.0;
//
//				// Ensure that f(T_INF)<0 and f(T_SUP)>0
//				f_inf = Psolidus(T_INF) - P[ir];
//				f_sup = Psolidus(T_SUP) - P[ir];
//
//				if (f_inf*f_sup > 0.0) {
//					printf("ExoCcycleGeo: Solidus not between Tsurf = %g K and Tmantle = %g K, cannot determine Tsolidus and melt fraction!\n", Tsurf, Tmantle);
//				}
//				else {
//					if (f_inf > 0.0) { // Swap INF and SUP if f_inf > 0 and f_sup < 0
//						T_TEMP = T_INF;
//						T_INF = T_SUP;
//						T_SUP = T_TEMP;
//					}
//					Tsolidus = 0.5*(T_INF+T_SUP); // Initialize the guess for the root,
//					dTold = fabs(T_INF-T_SUP); // Initialize the "stepsize before last"
//					dT = dTold;                // Initialize the last stepsize
//
//					f_x = Psolidus(T[ir]) - P[ir];
//					f_prime_x = dPsolidusdT(T[ir]);
//
//					// Loop over allowed iterations to find root of f
//					n_iter = 0;
//					while (fabs(f_x) > NewtRaphThresh) {
//
//						// Bisect if Newton is out of range, or if not decreasing fast enough
//						if ((((Tsolidus-T_SUP)*f_prime_x-f_x)*((Tsolidus-T_INF)*f_prime_x-f_x) > 0.0) || (fabs(2.0*f_x) > fabs(dTold*f_prime_x))) {
//							dTold = dT;
//							dT = 0.5*(T_SUP-T_INF);
//							Tsolidus = T_INF + dT;
//						}
//						else { // Do Newton-Raphson
//							dTold = dT;
//							dT = f_x/f_prime_x;
//							Tsolidus = Tsolidus - dT;
//						}
//						// Calculate updated f and f'
//						f_x = Psolidus(Tsolidus) - P[ir];
//						f_prime_x = dPsolidusdT(Tsolidus);
//
//						if (f_x < 0.0) T_INF = Tsolidus; // Maintain the bracket on the root
//						else T_SUP = Tsolidus;
//
//						n_iter++;
//						if (n_iter>=n_iter_max) {
//							printf("ExoCcycleGeo: could not find the brittle-ductile transition after %d iterations\n",n_iter_max);
//							break;
//						}
//					}
//				}
//
//				if (T[ir] < Tsolidus) Xmelt[i] = 0.0;
//				else {
//					// Determine Tliquidus
//					Tliquidus = Kelvin + 1736.2 + 4.343*P[ir]/1.0e9 + 180.0 * atan(P[ir]/1.0e9/2.2169);
//					if (T[ir] > Tliquidus) {
//						Xmelt[i] = 1.0;
//					}
//					else {
//						Tprime = (T[ir] - (Tsolidus+Tliquidus)/2.0) / (Tliquidus-Tsolidus); // K or C doesn't matter as the numerator and the denominator are both temp differences
//						Xmelt[i] = 0.5 + Tprime + (Tprime*Tprime - 0.25)*(0.4256+2.988*Tprime); // Again, K or C doesn't matter, Tprime is dimensionless
//					}
//				}
//				printf("%g \t %g \t %g\n", P[ir]/bar2Pa, Xmelt[i], T[ir]);
//			}
//		}

		// 2e. Convert to melt fraction vs. depth and determine amount of melt generated.
		// Assumes all melt generated reaches the surface and all ascending mantle parcels reach pressures < Psolidus (Kite et al. 2009)

//		if (imax > imin) {
			if(realtime > tstart) {
				title[0] = '\0';
				if (cmdline == 1) strncat(title,path,strlen(path)-20);
				else strncat(title,path,strlen(path)-18);
				strcat(title,"Outputs/Geotherm.txt");
				fout = fopen(title,"a");
				if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
				else {
					fprintf(fout, "Time (Gyr) \t Pressure (bar) \t Depth (km) \t Melt fraction \t Temp (K) \t Density (kg m-3)\n");
					for (i=iCMB;i<=NR;i++) {
						if (P[i] > PminMELTS && P[i] < PmaxMELTS)
							fprintf(fout, "%g \t %g \t %g \t %g \t %g \t %g\n", realtime/Gyr2sec, P[i]/bar2Pa, (r_p-r[i])/km2m, Meltfrac[i], T[i], rho[i]);
					}
				}
				fclose (fout);

				if (staglid) {
					for (i=imin-ir;i<=imax;i++) {
						meltmass += Meltfrac[i]*4.0/3.0*PI_greek*(pow(r[i+1],3)-pow(r[i],3))*rho[i];
					}
					if (Meltfrac[imin-ir]*4.0/3.0*PI_greek*(pow(r[imin-ir+1],3)-pow(r[imin-ir],3))*rho[imin-ir] > 0.5*meltmass)
						meltmass -= Meltfrac[imin-ir]*4.0/3.0*PI_greek*(pow(r[imin-ir+1],3)-pow(r[imin-ir],3))*rho[imin-ir]; // Spurious nonzero melt fraction only in the bottom grid zone of the MELTS calculation
				}
			}
//		}

        if (staglid) dMmelt_dt = meltmass / tConv;
        else dMmelt_dt = LplateRdg * WplateRdg * vConv * Meltfrac[0] * rhomelt;
		
		FCoutgas = dMmelt_dt * magmaCmassfrac / molmassC;

		if (rhomelt > 0.0 && r[iLith] > r_c) dzCrust_dt = dMmelt_dt/rhomelt / (4.0*PI_greek) / pow(r_p-zCrust,2.0); // Derived by expressing V = 4/3*pi*(r_p^3-(r_p-zCrust)^3), differentiating with time, and rearranging to express dzCrust_dt = f(dV/dt, zCrust, r_p)
		else dzCrust_dt = 0.0;
		zCrust = zCrust + dzCrust_dt*dtime;

//		FCoutgas = meltmass*magmaCmassfrac/0.044*vConv/(r[iLith]-r[iCMB])*0.4; // mol C s-1, 40% melt reaches the surface
//		if (realtime < 1.0*Gyr2sec) FCoutgas = FCoutgas * (realtime/(0.4*Gyr2sec)-1.5); // Ramp up outgassing from 0.6 to 1.0 Gyr to avoid step function

//      // Alternative: Kite et al. (2009) eq. 25 They didn't scale with mass: [sum melt fraction (depth)] * [mass (depth)] / [total mass between surf and Psolidus]. Also their typo: P0>Pf.
//		double Rmelt = 0.0;             // Rate of melt generation (m-2 s-1)
//		if (staglid) Rmelt = Asurf*vConv*rhoMagma; // *  (rhoCrust*gsurf*zCrust/(Psolidus-Psurf)); // kg s-1
//		// Enhancement of rate of melting due to plate tectonics. Focuses solely on mid-ocean ridges (their Sec. 2.3). Also a function of age. See Sasaki & Tajika 1995.
//		else Rmelt = vConv*MORlength*zCrust*rhoCrust; // MORlength presumably changes over time (Sasaki & Tajika 1995)
//
//		// Includes extrusive and intrusive melting because both degas (K09).
//		FCoutgas = deltaCvolcEarth * Rmelt/(RmeltEarth*mEarth); // mol C s-1. Accurate for Earth magma C concentrations, but for other mantle concentrations, should be calculated explicitly (MELTS?)

		// If MELTS calculation results are too erratic, stick to an approximate time envelope of the MELTS calculation.
		// Could improve this extrapolation by memorizing FCoutgas through time.
//		if (realtime > tstart) {
//			FCoutgas_mem[t][0] = log(FCoutgas);
//			FCoutgas_mem[t][1] = realtime/Gyr2sec;
//
//			if (t >= nslopeAvg) {
//				if (t >= t_reset+nslopeAvg) { // Enough reliable data to update the slope
//					slope = 0.0;
//					islope = 0;
//					for (i=t-nslopeAvg;i<t;i++) {
//						if (FCoutgas_mem[i][1]-FCoutgas_mem[i-1][1] == 0) printf("ExoCcycleGeo: Can't extrapolate outgassing flux curve with time because denominator is zero\n");
//						else {
//							slope = slope + (FCoutgas_mem[i][0]-FCoutgas_mem[i-1][0])/(FCoutgas_mem[i][1]-FCoutgas_mem[i-1][1]);
//							islope++;
//						}
//					}
//					slope = slope/(double) islope;
//					printf("slope = %g (log mol C s-1)/Gyr\n", slope);
//				}
//
//				if (FCoutgas <=0 || fabs(1.0 - FCoutgas / (FCoutgas_mem[t-1][0] + slope*(realtime/Gyr2sec-FCoutgas_mem[t-1][1]))) > 0.1) { // Threshold: 10% variation from trend
//					FCoutgas = exp(FCoutgas_mem[t-1][0] + slope * (realtime/Gyr2sec - FCoutgas_mem[t-1][1]));
//					t_reset = t;
//				}
//			}
//			t++;
//		}

		printf("\n");
		printf("Time: %g Gyr\n", realtime/Gyr2sec);
		printf("---------------------------------------------------------------------------\n");
		printf("Quantity                                | Earth benchmark   | Model result \n");
		printf("---------------------------------------------------------------------------\n");
		printf("New crust generation rate (m Myr-1)     | 40                | %.3g \n", dzCrust_dt*Myr2sec);
		printf("New crust density (kg m-3)              | 2800              | %.4g \n", rhomelt);
		printf("C and H2O outgassing rate (mol s-1)     | 129000-407000     | %.6g \n", FCoutgas);
		printf("Magma C mass fraction (ppm)             | 300-10000         | %.3g \n", magmaCmassfrac*1.0e6);
		printf("Convective velocity (cm yr-1)           | ~1                | %.2g \n", vConv*100.0*1.0e-6*Myr2sec);
		printf("Mantle convection timescale (Myr)       | ~50-200           | %.4g \n", tConv/Myr2sec);
		printf("Rayleigh number                         | ~1e6-1e7          | %.2g \n", Ra);
		printf("Half-mantle depth (km)                  | 1450              | %.4g \n", (r_p-r_c-zLith)/2.0/km2m);
		printf("Heat flux (mW m-2)                      | 86 (radiodecay 40)| %.3g \n", k*(Tmantle-Tsurf)/zLith*1000.0);
		printf("Thickness of lithosphere (km)           | 100               | %.3g \n", zLith/km2m);
		printf("Temperature at base of lithosphere (ºC) | 1400              | %.4g \n", Tmantle-Kelvin);
		printf("Adiabatic temperature gradient (K km-1) | 0.5               | %.2g \n", (T[iCMB+1]-Tmantle)/(r_p-zLith-r_c)*km2m);
		printf("---------------------------------------------------------------------------\n");
		printf("\n");

		// ------------------------------------
		// 3. Mantle thermal evolution
		// Instantaneous heating rate (Kite et al. 2009 Table 1): H = X_4.5 * W * exp(ln(1/2) / t1/2 * (t-4.5)) in W kg-1
		H =   x40K * 2.92e-5 * exp(log(0.5)/( 1.26 *Gyr2sec) * (realtime - 4.5*Gyr2sec))  //  40-K
		  + x232Th * 2.64e-5 * exp(log(0.5)/(14.0  *Gyr2sec) * (realtime - 4.5*Gyr2sec))  // 232-Th
		  +  x235U * 56.9e-5 * exp(log(0.5)/( 0.704*Gyr2sec) * (realtime - 4.5*Gyr2sec))  // 235-U
		  +  x238U * 9.46e-5 * exp(log(0.5)/( 4.47 *Gyr2sec) * (realtime - 4.5*Gyr2sec)); // 238-U
		H = H*radioMult; // Multiplier input (default is 1)
		// Correct for radionuclide loss from the mantle due to partitioning into the crust (Turcotte & Schubert 2002, Table 4.3, p. 248):
		// Depleted mantle has 1/30, continental crust 50x, and oceanic crust (tholeiitic basalt) 2.5x the radionuclide content of the reference mantle.
		// Mass of crust = 4/3*pi*(r_p^3-(r_p-zCrust)^3))*rhomelt. Cap crustal thickness at 50 km to account for material recycling into the mantle, even in the stagnant lid regime (Bedard 2018)
		// Mass of mantle = m_p-m_c
		H = H * (1.0 - (1.0-1.0/30.0) * 4.0/3.0*PI_greek*(pow(r_p,3.0)-pow(r_p-fmin(zCrust,50.0e3),3.0))*rhomelt / (m_p-m_c));
		// Effective thermal conductivity scaled with Nu
		Tmantle = Tmantle + dtime*(H/Cp - 3.0*kappa*Nu*(Tmantle-Tref)/(r_p-r_c)*r_p*r_p/(pow(r_p,3.0)-pow(r_c,3.0)));

		//-------------------------------------------------------------------
		// Calculate surface C flux from continental weathering
		//-------------------------------------------------------------------

		FCcontW = 0.0;
//		if (Tsurf > Tfreeze+0.01 && realtime > tstart && realtime <= tend) {
//			// Analytical calculation from Edson et al. (2012) Eq. 1; Abbot et al. (2012) Eq. 2
////			FCcontW = -L * 0.5*deltaCcontwEarth*Asurf * pow(xgas[0]/xCO2g0,0.3) * runoff/runoff_Earth * exp((Tsurf-TsurfEarth)/17.7);
//			L = Lnow*realtime/(4.56*Gyr2sec);
//			tResLand = tResLandNow*realtime/(4.56*Gyr2sec);
////			tResLand = tResLandNow;
//
//		    kintime = 1.01*tResLand; // Total time of kinetic simulation
//		    iResTime = floor(tResLand / kintime * (double)kinsteps); // Index in xriver corresponding to output at riverResTime;
//
//			for (i=0;i<kinsteps;i++) {
//				for (j=0;j<nvarKin;j++) xriver[i][j] = 0.0;
//			}
//
//			runoff = runoff0*pow(1.025,Tsurf-Tsurf0); // 2.5% increase in global mean precipitation per K of temperature increase (Allen and Ingram 2002, Trenberth et al. 2005, Pendergrass 2020)
//			WRcontW = 5000.0*runoff/runoff_Earth;
//
//			printf("Continental weathering... ");
//			AqueousChem(path, "io/ContWeather.txt", Tsurf, &Psurf, &Vatm, &nAir, &pH, &pe, &WRcontW, &xgas, &xaq, &xriver, 0, 1, kintime, kinsteps, nvarKin, 0.0, 0.0, &deltaCreac, staglid, dtime);
//
//			rainpH = xriver[1][3]; // xriver[0][3], the initial rain speciation, is returned as = 0, so this is as close as it gets (smallest reaction time)
//			massH2Oriver = xriver[iResTime][7];
//
//			if (hChazeFallout > 1.0) {
//				for (i=9;i<nvarKin;i++) xriver[iResTime][i] /= hChazeFallout; // Continental weathering decreasingly effective with increasing haze fallout thickness
//			}
//
//			// River abundances of cations (mol/kg)
////			xriver_Mg_evap = xriver[iResTime][11]/(1.0-fracEvap);
////			xriver_Ca_evap = xriver[iResTime][13]/(1.0-fracEvap);
////			xriver_Si_evap = xriver[iResTime][12]/(1.0-fracEvap);
////			xriver_Na_evap = xriver[iResTime][10]/(1.0-fracEvap);
////			xriver_Fe_evap = xriver[iResTime][14]/(1.0-fracEvap);
//
//			// Dissolved carbonates, sulfates, sulfides (mol)
//			Mg_carb_consumed = - xriver[iResTime][nvarKin-7] - xriver[iResTime][nvarKin-6] - xriver[iResTime][nvarKin-5]; // Dolomite-dis, -ord, and magnesite
//			Ca_carb_consumed = - xriver[iResTime][nvarKin-8] - xriver[iResTime][nvarKin-7] - xriver[iResTime][nvarKin-6]; // Dolomite-dis, -ord, and calcite
//			Ca_sulf_consumed = - xriver[iResTime][nvarKin-4] - xriver[iResTime][nvarKin-3]; // Anhydrite and gypsum
//			Fe_sulf_consumed = - xriver[iResTime][nvarKin-2] - xriver[iResTime][nvarKin-1]; // Pyrite and pyrrhotite
//
//			Mriver = Asurf*L*runoff/(1.0-fracEvap)*rhoH2O*dtime;
//
//			// Individual carbon fluxes (mol/s). Rainout = (river runoff)/(1-fracEvap).
//			FC_Mg_carb =     Mg_carb_consumed/massH2Oriver*Mriver/dtime; // Factor of 1 because despite divalent cation vs. monovalent bicarbonate, 1 biarbonate is released in dissolution, so net 1
//			FC_Mg_sil  = 2.0*(xriver[iResTime][11]        *Mriver/dtime - FC_Mg_carb); // Factor of 2 because of divalent cation and monovalent bicarbonate
//			if (FC_Mg_sil < 0.0) printf("time=%g Gyr, Continental weathering: Mg silicate net formation\n", realtime/Gyr2sec);
//			FC_Mg = FC_Mg_carb + FC_Mg_sil;
//
//			FC_Ca_carb =     Ca_carb_consumed/massH2Oriver*Mriver/dtime;
//			FC_Ca_sulf = 2.0*Ca_sulf_consumed/massH2Oriver*Mriver/dtime; // Factor of 2 because of divalent cation and monovalent bicarbonate
//			FC_Ca_sil  = 2.0*(xriver[iResTime][13]        *Mriver/dtime - FC_Ca_carb - 0.5*FC_Ca_sulf);
//			if (FC_Ca_sil < 0.0) printf("time=%g Gyr, Continental weathering: Ca silicate net formation\n", realtime/Gyr2sec);
//			FC_Ca = FC_Ca_carb + FC_Ca_sulf + FC_Ca_sil;
//
//			FC_Fe_sulf = 2.0*Fe_sulf_consumed/massH2Oriver*Mriver/dtime; // Factor of 2 because of divalent cation and monovalent bicarbonate
//			FC_Fe_sil  = 2.0*(xriver[iResTime][14]        *Mriver/dtime - 0.5*FC_Fe_sulf);
//			if (FC_Fe_sil < 0.0) printf("time=%g Gyr, Continental weathering: Fe silicate net formation\n", realtime/Gyr2sec);
//			FC_Fe = FC_Fe_sulf + FC_Fe_sil;
//
//			FCcontW = - FC_Mg - FC_Ca - FC_Fe; // Negative out of atmosphere
//
////			printf("%d %g %g %g %g %g %g\n", iResTime, xriver_Mg_evap, Mg_carb_consumed, massH2Oriver, FC_Mg_carb/1.0e12*Yr2sec, FC_Mg_sil/1.0e12*Yr2sec, FC_Mg/1.0e12*Yr2sec); // Tmol/yr
////			printf("%g %g %g %g %g %g %g %g\n", xriver_Ca_evap, Ca_carb_consumed, Ca_sulf_consumed, massH2Oriver, FC_Ca_carb/1.0e12*Yr2sec, FC_Ca_sulf/1.0e12*Yr2sec, FC_Ca_sil/1.0e12*Yr2sec, FC_Ca/1.0e12*Yr2sec); // Tmol/yr
////			printf("%g %g %g %g %g %g\n", xriver_Fe_evap, Fe_sulf_consumed, massH2Oriver, FC_Fe_sulf/1.0e12*Yr2sec, FC_Fe_sil/1.0e12*Yr2sec, FC_Fe/1.0e12*Yr2sec); // Tmol/yr
//		}
//
//		// Runoff removes part of the haze deposit at erosion rate scaled from modern Earth (assuming physical erosion here)
//		hChazeFallout -= 20.0e12/Yr2sec /2730.0 /0.29 /Asurf *dtime *runoff/runoff_Earth; // 20e12 kg: modern erosion rate on Earth (Borrelli et al. 2017, https://doi.org/10.1038/s41467-017-02142-7),
//		                                                                                  // scaled to global thickness using continental crust density of 2730 kg/m3 from Rudnick & Gao 2003 and modern Earth land areal fraction of 0.29
//		if (hChazeFallout < 0.0) hChazeFallout = 0.0;
//
//		//-------------------------------------------------------------------
//		// Calculate surface C flux from seafloor weathering and subduction
//		//-------------------------------------------------------------------
//
//		FCseafsubd = 0.0;
//		if (Tsurf > Tfreeze+0.01 && realtime > tstart && realtime <= tend) {
//
//			mix = Mocean/Mriver;
//
////			LplateRdg = 2.0*2.0*PI_greek*r_p; // Current length of mid-ocean ridges + continental rifts = length of subduction zones = 60000-65000 km (MOR) + 32000 km (rifts, Wong et al. 2022) = 87000 km.
			LplateRdg = pow(Ra/2.3e6,beta) * 2.0*2.0*PI_greek*r_p;
//			// 2.3e6 is canonical Ra for Earth today. Scaling with Ra^beta is consistent with scaling with Nu and also consistent with 3-5 times greater ridge length in Archean from Kadko et al. (1995).
//			// Today indicative size of 7 major plates is 70e6 km2, corresponding to diameter 2*sqrt(70e6/(4*pi)) = 4720 km > mantle depth, even though (isoviscous) convection cell should have aspect ratio of 1 (Bercovici et al. 2015).
//			// For Earth inputs it is = mantle depth (2900 km) for Ra = 1e7, i.e., 2 billion years ago (oldest evidence of plate tectonics 2-3 Ga; Brown et al. 2020)
//
//			zCrack = realtime/Myr2sec; // Reflects secular cooling below crust, which prevents cracks from healing as fast
////			zCrack = 6000.0;
//
////			volSeafCrust = (1.0-L)/(1.0-0.29) * LplateRdg*vConv*zCrack*dtime;
//			volSeafCrust = (1.0-L)/(1.0-0.29) * LplateRdg*(vConv*fmin(1.0,realtime/(1.5*Gyr2sec)))*zCrack*dtime; // vConv*fmin() term represents sluggish convection until full plate tectonics
//
//			Pseaf = rhoH2O * g[NR] * (r_p - pow(pow(r_p,3) - Mocean/rhoH2O/(4.0/3.0*PI_greek),1.0/3.0)) / bar2Pa; // Seafloor hydrostatic pressure: density*surface gravity*ocean depth TODO scale with land coverage
//			WRseafW = Mocean*dtime/tcirc / volSeafCrust; // Will be multiplied in AqueousChem() by (mass rock input to PHREEQC / seafloor crust density) = volume rock input to PHREEQC to get a mass of water for reaction with the mass of rock input to PHREEQC.
//
//			// Memorize aqueous C abundances
//			xaq0 = xaq[0];
//			xaq1 = xaq[1];
//
//			printf("Mixing river input into ocean...\n");
//			printf("Ocean will react with seafloor crust at W/R by vol. = %g\n", Mocean*dtime/tcirc / 1000.0 / volSeafCrust);
//			AqueousChem(path, "io/MixRiverOcean.txt", Tsurf, &Psurf, &Vatm, &nAir, &pH, &pe, &mix, &xgas, &xaq, &xriver, iResTime, 0, 0.0, 1, nvarEq, Pseaf, WRseafW, &deltaCreac, 0, dtime);
//			printf("deltaCreac = %g mol/kg\n", deltaCreac);
//
//			// Reset aqueous C abundances as they'll be adjusted by ocean-atmosphere equilibrium at the next time step, when netFC C is added to {atmosphere+ocean}
//			xaq[0] = xaq0;
//			xaq[1] = xaq1;
//
//			FCseafsubd = deltaCreac * Mocean / dtime; // Essentially independent of dtime because deltaCreac, for an ocean saturated in C before and after reaction, is essentially proportional to Mriver, i.e., dtime.
//		}                                             // There is a slight dependence on dtime if there was no saturation before reaction and because the solubility of C changes slightly with slight changes in ocean pressure, temperature, and composition, but it should be acceptable
//
//		//-------------------------------------------------------------------
//		// Compute net geo C fluxes
//		//-------------------------------------------------------------------
//		if (!staglid) farc = 0.25;
//		netFC = FCoutgas + (1.0-farc)*FCseafsubd; // FCcontw < 0, FCseafsubd < 0, ignore FCcontW since that's a flux to the ocean manifested in seafW

		dtime_old = dtime; // Memorize in case reservoirs need to be readjusted at the next time step so a negative netFC*dtime doesn't exceed xgas*nAir
		RCmantle_old = RCmantle;
		RCatmoc_old = RCatmoc;

		RCmantle = RCmantle - dtime*netFC; // Assuming plate = mantle (unlike Foley et al. 2015) and instantaneous mixing into the mantle once subducted (in practice could take 1 Gyr)
		if (RCmantle < 0.0) {
			printf("Mantle reservoir depleted\n");
			exit(0);
		}

		// Life. On Earth today, org C sequesters today 1e5 x more carbon than atmosphere. Org C has had same level since Archean, maybe was 2x lower back then (Martin et al. 2008 Fig. 2B)
//		if (realtime > tstart && realtime <= tend && RCorg < 1.0e21) { // Reservoir size is of order 1e21 mol today (Hayes and Waldbauer 2006 Table 2)
//			RCorg += 0.8*dtime*netFC;
//			netFC = 0.2*netFC;
//		}

		RCatmoc = RCatmoc + dtime*netFC; // Sum of atmospheric and ocean reservoirs, still needs partitioning
		if (RCatmoc < 0.0) {
			printf("Atmosphere-ocean reservoir depleted\n");
			exit(0);
		}

		// Write outputs
		if (realtime > tstart) {
			title[0] = '\0';
			if (cmdline == 1) strncat(title,path,strlen(path)-20);
			else strncat(title,path,strlen(path)-18);
			strcat(title,"Outputs/ReservoirsFluxes.txt");
			fout = fopen(title,"a");
			if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
			else {
				fprintf(fout, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n",
						realtime/Gyr2sec, Tmantle, RCmantle, RCatm, RCocean, RCatm+RCocean, RCatmoc, FCoutgas, FCcontW, FCseafsubd, netFC, xaq[3]*Mocean + xgas[3]*2.0*nAir, RCorg);
			}
			fclose (fout);

			title[0] = '\0';
			if (cmdline == 1) strncat(title,path,strlen(path)-20);
			else strncat(title,path,strlen(path)-18);
			strcat(title,"Outputs/CompoAtmosph.txt");
			fout = fopen(title,"a");
			if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
			else fprintf(fout, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n",
					realtime/Gyr2sec, xgas[0], xgas[1], xgas[2], xgas[3], xgas[4], Psurf, Tsurf, DeltaTghe, runoff/(1-fracEvap)*Yr2sec, nAir);
			fclose (fout);

			title[0] = '\0';
			if (cmdline == 1) strncat(title,path,strlen(path)-20);
			else strncat(title,path,strlen(path)-18);
			strcat(title,"Outputs/Outgassing.txt");
			fout = fopen(title,"a");
			if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
			else {
				fprintf(fout, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %d\n", realtime/Gyr2sec, Tmantle, Ra,
						nu*rho[(int)((NR+iCMB)/2.0)], k*(Tmantle-Tsurf)/zLith*1000.0, zLith/km2m, p, FCoutgas, vConv*1.0e-6*Myr2sec, tConv, zCrust, zCrust/dzCrust_dt, staglid);
			}
			fclose (fout);

			title[0] = '\0';
			if (cmdline == 1) strncat(title,path,strlen(path)-20);
			else strncat(title,path,strlen(path)-18);
			strcat(title,"Outputs/ContWeatherFluxes.txt");
			fout = fopen(title,"a");
			if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
			else {
				fprintf(fout, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n",
						realtime/Gyr2sec, L, tResLand/Yr2sec, FC_Mg, FC_Mg_sil, FC_Mg_carb, FC_Ca, FC_Ca_sil, FC_Ca_carb, FC_Ca_sulf, FC_Fe, FC_Fe_sil, FC_Fe_sulf, FC_Mg+FC_Ca+FC_Fe,
						xriver[iResTime][3], xriver[iResTime][17], xriver[iResTime][18], xriver[iResTime][19], xriver[iResTime][20], xriver[iResTime][21], xriver[iResTime][22], xriver[iResTime][23], xriver[iResTime][24], hChazeFallout);
			}
			fclose (fout);

			title[0] = '\0';
			if (cmdline == 1) strncat(title,path,strlen(path)-20);
			else strncat(title,path,strlen(path)-18);
			strcat(title,"Outputs/CompoOcean.txt");
			fout = fopen(title,"a");
			if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
			else {
				fprintf(fout, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n",
						realtime/Gyr2sec, Mocean, Mriver, volSeafCrust, pHout, pH, 4.0*(pe+pH)-logKO2H2O, rainpH, xaq[0], xaq[1], xaq[6], xaq[7], xaq[8], xaq[9], xaq[10], xaq[11], xaq[12]);
			}
			fclose (fout);
		}
	} // End time loop

	printf("\nExiting ExoCcycleGeo...\n");

	for (i=0;i<NR;i++) free (sys_tbl[i]);
	free (sys_tbl);
	for (i=0;i<kinsteps;i++) free (xriver[i]);
	free (xriver);
	free (input);
	free (title);
	free (xgas);
	free (xaq);
	free (r);
	free (rho);
	free (P);
	free (T);
	free (g);
	free (Meltfrac);

	Rf_endEmbeddedR(0);                                     // Close R and CHNOSZ

	kill(-1*pid, SIGKILL);
	return EXIT_SUCCESS;
}
