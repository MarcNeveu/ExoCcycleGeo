/*
 ============================================================================
 Name        : ExoCcycleGeo.c
 Author      : Marc Neveu
 Computes net C fluxes (in mol C m-2 s-1) at the surface-atmosphere interface
 of a terrestrial planet due to geophysical and geochemical processes.
 TODO:
 1- Lookup table based on PHREEQC calculations of kinetics of weathering for
 felsic and mafic rock, as a function of temperature and atmospheric carbon
 abundance
 2- Seafloor weathering on the timescale of emplacement of new crust. km3
 created per year * fraction getting within diffusional reach of seawater, if
 slow enough, let that determine the rate of consumption by assuming reaction
 has time to go to equilibrium. Confirm with PHREEQC simulations of rates.
 3- Redox of outgassing. MELTS?
 ============================================================================
 */

#include "ExoCcycleGeo.h"
#include "CHNOSZcmds.h"
#include "Compression.h"
#include "Geochem.h"
#include "Geodynamics.h"

//-------------------------------------------------------------------
// MAIN PROGRAM
//-------------------------------------------------------------------

int main(int argc, char *argv[]) {

	int i = 0;
	int ir = 0;
	int j = 0;
	int NR = 10000;                    // Number of radial grid zones used in determining planetary structure at setup. 10000 seems to achieve numerical convergence.
	int iLith = 0;                     // Index of grid zone that corresponds to base of lithosphere (geotherm inflexion to mantle adiabat)
	int imin = 0;                      // Index of innermost grid zone of melting in MELTS output
	int imax = 0;                      // Index of outermost grid zone of melting in MELTS output
	int iCMB = 0;                      // Index of grid zone that corresponds to core-mantle boundary
	int iBDT = 0;                      // Index of grid zone of the brittleo-ductile transition
	int itime = 0;                     // Time step counter
	int ntime = 0;                     // Total number of steps

	int iter = 0;                      // Iteration counter
	int niter = 10;                    // Number of iterations for ocean-atmosphere equilibrium

	int nslopeAvg = 5;                 // Target number of data points used for averaging melt fraction slope in the mantle
	int islope = 0;                    // Actual numbers of data points used for averaging melt fraction slope in the mantle

	double slope = 0.0;                // Average melt fraction slope in the mantle (bar-1)
	double meltfrac = 0.0;             // Melt fraction extrapolated at higher pressures than MELTS can handle

	double realtime = 0;               // Real time since birth of planetary system (s)

	// Planet parameters
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

	double *g = (double*) malloc((NR+1)*sizeof(double)); // Gravitational acceleration
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
	double RCatm = 0.0;                // Atmospheric C reservoir (mol)
	double RCocean = 0.0;              // Ocean C reservoir (mol)
	double RCatmoc = RCatm + RCocean;  // Combined atmospheric and ocean C reservoir (mol)

	// Fluxes
	double FCoutgas = 0.0;             // C flux from outgassing (subaerial+submarine) (mol s-1)
	double FCcontW = 0.0;              // C flux from continental weathering (mol s-1)
	double FCseafsubd = 0.0;           // C flux from seafloor weathering and subduction (mol s-1)
	double netFC = 0.0;                // Net surface C flux from all geological processes (mol s-1)
	double farc = 0.0;                 // Fraction of subducted C that makes it back out through arc volcanism (as opposed to into the mantle)

	// Geochem parameters
	double pH = 0.0;                   // pH of the surface ocean (default 8.22)
	double rainpH = 0.0;               // pH of rainwater
	double pe = 0.0;                   // pe (-log activity e-) corresponding to logfO2
	double logfO2 = 0.0;               // log O2 fugacity
	double logKO2H2O = 0.0;            // log K for reaction 4 H+ + 4 e- + O2 = 2 H2O, from CHNOSZ: subcrt(c("H+","e-","O2","H2O"),c(-4,-4,-1,2),c("aq","aq","g","liq"),T=25,P=1)
	double Tfreeze = 273.15;           // Temperature at which water freezes at surface pressure (K)

	// Atmosphere parameters
	double nAir = 0.0;                 // Number of mol in atmosphere (mol)
	double S = 0;                      // Stellar flux at planet (W m-2)
	double Teff = 0.0;                 // Effective temperature (K)
	double DeltaTghe = 0.0;            // Temperature increase over effective temperature due to greenhouse effect (K)
	double psi = 0.0;                  // log10(pCO2) (bar)
	double Vatm1 = 0.0;                // Initial atmospheric volume (m3), determined at second time step to avoid skew of starting values
	double Tsurf1 = 0.0;               // Surface temperature used to determine Vatm1 (K)
	double Vatm = 0.0;                 // Atmospheric volume (m3)

	// Kinetic parameters
	int kinsteps = 0;                  // Number of time steps of PHREEQC kinetic simulation
	double WRcontW = 0.0;              // Water:rock mass ratio for continental weathering reactions (kg/kg)
	double zContW = 0.0;               // Depth of continental weathering (m)
	double rhoContCrust = 0.0;         // Density of continental crust (kg m-3)
	double molmassContCrust = 0.0;     // Molar mass of continental crust (kg mol-1)
	double dCdtContW = 0.0;            // Carbon drawdown rate due to continental weathering (mol (kg H2O)-1 s-1)

	// Quantities to be computed by thermal/geodynamic model
	int staglid = 1;                   // 1 if stagnant-lid, 0 if mobile-lid
	int n_iter = 0, n_iter_max = 100;  // Iteration counter and max for the Bisection/Newton-Raphson loop to determine T_BDT
    double Tmantle = 0.0;              // Temperature at mid-mantle depth (K), initially set by accretion
	double H = 0.0;                    // Specific radiogenic heating rate (J s-1 kg-1)
	double x40K = 0.0;                 // Abundance of the radionuclide 40 K by mass
	double x232Th = 0.0;               // Abundance of the radionuclide 232 Th by mass
	double x235U = 0.0;                // Abundance of the radionuclide 235 U by mass
	double x238U = 0.0;                // Abundance of the radionuclide 238 U by mass
	double kappa = 0.0;                // Mantle thermal diffusivity (m2 s-1)
	double nu = 0.0;                   // Mantle kinematic viscosity (m2 s-1)
	double zCrust = 0.0;               // Crustal thickness (m), i.e. depth of layer crystallized from a melt
	double Ra = 0.0;                   // Rayleigh number for mantle convection (no dim)
	double Nu = 0.0;                   // Nusselt number (no dim)
	double Tref = 0.0;                 // Temperature at outer boundary of convective zone (surface or base of stagnant lid)
	double driveStress = 0.0;          // Driving stress at the base of lithosphere, used to switch between stagnant and mobile-lid modes
	double yieldStress = 0.0;          // Lithospheric yield stress, taken to be the strength at the brittle-ductile transition (brittle strength = ductile strength)
	double T_INF = 0.0, T_SUP = 0.0, T_TEMP = 0.0; // Intermediate temperatures used in determination of T_BDT by Bisection/Newton-Raphson loop
	double dT = 0.0, dTold = 0.0;      // Intermediate temperature changes used in determination of T_BDT by Bisection/Newton-Raphson loop
	double f_inf = 0.0, f_sup = 0.0;   // Differences between brittle and ductile strengths, used in determination of T_BDT by Bisection/Newton-Raphson loop (P)
	double f_x = 0.0, f_prime_x = 0.0; // Derivatives of the above (Pa K-1)
	double NewtRaphThresh = 1.0e5;     // Threshold for the Bisection/Newton-Raphson loop, here in Pa
	double T_BDT = 0.0;                // Temperature at brittle-ductile transition (K)
	double P_BDT = 0.0;                // Pressure at brittle-ductile transition (Pa)
	double meltmass = 0.0;             // Total mass of outgassing mantle melt (kg)
	double rhomelt = 0.0;              // Density of the melt (kg m-3)
	double zNewcrust = 0.0;            // Thickness of crust generated from mantle melt (km)
	double Tadiab = 0.0;               // Temperature in mantle adiabat (K)
	double Tupbnd = 0.0;               // Temperature in conductive upper boundary layer (K)
	double bndcoef = 0.0;              // Coefficient linking boundary layer thickness and (height of convection cell / Nusselt number)

	// Viscosity law constants (can't be #define'd because they are used as exponents)
	double KK08DryOlDiff[5];          // Dry diffusion creep flow law (Korenaga & Karato 2008)
	KK08DryOlDiff[0] = 261.0e3;       // Activation energy (J mol-1), default 261.0±28e3
	KK08DryOlDiff[1] = 6.0e-6;        // Activation volume (m3 mol-1), default 6±5e-6
	KK08DryOlDiff[2] = 5.25;          // Exponent of pre-exponential factor, default 5.25±0.03
	KK08DryOlDiff[3] = 1.0;           // Stress exponent
	KK08DryOlDiff[4] = 2.98;          // Grain size exponent, default 2.98±0.02

	double KK08WetOlDiff[5];		  // Wet diffusion creep flow law
	KK08WetOlDiff[0] = 387.0e3;       // Activation energy (J mol-1), default 387±53e3
	KK08WetOlDiff[1] = 25.0e-6;       // Activation volume (m3 mol-1), default 25±4e-6
	KK08WetOlDiff[2] = 4.32;          // Exponent of pre-exponential factor, default 4.32±0.38
	KK08WetOlDiff[3] = 1.0;           // Stress exponent
	KK08WetOlDiff[4] = 2.56;          // Grain size exponent, default 2.56±0.24

	double KK08DryOlDisl[5];          // Dry dislocation creep flow law (Korenaga & Karato 2008)
	KK08DryOlDisl[0] = 610.0e3;       // Activation energy (J mol-1), default 610±30e3
	KK08DryOlDisl[1] = 13.0e-6;       // Activation volume (m3 mol-1), default 13±8e-6
	KK08DryOlDisl[2] = 6.09;          // Exponent of pre-exponential factor, default 6.09±0.11
	KK08DryOlDisl[3] = 4.94;          // Stress exponent, default 4.94±0.05
	KK08DryOlDisl[4] = 0.0;           // Grain size exponent

	double KK08WetOlDisl[5];          // Wet dislocation creep flow law
	KK08WetOlDisl[0] = 523.0e3;       // Activation energy (J mol-1), default 523±100e3
	KK08WetOlDisl[1] = 4.0e-6;        // Activation volume (m3 mol-1), default 4±3e-6
	KK08WetOlDisl[2] = 0.6;           // Exponent of pre-exponential factor, default 0.6±0.5
	KK08WetOlDisl[3] = 3.60;          // Stress exponent, default 3.60±0.24
	KK08WetOlDisl[4] = 0.0;           // Grain size exponent

	double flowLawDiff[5];
	for (i=0;i<5;i++) flowLawDiff[i] = 0.0;

	double flowLawDisl[5];
	for (i=0;i<5;i++) flowLawDisl[i] = 0.0;

	const double grainSize = 1.0e3;	   // Grain size for computation of viscosity and creep laws (µm)

	// Quantities to be computed by melting model
	double vConv = 0.0;                // Convective velocity (m s-1)

	// Quantities to be computed by seafloor weathering model
	double zCrack = 0.0;               // Depth of fracturing below seafloor (m)
	double tcirc = 0.0;                // Time scale of hydrothermal circulation (s)
	double deltaCreac = 0.0;           // Net C leached/precipitated per kg of rock (mol kg-1)

	double *xgas = (double*) malloc(nAtmSpecies*sizeof(double));
	if (xgas == NULL) printf("ExoCcycleGeo: Not enough memory to create xgas[nAtmSpecies]\n"); // Mixing ratios by volume (or by mol since all gases are pretty much ideal and have the same molar volume) of atmospheric gases
    for (i=0;i<nAtmSpecies;i++) xgas[i] = 0.0;

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
	// Inputs
	//-------------------------------------------------------------------

	double dtime = 10.0*Myr2sec;     // Time step (s)

	ntime = (int) (5000.0*Myr2sec/dtime); // Number of time steps of simulation

	// Planet parameters
	double m_p = 1.0*mEarth;           // Planet mass (kg)
	double m_c = 0.325*m_p;            // Core mass (kg), default ≈0.3*m_p for Earth, Kite et al. (2009) use 0.325*m_p
	double L = 0.29;                   // Fraction of planet surface covered by land
	double Mocean = 1.4e21;            // Mass of ocean (kg, default: Earth=1.4e21)
    int redox = 2;                     // 1: current Earth surface, 2: hematite-magnetite, 3: fayalite-magnetite-quartz, 4: iron-wustite, code won't run with other values.
    double magmaCmassfrac = 0.002;     // Mass fraction of C in magmas. Default 0.004 = 0.4±0.25% H2O and CO2 in MORB and OIB parent magmas (Jones et al. 2018; HArtley et al. 2014; Hekinian et al. 2000; Gerlach & Graeber 1985; Anderson 1995)
	double zLith = 100.0*km2m;         // Depth of lithosphere (both thermal: geotherm inflexion to adiabat and mechanical: brittle-ductile transition) (m)
	double tConv = 10.0*Myr2sec;       // Convection timescale (s), default 200 Myr for deep Earth mantle today
	int rheology = 0;                  // 0 = dry olivine (KK08), 1 = wet olivine (KK08)
	int radionuclides = 1;             // 0 = Custom (LK07 lo), 1 = High (TS02), 2 = Intermediate (R91), 3 = Low (LK07), default = Intermediate (McDS95)

	// Atmospheric inputs
	double Tsurf0 = 288.0;             // Initial surface temperature (K)
	double S0 = 1368.0;                // Stellar flux for G2V star at 1 AU and 4.55 Gyr
	double albedo = 0.3;
	double Psurf = 1.0;                // Surface pressure (bar)
	double runoff = 0.7e-3/86400.0;    // Atmospheric runoff (m s-1), default runoff_0 = 0.665e-3 m day-1
	xgas[0] = 0.5;                     // CO2
    xgas[1] = 0.0;                     // CH4
    xgas[2] = 0.0;                     // O2
    xgas[3] = 0.5;                     // N2
    xgas[4] = 0.0;                     // H2O

    // Ocean inputs
	pH = 8.22;
    xaq[0] = 27.0/12.0/1000.0;         // 27 ppm C in today's oceans (scaled from 141 ppm Alk*M(C)/M(HCO3), M being molecular mass)
    xaq[1] = 0.0/12.0/1000.0;
    xaq[3] = 0.0/14.0/1000.0;          // ppm N in ocean

    // Continental weathering inputs
	double deplCcrit = 0.18;           // Chemical depletion fraction, i.e. how much carbon has been consumed from rain water by the time fresh rock is emplaced
	double kintime = 1.0e-6*Yr2sec;    // Total time of PHREEQC kinetic simulation (s), default: 1.0e-6 year
	kinsteps = 100;                    // 100 steps by default, could be more
	zContW = 10.0;
	rhoContCrust = 2740.0;
	molmassContCrust = 0.149;

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

	// Get pressure and density profiles with depth, accounting for compression and self-gravity
	compression(NR, m_p, m_c, Tsurf, 101, 107, 203, &r, &P, &rho, &g, &iCMB, path);

	r_p = r[NR];
	r_c = r[iCMB];
	kappa = k/(rho[(int)(NR+iCMB)/2]*Cp);
	rho_p = m_p / (4.0/3.0*PI_greek*pow(r_p,3));

	printf("\nPlanet radius = %g km, Core radius = %g km, Mid-mantle density = %g kg m-3, avg density = %g kg m-3\n", r_p/km2m, r_c/km2m, rho[(NR+iCMB)/2], m_p/(4.0/3.0*PI_greek*pow(r_p,3.0)));
	Tsurf = Tsurf0;
	Asurf = 4.0*PI_greek*r_p*r_p;
	nAir = Psurf*bar2Pa*Asurf/g[NR]/molmass_atm(xgas);

	vConv = (r_p-zLith-r[iCMB])/tConv;

	// Calculate Tmantle from accretion (Canup et al. 2020 Pluto chapter; volume averaging from leads to 3/5 term)
	Tmantle = Tsurf0 + 0.6*4.0*h_frac*chi*G*PI_greek*rho_p*r_p*r_p/(3.0*Cp);

	// Calculate temperatures (not computed in the core)
	for (i=iCMB;i<NR;i++) {
		Tadiab = Tmantle - alpha*g[(int)(0.5*(NR+iCMB))]*Tmantle/Cp * (r[i] - 0.5*(r_p+r_c)); // Mantle adiabat (Katsura et al. (2010) eq. 1 or TS02 eq. 4.254 at mid-mantle depth), technically stops at r[iCMB] + (zLith-zBDT) but that's OK.
		Tupbnd = Tmantle + (Tsurf - Tmantle)*erfc((r_p-r[i])/2.0*sqrt(vConv/(kappa*(r_p-r_c)/2.0))); // Upper mantle boundary, Turcotte & Schubert 2002 eq. 6.347 at midpoint between ascending and descending plumes
		if (Tadiab < Tupbnd) {
			T[i] = Tadiab;
			iLith = i+1;
			zLith = r_p-r[i+1];
		}
		else T[i] = Tupbnd;
	}

	// Initial guess: BDT = base of upper boundary layer
	T_BDT = T[iLith];
	P_BDT = P[iLith];

	// Rheology
	switch (rheology) {
	case 0: // Dry olivine rheology of Korenaga & Karato (2008)
		for (i=0;i<5;i++) {
			flowLawDiff[i] = KK08DryOlDiff[i];
			flowLawDisl[i] = KK08DryOlDisl[i];
		}
		break;
	case 1: // Wet olivine rheology of Korenaga & Karato (2008)
		for (i=0;i<5;i++) {
			flowLawDiff[i] = KK08WetOlDiff[i];
			flowLawDisl[i] = KK08WetOlDisl[i];
		}
		break;
	default:
		printf("ExoCcycleGeo: Choice of rheology is not recognized. Exiting.\n");
		exit(0);
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

	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("ExoCcycleGeo v20.5\n");
	printf("This code is in development and cannot be used for science yet.\n");
	if (cmdline == 1) printf("Command line mode\n");
	printf("-------------------------------------------------------------------\n");

	printf("\n");
	printf("Inputs:\n");
	printf("Planet mass \t \t \t \t %g M_Earth\n", m_p/mEarth);
	printf("Land coverage of planet surface \t %g%%\n", L*100.0);
	printf("Surface temperature \t \t \t %g K\n", Tsurf);
	printf("Mantle temperature \t \t \t %g K\n", Tmantle);
	printf("Surface pressure \t \t \t %g bar\n", Psurf);
	printf("Atmospheric CO2 molar mixing ratio \t %g\n", xgas[0]);
	printf("Atmospheric CH4 molar mixing ratio \t %g\n", xgas[1]);
	printf("Atmospheric O2 molar mixing ratio \t %g\n", xgas[2]);
	printf("Atmospheric N2 molar mixing ratio \t %g\n", xgas[3]);
	printf("Ocean pH \t \t \t \t %g \n", pH);
	printf("Surface runoff \t \t \t \t %g mm d-1\n", runoff);
	printf("-------------------------------------------------------------------\n");

	printf("\n");
	printf("Computing geo C fluxes through time...\n");

	//-------------------------------------------------------------------
	// Choose redox state (from most oxidized to most reduced)
	//-------------------------------------------------------------------

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
	// Initialize reservoirs
	//-------------------------------------------------------------------

	RCatm = (xgas[0]+xgas[1])*nAir;

	// Equilibrate ocean and atmosphere at input pressure and atmospheric C/N ratio
	if (Psurf < 0.01) {
		printf("ExoCcycleGeo: Pressure = %g bar too close to the triple point of H2O, oceans not stable at the surface. Exiting.\n", Psurf);
		exit(0);
	}

	printf("Equilibrating ocean and atmosphere at input pressure and atmospheric C/N ratio...\n");
	for (i=0;i<nAtmSpecies;i++) {
		if (xgas[i] > 0.0 && xaq[i] == 0.0) xaq[i] = xgas[i]; // xaq must be >0 otherwise PHREEQC ignores it, set to xgas (initial guess).
	}

	AqueousChem(path, "io/OceanStart.txt", itime, Tsurf, &Psurf, &Vatm, &nAir, &pH, &pe, &Mocean, &xgas, &xaq, NULL, 1, 0.0, 1, nvarEq);

	RCocean = (xaq[0]+xaq[1])*Mocean;
	RCatmoc = RCatm + RCocean;
	RCmantle = 10.0*(RCatmoc);

	// Print first line of outputs
	title[0] = '\0';
	if (cmdline == 1) strncat(title,path,strlen(path)-20);
	else strncat(title,path,strlen(path)-18);
	strcat(title,"Output.txt");
	fout = fopen(title,"w");
	if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
	else {
		fprintf(fout, "'Reservoirs in mol, fluxes in mol s-1'\n");
		fprintf(fout, "'Time (Myr)' \t 'Tmantle (K)' \t RCmantle \t RCatmoc \t RCocean \t FCoutgas \t FCcontw \t FCseafsubd \t 'Net C flux'\n");
		fprintf(fout, "Init \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \n", Tmantle, RCmantle, RCatmoc, RCocean, FCoutgas, FCcontW, FCseafsubd, netFC);
	}
	fclose (fout);

	title[0] = '\0';
	if (cmdline == 1) strncat(title,path,strlen(path)-20);
	else strncat(title,path,strlen(path)-18);
	strcat(title,"CompoOceanAtm.txt");
	fout = fopen(title,"w");
	if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
	else {
		fprintf(fout, "'Atmospheric species in mixing ratio (by mol, i.e. by volume for ideal gases), aqueous species in mol/(kg H2O), total C and N in mol'\n");
		fprintf(fout, "'Time (Myr)' \t CO2(g) \t CH4(g) \t O2(g) \t N2(g) \t H2O(g) \t 'P_surface (bar)' \t 'T_surface (K)' \t 'DeltaT GHE (K)' \t 'Ox C(aq)' \t 'Red C(aq)' \t 'Total N(aq)' \t "
				"'Ocean pH' \t 'Ocean log f(O2) at Tsurf0' \t 'Rain pH' \t 'Total C g+aq' \t 'Total N g+aq' \t 'Mocean (kg)' \t 'nAir (mol)'\n");
		fprintf(fout, "Init \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n",
				xgas[0], xgas[1], xgas[2], xgas[3], xgas[4], Psurf, Tsurf, DeltaTghe, xaq[0], xaq[1], xaq[3],
				pH, 4.0*(pe+pH)-logKO2H2O, rainpH, (xaq[0]+xaq[1])*Mocean + (xgas[0]+xgas[1])*nAir, xaq[3]*Mocean + xgas[3]*2.0*nAir, Mocean, nAir);
	}
	fclose (fout);

	//-------------------------------------------------------------------
	// Start time loop
	//-------------------------------------------------------------------

	double TotalN0 = xgas[3]*2.0*nAir + xaq[3]*Mocean;

	printf("Starting time loop...\n");
	for (itime = 0;itime<ntime;itime++) {

		realtime = (double)itime*dtime;                // Start at birth of planetary system
//		realtime = (double)itime*dtime + 4.55*Gyr2sec; // Start at present day

//		netFC = 0.0; // !! Debug

		// Print relative change in Ntot to ensure N is conserved. Also print net C flux relative to outgassing flux; if small the C cycle should be balanced
//		printf("Iteration %d/%d, Time: %g Myr, Tsurf: %g K, pCO2: %g bar, Delta_Ntot/Ntot: %.2f pct, NetFC/outgas: %.2f pct\n",
//				 itime, ntime, (double)itime*dtime/Myr2sec, Tsurf, xgas[0]*Psurf,
//				 (1.0-(xgas[3]*2.0*nAir + xaq[3]*Mocean)/TotalN0)*100.0, netFC/FCoutgas*100.0);

		//-------------------------------------------------------------------
		// Update surface temperature (unnecessary once coupled to Atmos)
		//-------------------------------------------------------------------

		// Parameterization of Caldeira & Kasting (1992) equation (4), valid between psi=-8 to -2 and T=0 to 100C
//		S = S0/(1.0-0.38*(realtime/Gyr2sec/4.55-1.0));
		S = S0;
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

//		Tsurf = Teff + DeltaTghe;

		if (Tsurf < Tfreeze+0.01) printf("ExoCcycleGeo: Surface temperature = %g K < 0.01 C. DeltaTghe=%g K.\n", Tsurf, DeltaTghe);

		//-------------------------------------------------------------------
		// Update atmosphere
		//-------------------------------------------------------------------

		if (itime < 2) { // First 2 time steps: re-compute Vatm to avoid it being too influenced by ballpark initial conditions
			Vatm1 = nAir*R_G*Tsurf/(Psurf*bar2Pa);
			Tsurf1 = Tsurf;
		}
		Vatm = Vatm1*Tsurf/Tsurf1; // After 2nd time step, keep Vatm constant

		// Redox of outgassing
		if (redox <= 4) {
			xgas[0] = (xgas[0]*nAir + 1.0*dtime*netFC)/(nAir + dtime*netFC); // CO2 dominates over CH4, assume 100% added gas is CO2 and let equilibration with ocean speciate accurately
			xgas[1] = xgas[1]*nAir/(nAir + dtime*netFC);                     // Dilute CH4
		}
		else {
			xgas[1] = (xgas[1]*nAir + 1.0*dtime*netFC)/(nAir + dtime*netFC); // CH4 dominates over CO2, assume 100% added gas is CH4
			xgas[0] = xgas[0]*nAir/(nAir + dtime*netFC);                     // Dilute CO2
		}

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

		if (Tsurf > Tfreeze) {

			AqueousChem(path, "io/OceanDiss.txt", itime, Tsurf, &Psurf, &Vatm, &nAir, &pH, &pe, &Mocean, &xgas, &xaq, NULL, 0, 0.0, 1, nvarEq);

			if (Psurf < 0.01) {
				printf("ExoCcycleGeo: Pressure = %g bar too close to the triple point of H2O, oceans not stable at the surface. Exiting.\n", Psurf);
				exit(0);
			}

			cleanup(path); // Remove PHREEQC selected output file
		}

		//-------------------------------------------------------------------
		// Calculate surface C flux from outgassing
		//-------------------------------------------------------------------

		// 1. Determine tectonic mode
		// If convective driving stress > lithospheric yield stress, mobile-lid (plate tectonics) regime. Otherwise, stagnant lid regime (O'Neill & Lenardic 2007).

		// 1a. Determine T and P at the brittle-ductile transition (base of lithosphere)
		// Solve for brittleStrength-ductileStrength = f(P,T) = f(T) = 0
		// Binary search polished with Newton-Raphson (see Numerical Recipes and Neveu et al. 2015 Icarus Appendix)

		// Initialize bounds for T_BDT
		T_INF = Tsurf;
		T_SUP = Tmantle;

		// Ensure that f(T_INF)<0 and f(T_SUP)>0
		f_inf = brittleDuctile(T_INF, P_BDT, flowLawDiff, flowLawDisl, grainSize, tConv);
		f_sup = brittleDuctile(T_SUP, P_BDT, flowLawDiff, flowLawDisl, grainSize, tConv);

		if (f_inf*f_sup > 0.0) {
			if (f_sup < 0.0) printf("ExoCcycleGeo: Brittle strength < ductile strength, i.e. brittle regime even at Tmantle=%g K.\n"
					                 "Brittle-ductile transition could not be determined.\n", Tmantle); // Brittle regime in mantle
			else P_BDT = Psurf*bar2Pa; // Ductile regime at surface
		}
		else {
			if (f_inf > 0.0) {         // Swap INF and SUP if f_inf > 0 and f_sup < 0
				T_TEMP = T_INF;
				T_INF = T_SUP;
				T_SUP = T_TEMP;
			}
			T_BDT = 0.5*(T_INF+T_SUP); // Initialize the guess for the root T_BDT,
			dTold = fabs(T_INF-T_SUP); // Initialize the "stepsize before last"
			dT = dTold;                // Initialize the last stepsize

			f_x = brittleDuctile(T_BDT, P_BDT, flowLawDiff, flowLawDisl, grainSize, tConv);
			f_prime_x = brittleDuctile_prime(T_BDT, P_BDT, Tsurf, Psurf*bar2Pa, flowLawDiff, flowLawDisl, grainSize, tConv);

			// Loop over allowed iterations to find T_BDT that is a root of f
			n_iter = 0;
			while (fabs(f_x) > NewtRaphThresh) {

				// Bisect if Newton is out of range, or if not decreasing fast enough
				if ((((T_BDT-T_SUP)*f_prime_x-f_x)*((T_BDT-T_INF)*f_prime_x-f_x) > 0.0) || (fabs(2.0*f_x) > fabs(dTold*f_prime_x))) {
					dTold = dT;
					dT = 0.5*(T_SUP-T_INF);
					T_BDT = T_INF + dT;
				}
				else { // Do Newton-Raphson
					dTold = dT;
					dT = f_x/f_prime_x;
					T_BDT = T_BDT - dT;
				}
				// Calculate updated f and f'
				f_x = brittleDuctile(T_BDT, P_BDT, flowLawDiff, flowLawDisl, grainSize, tConv);
				f_prime_x = brittleDuctile_prime(T_BDT, P_BDT, Tsurf, Psurf*bar2Pa, flowLawDiff, flowLawDisl, grainSize, tConv);

				if (f_x < 0.0) T_INF = T_BDT; // Maintain the bracket on the root
				else T_SUP = T_BDT;

				n_iter++;
				if (n_iter>=n_iter_max) {
					printf("ExoCcycleGeo: could not find the brittle-ductile transition after %d iterations\n",n_iter_max);
					break;
				}
			}

			for (i=NR;i>=0;i--) {
				if (T[i-1] > T_BDT && T[i] <= T_BDT) break;
			}
			iBDT = i;
			P_BDT = P[i];
		}

		// 1b. Determine yield stress, equated to brittle strength at BDT (equivalently, ductile strength)
		if (P_BDT < 200.0e6) yieldStress = 0.85*P_BDT;
		else yieldStress = 0.6*P_BDT + 50.0e6;

		// 1c. Determine convective drive stress
		// Stress = viscosity * strain rate (assumed Newtonian at base of crust - see Deschamps & Sotin 2000, confirmed because diffusion creep dominates over
		// dislocation creep at base of crust), i.e. stress = viscosity * vConv/δ, w/ δ: boundary layer thickness
		// δ = thickness of thermal gradient = 2.32*(kappa*δ/vConv)^0.5 (Turcotte & Schubert 2002, eq. 6.327 with δ/vConv = conduction timescale)
		// So δ^0.5 = 2.32*(kappa/vConv)^0.5 and δ ≈ 5.38*kappa/vConv

		// Mid-mantle:
		nu = combVisc(Tmantle, P[(int)((iBDT+iCMB)/2.0)], flowLawDiff, flowLawDisl, grainSize, tConv)/rho[(int)((iBDT+iCMB)/2)];
		// BDT:
//		nu = combVisc(Tmantle - alpha*g[iBDT]*Tmantle/Cp * (r[iBDT] - 0.5*(r_p+r_c)), P[iBDT], flowLawDiff, flowLawDisl, grainSize, tConv)/rho[iBDT];
		// Lower mantle:
//		nu = 1.0e16*exp((2.0e5 + P[(int)((iBDT+iCMB)/2.0)]*1.1e-6)/(R_G*Tmantle))/rho[(int)((iBDT+iCMB)/2)]; // Cízková et al. (2012)
		// CMB:
//		nu = 1.0e16*exp((2.0e5 + P[iCMB]*1.1e-6)/(R_G*(Tmantle - alpha*g[iCMB]*Tmantle/Cp * (r[iCMB] - 0.5*(r_p+r_c)))))/rho[iCMB]; // Cízková et al. (2012)

		// Compute vigor of convection, assuming whole-mantle convection (Kite et al. 2009)
		if (staglid) Tref = T_BDT; // Instead of scaling used by Kite et al. (2009, equation 8)
		else Tref = Tsurf;

		Ra = alpha*g[(int)((iBDT+iCMB)/2.0)]*(Tmantle-Tref)*pow(r[iBDT]-r_c,3.0)/(kappa*nu);
		Nu = pow(Ra/Ra_c, beta);

		// Convective velocity
//		vConv = 2.0*(Nu-1.0) * k / (rho[(int)((iBDT+iCMB)/2)]*Cp*(r[iBDT]-r_c)) * (Tmantle-Tsurf) / (Tmantle-Tref); // Kite et al. 2009 equation 24
		vConv = 0.354*kappa/(r[iBDT]-r[iCMB])*sqrt(Ra); // Eq. (6.379) of Turcotte & Schubert (2002), p. 514
		driveStress = nu * rho[(int)((iBDT+iCMB)/2)] * vConv*vConv / (5.38*kappa);
		tConv = (r[iBDT]-r[iCMB])/vConv; // Update convection timescale

		// 1d. Determine tectonic regime
		if (yieldStress < driveStress) staglid = 0; // Mobile-lid regime, i.e. plate tectonics, also requires surface liquid water to weaken the lithosphere.
		else staglid = 1;                           // Stagnant-lid regime

		// ------------------------------------
		// 2. Melting & outgassing model (Kite et al. 2009)

		// 2a. Calculate temperatures
		for (i=0    ;i<iCMB ;i++) T[i] = 0.0; // Temperatures not computed in the core
		for (i=iCMB ;i<NR;i++) {
			Tadiab = Tmantle - alpha*g[(int)(0.5*(iBDT+iCMB))]*Tmantle/Cp * (r[i] - 0.5*(r[iBDT]+r_c)); // Mantle adiabat, technically stops at r[iCMB] + (zLith-zBDT) but that's OK
			Tupbnd = Tmantle + (Tsurf - Tmantle)*erfc((r_p-r[i])/2.0*sqrt(vConv/(kappa*(r[iBDT]-r_c)/2.0))); // Upper mantle boundary, Turcotte & Schubert 2002 eq. 6.347 at midpoint between ascending and descending plumes
			if (Tadiab < Tupbnd) {
				T[i] = Tadiab;
				iLith = i+1;
				zLith = r_p-r[i+1];
			}
			else T[i] = Tupbnd;
		}
//		for (i=iBDT;i<NR;i++) T[i] = Tsurf + (T_BDT-Tsurf)*(P[i]-Psurf)/(P_BDT-Psurf); // Assumes linear relationship between P and T

		bndcoef = (zLith-r[iBDT]) / (pow(Ra/Ra_c,-beta)*(r_p-r[iBDT]-r_c)); // Turcotte & Schubert 2002 eq. 6.387 and Fig. 6.39; Shi et al. (2012) eq. 1
//		if (bndcoef < 0.15 || bndcoef > 0.5) printf("ExoCcycleGeo: convective boundary layer thickness coefficient %g outside common thickness bounds 0.15 to 0.5\n", bndcoef); TODO put back

		// 2b. Update alphaMELTS input files
		// Output P-T profile in lithosphere+boundary layer to be fed into alphaMELTS through PTexoC.txt
		// MELTS crashes below solidus for default composition, but that's enough to capture all depths above BDT at which melting occurs.

		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"alphaMELTS-1.9/ExoC/PTexoC.txt");

		fout = fopen (title,"w");
		if (fout == NULL) printf("ExoCcycleGeo: Missing PT file path: %s\n", title);

		ir = 0;
		for (i=iCMB;i<NR;i++) {
			if (P[i] < 1.0e10) { // 10 GPa is above MELTS upper limit.
				if (T[i]-Kelvin > 750.0) {
					fprintf(fout, "%g %g\n", P[i]/bar2Pa, T[i]-Kelvin); // Don't input below 950 C to avoid MELTS crashing.
					ir++;
				}
				else break;
			}
		}
		fclose(fout);

		// 2c. Run alphaMELTS to compute melt fraction along P-T profile in lithosphere and mantle
		// Reset sys_tbl
		for (i=0;i<NR;i++) {
			for (j=0;j<18;j++) sys_tbl[i][j] = 0.0;
		}

		// Call alphaMELTS
		if (realtime > 4.0*Gyr2sec && realtime <= 4.57*Gyr2sec) alphaMELTS(path, 0, ir, "ExoC/ExoC_env.txt", &sys_tbl);

		imin = 0;
		imax = 0;
		j = 0;
		for (i=iCMB+1;i<NR;i++) {
			Meltfrac[i] = 0.0;
			if (P[i] < 1.0e10) { // 10 GPa is above MELTS upper limit.
				if (sys_tbl[j][4] > 0.0 && (sys_tbl[j][13] < 1.5*sys_tbl[j][14] || sys_tbl[j][4] == 1.0)) { // Don't store if volume melt fraction negative or if density of melt > n*density of solid
					if (T[i]-Kelvin > 750.0) {
						if (fabs(1.0-P[i]/(sys_tbl[j][0]*bar2Pa)) > 1.0e-3) printf("ExoCcycleGeo: pressures from grid (%g bar) and MELTS (%g bar) misaligned at grid index %d\n", P[i]/bar2Pa, sys_tbl[j][0], i);
						if (fabs(1.0-T[i]/ sys_tbl[j][1]        ) > 1.0e-3) printf("ExoCcycleGeo: temperatures from grid (%g K) and MELTS (%g K) misaligned at grid index %d\n", T[i], sys_tbl[j][1], i);
						Meltfrac[i] = sys_tbl[j][3]; // Don't input below 950 C to avoid MELTS crashing.
						rhomelt = sys_tbl[j][13]*1000.0; // g cm-3 to kg m-3, taken at depth of highest melt fraction, for calculation of thickness of new crust generated
						if (imin == 0) imin = i;
					}
					imax = i;
				}
				j++;
			}
		}

		meltmass = 0.0;
		slope = 0.0;
		ir = 0;
		if (imax > imin) {
			if (Meltfrac[imin] > 0.1) {
				printf("ExoCcycleGeo: alphaMELTS could only calculate melting down to depth %g km (melt fraction %.2g > 0.1), extrapolating melting curve to 0 linearly with depth\n", (r_p-r[imin])/km2m, Meltfrac[imin]);
				// Extrapolate starting from last 5 indices, provided melt fraction of all is < 1 (otherwise, use less indices)
				if (imax-imin > 5) { // Only do this if there are > 5 grid zones of melt, otherwise chances are the slope will be skewed.
					for (i=imin;i<imin+nslopeAvg;i++) {
						if (P[i+1]-P[i] == 0) printf("ExoCcycleGeo: Can't extrapolate rock melting curve with pressure because denominator is zero\n");
						else {
							if (Meltfrac[i] < Meltfrac[i+1]) { // Ensure we're below depth of max melt
								slope = slope + (Meltfrac[i+1]-Meltfrac[i])/(P[i+1]-P[i]);
								islope++;
							}
						}
					}
					slope = slope/(double) islope;

					if (islope) {
						for (i=imin-1;i>iCMB;i--) {
							meltfrac = Meltfrac[i+1] + slope*(P[i]-P[i+1]);
							if (meltfrac > 0.0) Meltfrac[i] = meltfrac;
							else break;
							ir++;
						}
					}
				}
			}
		}

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
//			// Ensure that f(T_INF)<0 and f(T_SUP)>0 TODO Ductile strength depends on time step (deps/dt). Make time step-independent?
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
//				Tsolidus = 0.5*(T_INF+T_SUP); // Initialize the guess for the root T_BDT,
//				dTold = fabs(T_INF-T_SUP); // Initialize the "stepsize before last"
//				dT = dTold;                // Initialize the last stepsize
//
//				f_x = Psolidus(Tgeotherm[i]) - Pgeotherm[i];
//				f_prime_x = dPsolidusdT(Tgeotherm[i]);
//
//				// Loop over allowed iterations to find T_BDT that is a root of f
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
//			if (P[ir] > P_BDT && P[ir] < 100000.0*bar2Pa) {
//				// Adiabat
//				T[ir] = T_BDT + gsurf*alpha*Tmantle/Cp*(r_p-zLith-r[ir]);
//
//				// Determine Tsolidus
//
//				// Initialize bounds for Tsolidus
//				T_INF = Tsurf;
//				T_SUP = 10000.0;
//
//				// Ensure that f(T_INF)<0 and f(T_SUP)>0 TODO Ductile strength depends on time step (deps/dt). Make time step-independent?
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
//					Tsolidus = 0.5*(T_INF+T_SUP); // Initialize the guess for the root T_BDT,
//					dTold = fabs(T_INF-T_SUP); // Initialize the "stepsize before last"
//					dT = dTold;                // Initialize the last stepsize
//
//					f_x = Psolidus(T[ir]) - P[ir];
//					f_prime_x = dPsolidusdT(T[ir]);
//
//					// Loop over allowed iterations to find T_BDT that is a root of f
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

		if (imax > imin) {

			printf("Pressure (bar) \t Depth (km) \t Melt fraction \t Temp (K) \t Density (kg m-3)\n");
			for (i=imin-ir-100;i<=NR;i++) printf("%g \t %g \t %g \t %g \t %g\n", P[i]/bar2Pa, (r_p-r[i])/km2m, Meltfrac[i], T[i], rho[i]);

			for (i=imin-ir;i<=imax;i++) {
				meltmass = meltmass + Meltfrac[i]*4.0/3.0*PI_greek*(pow(r[i+1],3)-pow(r[i],3))*rho[i];
			}
		}

		if (rhomelt > 0.0 && r[iBDT] > r[iCMB]) zNewcrust = (r_p - pow(pow(r_p,3) - meltmass/rhomelt/(4.0/3.0*PI_greek),1.0/3.0)) * vConv/(r[iBDT]-r[iCMB]);
		else zNewcrust = 0.0;
		zCrust = zCrust + zNewcrust*dtime;

		FCoutgas = meltmass*magmaCmassfrac/0.044*vConv/(r[iBDT]-r[iCMB]); // mol C s-1

//      // Alternative: Kite et al. (2009) eq. 25 They didn't scale with mass: [sum melt fraction (depth)] * [mass (depth)] / [total mass between surf and Psolidus]. Also their typo: P0>Pf=P_BDT.
//		double Rmelt = 0.0;             // Rate of melt generation (m-2 s-1)
//		if (staglid) Rmelt = Asurf*vConv*rhoMagma; // *  (rhoCrust*gsurf*zCrust/(Psolidus-Psurf)); // kg s-1
//		// Enhancement of rate of melting due to plate tectonics. Focuses solely on mid-ocean ridges (their Sec. 2.3). Also a function of age. See Sasaki & Tajika 1995.
//		else Rmelt = vConv*MORlength*zCrust*rhoCrust; // MORlength presumably changes over time (Sasaki & Tajika 1995)
//
//		// Includes extrusive and intrusive melting because both degas (K09).
//		FCoutgas = deltaCvolcEarth * Rmelt/(RmeltEarth*mEarth); // mol C s-1. Accurate for Earth magma C concentrations, but for other mantle concentrations, should be calculated explicitly (MELTS?)

//		printf("\n");
//		printf("Time: %g Gyr\n", realtime/Gyr2sec);
//		printf("---------------------------------------------------------------------------\n");
//		printf("Quantity                                | Earth benchmark   | Model result \n");
//		printf("---------------------------------------------------------------------------\n");
//		printf("New crust generation rate (m Myr-1)     | 40                | %.3g \n", zNewcrust*Myr2sec);
//		printf("New crust density (kg m-3)              | 2800              | %.4g \n", rhomelt);
//		printf("C and H2O outgassing rate (mol s-1)     | 115000            | %.5g \n", FCoutgas);
//		printf("Magma C mass fraction                   | 0.1-0.65%%         | %.3g%% \n", magmaCmassfrac*100.0);
//		printf("Convective boundary layer thickness coef| 0.15 to 0.5       | %.2g \n", bndcoef);
//		printf("Convective velocity (cm yr-1)           | ~1                | %.2g \n", vConv*100.0*1.0e-6*Myr2sec);
//		printf("Mantle convection timescale (Myr)       | ~50-200           | %.4g \n", tConv/Myr2sec);
//		printf("Rayleigh number                         | ~1e6-1e7          | %.2g \n", Ra);
//		printf("Mantle adiabatic gradient (K km-1)      | 0.5               | %.3g \n", alpha*g[(int)((iBDT+iCMB)/2)]*Tmantle/Cp * km2m);
//		printf("Half-mantle depth (km)                  | 1450              | %.4g \n", (r_p-r_c-zLith)/2.0/km2m);
//		printf("Heat flux (mW m-2)                      | 86 (radiodecay 40)| %.2g \n", k*(Tmantle-Tsurf)/zLith*1000.0);
//		printf("---------------------------------------------------------------------------\n");
//		printf("Mid-mantle temperature: %g K \t Brittle-ductile transition temperature: %.3g K \t Temperature at base of lithosphere: %.4g K\n", Tmantle, T_BDT, T[iLith]);
//		printf("Depth of brittle-ductile transition: %.3g km \t Thickness of lithosphere: %.3g km \n", (r_p-r[iBDT])/km2m, zLith/km2m);
//		printf("Convective driving stress: %.4g MPa \t Lithospheric yield stress: %.4g MPa \t", driveStress/MPa2Pa, yieldStress/MPa2Pa);
//		if (staglid) printf("Stagnant lid\n");
//		else printf ("Plate tectonics\n");
//		printf("\n");

		// ------------------------------------
		// 3. Mantle thermal evolution
		// Instantaneous heating rate (Kite et al. 2009 Table 1): H = X_4.5 * W * exp(ln(1/2) / t1/2 * (t-4.5))
		H =   x40K * 2.92e-5 * exp(log(0.5)/( 1.26 *Gyr2sec) * (realtime - 4.5*Gyr2sec))  //  40-K
		  + x232Th * 2.64e-5 * exp(log(0.5)/(14.0  *Gyr2sec) * (realtime - 4.5*Gyr2sec))  // 232-Th
		  +  x235U * 56.9e-5 * exp(log(0.5)/( 0.704*Gyr2sec) * (realtime - 4.5*Gyr2sec))  // 235-U
		  +  x238U * 9.46e-5 * exp(log(0.5)/( 4.47 *Gyr2sec) * (realtime - 4.5*Gyr2sec)); // 238-U
		// Effective thermal conductivity scaled with Nu
		Tmantle = Tmantle + dtime*(H/Cp - kappa*Nu*(Tmantle-Tref)/pow((r[iBDT]-r_c)/2.0,2));

		printf("%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n", realtime/Gyr2sec, Tmantle, Ra, nu*rho[(int)((iBDT+iCMB)/2)], k*(Tmantle-Tsurf)/zLith*1000.0, zLith/km2m, bndcoef, FCoutgas, vConv*1.0e-6*Myr2sec);

		//-------------------------------------------------------------------
		// Calculate surface C flux from continental weathering
		//-------------------------------------------------------------------

		if (Tsurf > Tfreeze) {
			// Analytical calculation from Edson et al. (2012) Eq. 1; Abbot et al. (2012) Eq. 2
//			FCcontW = -L * 0.5*deltaCcontwEarth*Asurf * pow(xgas[0]/xCO2g0,0.3) * runoff/runoff_0 * exp((Tsurf-TsurfEarth)/17.7);

//			double **xrain = (double**) malloc(nAqSpecies*sizeof(double*));
//			if (xrain == NULL) printf("ExoCcycleGeo: Not enough memory to create xrain[nAqSpecies]\n"); // Molalities of aqueous species in rain at different times (mol (kg H2O)-1, rain[0][kinstep]: time in s)
//		    for (i=0;i<nAqSpecies;i++) {
//		    	xrain[i] = (double*) malloc(kinsteps*sizeof(double));
//		    	if (xrain[i] == NULL) printf("ExoCcycleGeo: Not enough memory to create xrain[nAqSpecies][kinsteps]\n");
//		    }
//
//		    kintime = 1.0e-6*Yr2sec; // Reset to very small time span
//		    double tol = 1.5; // Tolerance factor for deplCcrit
//			for (iter=0;iter<niter;iter++) {
//			    for (i=0;i<nAqSpecies;i++) {
//			    	for (j=0;j<kinsteps;j++) xrain[i][j] = 0.0;
//			    }
//
//				WRcontW = runoff*kintime/zContW*1000.0/rhoContCrust*molmassContCrust;
//
//				AqueousChem(path, "io/ContWeather.txt", itime, Tsurf, &Psurf, 0.0, &pH, &pe, WRcontW, &xgas, &xaq, &xrain, 1, kintime, kinsteps, nvarKin);
//
//				rainpH = xrain[1][1];
//
//				for (i=2;i<kinsteps-1;i++) {
//					if (xrain[2][i] == 0.0) break; // PHREEQC did not return a result (sim interrupted)
//					xrain[3][i] = 1.0-xrain[2][i]/xrain[2][1]; // Chemical depletion fraction
//					printf("%d\t Time: %g yr\t pH: %g\t C(aq): %g\t deplC: %g\n", i, xrain[0][i]/Yr2sec, xrain[1][i], xrain[2][i], xrain[3][i]);
//				}
//
//				if (xrain[3][2] > deplCcrit*tol) {
//					kintime = 0.1*kintime;
//					printf("Continental weathering: Time span too coarse, decreasing kinetic time span 10-fold to %g years...\n", kintime/Yr2sec);
//				}
//				else if (i > 2 && xrain[3][i-1] < deplCcrit/tol) { // Last step for which PHREEQC returned a result (kinsteps-1 unless sim interrupted)
//					// Roughly estimate drawdown rate anyway (scaled down from closest determined value), in case next iteration fails
//					xrain[4][i] = (xrain[2][i-1] - xrain[2][i-2])/kintime*(double)kinsteps;
//					dCdtContW = sqrt(xrain[4][i-1]*xrain[4][i-2])*xrain[3][i-1]/deplCcrit;
//					kintime = 10.0*kintime;
//					printf("Continental weathering: Didn't go far enough in time, increasing time span of simulation 10-fold to %g years...\n", kintime/Yr2sec);
//				}
//				else if (i == 2) {
//					kintime = 10.0*kintime; // Increase time step anyway
//					printf("PHREEQC kinetic simulation failed at first step. Increasing time span of simulation 10-fold to %g years...\n", kintime/Yr2sec);
//					kintime++;
//					// break;
//				}
//				else break;
//			}
//
//			if (iter == niter) printf("Could not get timescale of continental weathering after %d iterations. FCcontW=%g not accurately updated at this time step.\n", iter, FCcontW);
//			else {
//				dCdtContW = 0.0;
//				for (j=2;j<i-1;j++) {
//					xrain[4][j] = (xrain[2][j-1] - xrain[2][j-2])/kintime*(double)kinsteps; // Instantaneous C drawdown rate (mol (kg H2O)-1 s-1)
//					printf("%d \t Time: %g yr\t %g\t %g\t %g\t %g\n", j, xrain[0][j]/Yr2sec, xrain[1][j], xrain[2][j], xrain[3][j], xrain[4][j]);
//					dCdtContW = dCdtContW + xrain[4][j]; // Arithmetically average over time span to get rid of PHREEQC numerical noise
//				}                                        // (geom average would be more accurate but has 50% chance of yielding a NaN = sqrt(<0))
//				dCdtContW = dCdtContW/(double)(i-3);
//
//				printf("dCdtContW=%g mol (kg H2O)-1 s-1\n", dCdtContW);
//				FCcontW = dCdtContW*WRcontW*Asurf*zContW*L*rhoContCrust;
//
//				printf("Weathering time scale: %g years\n", kintime);
//			}
//
//			for (i=0;i<nAqSpecies;i++) free (xrain[i]);
//			free (xrain);
		}

		//-------------------------------------------------------------------
		// Calculate surface C flux from seafloor weathering and subduction
		//-------------------------------------------------------------------

		// 1. Call PHREEQC to get deltaCreac, net mol C leached/precipitated per kg of rock
		deltaCreac = 0.0;

//		AqueousChem(path, "io/OceanStart.txt", itime, Tsurf, &Psurf, &Vatm, &nAir, &pH, &pe, &Mocean, &xgas, &xaq, NULL, 1, 0.0, 1, nvarEq);

		tcirc = 1.0;      // TODO compute tcirc based on Nu(Ra(basal heat flux))
		FCseafsubd = -(1.0-L) * 4.0/3.0*PI_greek*(pow(r_p,3)-pow(r_p-zCrack,3))/tcirc*deltaCreac*rho[NR] / (4.0*PI_greek*r_p*r_p);

		//-------------------------------------------------------------------
		// Compute net geo C fluxes
		//-------------------------------------------------------------------

		// Usually, continental weathering fluxes /2 because 1/2 of weathered crust is carbonate that re-precipitates in the ocean, releasing CO2 (Gaillardet+ 1999; Foley+ 2015).
		// The only significant reactions are weathering of Ca and Mg *silicates* on continents, as is the case here.
		// The control of atm CO2 by weathering of Na & K silicates on land is less obvious: alkalis are involved in
		// ‘reverse weathering’ reactions in seawater, forming Na and K silicates and converting HCO3- back to atm CO2.
//		RCplate = RCplate + dtime*(FCcontw/2.0 + FCseafw - FCsubd);
//		RCmantle = RCmantle + dtime*((1.0-farc)*FCsubd - FCoutgas); // Assuming instantaneous mixing into the mantle once subducted (in practice could take 1 Gyr)
		netFC = FCoutgas + FCcontW + (1.0-farc)*FCseafsubd; // FCcontw < 0, FCseafw < 0, FCsubd < 0
		RCmantle = RCmantle - dtime*netFC; // Assuming plate = mantle (unlike Foley et al. 2015) and instantaneous mixing into the mantle once subducted (in practice could take 1 Gyr)
		RCatmoc = RCatmoc + dtime*netFC; // Sum of atmospheric and ocean reservoirs, still needs partitioning

		RCatm = (xgas[0]+xgas[1])*nAir;
		RCocean = (xaq[0]+xaq[1])*Mocean;

		// Write outputs
		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"Output.txt");
		fout = fopen(title,"a");
		if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
		else {
			fprintf(fout, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n", realtime/Gyr2sec, Tmantle, Ra, nu*rho[(int)((iBDT+iCMB)/2)], k*(Tmantle-Tsurf)/zLith*1000.0, zLith/km2m, bndcoef, FCoutgas, vConv*1.0e-6*Myr2sec);

//			fprintf(fout, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n",
//					(double)itime*dtime/Myr2sec, Tmantle, RCmantle, RCatmoc, RCocean, FCoutgas, FCcontW, FCseafsubd, netFC);
//
//			fprintf(fout, "\n");
//			fprintf(fout, "Time: %g Gyr\n", realtime/Gyr2sec);
//			fprintf(fout, "---------------------------------------------------------------------------\n");
//			fprintf(fout, "Quantity                                | Earth benchmark   | Model result \n");
//			fprintf(fout, "---------------------------------------------------------------------------\n");
//			fprintf(fout, "New crust generation rate (m Myr-1)     | 40                | %.3g \n", zNewcrust*Myr2sec);
//			fprintf(fout, "New crust density (kg m-3)              | 2800              | %.4g \n", rhomelt);
//			fprintf(fout, "C and H2O outgassing rate (mol s-1)     | 115000            | %.5g \n", FCoutgas);
//			fprintf(fout, "Magma C mass fraction                   | 0.1-0.65%%         | %.3g%% \n", magmaCmassfrac*100.0);
//			fprintf(fout, "Convective boundary layer thickness coef| 0.15 to 0.5       | %.2g \n", bndcoef);
//			fprintf(fout, "Convective velocity (cm yr-1)           | ~1                | %.2g \n", vConv*100.0*1.0e-6*Myr2sec);
//			fprintf(fout, "Mantle convection timescale (Myr)       | ~50-200 Myr       | %.4g \n", tConv/Myr2sec);
//			fprintf(fout, "Rayleigh number                         | ~1e6-1e7          | %.2g \n", Ra);
//			fprintf(fout, "Mantle adiabatic gradient (K km-1)      | 0.5               | %.3g \n", alpha*g[(int)((iBDT+iCMB)/2)]*Tmantle/Cp * km2m);
//			fprintf(fout, "Half-mantle depth (km)                  | 1450              | %.4g \n", (r_p-r_c-zLith)/2.0/km2m);
//			fprintf(fout, "Heat flux (mW m-2)                      | 86 (radiodecay 40)| %.2g \n", k*(Tmantle-Tsurf)/zLith*1000.0);
//			fprintf(fout, "---------------------------------------------------------------------------\n");
//			fprintf(fout, "Mid-mantle temperature: %g K \t Brittle-ductile transition temperature: %.3g K \t Temperature at base of lithosphere: %.4g K\n", Tmantle, T_BDT, T[iLith]);
//			fprintf(fout, "Depth of brittle-ductile transition: %.3g km \t Thickness of lithosphere: %.3g km \n", (r_p-r[iBDT])/km2m, zLith/km2m);
//			fprintf(fout, "Convective driving stress: %.4g MPa \t Lithospheric yield stress: %.4g MPa \t", driveStress/MPa2Pa, yieldStress/MPa2Pa);
//			if (staglid) fprintf(fout, "Stagnant lid\n");
//			else fprintf (fout, "Plate tectonics\n");
//			fprintf(fout, "\n\n");
		}
		fclose (fout);

		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"CompoOceanAtm.txt");
		fout = fopen(title,"a");
		if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
		else fprintf(fout, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n",
				(double)itime*dtime/Myr2sec, xgas[0], xgas[1], xgas[2], xgas[3], xgas[4], Psurf, Tsurf, DeltaTghe, xaq[0], xaq[1], xaq[3],
				pH, 4.0*(pe+pH)-logKO2H2O, rainpH, (xaq[0]+xaq[1])*Mocean + (xgas[0]+xgas[1])*nAir, xaq[3]*Mocean + xgas[3]*2.0*nAir, Mocean, nAir);
		fclose (fout);
	} // End time loop

	printf("\nExiting ExoCcycleGeo...\n");

	for (i=0;i<NR;i++) free (sys_tbl[i]);
	free (sys_tbl);

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

	return EXIT_SUCCESS;
}
