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
 2- Refine convection model with comparison of yield and driving stress to
 enable plate tectonics mode, enables subduction sink - then scale
 melting from Sasaki & Tajika 1995.
 3- Seafloor weathering on the timescale of emplacement of new crust. km3
 created per year * fraction getting within diffusional reach of seawater, if
 slow enough, let that determine the rate of consumption by assuming reaction
 has time to go to equilibrium. Confirm with PHREEQC simulations of rates.
 4- Redox of outgassing. MELTS?
 ============================================================================
 */

#include "ExoCcycleGeo.h"
#include "CHNOSZcmds.h"
#include "Compression.h"

//-------------------------------------------------------------------
// SUBROUTINE DECLARATIONS
//-------------------------------------------------------------------

int AqueousChem (char path[1024], char filename[64], int itime, double T, double *P, double *V, double *nAir, double *pH, double *pe,
		double *mass_w, double **xgas, double **xaq, double ***xrain, int forcedPP, double kintime, int kinsteps, int nvar);
int ExtractWrite(int instance, double*** data, int line, int nvar);
const char* ConCat (const char *str1, const char *str2);
int WritePHREEQCInput(const char *TemplateFile, int itime, double temp, double pressure, double gasvol, double pH, double pe, double mass_w,
		double *xgas, double *xaq, int forcedPP, double kintime, int kinsteps, char **tempinput);
int cleanup (char path[1024]);
double molmass_atm (double *xgas);
double brittleDuctile (double T, double rhoCrust, double zCrust, double gsurf, double Tmantle, double Tsurf, double flowLawDiff[5], double flowLawDisl[5], double grainSize,
		double dtime);
double brittleDuctile_prime (double T, double rhoCrust, double zCrust, double gsurf, double Tmantle, double Tsurf, double flowLawDiff[5], double flowLawDisl[5],
		double grainSize, double dtime);
double viscosity (double T, double P, double flowLaw[5], double grainSize, double dtime);
double combVisc (double T, double P, double flowLawDiff[5], double flowLawDisl[5], double grainSize, double dtime);
double dviscdT (double T, double dPdT, double flowLaw[5], double grainSize, double dtime, double Tsurf);

//-------------------------------------------------------------------
// MAIN PROGRAM
//-------------------------------------------------------------------

int main(int argc, char *argv[]) {

	int i = 0;
	int j = 0;
	int NR = 100;                   // Number of radial grid zones used in determining planetary structure at setup
	int cmbindex = 0;				// Index of core-mantle boundary
	int itime = 0;                  // Time step counter
	int ntime = 0;                  // Total number of steps

	int iter = 0;                   // Iteration counter
	int niter = 10;                 // Number of iterations for ocean-atmosphere equilibrium

	double realtime = 0;            // Real time since birth of planetary system (s)

	// Planet parameters
	double r_p = 0.0;               // Planet radius (m)
	double r_c = 0.0;				// Planet core radius (m)
	double gsurf = 0.0;             // Surface gravity (m s-2)
	double Asurf = 0.0;             // Planet surface area (m2)
	double Tsurf = 0.0;             // Surface temperature (K)

	double *r = (double*) malloc((NR+1)*sizeof(double)); // Radius
	if (r == NULL) printf("Compression: Not enough memory to create r[NR]\n");

	double *rho = (double*) malloc((NR+1)*sizeof(double)); // Density
	if (rho == NULL) printf("Compression: Not enough memory to create rho[NR]\n");

	double *P = (double*) malloc((NR+1)*sizeof(double)); // Pressure
	if (P == NULL) printf("Compression: Not enough memory to create P[NR]\n");

	// Reservoirs
//	double RCplate = 0.0;           // Plate/crust C reservoir (mol)
	double RCmantle = 0.0;          // Mantle C reservoir (mol)
	double RCatm = 0.0;             // Atmospheric C reservoir (mol)
	double RCocean = 0.0;           // Ocean C reservoir (mol)
	double RCatmoc = RCatm + RCocean; // Combined atmospheric and ocean C reservoir (mol)

	// Fluxes
	double FCoutgas = 0.0;          // C flux from outgassing (subaerial+submarine) (mol s-1)
	double FCcontW = 0.0;           // C flux from continental weathering (mol s-1)
	double FCseafW = 0.0;           // C flux from seafloor weathering (mol s-1)
	double FCsubd = 0.0;            // C flux from subduction (mol s-1)
	double netFC = 0.0;             // Net surface C flux from all geological processes (mol s-1)
	double farc = 0.0;              // Fraction of subducted C that makes it back out through arc volcanism (as opposed to into the mantle)

	// Geochem parameters
	double pH = 0.0;                // pH of the surface ocean (default 8.22)
	double rainpH = 0.0;            // pH of rainwater
	double pe = 0.0;                // pe (-log activity e-) corresponding to logfO2
	double logfO2 = 0.0;            // log O2 fugacity
	double logKO2H2O = 0.0;         // log K for reaction 4 H+ + 4 e- + O2 = 2 H2O, from CHNOSZ: subcrt(c("H+","e-","O2","H2O"),c(-4,-4,-1,2),c("aq","aq","g","liq"),T=25,P=1)
	double Tfreeze = 273.15;        // Temperature at which water freezes at surface pressure (K)

	// Atmosphere parameters
	double nAir = 0.0;              // Number of mol in atmosphere (mol)
	double S = 0;                   // Stellar flux at planet (W m-2)
	double Teff = 0.0;              // Effective temperature (K)
	double DeltaTghe = 0.0;         // Temperature increase over effective temperature due to greenhouse effect (K)
	double psi = 0.0;               // log10(pCO2) (bar)
	double Vatm1 = 0.0;             // Initial atmospheric volume (m3), determined at second time step to avoid skew of starting values
	double Tsurf1 = 0.0;            // Surface temperature used to determine Vatm1 (K)
	double Vatm = 0.0;              // Atmospheric volume (m3)

	// Kinetic parameters
	int kinsteps = 0;               // Number of time steps of PHREEQC kinetic simulation
	double WRcontW = 0.0;           // Water:rock mass ratio for continental weathering reactions (kg/kg)
	double zContW = 0.0;            // Depth of continental weathering (m)
	double rhoContCrust = 0.0;      // Density of continental crust (kg m-3)
	double molmassContCrust = 0.0;  // Molar mass of continental crust (kg mol-1)
	double dCdtContW = 0.0;         // Carbon drawdown rate due to continental weathering (mol (kg H2O)-1 s-1)

	// Quantities to be computed by thermal/geodynamic model
	int staglid = 1;                // 1 if stagnant-lid, 0 if mobile-lid
	int n_iter = 0, n_iter_max = 100; // Iteration counter and max for the Bisection/Newton-Raphson loop to determine T_BDT
	double d = 0.0;                 // Depth to core-mantle boundary (m)
	double Tmantle = 0.0;           // Mantle temperature (K)
	double H = 0.0;                 // Specific radiogenic heating rate (J s-1 kg-1)
	double kappa = k/(rhoMantle*Cp); // Mantle thermal diffusivity (m2 s-1)
	double nu = 0.0;                // Mantle kinematic viscosity (m2 s-1)
	double zCrust = 10.0*km2m;      // Crustal thickness (m)
	double Ra = 0.0;                // Rayleigh number for mantle convection (no dim)
	double Nu = 0.0;                // Nusselt number (no dim)
	double Tref = 0.0;              // Temperature at outer boundary of convective zone (surface or base of stagnant lid)
	double driveStress = 0.0;       // Driving stress at the base of lithosphere, used to switch between stagnant and mobile-lid modes
	double yieldStress = 0.0;       // Lithospheric yield stress, taken to be the strength at the brittle-ductile transition (brittle strength = ductile strength)
	double T_BDT = 0.0;             // Temperature at brittle-ductile transition (K)
	double T_INF = 0.0, T_SUP = 0.0, T_TEMP = 0.0; // Intermediate temperatures used in determination of T_BDT by Bisection/Newton-Raphson loop
	double dT = 0.0, dTold = 0.0;   // Intermediate temperature changes used in determination of T_BDT by Bisection/Newton-Raphson loop
	double f_inf = 0.0, f_sup = 0.0; // Differences between brittle and ductile strengths, used in determination of T_BDT by Bisection/Newton-Raphson loop (P)
	double f_x = 0.0, f_prime_x = 0.0; // Derivatives of the above (Pa K-1)
	double NewtRaphThresh = 1.0e5;  // Threshold for the Bisection/Newton-Raphson loop, here in Pa
	double P_BDT = 0.0;             // Pressure at brittle-ductile transition (Pa)

	// Viscosity law constants (can't be #define'd because they are used as exponents)
	double flowLawDryDiff[5]; // Dry diffusion creep flow law (Korenaga & Karato 2008)
	flowLawDryDiff[0] = 261.0e3;    // Activation energy (J)
	flowLawDryDiff[1] = 6.0e-6;     // Activation volume (m3)
	flowLawDryDiff[2] = 5.25;       // Exponent of pre-exponential factor
	flowLawDryDiff[3] = 1.0;        // Stress exponent
	flowLawDryDiff[4] = 2.98;       // Grain size exponent

	double flowLawDryDisl[5]; // Dry dislocation creep flow law (Korenaga & Karato 2008)
	flowLawDryDisl[0] = 610.0e3;    // Activation energy (J)
	flowLawDryDisl[1] = 13.0e-6;    // Activation volume (m3)
	flowLawDryDisl[2] = 6.09;       // Exponent of pre-exponential factor
	flowLawDryDisl[3] = 4.94;       // Stress exponent
	flowLawDryDisl[4] = 0.0;        // Grain size exponent

	const double grainSize = 1.0e3;	// Grain size for computation of viscosity and creep laws (µm)

	// Quantities to be computed by melting model
	double Rmelt = 0.0;             // Rate of melt generation (m-2 s-1)
	double vConv = 0.0;             // Convective velocity (m s-1)
	double P0 = 0.0;                // Pressure at which geotherm crosses mantle solidus (Pa)
	double Pf = 0.0;                // Pressure at base of crust (Pa)

	// Quantities to be computed by seafloor weathering model
	double zCrack = 0.0;            // Depth of fracturing below seafloor (m)
	double tcirc = 0.0;             // Time scale of hydrothermal circulation (s)
	double deltaCreac = 0.0;        // Net C leached/precipitated per kg of rock (mol kg-1)

	double *xgas = (double*) malloc(nAtmSpecies*sizeof(double));
	if (xgas == NULL) printf("ExoCcycleGeo: Not enough memory to create xgas[nAtmSpecies]\n"); // Mixing ratios by volume (or by mol since all gases are pretty much ideal and have the same molar volume) of atmospheric gases
    for (i=0;i<nAtmSpecies;i++) xgas[i] = 0.0;

	double *xaq = (double*) malloc(nAqSpecies*sizeof(double));
	if (xaq == NULL) printf("ExoCcycleGeo: Not enough memory to create xaq[nAqSpecies]\n"); // Molalities of aqueous species (mol (kg H2O)-1)
    for (i=0;i<nAqSpecies;i++) xaq[i] = 0.0;

	FILE *fout;
	char *title = (char*)malloc(1024*sizeof(char));

	//-------------------------------------------------------------------
	// Inputs
	//-------------------------------------------------------------------

	double dtime = 1.0e-3*Myr2sec;     // Time step (s)

	ntime = (int) (5000.0*Myr2sec/dtime); // Number of time steps of simulation

	// Planet parameters
	double m_p = 1.0*mEarth;    // Planet mass (kg)
	double m_c = 0.3*m_p;       // Core mass (kg), default ≈0.3*m_p for Earth
	double L = 0.29;            // Fraction of planet surface covered by land
	double Mocean = 1.4e21;     // Mass of ocean (kg, default: Earth=1.4e21)
	double rhoMagma = 3500.0;   // Magma density (kg m-3)
	double rhoCrust = 3500.0;   // Crustal density (kg m-3)
    int redox = 2; // 1: current Earth surface, 2: hematite-magnetite, 3: fayalite-magnetite-quartz, 4: iron-wustite, code won't run with other values.
    Tmantle = 1850.0;

	// Atmospheric inputs
	double Tsurf0 = 300.0;        // Initial surface temperature (K)
	double S0 = 1368.0;         // Stellar flux for G2V star at 1 AU and 4.55 Gyr
	double albedo = 0.3;
	double Psurf = 1.0;         // Surface pressure (bar)
	double runoff = 0.7e-3/86400.0; // Atmospheric runoff (m s-1), default runoff_0 = 0.665e-3 m day-1
	xgas[0] = 0.5;       // CO2
    xgas[1] = 0.0;       // CH4
    xgas[2] = 0.0;       // O2
    xgas[3] = 0.5;       // N2
    xgas[4] = 0.0;       // H2O

    // Ocean inputs
	pH = 8.22;
    xaq[0] = 27.0/12.0/1000.0;  // 27 ppm C in today's oceans (scaled from 141 ppm Alk*M(C)/M(HCO3), M being molecular mass)
    xaq[1] = 0.0/12.0/1000.0;
    xaq[3] = 0.0/14.0/1000.0;   // ppm N in ocean

    // Continental weathering inputs
	double deplCcrit = 0.18;    // Chemical depletion fraction, i.e. how much carbon has been consumed from rain water by the time fresh rock is emplaced
	double kintime = 1.0e-6*Yr2sec; // Total time of PHREEQC kinetic simulation (s), default: 1.0e-6 year
	kinsteps = 100;             // 100 steps by default, could be more
	zContW = 10.0;
	rhoContCrust = 2740.0;
	molmassContCrust = 0.149;

	//-------------------------------------------------------------------
	// Startup
	//-------------------------------------------------------------------

	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("ExoCcycleGeo v19.6\n");
	printf("This code is in development and cannot be used for science yet.\n");
	if (cmdline == 1) printf("Command line mode\n");
	printf("-------------------------------------------------------------------\n");

	printf("\n");
	printf("Inputs:\n");
	printf("Planet mass \t \t \t \t %g M_Earth\n", m_p/mEarth);
	printf("Land coverage of planet surface \t %g%%\n", L*100.0);
	printf("Surface temperature \t \t \t %g K\n", Tsurf);
	printf("Surface pressure \t \t \t %g K\n", Psurf);
	printf("Atmospheric CO2 molar mixing ratio \t %g\n", xgas[0]);
	printf("Atmospheric CH4 molar mixing ratio \t %g\n", xgas[1]);
	printf("Atmospheric O2 molar mixing ratio \t %g\n", xgas[2]);
	printf("Atmospheric N2 molar mixing ratio \t %g\n", xgas[3]);
	printf("Ocean pH \t \t \t \t %g \n", pH);
	printf("Surface runoff \t \t \t \t %g mm d-1\n", runoff);
	printf("-------------------------------------------------------------------\n");

	printf("\n");
	printf("Computing geo C fluxes through time...\n");

	// Initialize the R environment. We do it here, in the main loop, because this can be done only once.
	setenv("R_HOME","/Library/Frameworks/R.framework/Resources",1);     // Specify R home directory
	Rf_initEmbeddedR(argc, argv);                                       // Launch R
	CHNOSZ_init(1);                                                     // Launch CHNOSZ

	// Get current directory. Works for Mac only! To switch between platforms, see:
	// http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe
	char path[1024];
	unsigned int size = sizeof(path);
	path[0] = '\0';
	if (_NSGetExecutablePath(path, &size) == 0) printf("\n");
	else printf("ExoCcycleGeo: Couldn't retrieve executable directory. Buffer too small; need size %u\n", size);

	compression(NR, m_p, m_c, Tsurf, 101, 101, 202, &r, &P, &rho, &cmbindex, path);

	r_p = r[NR];
	r_c = r[cmbindex];
	printf("Planet radius = %g km, Core radius = %g km\n", r_p/km2m, r_c/km2m);
	d = r_p-r_c;
	gsurf = G*m_p/r_p/r_p;
	Asurf = 4.0*PI_greek*r_p*r_p;
	Tsurf = Tsurf0;
	nAir = Psurf*bar2Pa*Asurf/gsurf/molmass_atm(xgas);

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

	AqueousChem(path, "io/OceanStart", itime, Tsurf, &Psurf, &Vatm, &nAir, &pH, &pe, &Mocean, &xgas, &xaq, NULL, 1, 0.0, 1, nvarEq);

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
		fprintf(fout, "'Time (Myr)' \t 'Tmantle (K)' \t RCmantle \t RCatmoc \t RCocean \t FCoutgas \t FCcontw \t FCseafw \t FCsubd \t 'Net C flux'\n");
		fprintf(fout, "Init \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n", Tmantle, RCmantle, RCatmoc, RCocean, FCoutgas, FCcontW, FCseafW, FCsubd, netFC);
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

//		realtime = (double)itime*dtime;                // Start at birth of planetary system
		realtime = (double)itime*dtime + 4.55*Gyr2sec; // Start at present day

//		netFC = 0.0; // !! Debug

		// Print relative change in Ntot to ensure N is conserved. Also print net C flux relative to outgassing flux; if small the C cycle should be balanced
		printf("Iteration %d/%d, Time: %g Myr, Tsurf: %g K, pCO2: %g bar, Delta_Ntot/Ntot: %.2f pct, NetFC/outgas: %.2f pct\n",
				 itime, ntime, (double)itime*dtime/Myr2sec, Tsurf, xgas[0]*Psurf,
				 (1.0-(xgas[3]*2.0*nAir + xaq[3]*Mocean)/TotalN0)*100.0, netFC/FCoutgas*100.0);

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

		Tsurf = Teff + DeltaTghe;
//		Tsurf = Tsurf0;

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

		if (Psurf < 0.01) {
			printf("ExoCcycleGeo: Pressure = %g bar too close to the triple point of H2O, oceans not stable at the surface. Exiting.\n", Psurf);
			exit(0);
		}

		//-------------------------------------------------------------------
		// Equilibrate ocean and atmosphere
		//-------------------------------------------------------------------

		if (Tsurf > Tfreeze) {

			AqueousChem(path, "io/OceanDiss", itime, Tsurf, &Psurf, &Vatm, &nAir, &pH, &pe, &Mocean, &xgas, &xaq, NULL, 0, 0.0, 1, nvarEq);

			if (Psurf < 0.01) {
				printf("ExoCcycleGeo: Pressure = %g bar too close to the triple point of H2O, oceans not stable at the surface. Exiting.\n", Psurf);
				exit(0);
			}

			cleanup(path); // Remove PHREEQC selected output file
		}

		//-------------------------------------------------------------------
		// Calculate surface C flux from outgassing
		//-------------------------------------------------------------------

		// Kinematic viscosity
		// Upper mantle:
		nu = combVisc(Tmantle, rhoMantle*gsurf*d/2.0, flowLawDryDiff, flowLawDryDisl, grainSize, dtime)/rhoMantle;
		printf("nu KK08 = %g\n", nu);
		// Lower mantle:
//		nu = 1.0e16*exp((2.0e5 + rhoMantle*gsurf*d/2.0*1.1e-6)/(R_G*Tmantle))/rhoMantle; // Cízková et al. (2012)
//		printf("nu C12 = %g\n", nu);

		// Compute instantaneous heating rate (Kite et al. 2009 Table 1): H = X_4.5 * W * exp(ln(1/2) / t1/2 * (t-4.5))
		H =  36.9e-9  * 2.92e-5 * exp(log(0.5)/( 1.26 *Gyr2sec) * (realtime - 4.5*Gyr2sec))  //  40-K
		  + 124.0e-9  * 2.64e-5 * exp(log(0.5)/(14.0  *Gyr2sec) * (realtime - 4.5*Gyr2sec))  // 232-Th
		  +   0.22e-9 * 56.9e-5 * exp(log(0.5)/( 0.704*Gyr2sec) * (realtime - 4.5*Gyr2sec))  // 235-U
		  +  30.8e-9  * 9.46e-5 * exp(log(0.5)/( 4.47 *Gyr2sec) * (realtime - 4.5*Gyr2sec)); // 238-U

		// Compute effective thermal conductivity, assuming whole-mantle convection (Kite et al. 2009)
		if (!staglid) Tref = Tsurf; // TODO circular logic to then have value of !staglid depend on Tref below.
		else          Tref = Tmantle - 2.23*Tmantle*Tmantle/A0; // Tmantle - Tc in K09

		Ra = alpha*gsurf*(Tmantle-Tref)*pow(d,3)/(kappa*nu);
		Nu = pow(Ra/Ra_c, beta);

		// Convective velocity (K09 equation 24), replaced Tc (temp at base of stagnant lid) with Tref for compatibility with !staglid regime
		vConv = 2.0*(Nu-1.0) * (k/(rhoMagma*Cp*(r_p-r_c))) * (Tmantle-Tsurf) / (Tmantle-Tref); // Check against eq. (6.379) of Turcotte & Schubert (2002), p. 514

		// Determine if convective driving stress exceeds lithosphere yield stress, in which case switch to mobile-lid (plate tectonics) regime (O'Neill & Lenardic 2007)
		/* Convective drive stress:
		 * Viscosity = stress / strain rate (assumed Newtonian - see Deschamps & Sotin 2000), i.e. stress = viscosity * vConv/δ, w/ δ: boundary layer thickness
		 * δ = thickness of thermal gradient = 2.32*(kappa*δ/vConv)^0.5 (Turcotte & Schubert 2002, eq. 6.327 with δ/vConv = conduction timescale)
		 * So δ^0.5 = 2.32*(kappa/vConv)^0.5 and δ ≈ 5.38*kappa/vConv */
		driveStress = nu*rhoMantle*vConv*vConv/(5.38*kappa);

		/* Lithospheric yield stress = strength at brittle-ductile transition (BDT).
		 * Geotherm near surface is T prop to depth, with T=Tsurf at depth = 0 and T=Tmantle at depth = zCrust = P(BDT)/(rhoCrust*gsurf). [1]
		 * So T = (Tmantle-Tsurf)*depth/zCrust + Tsurf. [2]
		 * P ≈ rhoCrust*gsurf*depth = rhoCrust*gsurf*zCrust*(T-Tsurf)/(Tmantle-Tsurf) [3]
		 * Initiate zCrust. Solve for temperature of BDT: brittleStrength-ductileStrength = f(P,T) = f(T) = 0. Get updated zCrust from [2] (= depth). */

		// Initialize bounds for T_BDT
		T_INF = Tsurf;
		T_SUP = Tmantle;

		// Ensure that f(T_INF)<0 and f(T_SUP)>0
		f_inf = brittleDuctile(T_INF, rhoCrust, zCrust, gsurf, Tmantle, Tsurf, flowLawDryDiff, flowLawDryDisl, grainSize, dtime);
		f_sup = brittleDuctile(T_SUP, rhoCrust, zCrust, gsurf, Tmantle, Tsurf, flowLawDryDiff, flowLawDryDisl, grainSize, dtime);

		if (f_inf*f_sup > 0.0) {
			if (f_sup < 0.0) printf("ExoCcycleGeo: Brittle strength < ductile strength, i.e. brittle regime even at Tmantle=%g K.\n"
					                 "Brittle-ductile transition could not be determined.\n", Tmantle); // Brittle regime in mantle
			else zCrust = 0.0; // Ductile regime at surface
		}
		else {
			if (f_inf > 0.0) { // Swap INF and SUP if f_inf > 0 and f_sup < 0
				T_TEMP = T_INF;
				T_INF = T_SUP;
				T_SUP = T_TEMP;
			}
			T_BDT = 0.5*(T_INF+T_SUP); // Initialize the guess for the root T_BDT,
			dTold = fabs(T_INF-T_SUP); // Initialize the "stepsize before last"
			dT = dTold;                // Initialize the last stepsize

			f_x = brittleDuctile(T_BDT, rhoCrust, zCrust, gsurf, Tmantle, Tsurf, flowLawDryDiff, flowLawDryDisl, grainSize, dtime);
			f_prime_x = brittleDuctile_prime(T_BDT, rhoCrust, zCrust, gsurf, Tmantle, Tsurf, flowLawDryDiff, flowLawDryDisl, grainSize, dtime);

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
				f_x = brittleDuctile(T_BDT, rhoCrust, zCrust, gsurf, Tmantle, Tsurf, flowLawDryDiff, flowLawDryDisl, grainSize, dtime);
				f_prime_x = brittleDuctile_prime(T_BDT, rhoCrust, zCrust, gsurf, Tmantle, Tsurf, flowLawDryDiff, flowLawDryDisl, grainSize, dtime);

				if (f_x < 0.0) T_INF = T_BDT; // Maintain the bracket on the root
				else T_SUP = T_BDT;

				n_iter++;
				if (n_iter>=n_iter_max) {
					printf("ExoCcycleGeo: could not find the brittle-ductile transition after %d iterations\n",n_iter_max);
					break;
				}
			}

			zCrust = (T_BDT-Tsurf)/(Tmantle-Tsurf)*zCrust; // Crustal thickness is the depth of the brittle-ductile transition
		}

		P_BDT = rhoCrust*gsurf*zCrust; // Pressure at brittle-ductile transition (Pa)

		// Equate yield stress to brittle strength at brittle-ductie transition (equivalently, ductile strength but that expression is more complicated)
		if (P_BDT < 200.0e6) yieldStress = 0.85*P_BDT;
		else yieldStress = 0.6*P_BDT + 50.0e6;

		// Compare driving and yield stresses to assess tectonic regime
		if (yieldStress < driveStress) staglid = 0; // Mobile-lid regime, i.e. plate tectonics
		else staglid = 1;                           // Stagnant-lid regime

	// Debug: Check brittle and ductile yield stresses equal
		printf("T_BDT=%g K, P_BDT=%g MPa, driveStress=%g MPa, yieldStress=%g MPa, zCrust=%g km\n", T_BDT, P_BDT/MPa2Pa, driveStress/MPa2Pa, yieldStress/MPa2Pa, zCrust/km2m);

		if (yieldStress != combVisc(T_BDT, P_BDT, flowLawDryDiff, flowLawDryDisl, grainSize, dtime)/dtime)
			printf("Brittle yield stress=%g MPa different from ductile yield stress=%g MPa\n", yieldStress/MPa2Pa,
					combVisc(T_BDT, P_BDT, flowLawDryDiff, flowLawDryDisl, grainSize, dtime)/dtime/MPa2Pa);

		if (staglid) printf("stagnant lid\n\n");
		else printf ("plate tectonics\n\n");
	// End debug -----

		// Melting model of Kite et al. (2009)
		if (!staglid) Pf = 0.0; // K09 equations (10-12)
		else Pf = rhoCrust*gsurf*zCrust;
		P0 = 2.0*P_BDT; // Should be higher than Pf. Arbitrary for now, rhoCrust*gsurf*2.0*zCrust results in 50% melting rate

		double MORlength = 60000*km2m; // Unconstrained parameter, default 60000 km (Earth today), likely did not vary monotonically in the past

		if (!staglid) Rmelt = vConv*MORlength*zCrust*rhoCrust; // Enhancement of rate of melting due to plate tectonics. Also a function of age. See Sasaki & Tajika 1995.
		else Rmelt = Asurf*vConv*rhoMagma * (rhoCrust*gsurf*zCrust/(P0-Pf)); // kg s-1, Kite et al. (2009) eq. 25 not normalized by planet mass. Note typo: P0>Pf

		FCoutgas = deltaCvolcEarth * Rmelt/(RmeltEarth*mEarth);         // mol C s-1. Accurate for Earth magma C concentrations, but for other mantle concentrations, should be calculated explicitly (MELTS?)

		// Update mantle temperature using effective thermal conductivity scaled with Nu
		Tmantle = Tmantle + dtime*(H/Cp - kappa*Nu*(Tmantle-Tref)/(d*d));

		//-------------------------------------------------------------------
		// Calculate surface C flux from continental weathering
		//-------------------------------------------------------------------

		if (Tsurf > Tfreeze) {
			// Analytical calculation from Edson et al. (2012) Eq. 1; Abbot et al. (2012) Eq. 2
			FCcontW = -L * 0.5*deltaCcontwEarth*Asurf * pow(xgas[0]/xCO2g0,0.3) * runoff/runoff_0 * exp((Tsurf-TsurfEarth)/17.7);

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
//				AqueousChem(path, "io/ContWeather", itime, Tsurf, &Psurf, 0.0, &pH, &pe, WRcontW, &xgas, &xaq, &xrain, 1, kintime, kinsteps, nvarKin);
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
		// Calculate surface C flux from seafloor weathering TODO include kinetics, manage reservoir size
		//-------------------------------------------------------------------

		deltaCreac = 0.0; // TODO call PHREEQC to get deltaCreac, net mol C leached/precipitated per kg of rock
		tcirc = 1.0;      // TODO compute tcirc based on Nu(Ra(basal heat flux))
		FCseafW = -(1.0-L) * 4.0/3.0*PI_greek*(pow(r_p,3)-pow(r_p-zCrack,3))/tcirc*deltaCreac*rhoCrust / (4.0*PI_greek*r_p*r_p);

		//-------------------------------------------------------------------
		// Calculate surface C flux from subduction
		//-------------------------------------------------------------------

		// Function of mantle convective vigor (planet size and age) as for volcanism
		FCsubd = 0.0;

		//-------------------------------------------------------------------
		// Compute net geo C fluxes
		//-------------------------------------------------------------------

		// Usually, continental weathering fluxes /2 because 1/2 of weathered crust is carbonate that re-precipitates in the ocean, releasing CO2 (Gaillardet+ 1999; Foley+ 2015).
		// The only significant reactions are weathering of Ca and Mg *silicates* on continents, as is the case here.
		// The control of atm CO2 by weathering of Na & K silicates on land is less obvious: alkalis are involved in
		// ‘reverse weathering’ reactions in seawater, forming Na and K silicates and converting HCO3- back to atm CO2.
//		RCplate = RCplate + dtime*(FCcontw/2.0 + FCseafw - FCsubd);
//		RCmantle = RCmantle + dtime*((1.0-farc)*FCsubd - FCoutgas); // Assuming instantaneous mixing into the mantle once subducted (in practice could take 1 Gyr)
		netFC = FCoutgas + FCcontW + FCseafW + (1.0-farc)*FCsubd; // FCcontw < 0, FCseafw < 0, FCsubd < 0
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
		else fprintf(fout, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n",
				(double)itime*dtime/Myr2sec, Tmantle, RCmantle, RCatmoc, RCocean, FCoutgas, FCcontW, FCseafW, FCsubd, netFC);
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

	free (title);
	free (xgas);
	free (xaq);
	free (r);
	free (rho);
	free (P);

	Rf_endEmbeddedR(0);                                     // Close R and CHNOSZ

	return EXIT_SUCCESS;
}

/*--------------------------------------------------------------------
 *
 * Subroutine AqueousChem
 *
 * Calculates aqueous chemistry, including:
 * - ocean-atmospheric equilibrium
 * - continental weathering
 *
 * forced PP forces the aqueous phase to be in equilibrium with gas
 * partial pressures
 *
 *--------------------------------------------------------------------*/

int AqueousChem (char path[1024], char filename[64], int itime, double T, double *P, double *V, double *nAir, double *pH, double *pe,
		double *mass_w, double **xgas, double **xaq, double ***xrain, int forcedPP, double kintime, int kinsteps, int nvar) {

	int phreeqc = 0;
	int i = 0;
	int j = 0;

	char *dbase = (char*)malloc(1024);                           // Path to thermodynamic database
	char *infile = (char*)malloc(1024);                          // Path to initial (template) input file
	char *tempinput = (char*)malloc(1024);                       // Temporary PHREEQC input file, modified from template

	double nAir0 = *nAir; // Scaling factor for gas volume and water mass, so PHREEQC doesn't have to handle large numbers

	double **simdata = (double**) malloc(nvar*sizeof(double*));
	if (simdata == NULL) printf("AqueousChem: Not enough memory to create simdata[nvar][kinsteps]\n");
	for (i=0;i<nvar;i++) {
		simdata[i] = (double*) malloc(kinsteps*sizeof(double));
		if (simdata[i] == NULL) printf("AqueousChem: Not enough memory to create simdata[nvar][kinsteps]\n");
	}

	// Initializations
	dbase[0] = '\0';
	infile[0] = '\0';
	tempinput[0] = '\0';
	for (i=0;i<nvar;i++) {
		for (j=0;j<kinsteps;j++) simdata[i][j] = 0.0;
	}

	if (cmdline == 1) strncat(dbase,path,strlen(path)-20);
	else strncat(dbase,path,strlen(path)-18);
	strcat(dbase,"PHREEQC-3.1.2/core10_idealgas.dat");

	strncat(infile,dbase,strlen(dbase)-19);
	strcat(infile, filename);

	WritePHREEQCInput(infile, itime, T-Kelvin, *P, *V/nAir0, *pH, *pe, *mass_w/nAir0, *xgas, *xaq, forcedPP, kintime, kinsteps, &tempinput);

	phreeqc = CreateIPhreeqc(); // Run PHREEQC
	if (LoadDatabase(phreeqc,dbase) != 0) OutputErrorString(phreeqc);
	SetSelectedOutputFileOn(phreeqc,1);
	if (RunFile(phreeqc,tempinput) != 0) OutputErrorString(phreeqc);

	if      (strcmp(filename, "io/OceanStart") == 0) ExtractWrite(phreeqc, &simdata, 1, nvarEq); // 1st line of selected output has initial solution composition = equilibrium when there are no equilibrium phases
	else if (strcmp(filename, "io/OceanDiss") == 0) ExtractWrite(phreeqc, &simdata, 2, nvarEq);  // 2nd line of PHREEQC selected output solution and mineral+gas composition at equilibrium
	else if (strcmp(filename, "io/ContWeather") == 0) ExtractWrite(phreeqc, &simdata, kinsteps, nvarKin);
	else printf("AqueousChem: Cannot extract PHREEQC output data because input path %s is inaccurate. Exiting.\n", filename);

	if (DestroyIPhreeqc(phreeqc) != IPQ_OK) OutputErrorString(phreeqc);

//	if (strcmp(filename, "io/OceanDiss") == 0) {
//		for (i=0;i<nvar;i++) {
//			printf("%d\t", i);
//			for (j=0;j<kinsteps;j++) printf("%g\t", simdata[i][j]);
//			printf("\n");
//		}
//		exit(0);
//	}

	// Setting initial ocean chemistry
	if (strcmp(filename, "io/OceanStart") == 0) {
		(*pH) = simdata[7][0];
		(*pe) = simdata[8][0];

		(*xaq)[0] = simdata[45][0];             // C(4), i.e. dissolved CO2 and carbonate
		(*xaq)[1] = simdata[43][0];             // C(-4), i.e. dissolved methane
		(*xaq)[2] = simdata[72][0];             // O(0), i.e. dissolved O2
		(*xaq)[3] = simdata[28][0];             // Total dissolved N
		(*xaq)[4] = simdata[69][0];             // N(0), i.e. dissolved N2
		(*xaq)[5] = simdata[68][0];             // N(-3), i.e. dissolved NH3 and NH4+
	}

	// Ocean dissolution
	if (strcmp(filename, "io/OceanDiss") == 0) {
		(*pH) = simdata[7][0];
		(*pe) = simdata[8][0];
		(*nAir) = simdata[1007][0]*nAir0;
		(*P) = simdata[1006][0]*atm2bar;
		(*V) = simdata[1008][0]/1000.0*nAir0;
		(*mass_w) = simdata[11][0]*nAir0;

		(*xaq)[0] = simdata[45][0];             // C(4), i.e. dissolved CO2 and carbonate
		(*xaq)[1] = simdata[43][0];             // C(-4), i.e. dissolved methane
		(*xaq)[2] = simdata[72][0];             // O(0), i.e. dissolved O2
		(*xaq)[3] = simdata[28][0];             // Total dissolved N
		(*xaq)[4] = simdata[69][0];             // N(0), i.e. dissolved N2
		(*xaq)[5] = simdata[68][0];             // N(-3), i.e. dissolved NH3 and NH4+

		(*xgas)[0] = simdata[1014][0];          // CO2(g)
		(*xgas)[1] = simdata[1012][0];          // CH4(g)
		(*xgas)[2] = simdata[1022][0];          // O2(g)
		(*xgas)[3] = simdata[1018][0];          // N2(g)
		(*xgas)[4] = simdata[1016][0];          // H2O(g)

		for (i=0;i<nAtmSpecies;i++) (*xgas)[i] = (*xgas)[i]/simdata[1007][0]; // Divide by total mol gas to return mixing ratio
	}

	// Continental weathering
	if (strcmp(filename, "io/ContWeather") == 0) {

		(*xrain)[1][1] = simdata[3][1]; // Rainwater pH
		(*xrain)[2][1] = simdata[9][1]; // Initial total dissolved C
		for (i=2;i<kinsteps;i++) {
			(*xrain)[0][i] = simdata[1][i]; // Time (s)
			(*xrain)[1][i] = simdata[3][i]; // pH
			(*xrain)[2][i] = simdata[9][i]; // Total dissolved C
		}
	}

	for (i=0;i<nvar;i++) free (simdata[i]);
	free (dbase);
	free (infile);
	free (tempinput);
	free (simdata);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine ExtractWrite
 *
 * Write selected output from PHREEQC
 *
 *--------------------------------------------------------------------*/

int ExtractWrite(int instance, double*** data, int line, int nvar) {
	VAR v;
	int i = 0;
	int j = 0;
	VarInit(&v);

	if (nvar == nvarEq) { // Equilibrium simulation
		GetSelectedOutputValue(instance,1,3,&v);           // temp
		(*data)[0][0] = v.dVal;

		GetSelectedOutputValue(instance,1,1,&v);           // pH
		(*data)[2][0] = v.dVal;

		GetSelectedOutputValue(instance,1,2,&v);           // pe
		(*data)[3][0] = v.dVal;

		GetSelectedOutputValue(instance,1,5,&v);           // W:R
		(*data)[6][0] = v.dVal;

		for (i=1;i<nvar-6;i++) {                           // Rest of parameters
			GetSelectedOutputValue(instance,line,i,&v);
			if (fabs(v.dVal) < 1e-50) (*data)[i+6][0] = 0.0;
			else (*data)[i+6][0] = v.dVal;
		}
	}
	else if (nvar == nvarKin) { // Kinetic simulation
		for (i=0;i<nvar;i++) {
			for (j=0;j<line;j++) {
				GetSelectedOutputValue(instance,j,i,&v);
				if (fabs(v.dVal) < 1e-50) (*data)[i][j] = 0.0;
				else (*data)[i][j] = v.dVal;
			}
		}
	}
	else printf("ExtractWrite: nvar=%d different from nvarEq=%d and nvarKin=%d, must be equal to one of these. Exiting\n", nvar, nvarEq, nvarKin);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine ConCat
 *
 * Concatenation. Takes 2 strings and returns the concatenated string.
 *
 *--------------------------------------------------------------------*/

const char* ConCat(const char *str1, const char *str2) {
	char buffer[100];
	buffer[0] = '\0';

	strcpy(buffer,str1);
	return strcat(buffer,str2);
}

/*--------------------------------------------------------------------
 *
 * Subroutine WritePHREEQCInput
 *
 * Generate input file from a template.
 *
 *--------------------------------------------------------------------*/

int WritePHREEQCInput(const char *TemplateFile, int itime, double temp, double pressure, double gasvol, double pH, double pe, double mass_w,
		double *xgas, double *xaq, int forcedPP, double kintime, int kinsteps, char **tempinput) {

	int i = 0;
	int eqphases = 0; // Switch: is the EQUILIBRIUM_PHASES block being read?
	int gasphase = 0; // Switch: is the GAS_PHASE block being read?
	int sel = 0;      // Switch: is the SELECTED_OUTPUT block being read?

	FILE *fin;
	FILE *fout;
	char itime_str[32];
	char temp_str[32];
	char pressure_str[32];
	char vol_str[32];
	char pH_str[32];
	char pe_str[32];
	char mass_w_str[32];
	char steps_str1[64]; char steps_str2[64];

	char **gas_str1 = (char**)malloc(nAtmSpecies*sizeof(char*));
	for (i=0;i<nAtmSpecies;i++) gas_str1[i] = (char*)malloc(1024);

	char **gas_str2 = (char**)malloc(nAtmSpecies*sizeof(char*));
	for (i=0;i<nAtmSpecies;i++) gas_str2[i] = (char*)malloc(1024);

	char **gas_str3 = (char**)malloc(nAtmSpecies*sizeof(char*));
	for (i=0;i<nAtmSpecies;i++) gas_str3[i] = (char*)malloc(1024);

	char **aq_str = (char**) malloc(nAqSpecies*sizeof(char*));
	for (i=0;i<nAqSpecies;i++) aq_str[i] = (char*)malloc(1024);

	itime_str[0] = '\0';
	temp_str[0] = '\0';
	pressure_str[0] = '\0';
	vol_str[0] = '\0';
	pH_str[0] = '\0';
	pe_str[0] = '\0';
	mass_w_str[0] = '\0';
	steps_str1[0] = '\0'; steps_str2[0] = '\0';

	for (i=0;i<nAtmSpecies;i++) {
		gas_str1[i][0] = '\0';
		gas_str2[i][0] = '\0';
		gas_str3[i][0] = '\0';
	}
	for (i=0;i<nAqSpecies;i++) aq_str[i][0] = '\0';

	int line_length = 300;
	char line[line_length]; // Individual line

	// Convert to PHREEQC input units
	pressure = pressure/atm2bar; // Convert from bar to atm
	gasvol = gasvol*1000.0;      // Convert from m3 to L

	// Assemble file title
	sprintf(itime_str, "%d", itime);
	sprintf(temp_str, "%g", temp);
	sprintf(pressure_str, "%g", pressure);
	sprintf(vol_str, "%g", gasvol);
	sprintf(pH_str, "%g", pH);
	sprintf(pe_str, "%g", pe);
	sprintf(mass_w_str, "%g", mass_w);

	if (!forcedPP) {
		for (i=0;i<nAtmSpecies;i++) {
			if (xgas[i] > 0.0)
				sprintf(gas_str2[i], "%g", log(xgas[i]*pressure)/log(10.0)); // Convert from mixing ratio to -log partial pressure
			else
				sprintf(gas_str2[i], "%g", 0.0);
			strcat(gas_str1[i], gas_str2[i]);
			sprintf(gas_str2[i], "%g", xgas[i]);
			strcat(gas_str1[i], "\t");
			strcat(gas_str1[i], gas_str2[i]);
		}
	}
	else {
		for (i=0;i<nAtmSpecies;i++) {
			if (xgas[i] > 0.0)
				sprintf(gas_str2[i], "%g", log(xgas[i]*pressure)/log(10.0)); // Convert from mixing ratio to -log partial pressure
			else
				sprintf(gas_str2[i], "%g", 0.0);
		}
	}
	for (i=0;i<nAtmSpecies;i++) sprintf(gas_str3[i], "%g", xgas[i]*pressure); // Partial pressure

	strcpy(*tempinput,TemplateFile);
//	strcat(*tempinput,itime_str);
	strcat(*tempinput,"Exec"); // File title complete

	sprintf(aq_str[0],"%g mol/kgw", xaq[0]);  // mol per kg H2O of oxidized C
	sprintf(aq_str[1],"%g mol/kgw", xaq[1]);  // mol per kg H2O of reduced C
	sprintf(aq_str[3],"%g mol/kgw", xaq[3]);  // mol per kg H2O of N

	sprintf(steps_str1,"%g in ", kintime);    // Duration of kinetic simulation
	sprintf(steps_str2,"%d steps", kinsteps); // Number of time steps of kinetic simulation
	strcat(steps_str1, steps_str2);           // "kintime in kinsteps steps"

	fin = fopen (TemplateFile,"r"); 	// Open input file
	if (fin == NULL) printf("WritePHREEQCInput: Missing input file. Path: %s\n", TemplateFile);
	fout = fopen (*tempinput,"w");      // Open output file
	if (fout == NULL) printf("WritePHREEQCInput: Missing output file. Path: %s\n", *tempinput);

	while (fgets(line, line_length, fin)) {
		// Block switches
		if (line[0] == 'E' && line[1] == 'Q' && line[2] == 'U' && line[3] == 'I') eqphases = 1; // EQUILIBRIUM_PHASES block
		if (line[0] == 'G' && line[1] == 'A' && line[2] == 'S' && line[3] == '_') gasphase = 1; // EQUILIBRIUM_PHASES block
		if (line[0] == 'S' && line[1] == 'E' && line[2] == 'L' && line[3] == 'E') sel = 1;      // SELECTED_OUTPUT block
		// SOLUTION
		if (line[1] == 'p' && line[2] == 'H')                                                  fprintf(fout, "%s charge\n", ConCat("\tpH \t \t",pH_str));
		else if (line[1] == 't' && line[2] == 'e' && line[3] == 'm' && line[4] == 'p')         fprintf(fout, "%s\n", ConCat("\ttemp \t \t",temp_str));
		else if (line[1] == 'p' && line[2] == 'r' && line[3] == 'e' && line[4] == 's')         fprintf(fout, "%s\n", ConCat("\tpressure \t",pressure_str));
		else if (line[1] == 'p' && line[2] == 'e')                                             fprintf(fout, "%s\n", ConCat("\tpe \t \t",pe_str));
		else if (!sel && line[1] == '-' && line[2] == 'w' && line[3] == 'a' && line[4] == 't') fprintf(fout, "%s\n", ConCat("\t-water \t \t",mass_w_str));
		else if (!eqphases && line[1] == 'C' && line[2] == '(' && line[3] == '4' && line[4] == ')') {
			if (forcedPP && xgas[0] > 0.0) {
				strcat(aq_str[0], "\tCO2(g) \t");
				strcat(aq_str[0], gas_str2[0]);
				fprintf(fout, "%s\n", ConCat("\tC(4) \t\t", aq_str[0]));
			}
			else          fprintf(fout, "%s\n", ConCat("\tC(4) \t\t",aq_str[0]));
		}
		else if (!eqphases && line[1] == 'C' && line[2] == '(' && line[3] == '-' && line[4] == '4') {
			if (forcedPP && xgas[1] > 0.0) {
				strcat(aq_str[1], "\tCH4(g) \t");
				strcat(aq_str[1], gas_str2[1]);
				fprintf(fout, "%s\n", ConCat("\tC(-4) \t\t", aq_str[1]));
			}
			else          fprintf(fout, "%s\n", ConCat("\tC(-4) \t\t",aq_str[1]));
		}
		else if (!eqphases && line[1] == 'N' && line[2] == '\t') {
			if (forcedPP && xgas[3] > 0.0) {
				strcat(aq_str[3], "\tN2(g) \t");
				strcat(aq_str[3], gas_str2[3]);
				fprintf(fout, "%s\n", ConCat("\tN \t\t", aq_str[3]));
			}
			else          fprintf(fout, "%s\n", ConCat("\tN \t\t",aq_str[3]));
		}
		// KINETICS
		else if (line[0] == '-' && line[1] == 's' && line[2] == 't' && line[3] == 'e')         fprintf(fout, "%s\n", ConCat("-steps\t",steps_str1));
		// GAS_PHASE
		else if (gasphase && line[1] == 'v' && line[2] == 'o' && line[3] == 'l') fprintf(fout, "%s\n", ConCat("\tvolume \t\t",vol_str));
		else if (gasphase && line[1] == 'C' && line[2] == 'O' && line[3] == '2' && line[4] == '(') fprintf(fout, "%s\n", ConCat("\tCO2(g) \t\t",gas_str3[0]));
		else if (gasphase && line[1] == 'C' && line[2] == 'H' && line[3] == '4' && line[4] == '(') fprintf(fout, "%s\n", ConCat("\tCH4(g) \t\t",gas_str3[1]));
		else if (gasphase && line[1] == 'O' && line[2] == '2' && line[3] == '(' && line[4] == 'g') fprintf(fout, "%s\n", ConCat("\tO2(g) \t\t",gas_str3[2]));
		else if (gasphase && line[1] == 'N' && line[2] == '2' && line[3] == '(' && line[4] == 'g') fprintf(fout, "%s\n", ConCat("\tN2(g) \t\t",gas_str3[3]));
		else if (gasphase && line[1] == 'H' && line[2] == '2' && line[3] == 'O' && line[4] == '(') fprintf(fout, "%s\n", ConCat("\tH2O(g) \t\t",gas_str3[4]));
		// EQUILIBRIUM_PHASES
		else if (eqphases && line[2] == 'C' && line[3] == 'O' && line[4] == '2' && line[5] == '(' && line[6] == 'g') fprintf(fout, "%s\n", ConCat("\tCO2(g) \t\t",gas_str1[0]));
		else if (eqphases && line[2] == 'C' && line[3] == 'H' && line[4] == '4' && line[5] == '(' && line[6] == 'g') fprintf(fout, "%s\n", ConCat("\tCH4(g) \t\t",gas_str1[1]));
		else if (eqphases && line[2] == 'O' && line[3] == '2' && line[4] == '(' && line[5] == 'g' && line[6] == ')') fprintf(fout, "%s\n", ConCat("\tO2(g) \t\t",gas_str1[2]));
		else if (eqphases && line[2] == 'N' && line[3] == '2' && line[4] == '(' && line[5] == 'g' && line[6] == ')') fprintf(fout, "%s\n", ConCat("\tN2(g) \t\t",gas_str1[3]));
		else if (eqphases && line[2] == 'H' && line[3] == '2' && line[4] == 'O' && line[5] == '(' && line[6] == 'g') fprintf(fout, "%s\n", ConCat("\tH2O(g) \t\t",gas_str1[4]));

		else fputs(line,fout);
	}
	if (ferror(fin)) {
		printf("WritePHREEQCInput: Error reading template input file %s\n",TemplateFile);
		return 1;
	}

	fclose(fin);
	fclose(fout);

	for (i=0;i<nAtmSpecies;i++) {
		free(gas_str1[i]);
		free(gas_str2[i]);
		free(gas_str3[i]);
	}
	for (i=0;i<nAqSpecies;i++) free(aq_str[i]);
	free(gas_str1);
	free(gas_str2);
	free(gas_str3);
	free(aq_str);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine cleanup
 *
 * Remove PHREEQC selected output files so as not to clutter the
 * ExoCcycleGeo folder.
 *
 *--------------------------------------------------------------------*/

int cleanup (char path[1024]) {

	char *cmd = (char*)malloc(1024*sizeof(char)); // Don't forget to free!

	cmd[0] = '\0';
	strcat(cmd,"rm ");
	if (cmdline == 1) strncat(cmd,path,strlen(path)-20);
	else strncat(cmd,path,strlen(path)-18);
	strcat(cmd,"sel*");
	system(cmd);

	free(cmd);

	return(0);
}

/*--------------------------------------------------------------------
 *
 * Subroutine molmass_atm
 *
 * Compute the average molecular mass of the atmosphere in kg mol-1.
 *
 *--------------------------------------------------------------------*/

double molmass_atm (double *xgas) {

	double molmass = 0.0;

	// Masses of core10.dat
	double M_H = 1.0079;
	double M_C = 12.011;
	double M_N = 14.0067;
	double M_O = 15.994;

	molmass = (xgas[0]*(M_C + 2.0*M_O) // CO2
			 + xgas[1]*(M_C + 4.0*M_H) // CH4
			 + xgas[2]*2.0*M_O         // O2
			 + xgas[3]*2.0*M_N         // N2
		     + xgas[4]*(2.0*M_H+M_O))  // H2O
					 /(xgas[0]+xgas[1]+xgas[2]+xgas[3]+xgas[4])*0.001;

	return molmass;
}

/*--------------------------------------------------------------------
 *
 * Subroutine brittleductile
 *
 * Calculates the difference between the brittle strength and ductile
 * strength of the interior, both in Pa.
 * Brittle-ductile and brittle-plastic transitions are
 * mixed up, although they shouldn't (Kohlstedt et al. 1995).
 *
 * The brittle strength is given by a friction/low-P Byerlee type law:
 * strength = mu*P, assuming negligible water pressure since in practice
 * the brittle-ductile transition occurs in dehydrated rock (T>700 K)
 * even over long time scales.
 *
 * The ductile strength sigma is given by a flow law:
 * d epsilon/dt = A*sigma^n*d^-p*exp[(-Ea+P*V)/RT].
 * This law requires choosing timescale for mantle convection =
 * time step (3e-11 s-1 for dtime = 1e3 years, relatively insensitive
 * to that time step because viscosity varies by orders of mag with
 * temperature).
 *
 *--------------------------------------------------------------------*/

double brittleDuctile (double T, double rhoCrust, double zCrust, double gsurf, double Tmantle, double Tsurf, double flowLawDiff[5], double flowLawDisl[5], double grainSize,
		double dtime) {

	double f = 0.0;

	double P = rhoCrust*gsurf*zCrust*(T-Tsurf)/(Tmantle-Tsurf); // Pressure (Pa)
	double brittleStrength = 0.0;                               // Brittle strength (Pa)
	double ductileStrength = 0.0;                               // Ductile strength (Pa)

	if (P < 200.0e6) brittleStrength = 0.85*P;
	else brittleStrength = 0.6*P + 50.0e6;

	ductileStrength = combVisc(T, P, flowLawDiff, flowLawDisl, grainSize, dtime)/dtime;

	f = brittleStrength - ductileStrength;

	return f;
}

/*--------------------------------------------------------------------
 *
 * Subroutine brittleductile_prime
 *
 * Calculates the temperature derivative of brittleDuctile.
 *
 *--------------------------------------------------------------------*/

double brittleDuctile_prime (double T, double rhoCrust, double zCrust, double gsurf, double Tmantle, double Tsurf, double flowLawDiff[5], double flowLawDisl[5],
		double grainSize, double dtime) {

	double f_prime = 0.0;

	double P = rhoCrust*gsurf*zCrust*(T-Tsurf)/(Tmantle-Tsurf); // Pressure (Pa)
	double dPdT = rhoCrust*gsurf*zCrust/(Tmantle-Tsurf);        // Geotherm (Pa K-1), independent of T
	double brittle_prime = 0.0;                                 // Brittle strength (Pa)
	double ductile_prime = 0.0;                                 // Ductile strength (Pa)
	double ductileDiff = 0.0;                                   // Dry silicate diffusion (Pa)
	double ductileDisl = 0.0;                                   // Dry silicate dislocation (Pa)
	double ductileDiff_prime = 0.0;                             // Dry silicate diffusion (Pa)
	double ductileDisl_prime = 0.0;                             // Dry silicate dislocation (Pa)

	// Temperature derivative of brittle strength
	if (P < 200.0e6) brittle_prime = 0.85*dPdT;
	else brittle_prime = 0.6*dPdT;

	// Temperature derivative of ductile strength
	ductileDiff = viscosity(T, P, flowLawDiff, grainSize, dtime); // Korenaga & Karato (2008), dry diffusion
	ductileDisl = viscosity(T, P, flowLawDisl, grainSize, dtime); // Korenaga & Karato (2008), dry dislocation

	ductileDiff_prime = dviscdT(T, dPdT, flowLawDiff, grainSize,  dtime, Tsurf)/dtime;
	ductileDisl_prime = dviscdT(T, dPdT, flowLawDisl, grainSize,  dtime, Tsurf)/dtime;

	// Easiest to take the derivative of the parallel combination as d/dT [Diff*Disl/(Diff+Disl)]
	ductile_prime = ( (ductileDiff_prime*ductileDisl + ductileDiff*ductileDisl_prime) * (ductileDiff+ductileDisl)
			         -(ductileDiff_prime+ductileDisl_prime) * ductileDiff*ductileDisl ) / pow(ductileDiff_prime + ductileDisl_prime,2);

	f_prime = brittle_prime - ductile_prime;

	return f_prime;
}

/*--------------------------------------------------------------------
 *
 * Subroutine viscosity
 *
 * Calculates non-Newtonian viscosity in Pa s.
 *
 *--------------------------------------------------------------------*/

double viscosity (double T, double P, double flowLaw[5], double grainSize, double dtime) {

	double visc = 0.0;

	visc = MPa2Pa*pow(pow(10.0,-flowLaw[2]) * pow(grainSize,flowLaw[4]) * exp((flowLaw[0] + P*flowLaw[1])/(R_G*T)), 1.0/flowLaw[3]) * pow(1.0/dtime, 1.0/flowLaw[3]-1.0);

	return visc;
}

/*--------------------------------------------------------------------
 *
 * Subroutine comboVisc
 *
 * Calculates non-Newtonian viscosity in Pa s combined from parallel
 * diffusion and dislocation flow laws (e.g. Cizkova et al. 2012 eq. 1).
 *
 *--------------------------------------------------------------------*/

double combVisc (double T, double P, double flowLawDiff[5], double flowLawDisl[5], double grainSize, double dtime) {

	double visc = 0.0;
	double viscDryDiff = 0.0;
	double viscDryDisl = 0.0;

	viscDryDiff = viscosity(T, P, flowLawDiff, grainSize, dtime); // Korenaga & Karato (2008), dry diffusion
	viscDryDisl = viscosity(T, P, flowLawDisl, grainSize, dtime); // Korenaga & Karato (2008), dry dislocation
	visc = 1.0/(1.0/viscDryDiff + 1.0/viscDryDisl); // Parallel combination

	return visc;
}

/*--------------------------------------------------------------------
 *
 * Subroutine dviscdT
 *
 * Calculates temperature derivative of non-Newtonian viscosity along
 * linear geotherm in Pa s K-1.
 *
 *--------------------------------------------------------------------*/

double dviscdT (double T, double dPdT, double flowLaw[5], double grainSize, double dtime, double Tsurf) {

	double dviscdT = 0.0;

	dviscdT = MPa2Pa*pow(pow(10.0,-flowLaw[2]) * pow(grainSize,flowLaw[4]) * exp(dPdT*flowLaw[1]/R_G) * exp((flowLaw[0] - dPdT*flowLaw[1]*Tsurf)/(R_G*T)), 1.0/flowLaw[3])
	          * (dPdT*flowLaw[1]*Tsurf-flowLaw[0])/(flowLaw[3]*R_G*T*T) * pow(1.0/dtime, 1.0/flowLaw[3]-1.0);

	return dviscdT;
}
