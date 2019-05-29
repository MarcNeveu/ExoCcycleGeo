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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <R.h>                          // To use the external R software package
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rembedded.h>
#include <IPhreeqc.h>                   // To use the external PHREEQC geochemical code
#include "CHNOSZcmds.h"
#include <unistd.h>
#include <fcntl.h>

#define cmdline 0						// If execution from terminal as "./ExoCcycleGeo"

// Constants and conversion factors
#define PI_greek 3.14159265358979323846
#define G 6.67e-11                      // Gravitational constant (SI)
#define R_G 8.314                       // Universal gas constant, same number of decimal places as in PHREEQC (J mol-1 K-1)
#define sigStefBoltz 5.6704e-8          // Stefan-Boltzmann constant (W m-2 K-4)
#define rEarth 6378000.0                // Earth radius (m)
#define mEarth 6.0e24                   // Earth mass (kg)
#define km2m 1000.0                     // Convert from km to m
#define m2cm 100.0                      // Convert m to cm
#define bar2Pa 1.0e5                    // Convert bar to Pa
#define atm2bar 1.01325                 // Convert atm to bar, same conversion factor as for PHREEQC
#define Kelvin 273.15                   // Convert between Celsius and Kelvin, same conversion factor as for PHREEQC
#define Gyr2sec 3.15576e16              // Convert 1 billion years to s
#define Myr2sec 3.15576e13              // Convert 1 million years to s
#define Yr2sec 3.15576e7                // Convert 1 year to s

// Earth scalings
#define TsurfEarth 288.0                // Mean equilibrium surface temperature on Earth (K)
#define RmeltEarth 3.8e-19              // Rate of melt generation on Earth (s-1) Kite et al. 2009 Fig. 15; http://dx.doi.org/10.1088/0004-637X/700/2/1732
#define deltaCvolcEarth 2.2e5           // Surface C flux from subaerial+submarine volcanic outgassing (mol C s-1) Donnadieu et al. 2006; http://dx.doi.org/10.1029/2006GC001278
#define deltaCcontwEarth 8.4543e-10     // Surface C flux from continental weathering on Earth (mol C m-2 s-1)

// Mantle parameters TODO make function of composition?
#define k 4.18                          // Mantle thermal conductivity (W m-1 K-1)
#define alpha 3.0e-5                    // Mantle thermal expansivity (K-1)
#define Cp 914.0                        // Mantle heat capacity (J kg-1 K-1)
#define rhoMantle 4000.0                // Mantle density of the mantle TODO get through compression() routine of IcyDwarf?

// Thermal model parameters
#define beta 0.3                        // Exponent for scaling Nusselt number to Rayleigh number, 1/4 to 1/3 (Schubert et al. 2001)
#define Ra_c 1707.762                   // Critical Rayleigh number for convection, http://home.iitk.ac.in/~sghorai/NOTES/benard/node15.html, see also Koschmieder EL (1993) Benard cells and Taylor vortices, Cambridge U Press, p. 20.
#define A0 70000.0                      // Activation temperature for thermal model (K)

// Atmosphere parameters
#define xCO2g0 355.0e-6                 // Reference atmospheric CO2 mixing ratio (ppmv)
#define runoff_0 7.70e-9                // Reference runoff (m s-1) = 0.665e-3 m day-1 (Edson et al. 2012, http://dx.doi.org/10.1089/ast.2011.0762)

// Geochem parameters
#define nAtmSpecies 5                   // Number of atmospheric species whose abundances are passed between the physical and chemical models
#define nAqSpecies 6                    // Number of aqueous species whose concentrations are passed between the physical and chemical models
#define nvarEq 1036                     // Number of geochemical variables stored in each PHREEQC equilibrium simulation
#define nvarKin 60                      // Number of geochemical variables stored in each PHREEQC kinetic simulation
#define naq 258                         // Number of aqueous species (+ physical parameters)
#define nmingas 389                     // Number of minerals and gases
#define nelts 31                        // 30 elements + 1 extra column in WaterRock/Molar_masses.txt

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
int strain (double Pressure, double Xhydr, double T, double *strain_rate, double *Brittle_strength, double porosity);

//-------------------------------------------------------------------
// MAIN PROGRAM
//-------------------------------------------------------------------

int main(int argc, char *argv[]) {

	int i = 0;
	int j = 0;
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
	double d = 0.0;                 // Depth to core-mantle boundary (m)
	double Tmantle = 0.0;           // Mantle temperature (K)
	double H = 0.0;                 // Specific radiogenic heating rate (J s-1 kg-1)
	double kappa = k/(rhoMantle*Cp); // Mantle thermal diffusivity (m2 s-1)
	double nu = 0.0;                // Mantle kinematic viscosity (m2 s-1)
	double zCrust = 0.0;            // Crustal thickness (m)
	double Ra = 0.0;                // Rayleigh number for mantle convection (no dim)
	double Nu = 0.0;                // Nusselt number (no dim)
	double Tref = 0.0;              // Temperature at outer boundary of convective zone (surface or base of stagnant lid)

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
	double m_p = 0.5*mEarth;    // Planet mass (kg)
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
	printf("ExoCcycleGeo v18.12\n");
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

	r_p = pow(m_p/mEarth,0.27)*rEarth; // Planet radius, determined from mass scaling (m) // TODO implement IcyDwarf compression model or ExoPLEX for more accurate mass-radius relationship
	r_c = pow(0.007*pow(r_p,3.0),1.0/3.0); // Radius of planet core (m) such that the volume of the core is 0.7% of the overall volume as for Earth (r_c = 1220 km)
	printf("Core radius = %g km\n", r_c/1000.0);
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

		/* The geophysical scaling laws below for heat flux, outgassing rate, and subduction rate are solely valid for limited and very
		 * specific cases and are not generalizable.
		 *
		 * Vary tectonic mode with planet mass. More massive planets could be less effective at subducting.
		 */

		// Kinematic viscosity
		nu = 1.0e16*exp((2.0e5 + (rhoMantle*gsurf*d/2.0) *1.1e-6)/(R_G*Tmantle))/rhoMantle; // Cízková et al. (2012)

		// Compute instantaneous heating rate (Kite et al. 2009 Table 1): H = X_4.5 * W * exp(ln(1/2) / t1/2 * (t-4.5))
		H =  36.9e-9  * 2.92e-5 * exp(log(0.5)/( 1.26 *Gyr2sec) * (realtime - 4.5*Gyr2sec))  //  40-K
		  + 124.0e-9  * 2.64e-5 * exp(log(0.5)/(14.0  *Gyr2sec) * (realtime - 4.5*Gyr2sec))  // 232-Th
		  +   0.22e-9 * 56.9e-5 * exp(log(0.5)/( 0.704*Gyr2sec) * (realtime - 4.5*Gyr2sec))  // 235-U
		  +  30.8e-9  * 9.46e-5 * exp(log(0.5)/( 4.47 *Gyr2sec) * (realtime - 4.5*Gyr2sec)); // 238-U

		// Compute effective thermal conductivity
		if (!staglid) Tref = Tsurf; // TODO circular logic to then have value of !staglid depend on Tref below.
		else          Tref = Tmantle - 2.23*Tmantle*Tmantle/A0; // Tmantle - Tc in K09

		Ra = alpha*gsurf*(Tmantle-Tref)*pow(d,3)/(kappa*nu);
		Nu = pow(Ra/Ra_c, beta);

		// Convective velocity (K09 equation 24), replaced Tc (temp at base of stagnant lid) with Tref for compatibility with !staglid regime
		vConv = 2.0*(Nu-1.0) * (k/(rhoMagma*Cp*(r_p-r_c))) * (Tmantle-Tsurf) / (Tmantle-Tref); // Check against eq. (6.379) of Turcotte & Schubert (2002), p. 514

		// Determine if convective driving stress exceeds lithosphere yield stress, in which case switch to mobile-lid (plate tectonics) regime
		// Convective drive stress:
		// Viscosity = stress / strain rate (assumed Newtonian - see Deschamps & Sotin 2000), i.e. stress = viscosity * vConv/δ, w/ δ: boundary layer thickness
		// δ = thickness of thermal gradient = 2.32*(kappa*δ/vConv)^0.5 (Turcotte & Schubert 2002, eq. 6.327 with δ/vConv = conduction timescale)
		// So δ^0.5 = 2.32*(kappa/vConv)^0.5 and δ ≈ 5.38*kappa/vConv
		double driveStress = nu*rhoMantle*vConv*vConv/(5.38*kappa);

		// Lithospheric yield stress = strength at brittle-ductile transition TODO

		// Melting model of Kite et al. (2009)
		zCrust = 10.0*km2m;
		if (!staglid) Pf = 0.0; // K09 equations (10-12)
		else Pf = rhoCrust*gsurf*zCrust;
		P0 = rhoCrust*gsurf*(10.0*zCrust); // Should be higher than Pf. Arbitrary for now,  TODO calculate based on brittle-ductile transition?, 2.0*zCrust results in 50% melting rate

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
 * Subroutine strain
 *
 * Calculates the brittle strength in Pa and corresponding ductile
 * strain rate in s-1 of silicate rock.
 * Brittle-ductile and brittle-plastic transitions are
 * mixed up, although they shouldn't (Kohlstedt et al. 1995).
 * The brittle strength is given by a friction/low-P Byerlee type law:
 * stress = mu*P, assuming negligible water pressure since in practice
 * the brittle-ductile transition occurs in dehydrated rock (T>700 K)
 * even over long time scales.
 * The ductile strength is given by a flow law:
 * d epsilon/dt = A*sigma^n*d^-p*exp[(-Ea+P*V)/RT].
 *
 *--------------------------------------------------------------------*/

int strain (double Pressure, double Xhydr, double T, double *strain_rate, double *Brittle_strength, double porosity) {

	double Hydr_strength = 0.0;
	double Dry_strength = 0.0;

//	Hydr_strength = mu_f_serp*Pressure;
//	if (Pressure <= 200.0e6) Dry_strength = mu_f_Byerlee_loP*Pressure;
//	else Dry_strength = mu_f_Byerlee_hiP*Pressure + C_f_Byerlee_hiP;
//	(*Brittle_strength) = Xhydr*Hydr_strength + (1.0-Xhydr)*Dry_strength;
//	(*Brittle_strength) = (*Brittle_strength)/(1.0-porosity);
//
//	if (T > 140.0)
//		(*strain_rate) = pow(10.0,5.62)*pow((*Brittle_strength)/MPa,1.0)*pow(d_flow_law,-3.0)*exp((-240.0e3 + Pressure*0.0)/(1.0*R_G*T));
//	else // Set T at 140 K to calculate ductile strength so that it doesn't yield numbers too high to handle
//		(*strain_rate) = pow(10.0,5.62)*pow((*Brittle_strength)/MPa,1.0)*pow(d_flow_law,-3.0)*exp((-240.0e3 + Pressure*0.0)/(1.0*R_G*140.0));

	return 0;
}
