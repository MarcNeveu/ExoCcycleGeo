/*
 ============================================================================
 Name        : ExoCcycleGeo.c
 Author      : Marc Neveu

 Computes net C fluxes (in mol C m-2 s-1) at the surface-atmosphere interface
 of a terrestrial planet due to geophysical and geochemical processes.
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

#define PI_greek 3.14159265358979323846
#define G 6.67e-11                      // Gravitational constant (SI)
#define rEarth 6378000.0                // Earth radius (m)
#define mEarth 6.0e24                   // Earth mass (kg)
#define km2m 1000.0                     // Convert from km to m
#define m2cm 100.0                      // Convert m to cm
#define bar2Pa 101325.0                 // Convert bar to Pa
#define Kelvin 273.15                   // Convert between Celsius and Kelvin
#define Myr2sec 31.5576e13              // Convert 1 million years to s

#define TsurfEarth 288.0                // Mean equilibrium surface temperature on Earth (K)
#define RmeltEarth 3.8e-19              // Rate of melt generation on Earth (s-1) Kite et al. 2009 Fig. 15; http://dx.doi.org/10.1088/0004-637X/700/2/1732
#define deltaCvolcEarth 2.2e5           // Surface C flux from subaerial+submarine volcanic outgassing (mol C s-1) Donnadieu et al. 2006; http://dx.doi.org/10.1029/2006GC001278
#define deltaCcontwEarth 8.4543e-10     // Surface C flux from continental weathering on Earth (mol C m-2 s-1)

#define Ra_c 1707.762                   // Critical Rayleigh number for convection, http://home.iitk.ac.in/~sghorai/NOTES/benard/node15.html, see also Koschmieder EL (1993) Benard cells and Taylor vortices, Cambridge U Press, p. 20.
#define Cp 914.0                        // Specific heat capacity of mantle rock (J kg-1 K-1)
#define kcrust 4.18                     // Thermal conductivity of crust (W m-1 K-1)
#define A0 70000.0                      // Activation temperature for thermal model (K)

#define xCO2g0 355.0e-6                 // Reference atmospheric CO2 mixing ratio (ppmv)
#define runoff_0 0.665                  // Reference runoff (mm/day) Edson et al. 2012, http://dx.doi.org/10.1089/ast.2011.0762
#define nAtmSpecies 4                   // Number of atmospheric species whose abundances are passed between the physical and chemical models
#define nAqSpecies 6                    // Number of aqueous species whose concentrations are passed between the physical and chemical models

//-------------------------------------------------------------------
// WATER-ROCK PARAMETERS
//-------------------------------------------------------------------

#define nvar 1036                       // Number of geochemical variables stored in each PHREEQC simulation
#define naq 258                         // Number of aqueous species (+ physical parameters)
#define nmingas 389                     // Number of minerals and gases
#define nelts 31                        // 30 elements + 1 extra column in WaterRock/Molar_masses.txt

//-------------------------------------------------------------------
// SUBROUTINE DECLARATIONS
//-------------------------------------------------------------------

int AqueousChem (char path[1024], char filename[64], int itime, double T, double P, double *pH, double *pe, double mass_w,
		double **xgas, double **xaq, int forcedPP, double kintime, int kinsteps);
int ExtractWrite(int instance, double** data, int line);
const char* ConCat (const char *str1, const char *str2);
int WritePHREEQCInput(const char *TemplateFile, int itime, double temp, double pressure, double pH, double pe, double mass_w,
		double *xgas, double *xaq, int forcedPP, double kintime, int kinsteps, char **tempinput);
int cleanup (char path[1024]);
double molmass_atm (double *xgas);

//-------------------------------------------------------------------
// MAIN PROGRAM
//-------------------------------------------------------------------

int main(int argc, char *argv[]) {

	int i = 0;
	int itime = 0;                  // Time step counter
	int ntime = 0;                  // Total number of steps

	int iter = 0;                   // Iteration counter
	int subiter = 0;                // Secondary iteration counter
	int niter = 20;                 // Number of iterations for ocean-atmosphere equilibrium
	int moveon = 0;                 // Switch to break for loop
	double threshold = 0.1;         // Threshold for convergence of ocean composition in equilibrium with prescribed atmosphere

	double r_p = 0.0;               // Planet radius (m)

//	double RCplate = 0.0;           // Plate/crust C reservoir (mol)
	double RCmantle = 0.0;          // Mantle C reservoir (mol)
	double RCatm = 0.0;             // Atmospheric C reservoir (mol)
	double RCocean = 0.0;           // Ocean C reservoir (mol)
	double RCatmoc = RCatm + RCocean; // Combined atmospheric and ocean C reservoir (mol)

//	double nCatm = 0.0;             // Atmospheric C reservoir (mol), = RCatm but for scaling
//	double nCocean = 0.0;           // Ocean C reservoir (mol), = RCocean but for scaling

	double FCoutgas = 0.0;          // C flux from outgassing (subaerial+submarine) (mol s-1)
	double FCcontw = 0.0;           // C flux from continental weathering (mol s-1)
	double FCseafw = 0.0;           // C flux from seafloor weathering (mol s-1)
	double FCsubd = 0.0;            // C flux from subduction (mol s-1)
	double netFC = 0.0;             // Net surface C flux from all geological processes (mol s-1)

	double farc = 0.0;              // Fraction of subducted C that makes it back out through arc volcanism (as opposed to into the mantle)

	double pH = 0.0;                // pH of the surface ocean (default 8.22)
	double pe = 0.0;                // pe corresponding to logfO2
	double logfO2 = 0.0;            // O2 fugacity
	double logKO2H2O = 0.0;         // log K for reaction 4 H+ + 4 e- + O2 = 2 H2O, from CHNOSZ: subcrt(c("H+","e-","O2","H2O"),c(-4,-4,-1,2),c("aq","aq","g","liq"),T=25,P=1)

	double Rmelt = 0.0;             // Rate of melt generation (m-2 s-1)
	double vConv = 0.0;             // Convective velocity (m s-1)
	double gsurf = 0.0;             // Surface gravity (m s-2)
	double P0 = 0.0;                // Surface pressure (Pa)
	double Pf = 0.0;                // Pressure at base of crust (Pa)
	double TbaseLid = 0.0;          // Temperature at base of stagnant lid (K)
	double Nu = 0.0;                // Nusselt number (no dim)
	double zCrack = 0.0;            // Depth of fracturing below seafloor (m)
	double tcirc = 0.0;             // Time scale of hydrothermal circulation (s)
	double deltaCreac = 0.0;        // Net C leached/precipitated per kg of rock (mol kg-1)

	double Asurf = 0.0;             // Planet surface area (m2)
	double nAir = 0.0;              // Number of mol in atmosphere (mol)
	double sumPP = 0.0;             // Sum of partial pressures xgas*Psurf (ideally, should be equal to Psurf)
	double sumMR = 0.0;             // Sum of mixing ratios xgas (ideally, should be equal to 1)

	double nCatmoc = 0.0;           // Total C in {atmosphere+ocean}
	double nNatmoc = 0.0;           // Total N in {atmosphere+ocean}
	double nCatmoc_old = 0.0;
	double nNatmoc_old = 0.0;
	double gasdiff = 0.0;           // Difference in total gas mixing ratios before and after atmosphere-ocean equilibration

	// Kinetic parameters
	int kinsteps = 0;               // Number of time steps of PHREEQC kinetic simulation
	double kintime = 0.0;           // Total time of PHREEQC kinetic simulation

	// Quantities to be computed by thermal/geodynamic model
	double zCrust = 0.0;        // Crustal thickness (m)
	double Tmantle = 0.0;       // Mantle temperature (K)
	double Ra = 0.0;            // Rayleigh number for mantle convection (no dim)

	double *xgas = (double*) malloc(nAtmSpecies*sizeof(double));
	if (xgas == NULL) printf("ExoCcycleGeo: Not enough memory to create xgas[nAtmSpecies]\n"); // Mixing ratios by volume (or by mol since all gases are pretty much ideal and have the same molar volume) of atmospheric gases
    for (i=0;i<nAtmSpecies;i++) xgas[i] = 0.0;

	double *xaq = (double*) malloc(nAqSpecies*sizeof(double));
	if (xaq == NULL) printf("ExoCcycleGeo: Not enough memory to create xaq[nAqSpecies]\n"); // Molalities of aqueous species (mol (kg H2O)-1)
    for (i=0;i<nAqSpecies;i++) xaq[i] = 0.0;

	double *xgas_old = (double*) malloc(nAtmSpecies*sizeof(double));
	if (xgas_old == NULL) printf("ExoCcycleGeo: Not enough memory to create xgas_old[nAtmSpecies]\n"); // Old mixing ratios by volume (or by mol since all gases are pretty much ideal and have the same molar volume) of atmospheric gases
    for (i=0;i<nAtmSpecies;i++) xgas_old[i] = 0.0;

	double *xaq_old = (double*) malloc(nAqSpecies*sizeof(double));
	if (xaq_old == NULL) printf("ExoCcycleGeo: Not enough memory to create xaq_old[nAqSpecies]\n"); // Old molalities of aqueous species (mol (kg H2O)-1)
    for (i=0;i<nAqSpecies;i++) xaq_old[i] = 0.0;

	FILE *fout;
	char *title = (char*)malloc(1024*sizeof(char)); // Don't forget to free!

	//-------------------------------------------------------------------
	// Inputs
	//-------------------------------------------------------------------

	double dtime = 1.0*Myr2sec;     // Time step (s)

	ntime = (int) (50.0*Myr2sec/dtime); // Number of time steps of simulation

	// Planet parameters
	double m_p = 1.0*mEarth;    // Planet mass (kg)
	double L = 0.15;            // Fraction of planet surface covered by land
	double Mocean = 1.4e21;     // Mass of ocean (kg, default: Earth=1.4e21)
	double r_c = 1220000.0;     // Radius of planet core (m)
	double rhoMagma = 3500.0;   // Magma density (kg m-3)
	double rhoCrust = 3500.0;   // Crustal density (kg m-3)

	// Atmospheric inputs
	double Tsurf = 288.15;      // Surface temperature (K)
	double Psurf = 1.0;         // Surface pressure (bar)
	double runoff = 0.7;        // Atmospheric runoff (mm/day)
	xgas[0] = 1.0e-6;    // CO2
    xgas[1] = 0.0;       // CH4
    xgas[2] = 0.0;       // O2
    xgas[3] = 1.0;       // N2

    // Ocean inputs
	pH = 8.22;
    int redox = 1; // 1: current Earth surface, 2: hematite-magnetite, 3: fayalite-magnetite-quartz, code won't run with other values.
    xaq[0] = 27.0/12.0/1000.0;   // 27 ppm C in today's oceans (scaled from 141 ppm Alk*M(C)/M(HCO3), M being molecular mass)
    xaq[1] = 0.0/12.0/1000.0;
    xaq[3] = 0.0/14.0/1000.0;   // ppm N in ocean

	//-------------------------------------------------------------------
	// Startup
	//-------------------------------------------------------------------

	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("ExoCcycleGeo v18.10\n");
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
	gsurf = G*m_p/r_p/r_p;
	Asurf = 4.0*PI_greek*r_p*r_p;
	for(i=0;i<nAtmSpecies;i++) sumPP = sumPP + xgas[i]*Psurf;

	//-------------------------------------------------------------------
	// Choose redox state (from most oxidized to most reduced)
	//-------------------------------------------------------------------

	switch(redox) {
	case 1: // Present-day Earth surface
		printf("Redox set to present-day Earth surface\n");
		logfO2 = log(0.2)/log(10.0);
		printf("log f(O2) = %g\n", logfO2);
		break;
	case 2: // Use CHNOSZ to get log fO2 for hematite-magnetite (HM) buffer at given T and P.
		printf("Redox set to HM buffer\n");
		logfO2 = -6.0*CHNOSZ_logK("hematite", "cr", Tsurf, Psurf, "SUPCRT92")
			     +4.0*CHNOSZ_logK("magnetite", "cr", Tsurf, Psurf, "SUPCRT92")
			     +1.0*CHNOSZ_logK("O2", "g", Tsurf, Psurf, "SUPCRT92");
		logKO2H2O = -4.0*CHNOSZ_logK("H+", "aq", Tsurf, Psurf, "SUPCRT92")
					-4.0*CHNOSZ_logK("e-", "aq", Tsurf, Psurf, "SUPCRT92")
					-1.0*CHNOSZ_logK("O2", "g", Tsurf, Psurf, "SUPCRT92")
					+2.0*CHNOSZ_logK("H2O", "liq", Tsurf, Psurf, "SUPCRT92");
		printf("log f(O2) = %g\n", logfO2);
		break;
	case 3: // Use CHNOSZ to get log fO2 for fayalite-magnetite-quartz (FMQ) buffer at given T and P.
		printf("Redox set to FMQ buffer\n");
		logfO2 = -3.0*CHNOSZ_logK("quartz", "cr", Tsurf, Psurf, "SUPCRT92")
			     -2.0*CHNOSZ_logK("magnetite", "cr", Tsurf, Psurf, "SUPCRT92")
		         +3.0*CHNOSZ_logK("fayalite", "cr", Tsurf, Psurf, "SUPCRT92")
			     +1.0*CHNOSZ_logK("O2", "g", Tsurf, Psurf, "SUPCRT92");
		logKO2H2O = -4.0*CHNOSZ_logK("H+", "aq", Tsurf, Psurf, "SUPCRT92")
					-4.0*CHNOSZ_logK("e-", "aq", Tsurf, Psurf, "SUPCRT92")
					-1.0*CHNOSZ_logK("O2", "g", Tsurf, Psurf, "SUPCRT92")
					+2.0*CHNOSZ_logK("H2O", "liq", Tsurf, Psurf, "SUPCRT92");
		printf("log f(O2) = %g\n", logfO2);
		break;
	default:
		printf("ExoCcycleGeo: Redox switch incorrectly specified, should be 1 (current Earth surface), 2 (hematite-magnetite), or 3 (fayalite-magnetite-quartz). Exiting.\n");
		exit(0);
	}

	pe = -pH + 0.25*(logfO2+logKO2H2O);

	//-------------------------------------------------------------------
	// Initialize reservoirs
	//-------------------------------------------------------------------

	nAir = Psurf*bar2Pa*Asurf/gsurf/molmass_atm(xgas);
	RCatm = (xgas[0]+xgas[1])*nAir;

	// Equilibrate ocean and atmosphere at input pressure and atmospheric C/N ratio
	if (Psurf < 0.01) {
		printf("ExoCcycleGeo: Pressure = %g bar too close to the triple point of H2O, oceans not stable at the surface. Exiting.\n", Psurf);
		exit(0);
	}
	if (Mocean/nAir < 0.1) {
		printf("ExoCcycleGeo: Mass Ocean / moles of air = %g too low to compute ocean-atmosphere equilibrium. Exiting.\n", Mocean/nAir);
//			exit(0);
	}

	printf("Equilibrating ocean and atmosphere at input pressure and atmospheric C/N ratio...\n");
	for (i=0;i<nAtmSpecies;i++) {
		if (xgas[i] > 0.0 && xaq[i] == 0.0) xaq[i] = xgas[i]; // xaq must be >0 otherwise PHREEQC ignores, set to xgas (initial guess).
	}
	AqueousChem(path, "io/OceanStart", itime, Tsurf, Psurf, &pH, &pe, Mocean/nAir, &xgas, &xaq, 1, 0.0, 0);

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
		fprintf(fout, "'Time (Myr)' \t RCmantle \t RCatmoc \t RCocean \t FCoutgas \t FCcontw \t FCseafw \t FCsubd \t 'Net C flux'\n");
		fprintf(fout, "Init \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n", RCmantle, RCatmoc, RCocean, FCoutgas, FCcontw, FCseafw, FCsubd, netFC);
	}
	fclose (fout);

	title[0] = '\0';
	if (cmdline == 1) strncat(title,path,strlen(path)-20);
	else strncat(title,path,strlen(path)-18);
	strcat(title,"CompoOceanAtm.txt");
	fout = fopen(title,"w");
	if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
	else {
		fprintf(fout, "'Atmospheric species in mixing ratio (by mol, i.e. by volume for ideal gases), aqueous species in mol/(kg H2O)'\n");
		fprintf(fout, "'Time (Myr)' \t CO2(g) \t CH4(g) \t N2(g) \t 'P_surface (bar)' \t 'Ox C(aq)' \t 'Red C(aq)' \t 'Total N(aq)' \t pH \n");
		fprintf(fout, "Init \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n", xgas[0], xgas[1], xgas[3], Psurf, xaq[0], xaq[1], xaq[3], pH);
	}
	fclose (fout);

	//-------------------------------------------------------------------
	// Start time loop
	//-------------------------------------------------------------------

	printf("Starting time loop...\n");
	for (itime = 0;itime<ntime;itime++) {
		printf("Time: %g Myr, iteration %d/%d. Total N = %g mol, Total C = %g mol, Added C = %g mol\n", (double)itime*dtime/Myr2sec, itime, ntime,
				xaq[3]*Mocean + xgas[3]*2.0*nAir, (xaq[0]+xaq[1])*Mocean + (xgas[0]+xgas[1])*nAir + netFC*dtime, netFC*dtime*(double)itime);

		//-------------------------------------------------------------------
		// Update atmosphere
		//-------------------------------------------------------------------
//		Psurf_old = Psurf;

		for (i=2;i<nAtmSpecies;i++) {
			if (nAir > fabs(dtime*netFC)) xgas[i] = xgas[i]*nAir/(nAir+dtime*netFC); // Dilute other gases accordingly
			else {
				printf("ExoCcycleGeo: Too much atmospheric change in one timestep. Reduce timestep.\n"); exit(0);
			}
		}

		// Redox of outgassing
		if      (xgas[0] > 0.0 && xgas[0] > 10.0*xgas[1]) xgas[0] = (xgas[0]*nAir + dtime*netFC)/(nAir + dtime*netFC); // CO2 dominates over CH4, assume all added gas is CO2 and let equilibration with ocean speciate accurately
		else if (xgas[1] > 0.0 && xgas[1] > 10.0*xgas[0]) xgas[1] = (xgas[1]*nAir + dtime*netFC)/(nAir + dtime*netFC); // CH4 dominates over CO2, assume all added gas is CH4
		else {
			xgas[0] = (xgas[0]*nAir + dtime*netFC/2.0)/(nAir + dtime*netFC); // CH4 and CO2 about equal, split added gas equally
			xgas[1] = (xgas[1]*nAir + dtime*netFC/2.0)/(nAir + dtime*netFC);
		}

		nAir = nAir + dtime*netFC;
		Psurf = nAir/(bar2Pa*Asurf/gsurf/molmass_atm(xgas));

		if (Psurf < 0.01) {
			printf("ExoCcycleGeo: Pressure = %g bar too close to the triple point of H2O, oceans not stable at the surface. Exiting.\n", Psurf);
			exit(0);
		}
		if (Mocean/nAir < 0.1) {
			printf("ExoCcycleGeo: Mass Ocean / moles of air = %g too low to compute ocean-atmosphere equilibrium. Exiting.\n", Mocean/nAir);
//			exit(0);
		}

		//-------------------------------------------------------------------
		// Calculate partitioning of C between ocean and atmosphere
		//-------------------------------------------------------------------

		for (iter = 0;iter < niter;iter++) {
			moveon = 1;

			for (i=0;i<nAtmSpecies;i++) xgas_old[i] = xgas[i];
			for (i=0;i<nAqSpecies;i++) xaq_old[i] = xaq[i];

			nCatmoc_old = (xaq[0]+xaq[1])*Mocean + (xgas[0]+xgas[1])*nAir;
			nNatmoc_old = xaq[3]*Mocean + xgas[3]*2.0*nAir;

			// Equilibrate ocean and atmosphere
			for (subiter=0;subiter<niter;subiter++) AqueousChem(path, "io/OceanDissInput", itime, Tsurf, Psurf, &pH, &pe, Mocean/nAir, &xgas, &xaq, 0, 0.0, 0);

			// Force conservation of Total C (atmosphere+ocean) and Total N (atmosphere+ocean)
			nCatmoc = (xaq[0]+xaq[1])*Mocean + (xgas[0]+xgas[1])*nAir;
			nNatmoc = xaq[3]*Mocean + xgas[3]*2.0*nAir;

			xgas[0] = xgas[0]*nCatmoc_old/nCatmoc;
			xgas[1] = xgas[1]*nCatmoc_old/nCatmoc;
			xgas[3] = xgas[3]*nNatmoc_old/nNatmoc;

			xaq[0] = xaq[0]*nCatmoc_old/nCatmoc;
			xaq[1] = xaq[1]*nCatmoc_old/nCatmoc;
			xaq[3] = xaq[3]*nNatmoc_old/nNatmoc;

			// Update air quantity and surface pressure
			gasdiff = 0.0;
			for (i=0;i<nAtmSpecies;i++) gasdiff = gasdiff + xgas[i] - xgas_old[i];
			nAir = nAir*(1.0 + gasdiff);
			Psurf = nAir/(bar2Pa*Asurf/gsurf/molmass_atm(xgas));

			if (Psurf < 0.01) {
				printf("ExoCcycleGeo: Pressure = %g bar too close to the triple point of H2O, oceans not stable at the surface. Exiting.\n", Psurf);
				exit(0);
			}
			if (Mocean/nAir < 0.1) {
				printf("ExoCcycleGeo: Mass Ocean / moles of air = %g too low to compute ocean-atmosphere equilibrium. Exiting.\n", Mocean/nAir);
	//			exit(0);
			}

			// Rescale mixing ratios
			sumMR = 0.0;
			for (i=0;i<nAtmSpecies;i++) sumMR = sumMR + xgas[i];
			for (i=0;i<nAtmSpecies;i++) xgas[i] = xgas[i]/sumMR;

			for(i=0;i<nAtmSpecies;i++) {
				if (i!=2 && xgas[i] > 1.0e-10) { // Bottom threshold on abundance to avoid spending too much time converging on negligible species
					if (fabs(xgas[i]/xgas_old[i]-1.0) > threshold) moveon = 0;
				}
			}
			for(i=0;i<nAqSpecies;i++) {
				if (i!=2 && xaq[i] > 1.0e-10) { // Bottom threshold on abundance to avoid spending too much time converging on negligible species
					if (fabs(xaq[i]/xaq_old[i]-1.0) > threshold) moveon = 0;
				}
			}

			cleanup(path); // Remove PHREEQC selected output files

			if (moveon) break;
		}

		//-------------------------------------------------------------------
		// Calculate surface C flux from outgassing
		//-------------------------------------------------------------------

		/* TODO compute the following quantities instead of assuming them. See Kite et al. Section 2.2.
		 * i.e., compute how geodynamic processes vary with planet mass, structure, composition, and age.
		 * The geophysical scaling laws below for heat flux, outgassing rate, and subduction rate are solely valid for limited and very
		 * specific cases and are not generalizable.
		 * Vary tectonic mode with planet mass.
		 *
		 * Subduction and volcanism may not be more effective on more massive planets, as they could be stagnant-lid planets or have
		 * intrusive volcanism.
		 *
		 * See Stein et al. (2011), Deschamps & Sotin (2000), Stamenkovic et al. (2012), Wong & Solomatov (2016), Burgisser & Scaillet
		 * (2007), and Dasgupta & Hirschmann (2010).
		 *
		 * The mantle reservoir of carbon is large enough that CO2 outgassing does not depend on the amount of carbon subducted into the mantle (Abbot et al. 2012).
		 */

		// Model of Kite et al. (2009)
		zCrust = 10.0*km2m;
		Tmantle = 3000.0;
		Ra = 3000.0;

		Pf = rhoCrust*gsurf*zCrust;
		TbaseLid = Tmantle - 2.23*Tmantle*Tmantle/A0;
		Nu = pow(Ra/Ra_c,0.25);

		// Assumes that all melt generated reaches the surface
		vConv = 2.0*(Nu-1.0) * (kcrust/(rhoMagma*Cp*(r_p-r_c))) * (Tmantle-Tsurf) / (Tmantle-TbaseLid);
		Rmelt = vConv*rhoMagma/m_p * (rhoCrust*gsurf*zCrust/(Pf-P0)); // m-2 s-1, Kite et al. (2009) multiplied by planet's surface area
		FCoutgas = deltaCvolcEarth * Rmelt/RmeltEarth * Asurf;              // (mol C s-1) * (m-2 s-1) / (s-1) * m2 = mol C s-1

		//-------------------------------------------------------------------
		// Calculate surface C flux from continental weathering
		//-------------------------------------------------------------------

		kintime = 1.0*365.25*86400.0;
		kinsteps = 20;
		FCcontw = -L * 0.5*deltaCcontwEarth*Asurf * pow(xgas[0]/xCO2g0,0.3) * runoff/runoff_0 * exp((Tsurf-TsurfEarth)/17.7); // Edson et al. (2012) Eq. 1; Abbot et al. (2012) Eq. 2
		AqueousChem(path, "io/ContWeather", itime, Tsurf, Psurf, &pH, &pe, Mocean/nAir, &xgas, &xaq, 1, kintime, kinsteps);

		//-------------------------------------------------------------------
		// Calculate surface C flux from seafloor weathering TODO include kinetics, manage reservoir size
		//-------------------------------------------------------------------

		deltaCreac = 0.0; // TODO call PHREEQC to get deltaCreac, net mol C leached/precipitated per kg of rock
		tcirc = 1.0;      // TODO compute tcirc based on Nu(Ra(basal heat flux))
		FCseafw = -(1.0-L) * 4.0/3.0*PI_greek*(pow(r_p,3)-pow(r_p-zCrack,3))/tcirc*deltaCreac*rhoCrust / (4.0*PI_greek*r_p*r_p);

		//-------------------------------------------------------------------
		// Calculate surface C flux from subduction
		//-------------------------------------------------------------------

		// Function of mantle convective vigor (planet size and age) as for volcanism
		FCsubd = 0.0;

		//-------------------------------------------------------------------
		// Compute net geo C fluxes
		//-------------------------------------------------------------------

		// Continental weathering fluxes halved because of re-precipitation of carbonates in the ocean dissolved during continental weathering
//		RCplate = RCplate + dtime*(FCcontw/2.0 + FCseafw - FCsubd);
//		RCmantle = RCmantle + dtime*((1.0-farc)*FCsubd - FCoutgas); // Assuming instantaneous mixing into the mantle once subducted (in practice could take 1 Gyr)
		netFC = FCoutgas - FCcontw/2.0 - FCseafw - (1.0-farc)*FCsubd;
		RCmantle = RCmantle - dtime*netFC; // Assuming plate = mantle (unlike Foley et al. 2015) and instantaneous mixing into the mantle once subducted (in practice could take 1 Gyr)
		RCatmoc = RCatmoc + dtime*netFC; // Sum of atmospheric and ocean reservoirs, still needs partitioning

//		// Split RCatmoc between RCatm and RCocean
//		molmass_atm = (xgas[0]*44.0 + xgas[1]*16.0 + xgas[2]*32.0 + xgas[3]*28.0)/(xgas[0]+xgas[1]+xgas[2]+xgas[3])*0.001;
//		nCatm = xgas[0]*Psurf*bar2Pa*Asurf/gsurf/molmass_atm;
//		nCocean = xaq[0]*Mocean;
//		fatm = nCatm / (nCocean + nCatm); // RCatm / Rocean+RCatm. P*A/g=mass of atmosphere
//
//		RCatm = fatm * RCatmoc;
//		RCocean = RCatmoc - RCatm;

		// Write outputs
		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"Output.txt");
		fout = fopen(title,"a");
		if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
		else fprintf(fout, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n", (double)itime*dtime/Myr2sec, RCmantle, RCatmoc, RCocean, FCoutgas, FCcontw, FCseafw, FCsubd, netFC);
		fclose (fout);

		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"CompoOceanAtm.txt");
		fout = fopen(title,"a");
		if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
		else fprintf(fout, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n", (double)itime*dtime/Myr2sec, xgas[0], xgas[1], xgas[3], Psurf, xaq[0], xaq[1], xaq[3], pH);
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

int AqueousChem (char path[1024], char filename[64], int itime, double T, double P, double *pH, double *pe, double mass_w,
		double **xgas, double **xaq, int forcedPP, double kintime, int kinsteps) {

	int phreeqc = 0;
//	int i = 0;

	char *dbase = (char*)malloc(1024);                           // Path to thermodynamic database
	char *infile = (char*)malloc(1024);                          // Path to initial (template) input file
	char *tempinput = (char*)malloc(1024);                       // Temporary PHREEQC input file, modified from template

	double *simdata = (double*) malloc(nvar*sizeof(double));
	if (simdata == NULL) printf("AqueousChem: Not enough memory to create simdata[nvar]\n");

	// Initializations
	dbase[0] = '\0';
	infile[0] = '\0';
	tempinput[0] = '\0';

	if (cmdline == 1) strncat(dbase,path,strlen(path)-20);
	else strncat(dbase,path,strlen(path)-18);
	strcat(dbase,"PHREEQC-3.1.2/core10.dat");

	strncat(infile,dbase,strlen(dbase)-10);
	strcat(infile, filename);

	WritePHREEQCInput(infile, itime, T-Kelvin, P, *pH, *pe, mass_w, *xgas, *xaq, forcedPP, kintime, kinsteps, &tempinput);

	phreeqc = CreateIPhreeqc(); // Run PHREEQC
	if (LoadDatabase(phreeqc,dbase) != 0) OutputErrorString(phreeqc);
	SetSelectedOutputFileOn(phreeqc,1);
	if (RunFile(phreeqc,tempinput) != 0) OutputErrorString(phreeqc);

	if (strcmp(filename, "io/OceanStart") == 0) ExtractWrite(phreeqc, &simdata, 1); // 1st line of selected output has initial solution composition = equilibrium when there are no equilibrium phases
	else ExtractWrite(phreeqc, &simdata, 2);                                        // 2nd line of PHREEQC selected output solution and mineral+gas composition at equilibrium

	if (DestroyIPhreeqc(phreeqc) != IPQ_OK) OutputErrorString(phreeqc);

//	for (i=0;i<nvar;i++) printf("%d %g\n", i, simdata[i]);

	// Setting initial ocean chemistry
	if (strcmp(filename, "io/OceanStart") == 0) {
		(*pH) = simdata[7];
		(*pe) = simdata[8];

		(*xaq)[0] = simdata[45];             // C(4), i.e. dissolved CO2 and carbonate
		(*xaq)[1] = simdata[43];             // C(-4), i.e. dissolved methane
		(*xaq)[2] = simdata[72];             // O(0), i.e. dissolved O2
		(*xaq)[3] = simdata[28];             // Total dissolved N
		(*xaq)[4] = simdata[69];             // N(0), i.e. dissolved N2
		(*xaq)[5] = simdata[68];             // N(-3), i.e. dissolved NH3 and NH4+
	}

	// Ocean dissolution
	if (strcmp(filename, "io/OceanDissInput") == 0) {
		(*pH) = simdata[7];
		(*pe) = simdata[8];

		(*xaq)[0] = simdata[45];             // C(4), i.e. dissolved CO2 and carbonate
		(*xaq)[1] = simdata[43];             // C(-4), i.e. dissolved methane
		(*xaq)[2] = simdata[72];             // O(0), i.e. dissolved O2
		(*xaq)[3] = simdata[28];             // Total dissolved N
		(*xaq)[4] = simdata[69];             // N(0), i.e. dissolved N2
		(*xaq)[5] = simdata[68];             // N(-3), i.e. dissolved NH3 and NH4+

		(*xgas)[0] = simdata[1016];          // CO2(g)
		(*xgas)[1] = simdata[1012];          // CH4(g)
		(*xgas)[2] = simdata[1032];          // O2(g)
		(*xgas)[3] = simdata[1024];          // N2(g)
	}

	// Continental weathering
	if (strcmp(filename, "io/ContWeather") == 0) {

	}

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

int ExtractWrite(int instance, double** data, int line) {
	VAR v;
	int i = 0;
	VarInit(&v);

	GetSelectedOutputValue(instance,1,3,&v);           // temp
	(*data)[0] = v.dVal;

	GetSelectedOutputValue(instance,1,1,&v);           // pH
	(*data)[2] = v.dVal;

	GetSelectedOutputValue(instance,1,2,&v);           // pe
	(*data)[3] = v.dVal;

	GetSelectedOutputValue(instance,1,5,&v);           // W:R
	(*data)[6] = v.dVal;

	for (i=1;i<nvar-6;i++) {                           // Rest of parameters
		GetSelectedOutputValue(instance,line,i,&v);
		if (fabs(v.dVal) < 1e-50) (*data)[i+6] = 0.0;
		else (*data)[i+6] = v.dVal;
	}
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

int WritePHREEQCInput(const char *TemplateFile, int itime, double temp, double pressure, double pH, double pe, double mass_w,
		double *xgas, double *xaq, int forcedPP, double kintime, int kinsteps, char **tempinput) {

	int i = 0;

	// Open input file
	FILE *fin;
	FILE *fout;
	char itime_str[5];
	char temp_str[10];
	char pressure_str[10];
	char pH_str[10];
	char pe_str[10];
	char mass_w_str[10];
	char steps_str1[64]; char steps_str2[64];

	char **gas_str1 = (char**)malloc(nAtmSpecies*sizeof(char*));
	for (i=0;i<nAtmSpecies;i++) gas_str1[i] = (char*)malloc(1024);

	char **gas_str2 = (char**)malloc(nAtmSpecies*sizeof(char*));
	for (i=0;i<nAtmSpecies;i++) gas_str2[i] = (char*)malloc(1024);

	char **aq_str = (char**) malloc(nAqSpecies*sizeof(char*));
	for (i=0;i<nAqSpecies;i++) aq_str[i] = (char*)malloc(1024);

	itime_str[0] = '\0';
	temp_str[0] = '\0';
	pressure_str[0] = '\0';
	pH_str[0] = '\0';
	pe_str[0] = '\0';
	mass_w_str[0] = '\0';
	steps_str1[0] = '\0'; steps_str2[0] = '\0';

	for (i=0;i<nAtmSpecies;i++) {
		gas_str1[i][0] = '\0';
		gas_str2[i][0] = '\0';
	}
	for (i=0;i<nAqSpecies;i++) aq_str[i][0] = '\0';

	int line_length = 300;
	char line[line_length]; // Individual line
	int eqphases = 0; // Switch to determine if the EQUILIBRIUM_PHASES block is being read
	int sel = 0;      // Switch to determine if the SELECTED_OUTPUT block is being read

	// Assemble file title
	sprintf(itime_str, "%d", itime);
	sprintf(temp_str, "%g", temp);
	sprintf(pressure_str, "%g", pressure);
	sprintf(pH_str, "%g", pH);
	sprintf(pe_str, "%g", pe);
	sprintf(mass_w_str, "%g", mass_w);

//	for (i=0;i<nAtmSpecies;i++) {
//		if (xgas[i] > 0.0)
//			sprintf(gas_str2[i], "%g", log(xgas[i]*pressure)/log(10.0)); // Convert from mixing ratio to -log partial pressure
//		else
//			sprintf(gas_str2[i], "%g", 0.0);
//		strcat(gas_str1[i], gas_str2[i]);
//		sprintf(gas_str2[i], "%g", xgas[i]);
//		strcat(gas_str1[i], "\t");
//		strcat(gas_str1[i], gas_str2[i]);
//	}
	for (i=0;i<nAtmSpecies;i++) {
		if (xgas[i] > 0.0)
			sprintf(gas_str2[i], "%g", log(xgas[i]*pressure)/log(10.0)); // Convert from mixing ratio to -log partial pressure
		else
			sprintf(gas_str2[i], "%g", 0.0);
	}

	strcpy(*tempinput,TemplateFile);
//	strcat(*tempinput,itime_str);
//	strcat(*tempinput,"T");
	strcat(*tempinput,temp_str);
//	strcat(*tempinput,"P");
//	strcat(*tempinput,pressure_str);
//	strcat(*tempinput,"pH");
//	strcat(*tempinput,pH_str);
//	strcat(*tempinput,"pe");
//	strcat(*tempinput,pe_str);
//	strcat(*tempinput,"xCO2");
//	strcat(*tempinput,gas_str2[0]);
//	strcat(*tempinput,"xCH4");
//	strcat(*tempinput,gas_str2[1]); // File title complete

	sprintf(aq_str[0],"%g mol/kgw", xaq[0]);  // mol per kg H2O of oxidized C
	sprintf(aq_str[1],"%g mol/kgw", xaq[1]);  // mol per kg H2O of reduced C
	sprintf(aq_str[3],"%g mol/kgw", xaq[3]);  // mol per kg H2O of N

	sprintf(steps_str1,"%g in ", kintime);    // Duration of kinetic simulation
	sprintf(steps_str2,"%d steps", kinsteps); // Number of time steps of kinetic simulation
	strcat(steps_str1, steps_str2);           // "kintime in kinsteps steps"

	fin = fopen (TemplateFile,"r");
	if (fin == NULL) printf("WritePHREEQCInput: Missing input file. Path: %s\n", TemplateFile);
	fout = fopen (*tempinput,"w");
	if (fout == NULL) printf("WritePHREEQCInput: Missing output file. Path: %s\n", *tempinput);

	while (fgets(line, line_length, fin)) {
		if (line[0] == 'E' && line[1] == 'Q' && line[2] == 'U' && line[3] == 'I') eqphases = 1; // EQUILIBRIUM_PHASES block
		if (line[0] == 'S' && line[1] == 'E' && line[2] == 'L' && line[3] == 'E') sel = 1;      // SELECTED_OUTPUT block
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
		else if (line[0] == '-' && line[1] == 's' && line[2] == 't' && line[3] == 'e')         fprintf(fout, "%s\n", ConCat("-steps\t",steps_str1));
		else if (eqphases && line[2] == 'C' && line[3] == 'O' && line[4] == '2' && line[5] == '(' && line[6] == 'g') fprintf(fout, "%s\n", ConCat("\tCO2(g) \t\t",gas_str1[0]));
		else if (eqphases && line[2] == 'C' && line[3] == 'H' && line[4] == '4' && line[5] == '(' && line[6] == 'g') fprintf(fout, "%s\n", ConCat("\tCH4(g) \t\t",gas_str1[1]));
		else if (eqphases && line[2] == 'O' && line[3] == '2' && line[4] == '(' && line[5] == 'g' && line[6] == ')') fprintf(fout, "%s\n", ConCat("\tO2(g) \t\t",gas_str1[2]));
		else if (eqphases && line[2] == 'N' && line[3] == '2' && line[4] == '(' && line[5] == 'g' && line[6] == ')') fprintf(fout, "%s\n", ConCat("\tN2(g) \t\t",gas_str1[3]));
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
	}
	for (i=0;i<nAqSpecies;i++) free(aq_str[i]);
	free(gas_str1);
	free(gas_str2);
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
	molmass = (xgas[0]*44.0 + xgas[1]*16.0 + xgas[2]*32.0 + xgas[3]*28.0)/(xgas[0]+xgas[1]+xgas[2]+xgas[3])*0.001;

	return molmass;
}
