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
#define nAqSpecies 4                    // Number of aqueous species whose concentrations are passed between the physical and chemical models

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

int OceanDiss (char path[1024], int itime, double T, double P, double *pH, double *pe, double mass_w, double **xO2g, double **aq);
// int LoadMolMass (char path[1024], double ***molmass);
int EHandler(int phreeqc);
int ExtractWrite(int instance, double** data);
const char* ConCat(const char *str1, const char *str2);
int WritePHREEQCInput(const char *TemplateFile, int itime, double temp, double pressure, double pH, double pe, double mass_w, double *xgas, double *xaq, char *tempinput[1024]);
double **read_input (int H, int L, double **Input, char path[1024], const char filename[1024]);

//-------------------------------------------------------------------
// MAIN PROGRAM
//-------------------------------------------------------------------

int main(int argc, char *argv[]) {

	int i = 0;
	int itime = 0;                  // Time step counter
	int ntime = 0;                  // Total number of steps

	double dtime = 1.0*Myr2sec;     // Time step (s)

	double r_p = 0.0;               // Planet radius (m)

	double RCplate = 0.0;           // Plate/crust C reservoir (mol)
	double RCmantle = 0.0;          // Mantle C reservoir (mol)
	double RCatm = 0.0;             // Atmospheric C reservoir (mol)
	double RCocean = 0.0;           // Ocean C reservoir (mol)
	double RCatmoc = RCatm + RCocean; // Combined atmospheric and ocean C reservoir (mol)

	double nCatm = 0.0;             // Atmospheric C reservoir (mol), = RCatm but for scaling
	double nCocean = 0.0;           // Ocean C reservoir (mol), = RCocean but for scaling

	double FCoutgas = 0.0;          // C flux from outgassing (subaerial+submarine) (mol s-1)
	double FCcontw = 0.0;           // C flux from continental weathering (mol s-1)
	double FCseafw = 0.0;           // C flux from seafloor weathering (mol s-1)
	double FCsubd = 0.0;            // C flux from subduction (mol s-1)
	double netFC = 0.0;             // Net surface C flux from all geological processes (mol s-1)

	double fatm = 0.0;              // Fraction of C in {atm+ocean} reservoir that is in the atmosphere (as opposed to in the ocean)
	double farc = 0.0;              // Fraction of subducted C that makes it back out through arc volcanism (as opposed to into the mantle)

	double pH = 8.22;               // pH of the surface ocean (default 8.22)
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
	double molmass_atm = 0.0;		// Mean molecular mass of the atmosphere (kg mol-1)
	double nAir = 0.0;              // Number of mol in atmosphere (mol)

	//-------------------------------------------------------------------
	// Inputs
	//-------------------------------------------------------------------

	ntime = (int) (5000.0*Myr2sec/dtime); // Run for 5 Gyr

	double m_p = 1.0*mEarth;    // Planet mass (kg)
	double L = 0.15;            // Fraction of planet surface covered by land
	double Mocean = 1.4e21;     // Mass of ocean (kg, default: Earth=1.4e21)

	// Atmospheric inputs
	double Tsurf = 288.0;       // Surface temperature (K)
	double Psurf = 2.0;         // Surface pressure (bar)
	double runoff = 0.7;        // Atmospheric runoff (mm/day)

	double *xgas = (double*) malloc(nAtmSpecies*sizeof(double));
	if (xgas == NULL) printf("WaterRock: Not enough memory to create xgas[nAtmSpecies]\n"); // Mixing ratios by volume (or by mol since all gases are pretty much ideal and have the same molar volume) of atmospheric gases
    xgas[0] = 0.22;      // 400.0e-6;  // CO2
    xgas[1] = 0.0;       // 1.8e-6;    // CH4
    xgas[2] = 0.0;       // 0.2;       // O2
    xgas[3] = 0.78;                    // N2
    for (i=4;i<nAtmSpecies;i++) xgas[i] = 0.0;

	double *xaq = (double*) malloc(nAtmSpecies*sizeof(double));
	if (xaq == NULL) printf("WaterRock: Not enough memory to create xaq[nAtmSpecies]\n"); // Molalities of aqueous species (mol (kg H2O)-1)
    for (i=0;i<nAtmSpecies;i++) xaq[i] = 0.0;
    xaq[0] = 27.0/12.0/1000.0;  // 27 ppm C in today's oceans (scaled from 141 ppm Alk*M(C)/M(HCO3), M being molecular mass)

	// Interior inputs
	double r_c = 1220000.0;     // Radius of planet core (m)
	double rhoMagma = 3500.0;   // Magma density (kg m-3)
	double rhoCrust = 3500.0;   // Crustal density (kg m-3)

	// Quantities to be computed by thermal/geodynamic model
	double zCrust = 0.0;        // Crustal thickness (m)
	double Tmantle = 0.0;       // Mantle temperature (K)
	double Ra = 0.0;            // Rayleigh number for mantle convection (no dim)

	//-------------------------------------------------------------------
	// Startup
	//-------------------------------------------------------------------

	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("ExoCcycleGeo v18.9\n");
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

	r_p = pow(m_p/mEarth,0.27)*rEarth; // Planet radius, determined from mass scaling (m) // TODO implement IcyDwarf compression model for more accurate mass-radius relationship
	gsurf = G*m_p/r_p/r_p;
	Asurf = 4.0*PI_greek*r_p*r_p;

	// Use CHNOSZ to get log fO2 for fayalite-magnetite-quartz (FMQ) buffer at given T and P.
	logfO2 = -3.0*CHNOSZ_logK("quartz", "cr", Tsurf, Psurf, "SUPCRT92")
		     -2.0*CHNOSZ_logK("magnetite", "cr", Tsurf, Psurf, "SUPCRT92")
	         +3.0*CHNOSZ_logK("fayalite", "cr", Tsurf, Psurf, "SUPCRT92")
		     +1.0*CHNOSZ_logK("O2", "g", Tsurf, Psurf, "SUPCRT92");
	logKO2H2O = -4.0*CHNOSZ_logK("H+", "aq", Tsurf, Psurf, "SUPCRT92")
				-4.0*CHNOSZ_logK("e-", "aq", Tsurf, Psurf, "SUPCRT92")
				-1.0*CHNOSZ_logK("O2", "g", Tsurf, Psurf, "SUPCRT92")
				+2.0*CHNOSZ_logK("H2O", "liq", Tsurf, Psurf, "SUPCRT92");

	pe = -pH + 0.25*(logfO2+logKO2H2O);

	//-------------------------------------------------------------------
	// Initialize reservoirs
	//-------------------------------------------------------------------

	RCmantle = 0.0;

	molmass_atm = (xgas[0]*44.0 + xgas[1]*16.0 + xgas[2]*32.0 + xgas[3]*28.0)/(xgas[0]+xgas[1]+xgas[2]+xgas[3])*0.001; // TODO make a function
	nAir = Psurf*bar2Pa*Asurf/gsurf/molmass_atm;
	RCatm = (xgas[0]+xgas[1])*nAir;
	RCocean = 27.0e-6*Mocean; // Scaled from alkalinity of about 140*M(C)/M(HCO3) where M is molar mass
	RCatmoc = RCatm + RCocean;
	RCmantle = 10.0*(RCatmoc);

	// Print first line of output
	FILE *fout;
	char *title = (char*)malloc(1024*sizeof(char)); // Don't forget to free!
	title[0] = '\0';
	if (cmdline == 1) strncat(title,path,strlen(path)-20);
	else strncat(title,path,strlen(path)-18);
	strcat(title,"Output.txt");
	fout = fopen(title,"w");
	if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
	else {
		fprintf(fout, "'Reservoirs in mol, fluxes in mol s-1'\n");
		fprintf(fout, "'Time (Myr)' \t RCmantle \t RCatm \t RCocean \t FCoutgas \t FCcontw \t FCseafw \t FCsubd \t 'Net C flux'\n");
		fprintf(fout, "Init \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n", RCmantle, RCatm, RCocean, FCoutgas, FCcontw, FCseafw, FCsubd, netFC);
	}
	fclose (fout);
	free (title);

	//-------------------------------------------------------------------
	// Start time loop
	//-------------------------------------------------------------------

	for (itime = 0;itime<ntime;itime++) {

		double sumPP = 0.0;
		for(i=0;i<nAtmSpecies;i++) sumPP = sumPP + xgas[i]*Psurf;
		Psurf = sumPP;

		//-------------------------------------------------------------------
		// Calculate partitioning of C between ocean and atmosphere
		//-------------------------------------------------------------------

		/* Assumes that mixing gas and liquid at the surface is as inefficient as on Earth:
		 * - same wind agitating the liquid surface;
		 * - same mixing of the upper layers of the ocean with the next lower layers (takes many years on Earth, the time scale of
		 *   energy transfer in the upper layers);
		 * - as on Earth, movement of surface layers to the deep ocean takes centuries.
		 */

//		printf("xCO2(g) = %g ppm, xCH4(g) = %g ppm, xO2(g) = %g, xN2(g) = %g\n", xgas[0]/1.0e-6, xgas[1]/1.0e-6, xgas[2], xgas[3]);
//		printf("Calculating ocean dissolution...\n");
		OceanDiss(path, itime, Tsurf, Psurf, &pH, &pe, Mocean/nAir, &xgas, &xaq);
//		printf("xCO2(g) = %g ppm, xCH4(g) = %g ppm, xO2(g) = %g, xN2(g) = %g\n", xgas[0]/1.0e-6, xgas[1]/1.0e-6, xgas[2], xgas[3]);
//		printf("pH = %g\n", pH);
//		printf("\nx(aq) \t mol/(kg H2O)\n");
//		printf("----------------------\n");
//		printf("C +IV \t %g\n", xaq[0]);
//		printf("C -IV \t %g\n", xaq[1]);
//		printf("O2 \t %g\n", xaq[2]);
//		printf("N2 \t %g\n", xaq[3]);

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
		FCoutgas = deltaCvolcEarth * Rmelt/RmeltEarth;              // (mol C s-1) * (m-2 s-1) / (s-1) = mol C m-2 s-1

		//-------------------------------------------------------------------
		// Calculate surface C flux from continental weathering
		//-------------------------------------------------------------------

		FCcontw = -L * 0.5*deltaCcontwEarth * pow(xgas[0]/xCO2g0,0.3) * runoff/runoff_0 * exp((Tsurf-TsurfEarth)/17.7); // Edson et al. (2012) Eq. 1; Abbot et al. (2012) Eq. 2

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
//		RCatmoc = RCatmoc + dtime*(farc*FCsubd + FCoutgas - FCcontw/2.0 - FCseafw); // Sum of atmospheric and ocean reservoirs, still needs partitioning

		// Split RCatmoc between RCatm and RCocean
		molmass_atm = (xgas[0]*44.0 + xgas[1]*16.0 + xgas[2]*32.0 + xgas[3]*28.0)/(xgas[0]+xgas[1]+xgas[2]+xgas[3])*0.001;
		nCatm = xgas[0]*Psurf*bar2Pa*Asurf/gsurf/molmass_atm;
		nCocean = xaq[0]*Mocean;
		fatm = nCatm / (nCocean + nCatm); // RCatm / Rocean+RCatm. P*A/g=mass of atmosphere

		RCatm = fatm * RCatmoc;
		RCocean = RCatmoc - RCatm;

		netFC = FCoutgas - FCcontw - FCseafw + farc*FCsubd - (1.0-farc)*FCsubd;

		// Write output
		char *title = (char*)malloc(1024*sizeof(char)); // Don't forget to free!
		title[0] = '\0';
		if (cmdline == 1) strncat(title,path,strlen(path)-20);
		else strncat(title,path,strlen(path)-18);
		strcat(title,"Output.txt");
		fout = fopen(title,"a");
		if (fout == NULL) printf("ExoCcycleGeo: Error opening %s output file.\n",title);
		else fprintf(fout, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n", (double)itime*dtime/Myr2sec, RCmantle, RCatm, RCocean, FCoutgas, FCcontw, FCseafw, FCsubd, netFC);
		fclose (fout);
		free (title);

		if (itime == 2) exit(0);
	} // End time loop

	printf("\nExiting ExoCcycleGeo...\n");

	free (xgas);
	free (xaq);

	Rf_endEmbeddedR(0);                                     // Close R and CHNOSZ

	return EXIT_SUCCESS;
}

/*--------------------------------------------------------------------
 *
 * Subroutine OceanDiss
 *
 * Calculates C flux from ocean dissolution // TODO include ocean mass
 *
 *--------------------------------------------------------------------*/

int OceanDiss (char path[1024], int itime, double T, double P, double *pH, double *pe, double mass_w, double **xgas, double **xaq) {

	int phreeqc = 0;
	int i = 0;

	char *dbase = (char*)malloc(1024);                           // Path to thermodynamic database
	char *infile = (char*)malloc(1024);                          // Path to initial (template) input file
	char *tempinput = (char*)malloc(1024);                       // Temporary PHREEQC input file, modified from template

	double *simdata = (double*) malloc(nvar*sizeof(double));
	if (simdata == NULL) printf("WaterRock: Not enough memory to create simdata[nvar]\n");

	// Initializations
	dbase[0] = '\0';
	infile[0] = '\0';
	tempinput[0] = '\0';

	if (cmdline == 1) strncat(dbase,path,strlen(path)-20);
	else strncat(dbase,path,strlen(path)-18);
	strcat(dbase,"PHREEQC-3.1.2/core10.dat");

	strncat(infile,dbase,strlen(dbase)-10);
	strcat(infile,"io/OceanDissInput");

//	LoadMolMass (path, &molmass);

	printf("xCO2(aq) \t xCH4(aq) \t xgas\n");
	printf("%g %g\t", (*xaq)[0]*12.0*1000.0, (*xaq)[1]*12.0*1000.0);
	printf("%g\n", (*xgas)[0]);

	WritePHREEQCInput(infile, itime, T-Kelvin, P, *pH, *pe, mass_w, *xgas, *xaq, &tempinput);

	phreeqc = CreateIPhreeqc(); // Run PHREEQC
	if (LoadDatabase(phreeqc,dbase) != 0) OutputErrorString(phreeqc);
	SetSelectedOutputFileOn(phreeqc,1);
	if (RunFile(phreeqc,tempinput) != 0) OutputErrorString(phreeqc);

	ExtractWrite(phreeqc, &simdata);

	if (DestroyIPhreeqc(phreeqc) != IPQ_OK) OutputErrorString(phreeqc);

//	for (i=0;i<nvar;i++) printf("%d %g\n", i, simdata[i]);

	(*xaq)[0] = simdata[45];             // C(4), i.e. dissolved CO2 and carbonate
	(*xaq)[1] = simdata[43];             // C(-4), i.e. dissolved methane
	(*xaq)[2] = simdata[72];             // O(0), i.e. dissolved O2
	(*xaq)[3] = simdata[69];             // N(0), i.e. dissolved N2

	(*xgas)[0] = simdata[1016];          // CO2(g)
	(*xgas)[1] = simdata[1012];          // CH4(g)
	(*xgas)[2] = simdata[1032];          // O2(g)
	(*xgas)[3] = simdata[1024];          // N2(g)

	(*pH) = simdata[7];
	(*pe) = simdata[8];

	free (simdata);

	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine EHandler
 *
 * Error handler from Charlton & Parkhurst (2011), Computer &
 * Geosciences 37, 1653-1663.
 *
 *--------------------------------------------------------------------*/
int EHandler(int instance) {
	OutputErrorString(instance);
	// exit(EXIT_FAILURE);
	return 0;
}

/*--------------------------------------------------------------------
 *
 * Subroutine ExtractWrite
 *
 * Write selected output from PHREEQC
 *
 *--------------------------------------------------------------------*/
int ExtractWrite(int instance, double** data) {
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
		GetSelectedOutputValue(instance,2,i,&v);
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
int WritePHREEQCInput(const char *TemplateFile, int itime, double temp, double pressure, double pH, double pe, double mass_w, double *xgas, double *xaq, char **tempinput) {

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
	char **gas_str1 = (char**)malloc(nAtmSpecies*sizeof(char*));
	for (i=0;i<nAtmSpecies;i++) gas_str1[i] = (char*)malloc(1024);

	char **gas_str2 = (char**)malloc(nAtmSpecies*sizeof(char*));
	for (i=0;i<nAtmSpecies;i++) gas_str2[i] = (char*)malloc(1024);

	char *aqC_str = (char*) malloc(1024);

	itime_str[0] = '\0';
	temp_str[0] = '\0';
	pressure_str[0] = '\0';
	pH_str[0] = '\0';
	pe_str[0] = '\0';
	mass_w_str[0] = '\0';
	int line_length = 300;
	char line[line_length]; // Individual line
	int line_no = 0;  // Line number
	int eqphases = 0; // Switch to determine if the EQUILIBRIUM_PHASES block is being read

	// Assemble file title
	sprintf(itime_str, "%d", itime);
	sprintf(temp_str, "%g", temp);
	sprintf(pressure_str, "%g", pressure);
	sprintf(pH_str, "%g", pH);
	sprintf(pe_str, "%g", pe);
	sprintf(mass_w_str, "%g", mass_w);
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

	strcpy(*tempinput,TemplateFile);
	strcat(*tempinput,itime_str);
	strcat(*tempinput,"T");
	strcat(*tempinput,temp_str);
	strcat(*tempinput,"P");
	strcat(*tempinput,pressure_str);
	strcat(*tempinput,"pH");
	strcat(*tempinput,pH_str);
	strcat(*tempinput,"pe");
	strcat(*tempinput,pe_str);
	strcat(*tempinput,"xCO2");
	strcat(*tempinput,gas_str2[0]);
	strcat(*tempinput,"xCH4");
	strcat(*tempinput,gas_str2[1]); // File title complete

	sprintf(aqC_str,"%g mol/kgw", xaq[0]+xaq[1]); // mol per kg of oxidized and reduced C

	fin = fopen (TemplateFile,"r");
	if (fin == NULL) printf("OceanDiss: Missing input file. Path: %s\n", TemplateFile);
	fout = fopen (*tempinput,"w");
	if (fout == NULL) printf("OceanDiss: Missing output file. Path: %s\n", *tempinput);

	while (fgets(line, line_length, fin)) {
		line_no++;
		if (line_no == 5) {
			fprintf(fout, "%s charge\n", ConCat("\tpH \t \t",pH_str));
		}
		else if (line_no == 6) {
			fprintf(fout, "%s\n", ConCat("\ttemp \t \t",temp_str));
		}
		else if (line_no == 7) {
			fprintf(fout, "%s\n", ConCat("\tpressure \t",pressure_str));
		}
		else if (line_no == 8) {
			fprintf(fout, "%s\n", ConCat("\tpe \t \t",pe_str));
		}
		else if (line[1] == '-' && line[2] == 'w' && line[3] == 'a' && line[4] == 't' && line[5] == 'e') {
			fprintf(fout, "%s\n", ConCat("\t-water \t \t",mass_w_str));
		}
		else if (line[1] == '-' && line[2] == 'p' && line[3] == 'r' && line[4] == 'e' && line[5] == 's') {
			fprintf(fout, "%s\n", ConCat("\t-pressure \t",pressure_str));
		}
		else if (line[1] == '-' && line[2] == 't' && line[3] == 'e' && line[4] == 'm' && line[5] == 'p') {
			fprintf(fout, "%s\n", ConCat("\t-temperature \t",temp_str));
		}
		else if (eqphases == 1 && line[2] == 'C' && line[3] == 'O' && line[4] == '2' && line[5] == '(' && line[6] == 'g') {
			fprintf(fout, "%s\n", ConCat("\tCO2(g) \t\t",gas_str1[0]));
		}
		else if (eqphases == 1 && line[2] == 'C' && line[3] == 'H' && line[4] == '4' && line[5] == '(' && line[6] == 'g') {
			fprintf(fout, "%s\n", ConCat("\tCH4(g) \t\t",gas_str1[1]));
		}
		else if (eqphases == 1 && line[2] == 'O' && line[3] == '2' && line[4] == '(' && line[5] == 'g' && line[6] == ')') {
			fprintf(fout, "%s\n", ConCat("\tO2(g) \t\t",gas_str1[2]));
		}
		else if (eqphases == 1 && line[2] == 'N' && line[3] == '2' && line[4] == '(' && line[5] == 'g' && line[6] == ')') {
			fprintf(fout, "%s\n", ConCat("\tN2(g) \t\t",gas_str1[3]));
		}
		else if (!eqphases && line[1] == 'C' && line[2] == '\t') {
			fprintf(fout, "%s\n", ConCat("\tC \t\t",aqC_str));
		}
		else fputs(line,fout);
		if (line[0] == 'E' && line[1] == 'Q' && line[2] == 'U' && line[3] == 'I') eqphases = 1;
	}
	if (ferror(fin)) {
		printf("ParamExploration: Error reading template input file %s\n",TemplateFile);
		return 1;
	}

	fclose(fin);
	fclose(fout);

	for (i=0;i<nAtmSpecies;i++) {
		free(gas_str1[i]);
		free(gas_str2[i]);
	}
	free(gas_str1);
	free(gas_str2);
	free(aqC_str);

	return 0;
}

double **read_input (int H, int L, double **Input, char path[1024], const char filename[1024]) {

	FILE *fin;
	int l = 0;
	int h = 0;

	// Turn working directory into full file path by moving up two directories
	// to IcyDwarf (i.e., removing "Release/IcyDwarf" characters) and specifying
	// the right path end.

	char *title = (char*)malloc(1024);       // Don't forget to free!
	title[0] = '\0';
	if (cmdline == 1) strncat(title,path,strlen(path)-20);
	else strncat(title,path,strlen(path)-18);
	strcat(title,filename);

	fin = fopen (title,"r");
	if (fin == NULL) {
		printf("IcyDwarf: Error opening %s input file.\n",title);
	}
	else {
		for (l=0;l<L;l++) {
			for (h=0;h<H;h++) {
				int scan = fscanf(fin,"%lg",&Input[l][h]);
				if (scan != 1)
					printf("IcyDwarf: Error scanning %s file at l=%d, h=%d.\n",title,l,h);
			}
		}
	}

	fclose (fin);
	free (title);

	return Input;
}
