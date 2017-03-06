/*
 ============================================================================
 Name        : ExoCcycleGeo.c
 Author      : Marc Neveu
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <IPhreeqc.h>                   // To use the external PHREEQC geochemical code

#define cmdline 0						// If execution from terminal as "./ExoCcycleGeo"

#define PI_greek 3.14159265358979323846
#define G 6.67e-11                      // Gravitational constant (SI)
#define rEarth 6378000.0                // Earth radius (m)
#define mEarth 6.0e24                   // Earth mass (kg)
#define km2m 1000.0                     // Convert from km to m
#define m2cm 100.0                      // Convert m to cm
#define Kelvin 273.15                   // Convert between Celsius and Kelvin

#define TsurfEarth 288                  // Mean equilibrium surface temperature on Earth (K)
#define RmeltEarth 3.8e-19              // Rate of melt generation on Earth (s-1) Kite et al. 2009 Fig. 15; http://dx.doi.org/10.1088/0004-637X/700/2/1732
#define deltaCvolcEarth 2.2e5           // Surface C flux from subaerial+submarine volcanic outgassing (mol C s-1) Donnadieu et al. 2006; http://dx.doi.org/10.1029/2006GC001278
#define deltaCcontwEarth 8.4543e-10*4.0*PI_greek*rEarth*rEarth // Surface C flux from continental weathering on Earth (mol C s-1)

#define Ra_c 1707.762                   // Critical Rayleigh number for convection, http://home.iitk.ac.in/~sghorai/NOTES/benard/node15.html, see also Koschmieder EL (1993) Benard cells and Taylor vortices, Cambridge U Press, p. 20.
#define Cp 914.0                        // Specific heat capacity of mantle rock (J kg-1 K-1)
#define kcrust 4.18                     // Thermal conductivity of crust (W m-1 K-1)
#define A0 70000.0                      // Activation temperature for thermal model (K)

#define xCO2g0 355.0e-6                 // Reference atmospheric CO2 mixing ratio (ppmv)
#define runoff_0 0.665                  // Reference runoff (mm/day) Edson et al. 2012, http://dx.doi.org/10.1089/ast.2011.0762
#define nAtmSpecies 4                   // Number of atmospheric species in the model

//-------------------------------------------------------------------
// WATER-ROCK PARAMETERS
//-------------------------------------------------------------------

#define nvar 1016                                          // Number of geochemical variables stored in each PHREEQC simulation
#define naq 258                                            // Number of aqueous species (+ physical parameters)
#define ngases 15                                          // Number of gaseous species
#define nmingas 385                                        // Number of minerals and gases
#define nelts 31                                           // 30 elements + 1 extra column in WaterRock/Molar_masses.txt

//-------------------------------------------------------------------
// SUBROUTINE DECLARATIONS
//-------------------------------------------------------------------

int OceanDiss (char path[1024], double T, double P, double *pH, double **xO2g, double **aq);
// int LoadMolMass (char path[1024], double ***molmass);
int EHandler(int phreeqc);
int ExtractWrite(int instance, double** data);
const char* ConCat(const char *str1, const char *str2);
int WritePHREEQCInput(const char *TemplateFile, double temp, double pressure, double pH, double pe, double *xgas, char *tempinput[1024]);
double **read_input (int H, int L, double **Input, char path[1024], const char filename[1024]);

//-------------------------------------------------------------------
// MAIN PROGRAM
//-------------------------------------------------------------------

int main(int argc, char *argv[]) {

	int i = 0;

	double r_p = 0.0;               // Planet radius (m)

	double deltaCvolc = 0.0;        // Surface C flux from volcanic outgassing (subaerial+submarine) (mol C m-2 s-1)
	double deltaCocean = 0.0;       // Surface C flux from surface ocean dissolution (mol C m-2 s-1)
	double deltaCcontw = 0.0;       // Surface C flux from continental weathering (mol C m-2 s-1)
	double deltaCseafw = 0.0;       // Surface C flux from seafloor weathering (mol C m-2 s-1)
	double deltaC = 0.0;            // Net surface C flux from all geological processes (mol C m-2 s-1)

	double pH = 8.22;               // pH of the surface ocean (default 8.22)

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

	//-------------------------------------------------------------------
	// Inputs
	//-------------------------------------------------------------------

	double m_p = 2.0*mEarth;    // Planet mass (kg)
	double L = 0.3;             // Fraction of planet surface covered by land

	// Atmospheric inputs
	double Tsurf = 288.0;       // Surface temperature (K)
	double Psurf = 1.0;         // Surface pressure (bar)
	double runoff = 0.7;        // Atmospheric runoff (mm/day)

	double *xgas = (double*) malloc(nAtmSpecies*sizeof(double));
	if (xgas == NULL) printf("WaterRock: Not enough memory to create xgas[nAtmSpecies]\n"); // Mixing ratios (partial pressures in bar) of atmospheric gases
    xgas[0] = 400.0e-6;  // CO2
    xgas[1] = 1.8e-6;    // CH4
    xgas[2] = 0.2;       // O2
    xgas[3] = 0.78;      // N2

	double *xaq = (double*) malloc(nAtmSpecies*sizeof(double));
	if (xaq == NULL) printf("WaterRock: Not enough memory to create xaq[nAtmSpecies]\n"); // Molalities of aqueous species (mol (kg H2O)-1)
    for (i=0;i<nAtmSpecies;i++) xaq[i] = 0.0;

	// Interior inputs
	double r_c = 1220000.0;     // Radius of planet core (m)
	double rhoMagma = 3500.0;   // Magma density (kg m-3)
	double rhoCrust = 3500.0;   // Crustal density (kg m-3)

	// TODO Quantities to be computed by thermal/geodynamic model
	double zCrust = 10.0*km2m;  // Crustal thickness (m)
	double Tmant = 3000.0;      // Mantle temperature (K)
	double Ra = 3000.0;         // Rayleigh number for mantle convection (no dim)

	//-------------------------------------------------------------------
	// Startup
	//-------------------------------------------------------------------

	printf("\n");
	printf("-------------------------------------------------------------------\n");
	printf("ExoCcycleGeo v17.2\n");
	printf("This code is in development and cannot be used for science yet.\n");
	if (cmdline == 1) printf("Command line mode\n");
	printf("-------------------------------------------------------------------\n");

	printf("\n");
	printf("Inputs:\n");
	printf("Planet mass \t \t \t \t %g M_Earth\n", m_p/mEarth);
	printf("Land coverage of planet surface \t %g%%\n", L*100.0);
	printf("Surface temperature \t \t \t %g K\n", Tsurf);
	printf("Atmospheric CO2 molar mixing ratio \t %g\n", xgas[0]);
	printf("Atmospheric CH4 molar mixing ratio \t %g\n", xgas[1]);
	printf("Atmospheric O2 molar mixing ratio \t %g\n", xgas[2]);
	printf("Atmospheric N2 molar mixing ratio \t %g\n", xgas[3]);
	printf("Ocean pH \t \t \t \t %g \n", pH);
	printf("Surface runoff \t \t \t \t %g mm d-1\n", runoff);
	printf("-------------------------------------------------------------------\n");

	printf("\n");
	printf("Computing geo C fluxes...\n");

	// Get current directory. Works for Mac only! To switch between platforms, see:
	// http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe
	char path[1024];
	unsigned int size = sizeof(path);
	path[0] = '\0';
	if (_NSGetExecutablePath(path, &size) == 0) printf("\n");
	else printf("IcyDwarf: Couldn't retrieve executable directory. Buffer too small; need size %u\n", size);

	r_p = pow(m_p/mEarth,0.27)*rEarth;  // Planet radius (m)
	gsurf = G*m_p/r_p/r_p;
	Pf = rhoCrust*gsurf*zCrust;
	TbaseLid = Tmant - 2.23*Tmant*Tmant/A0;
	Nu = pow(Ra/Ra_c,0.25);

	//-------------------------------------------------------------------
	// Calculate surface C flux from volcanic outgassing (Kite et al. 2009) in mol C m-2 s-1
	//-------------------------------------------------------------------

	// Assumes that all melt generated reaches the surface
	vConv = 2.0*(Nu-1.0) * (kcrust/(rhoMagma*Cp*(r_p-r_c))) * (Tmant-Tsurf) / (Tmant-TbaseLid);
	Rmelt = vConv*rhoMagma/m_p * (rhoCrust*gsurf*zCrust/(Pf-P0)); // m-2 s-1
	deltaCvolc = deltaCvolcEarth * Rmelt/RmeltEarth; // (mol C s-1) * (m-2 s-1) / (s-1) = mol C m-2 s-1

	//-------------------------------------------------------------------
	// Calculate surface C flux from continental weathering (Edson et al. 2012) in mol C s-1 // TODO scale by surface
	//-------------------------------------------------------------------

	deltaCcontw = -L * 0.5*deltaCcontwEarth*r_p*r_p/rEarth/rEarth * pow(xgas[0]/xCO2g0,0.3) * runoff/runoff_0 * exp((Tsurf-TsurfEarth)/17.7);

	//-------------------------------------------------------------------
	// Calculate surface C flux from surface ocean dissolution in top m  TODO scale by surface, include kinetics, manage reservoir size
	//-------------------------------------------------------------------

	printf("xCO2(g) = %g ppm, xCH4(g) = %g ppm, xO2(g) = %g, xN2(g) = %g\n", xgas[0]/1.0e-6, xgas[1]/1.0e-6, xgas[2], xgas[3]);
	deltaCocean = xgas[0];
	printf("Calculating ocean dissolution...\n");
	OceanDiss(path, Tsurf, Psurf, &pH, &xgas, &xaq);
	printf("xCO2(g) = %g ppm, xCH4(g) = %g ppm, xO2(g) = %g, xN2(g) = %g\n", xgas[0]/1.0e-6, xgas[1]/1.0e-6, xgas[2], xgas[3]);
	printf("pH = %g\n", pH);
	printf("\nx(aq) \t mol/(kg H2O)\n");
	printf("----------------------\n");
	printf("C +IV \t %g\n", xaq[0]);
	printf("C -IV \t %g\n", xaq[1]);
	printf("O2 \t %g\n", xaq[2]);
	printf("N2 \t %g\n", xaq[3]);
	deltaCocean = -(deltaCocean + xgas[0])/Psurf;

	//-------------------------------------------------------------------
	// Calculate surface C flux from seafloor weathering TODO scale by surface, include kinetics, manage reservoir size
	//-------------------------------------------------------------------

	deltaCreac = 0.0; // TODO call PHREEQC to get deltaCreac, net mol C leached/precipitated per kg of rock
	tcirc = 1.0; // TODO compute tcirc based on Nu(Ra(basal heat flux))
	deltaCseafw = -(1.0-L) * 4.0/3.0*PI_greek*(pow(r_p,3)-pow(r_p-zCrack,3))/tcirc*deltaCreac*rhoCrust;

	//-------------------------------------------------------------------
	// Compute net geo C flux
	//-------------------------------------------------------------------

	deltaC = deltaCvolc + deltaCcontw + deltaCseafw;

	printf("\n");
	printf("Net surface C flux from... \t mol C cm-2 s-1\n");
	printf("-------------------------------------------\n");
	printf("volcanic outgassing \t \t %g\n", deltaCvolc);
	printf("continental weathering \t \t %g\n", deltaCcontw);
	printf("seafloor weathering \t \t %g\n", deltaCseafw);
	printf("-------------------------------------------\n");
	printf("all geo processes \t \t %g\n", deltaC);

	printf("\nExiting ExoCcycleGeo...\n");

	free (xgas);
	free (xaq);

	return EXIT_SUCCESS;
}

int OceanDiss (char path[1024], double T, double P, double *pH, double **xgas, double **xaq) {

	int phreeqc = 0;
	int i = 0;
	int j = 0;

	double pe = 0.0;                                             // pe corresponding to logfO2
	double logfO2 = log((*xgas)[2])/log(10.0);                   // O2 fugacity
	double logKO2H2O = 83.10494;                                 // log K for reaction 4 H+ + 4 e- + O2 = 2 H2O, , from CHNOSZ: subcrt(c("H+","e-","O2","H2O"),c(-4,-4,-1,2),c("aq","aq","g","liq"),T=25,P=1)
	double mass_water = 0.0;                                     // Mass of water at the end of the PHREEQC simulation
	double total_mol_gas = 0.0;								     // Moles of gas at the end of the PHREEQC simulation

	char *dbase = (char*)malloc(1024);                           // Path to thermodynamic database
	char *infile = (char*)malloc(1024);                          // Path to initial (template) input file
	char *tempinput = (char*)malloc(1024);                       // Temporary PHREEQC input file, modified from template

	double *simdata = (double*) malloc(nvar*sizeof(double));
	if (simdata == NULL) printf("WaterRock: Not enough memory to create simdata[nvar]\n");

	double **molmass = (double**) malloc(nvar*sizeof(double*));  // Rearranged data from Molar_masses.txt
	if (molmass == NULL) printf("ParamExploration_plot: Not enough memory to create molmass[nvar]\n");
	for (i=0;i<nvar;i++) {
		molmass[i] = (double*) malloc(nelts*sizeof(double));
		if (molmass[i] == NULL) printf("ParamExploration_plot: Not enough memory to create molmass[nvar][nelts]\n");
	}
	for (i=0;i<nvar;i++) {
		for (j=0;j<nelts;j++) {
			molmass[i][j] = 0.0;
		}
	}

	// Initializations
	dbase[0] = '\0';
	infile[0] = '\0';
	tempinput[0] = '\0';

	if (cmdline == 1) strncat(dbase,path,strlen(path)-20);
	else strncat(dbase,path,strlen(path)-18);
	strcat(dbase,"PHREEQC-3.1.2/core9.dat");

	strncat(infile,dbase,strlen(dbase)-9);
	strcat(infile,"io/PHREEQCOceanDissInput");

//	LoadMolMass (path, &molmass);

// Use CHNOSZ to get log fO2 for fayalite-magnetite-quartz (FMQ) buffer at given T and P
//	logfO2 = -3.0*CHNOSZ_logK("quartz", "cr", T, P, "SUPCRT92")
//		     -2.0*CHNOSZ_logK("magnetite", "cr", T, P, "SUPCRT92")
//	         +3.0*CHNOSZ_logK("fayalite", "cr", T, P, "SUPCRT92")
//		     +1.0*CHNOSZ_logK("O2", "g", T, P, "SUPCRT92");
//	logKO2H2O = -4.0*CHNOSZ_logK("H+", "aq", T, P, "SUPCRT92")
//				-4.0*CHNOSZ_logK("e-", "aq", T, P, "SUPCRT92")
//				-1.0*CHNOSZ_logK("O2", "g", T, P, "SUPCRT92")
//				+2.0*CHNOSZ_logK("H2O", "liq", T, P, "SUPCRT92");

	pe = -(*pH) + 0.25*(logfO2+logKO2H2O);

	WritePHREEQCInput(infile, T-Kelvin, P, *pH, pe, *xgas, &tempinput);

	phreeqc = CreateIPhreeqc(); // Run PHREEQC
	if (LoadDatabase(phreeqc,dbase) != 0) OutputErrorString(phreeqc);
	SetSelectedOutputFileOn(phreeqc,1);
	if (RunFile(phreeqc,tempinput) != 0) OutputErrorString(phreeqc);

	ExtractWrite(phreeqc, &simdata);

	if (DestroyIPhreeqc(phreeqc) != IPQ_OK) OutputErrorString(phreeqc);

	mass_water = simdata[11];
	if (simdata[999] > 0.0) total_mol_gas = simdata[999];
	else {
		printf("Gas dissolved completely\n");
		total_mol_gas = 1.0; // To avoid dividing by zero later
	}

	(*xaq)[0] = simdata[45];             // C(4), i.e. dissolved CO2 and carbonate
	(*xaq)[1] = simdata[43];             // C(-4), i.e. dissolved methane
	(*xaq)[2] = simdata[72];             // O(0), i.e. dissolved O2
	(*xaq)[3] = simdata[69];             // N(0), i.e. dissolved N2

	(*xgas)[0] = simdata[1006]/total_mol_gas;          // CO2(g)
	(*xgas)[1] = simdata[1004]/total_mol_gas;          // CH4(g)
	(*xgas)[2] = simdata[1014]/total_mol_gas;          // O2(g)
	(*xgas)[3] = simdata[1010]/total_mol_gas;          // N2(g)

	(*pH) = simdata[7];

	free (simdata);
	for (i=0;i<nvar;i++) free (molmass[i]);
	free (molmass);

	return 0;
}

//int LoadMolMass (char path[1024], double ***molmass) {
//
//	int i = 0;
//	int j = 0;
//	int k = 0;
//
//	double **molmass_read = (double**) malloc(nmingas*sizeof(double*));  // Data from Molar_masses.txt
//	if (molmass_read == NULL) printf("ParamExploration_plot: Not enough memory to create molmass_read[nmingas]\n");
//	for (i=0;i<nmingas;i++) {
//		molmass_read[i] = (double*) malloc(nelts*sizeof(double));
//		if (molmass_read[i] == NULL) printf("ParamExploration_plot: Not enough memory to create molmass_read[nmingas][nelts]\n");
//	}
//	for (i=0;i<nmingas;i++) {
//		for (j=0;j<nelts;j++) {
//			molmass_read[i][j] = 0.0;
//		}
//	}
//
//	// Read molmass database
//	molmass_read = read_input(nelts, nmingas, molmass_read, path, "Data/Molar_masses.txt");
//
//	// Shift to positions corresponding to simdata
//	// Gas species
//	for (i=0;i<ngases;i++) {
//		for (j=0;j<nelts;j++) {
//			(*molmass)[naq+2*(nmingas-ngases)+5-1+i][j] = molmass_read[nmingas-ngases+i][j];
//		}
//	}
//	// Solid species
//	k = naq-1;
//	for (i=0;i<nmingas-ngases;i++) {
//		for (j=0;j<nelts;j++) {
//			(*molmass)[k][j] = molmass_read[i][j];
//			(*molmass)[k+1][j] = (*molmass)[k][j];
//		}
//		k = k+2;
//	}
//	// First line with molar masses of elements
//	for (j=0;j<nelts;j++) (*molmass)[0][j] = molmass_read[0][j];
//
//	for (i=0;i<nmingas;i++) free (molmass_read[i]);
//	free (molmass_read);
//
//	return 0;
//}

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
 * Modifies P, T, pH, pe, W:R
 *
 *--------------------------------------------------------------------*/
int WritePHREEQCInput(const char *TemplateFile, double temp, double pressure, double pH, double pe, double *xgas, char **tempinput) {

	// Open input file
	FILE *fin;
	FILE *fout;
	char temp_str[10];
	char pressure_str[10];
	char pH_str[10];
	char pe_str[10];
	char CO2_str[10];
	char CH4_str[10];
	char O2_str[10];
	char N2_str[10];
	temp_str[0] = '\0';
	pressure_str[0] = '\0';
	pH_str[0] = '\0';
	pe_str[0] = '\0';
	CO2_str[0] = '\0';
	CH4_str[0] = '\0';
	int line_length = 300;
	char line[line_length]; // Individual line
	int line_no = 0; // Line number

	sprintf(temp_str,"%g",temp);
	sprintf(pressure_str,"%g",pressure);
	sprintf(pH_str,"%g",pH);
	sprintf(pe_str,"%g",pe);
	sprintf(CO2_str,"%g",xgas[0]/pressure); // Convert from mixing ratio to partial pressure
	sprintf(CH4_str,"%g",xgas[1]/pressure); // Convert from mixing ratio to partial pressure
	sprintf(O2_str,"%g",xgas[2]/pressure); // Convert from mixing ratio to partial pressure
	sprintf(N2_str,"%g",xgas[3]/pressure); // Convert from mixing ratio to partial pressure

	strcpy(*tempinput,TemplateFile);
	strcat(*tempinput,"T");
	strcat(*tempinput,temp_str);
	strcat(*tempinput,"P");
	strcat(*tempinput,pressure_str);
	strcat(*tempinput,"pH");
	strcat(*tempinput,pH_str);
	strcat(*tempinput,"pe");
	strcat(*tempinput,pe_str);
	strcat(*tempinput,"xCO2");
	strcat(*tempinput,CO2_str);
	strcat(*tempinput,"xCH4");
	strcat(*tempinput,CH4_str);

	fin = fopen (TemplateFile,"r");
	if (fin == NULL) printf("OceanDiss: Missing input file. Path: %s\n", TemplateFile);
	fout = fopen (*tempinput,"w");
	if (fout == NULL) printf("OceanDiss: Missing output file. Path: %s\n", *tempinput);

	while (fgets(line, line_length, fin)) {
		line_no++;
		if (line_no == 5) {
			fprintf(fout, "%s\n", ConCat("\tpH \t \t",pH_str));
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
		else if (line[1] == '-' && line[2] == 'p' && line[3] == 'r' && line[4] == 'e' && line[5] == 's') {
			fprintf(fout, "%s\n", ConCat("\t-pressure \t",pressure_str));
		}
		else if (line[1] == '-' && line[2] == 't' && line[3] == 'e' && line[4] == 'm' && line[5] == 'p') {
			fprintf(fout, "%s\n", ConCat("\t-temperature \t",temp_str));
		}
		else if (line[1] == 'C' && line[2] == 'O' && line[3] == '2' && line[4] == '(' && line[5] == 'g') {
			fprintf(fout, "%s\n", ConCat("\tCO2(g) \t\t",CO2_str));
		}
		else if (line[1] == 'C' && line[2] == 'H' && line[3] == '4' && line[4] == '(' && line[5] == 'g') {
			fprintf(fout, "%s\n", ConCat("\tCH4(g) \t\t",CH4_str));
		}
		else if (line[1] == 'O' && line[2] == '2' && line[3] == '(' && line[4] == 'g' && line[5] == ')') {
			fprintf(fout, "%s\n", ConCat("\tO2(g) \t\t",O2_str));
		}
		else if (line[1] == 'N' && line[2] == '2' && line[3] == '(' && line[4] == 'g' && line[5] == ')') {
			fprintf(fout, "%s\n", ConCat("\tN2(g) \t\t",N2_str));
		}
		else fputs(line,fout);
	}
	if (ferror(fin)) {
		printf("ParamExploration: Error reading template input file %s\n",TemplateFile);
		return 1;
	}

	fclose(fin);
	fclose(fout);

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
