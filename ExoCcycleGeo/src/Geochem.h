/*
 * Geochem.h
 *
 *  Created on: Jun 14, 2019
 *      Author: Marc Neveu (marc.f.neveu@nasa.gov)
 */

#include "ExoCcycleGeo.h"

#ifndef GEOCHEM_H_
#define GEOCHEM_H_

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
int alphaMELTS (char *path, int nPTstart, int nPTend, char *aMELTS_setfile, double ***sys_tbl);

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

	if      (strcmp(filename, "io/OceanStart.txt") == 0) ExtractWrite(phreeqc, &simdata, 1, nvarEq); // 1st line of selected output has initial solution composition = equilibrium when there are no equilibrium phases
	else if (strcmp(filename, "io/OceanDiss.txt") == 0) ExtractWrite(phreeqc, &simdata, 2, nvarEq);  // 2nd line of PHREEQC selected output solution and mineral+gas composition at equilibrium
	else if (strcmp(filename, "io/ContWeather.txt") == 0) ExtractWrite(phreeqc, &simdata, kinsteps, nvarKin);
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
	if (strcmp(filename, "io/OceanStart.txt") == 0) {
		(*pH) = simdata[7][0];
		(*pe) = simdata[8][0];

		(*xaq)[0] = simdata[45][0];             // C(4), i.e. dissolved CO2 and carbonate
		(*xaq)[1] = simdata[43][0];             // C(-4), i.e. dissolved methane
//		(*xaq)[2] = simdata[72][0];             // O(0), i.e. dissolved O2
		(*xaq)[3] = simdata[28][0];             // Total dissolved N
		(*xaq)[4] = simdata[69][0];             // N(0), i.e. dissolved N2
		(*xaq)[5] = simdata[68][0];             // N(-3), i.e. dissolved NH3 and NH4+
	}

	// Ocean dissolution
	if (strcmp(filename, "io/OceanDiss.txt") == 0) {
		(*pH) = simdata[7][0];
		(*pe) = simdata[8][0];
		(*nAir) = simdata[1007][0]*nAir0;
		(*P) = simdata[1006][0]*atm2bar;
		(*V) = simdata[1008][0]/1000.0*nAir0;
		(*mass_w) = simdata[11][0]*nAir0;

		(*xaq)[0] = simdata[45][0];             // C(4), i.e. dissolved CO2 and carbonate
		(*xaq)[1] = simdata[43][0];             // C(-4), i.e. dissolved methane
//		(*xaq)[2] = simdata[72][0];             // O(0), i.e. dissolved O2
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
	if (strcmp(filename, "io/ContWeather.txt") == 0) {

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
	strcat(*tempinput,"Exec.txt"); // File title complete

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

	char *temp = (char*)malloc(1024*sizeof(char));
	char *cmd = (char*)malloc(1024*sizeof(char));

	temp[0] = '\0';
	cmd[0] = '\0';

	// Store path
	strcat(temp,path);

	// Do cleanup (rm command)
	strcat(cmd,"rm ");
	if (cmdline == 1) strncat(cmd,path,strlen(path)-20);
	else strncat(cmd,path,strlen(path)-18);
	strcat(cmd,"sel*");
	system(cmd);

	// Restore path
	path[0] = '\0';
	strcat(path,temp);

	free(cmd);
	free(temp);

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
 * Subroutine alphaMELTS
 *
 * Compute the amount of melting in the mantle and crust using the
 * alphaMELTS implementation of the MELTS model
 * (https://magmasource.caltech.edu/alphamelts). See:
 *
 * Smith & Asimow (2005) doi: 10.1029/2004GC000816
 * Ghiorso & Sack (1995) Cont Min Petr 119, 197-212
 * Asimow & Ghiorso (1998) Am Min 83, 1127-1132
 *
 *--------------------------------------------------------------------*/

int alphaMELTS (char *path, int nPTstart, int nPTend, char *aMELTS_setfile, double ***sys_tbl) {

	int i = 0;
	int j = 0;
	int scan = 0;
	int lineno = 0;
	char tmp = '\0';
	FILE *f;

	// --- Build system command ---
	// /.../ExoCcycleGeo/alphaMELTS-1.9/run_alphameltsExoC.command -f /.../ExoCcycleGeo/alphaMELTS-1.9/ExoC/ExoC_env.txt -b /.../ExoCcycleGeo/alphaMELTS-1.9/ExoC/ExoCbatch.txt -p /.../ExoCcycleGeo/alphaMELTS-1.9/
	char *aMELTSsys = (char*)malloc(65536); // System command
	aMELTSsys[0] = '\0';
	char *aMELTStmp = (char*)malloc(1024);  // Temporary path
	aMELTStmp[0] = '\0';

	// Executable
	if (cmdline == 1) strncat(aMELTStmp,path,strlen(path)-20);
	else strncat(aMELTStmp,path,strlen(path)-18);
	strcat(aMELTStmp,"alphaMELTS-1.9/run_alphameltsExoC.command -f ");
	strcat(aMELTSsys, aMELTStmp);

	// Settings file
	aMELTStmp[0] = '\0';
	if (cmdline == 1) strncat(aMELTStmp,path,strlen(path)-20);
	else strncat(aMELTStmp,path,strlen(path)-18);
	strcat(aMELTStmp,"alphaMELTS-1.9/");
	strcat(aMELTStmp,aMELTS_setfile);
	strcat(aMELTStmp," -b ");
	strcat(aMELTSsys, aMELTStmp);

	// Batch file
	aMELTStmp[0] = '\0';
	if (cmdline == 1) strncat(aMELTStmp,path,strlen(path)-20);
	else strncat(aMELTStmp,path,strlen(path)-18);
	strcat(aMELTStmp,"alphaMELTS-1.9/ExoC/ExoCbatch.txt -p "); // Uses ExoCcycleGeo.melts
	strcat(aMELTSsys, aMELTStmp);

	// Output path
	aMELTStmp[0] = '\0';
	if (cmdline == 1) strncat(aMELTStmp,path,strlen(path)-20);
	else strncat(aMELTStmp,path,strlen(path)-18);
	strcat(aMELTStmp,"alphaMELTS-1.9/");
	strcat(aMELTSsys, aMELTStmp);

	// --- Run alphaMELTS ---
	system(aMELTSsys);
	// Alternative if ever needed: echo [interactive inputs] |.
	// system("echo \"1 /Users/mneveu/eclipse-workspace/ExoCcycleGeo/ExoCcycleGeo/alphaMELTS-1.9/ExoC/ExoCcycleGeo.melts 4 1 0\" | /Users/mneveu/eclipse-workspace/ExoCcycleGeo/ExoCcycleGeo/alphaMELTS-1.9/run_alphameltsExoC.command -f /Users/mneveu/eclipse-workspace/ExoCcycleGeo/ExoCcycleGeo/alphaMELTS-1.9/ExoC/Mantle_env.txt");

	// --- Extract melt fraction vs. P from System_main_tbl_geotherm.txt ---

	char *f_sys_tbl = (char*)malloc(1024); // Path to PT file input into alphaMELTS
	f_sys_tbl[0] = '\0';

	if (cmdline == 1) strncat(f_sys_tbl,path,strlen(path)-20);
	else strncat(f_sys_tbl,path,strlen(path)-18);
	strcat(f_sys_tbl,"alphaMELTS-1.9/System_main_tbl.txt");

	f = fopen (f_sys_tbl,"r");
	if (f == NULL) printf("ExoCcycleGeo: Missing f_sys_tbl file path: %s\n", f_sys_tbl);

	// Skip first 4 lines, which are text
	while (lineno < 4) {
		scan = fscanf(f, "%c", &tmp);
		if (!scan) printf("Error scanning alphaMELTS System_main_tbl.txt at line number %d\n",lineno);
		if (tmp == '\n') lineno++;
	}
	// Store the rest, won't fill full table if alphaMELTS run crashed but that's OK
	for(i=nPTstart;i<nPTend;i++) {
		for(j=0;j<18;j++) {
			scan = fscanf(f, "%lg", &(*sys_tbl)[i][j]);
			if (!scan) printf("Error scanning alphaMELTS System_main_tbl.txt at line number %d, column number %d\n", lineno+i+1,j+1);
		}
	}

	fclose(f);
	free(f_sys_tbl);

	free(aMELTSsys);
	free(aMELTStmp);

	return 0;
}

#endif /* GEOCHEM_H_ */
