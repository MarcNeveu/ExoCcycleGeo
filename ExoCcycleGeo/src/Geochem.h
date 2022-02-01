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
		double *mass_w, double **xgas, double **xaq, double ***xriver, int iResTime, int forcedPP, double kintime, int kinsteps, int nvar);
int ExtractWrite(int instance, double*** data, int line, int nvar);
const char* ConCat (const char *str1, const char *str2);
int WritePHREEQCInput(const char *TemplateFile, int itime, double temp, double pressure, double gasvol, double pH, double pe, double mass_w,
		double *xgas, double *xaq, double *xriver, int forcedPP, double kintime, int kinsteps, char **tempinput, int readRockVol);
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
		double *mass_w, double **xgas, double **xaq, double ***xriver, int iResTime, int forcedPP, double kintime, int kinsteps, int nvar) {

	int phreeqc = 0;
	int i = 0;
	int j = 0;
	int readRockVol = 0;                                         // Switch to read rock volume from TITLE row of PHREEQC input file in WritePHREEQCInput

	char *dbase = (char*)malloc(1024);                           // Path to thermodynamic database
	char *infile = (char*)malloc(1024);                          // Path to initial (template) input file
	char *tempinput = (char*)malloc(1024);                       // Temporary PHREEQC input file, modified from template
	double nAir0 = *nAir; // Scaling factor for gas volume and water mass, so PHREEQC doesn't have to handle large numbers

	double **simdata = (double**) malloc(kinsteps*sizeof(double*));
	if (simdata == NULL) printf("AqueousChem: Not enough memory to create simdata[kinsteps][kinsteps]\n");
	for (i=0;i<kinsteps;i++) {
		simdata[i] = (double*) malloc(nvar*sizeof(double));
		if (simdata[i] == NULL) printf("AqueousChem: Not enough memory to create simdata[kinsteps][nvar]\n");
	}

	// Initializations
	dbase[0] = '\0';
	infile[0] = '\0';
	tempinput[0] = '\0';
	for (i=0;i<kinsteps;i++) {
		for (j=0;j<nvar;j++) simdata[i][j] = 0.0;
	}

	if (cmdline == 1) strncat(dbase,path,strlen(path)-20);
	else strncat(dbase,path,strlen(path)-18);
	strcat(dbase,"PHREEQC-3.1.2/core11_idealgas.dat");

	strncat(infile,dbase,strlen(dbase)-19);
	strcat(infile, filename);

	if (kintime || strcmp(filename, "io/MixRiverOcean.txt") == 0  || strcmp(filename, "io/SeafWeather.txt") == 0) *mass_w *= nAir0; // Don't scale water mass for weathering, which does not require equilibrating ocean and atmosphere i.e. handling large numbers
	if (strcmp(filename, "io/SeafWeather.txt") == 0) readRockVol = 1;

	WritePHREEQCInput(infile, itime, T-Kelvin, *P, *V/nAir0, *pH, *pe, *mass_w/nAir0, *xgas, *xaq, (*xriver)[iResTime], forcedPP, kintime, kinsteps, &tempinput, readRockVol);

	phreeqc = CreateIPhreeqc(); // Run PHREEQC
	if (LoadDatabase(phreeqc,dbase) != 0) OutputErrorString(phreeqc);
	SetSelectedOutputFileOn(phreeqc,1);
	printf("Running PHREEQC\n");
	if (RunFile(phreeqc,tempinput) != 0) OutputErrorString(phreeqc);
	printf("PHREEQC ran successfully\n");

	if      (strcmp(filename, "io/OceanStart.txt") == 0) ExtractWrite(phreeqc, &simdata, 1, nvarEq); // 1st line of selected output has initial solution composition = equilibrium when there are no equilibrium phases
	else if (strcmp(filename, "io/OceanDiss.txt") == 0) ExtractWrite(phreeqc, &simdata, 2, nvarEq);  // 2nd line of PHREEQC selected output: solution and mineral+gas composition at equilibrium
	else if (strcmp(filename, "io/ContWeather.txt") == 0) ExtractWrite(phreeqc, &simdata, kinsteps, nvarKin);
	else if (strcmp(filename, "io/MixRiverOcean.txt") == 0) ExtractWrite(phreeqc, &simdata, 3, nvarEq); // 3rd line of PHREEQC selected output: mixed solution
	else if (strcmp(filename, "io/SeafWeather.txt") == 0) ExtractWrite(phreeqc, &simdata, 2, nvarEq);
	else printf("AqueousChem: Cannot extract PHREEQC output data because input path %s is inaccurate.\n", filename);

	if (DestroyIPhreeqc(phreeqc) != IPQ_OK) OutputErrorString(phreeqc);

	// Setting initial ocean chemistry
	if (strcmp(filename, "io/OceanStart.txt") == 0) {
		(*pH) = simdata[0][1];
		(*pe) = simdata[0][2];

		(*xaq)[0] = simdata[0][39];             // C(4), i.e. dissolved CO2 and carbonate
		(*xaq)[1] = simdata[0][37];             // C(-4), i.e. dissolved methane
		(*xaq)[2] = simdata[0][66];             // O(0), i.e. dissolved O2
		(*xaq)[3] = simdata[0][84];             // Ntg
//		(*xaq)[4] = simdata[0][22];             // Total dissolved N
//		(*xaq)[5] = simdata[0][62];             // N(-3), i.e. dissolved NH3 and NH4+
	}

	// Ocean dissolution
	if (strcmp(filename, "io/OceanDiss.txt") == 0) {
		(*pH) = simdata[0][1];
		(*pe) = simdata[0][2];
		(*nAir) = simdata[0][1005]*nAir0;
		(*P) = simdata[0][1004]*atm2bar;
		(*V) = simdata[0][1006]/1000.0*nAir0;
		(*mass_w) = simdata[0][5]*nAir0;
		(*xaq)[0] = simdata[0][39];             // C(4), i.e. dissolved CO2 and carbonate
		(*xaq)[1] = simdata[0][37];             // C(-4), i.e. dissolved methane
		(*xaq)[2] = simdata[0][66];             // O(0), i.e. dissolved O2
		(*xaq)[3] = simdata[0][84];             // Ntg
//		(*xaq)[4] = simdata[0][22];             // N excluding Ntg
//		(*xaq)[5] = simdata[0][62];             // N(-3), i.e. dissolved NH3 and NH4+
		(*xaq)[6] = simdata[0][19]; 			// Mg
		(*xaq)[7] = simdata[0][9]; 				// Ca
		(*xaq)[8] = simdata[0][15]; 			// Fe
		(*xaq)[9] = simdata[0][28]; 			// Si
		(*xaq)[10]= simdata[0][23]; 			// Na
		(*xaq)[11]= simdata[0][26]; 			// S
		(*xaq)[12]= simdata[0][10]; 			// Cl

		(*xgas)[0] = simdata[0][1012];          // CO2(g)
		(*xgas)[1] = simdata[0][1010];          // CH4(g)
		(*xgas)[2] = simdata[0][1020];          // O2(g)
		(*xgas)[3] = simdata[0][1024];          // Ntg(g)
		(*xgas)[4] = simdata[0][1014];          // H2O(g)

		for (i=0;i<nAtmSpecies;i++) (*xgas)[i] = (*xgas)[i]/simdata[0][1005]; // Divide by total mol gas to return mixing ratio
	}

	// Continental weathering
	if (strcmp(filename, "io/ContWeather.txt") == 0) {
		for (i=0;i<kinsteps;i++) {
			for (j=0;j<nvarKin-8;j++) (*xriver)[i][j] = simdata[i][j];
			if (i > 0) {
				(*xriver)[i][nvarKin-8] = (*xriver)[i-1][nvarKin-8] + (*xriver)[i][107]; // Cumulative calcite consumed (mol)
				(*xriver)[i][nvarKin-7] = (*xriver)[i-1][nvarKin-7] + (*xriver)[i][109]; // Cumulative dolomite-disordered consumed (mol)
				(*xriver)[i][nvarKin-6] = (*xriver)[i-1][nvarKin-6] + (*xriver)[i][111]; // Cumulative dolomite-ordered consumed (mol)
				(*xriver)[i][nvarKin-5] = (*xriver)[i-1][nvarKin-5] + (*xriver)[i][113]; // Cumulative magnesite consumed (mol)
				(*xriver)[i][nvarKin-4] = (*xriver)[i-1][nvarKin-4] + (*xriver)[i][115]; // Cumulative anhydrite consumed (mol)
				(*xriver)[i][nvarKin-3] = (*xriver)[i-1][nvarKin-3] + (*xriver)[i][117]; // Cumulative gypsum consumed (mol)
				(*xriver)[i][nvarKin-2] = (*xriver)[i-1][nvarKin-2] + (*xriver)[i][119]; // Cumulative pyrite consumed (mol)
				(*xriver)[i][nvarKin-1] = (*xriver)[i-1][nvarKin-1] + (*xriver)[i][121]; // Cumulative pyrrhotite consumed (mol)
			}
		}
	}

	// Mixing river and ocean. Ocean mass is not updated as it is assumed constant.
	if (strcmp(filename, "io/MixRiverOcean.txt") == 0) {
		(*pH) = simdata[0][1]; // 1st index is 2, i.e., 3rd row of selected output file for mixing calculation
		(*pe) = simdata[0][2];
		(*xaq)[0] = simdata[0][23]; // C(4), i.e. dissolved CO2 and carbonate
		(*xaq)[1] = simdata[0][21]; // C(-4), i.e. dissolved methane
//		(*xaq)[2] = simdata[0][36]; // O(0), i.e. dissolved O2
		(*xaq)[3] = simdata[0][46]; // Ntg
//		(*xaq)[4] = 0.0;                                              // N excluding Ntg
//		(*xaq)[5] = 0.0;                                              // N(-3), i.e. dissolved NH3 and NH4+
		(*xaq)[6] = simdata[0][12]; // Mg
		(*xaq)[7] = simdata[0][8] ; // Ca
		(*xaq)[8] = simdata[0][10]; // Fe
		(*xaq)[9] = simdata[0][17]; // Si
		(*xaq)[10]= simdata[0][14]; // Na
		(*xaq)[11]= simdata[0][16]; // S
		(*xaq)[12]= simdata[0][14]; // Cl
	}

	// Seafloor weathering
	if (strcmp(filename, "io/SeafWeather.txt") == 0) {
		(*pH) = simdata[0][1];
		(*pe) = simdata[0][2];
		// Scale by (original water mass)/(new water mass) under assumption that hydration and dehydration (but not carbonation/decarbonation) of ocean crust are balanced
		(*xaq)[0] = simdata[0][23] * (*mass_w)/nAir0 / simdata[0][5]; // C(4), i.e. dissolved CO2 and carbonate
		(*xaq)[1] = simdata[0][21] * (*mass_w)/nAir0 / simdata[0][5]; // C(-4), i.e. dissolved methane
//		(*xaq)[2] = simdata[0][36] * (*mass_w)/nAir0 / simdata[0][5]; // O(0), i.e. dissolved O2
		(*xaq)[3] = simdata[0][46] * (*mass_w)/nAir0 / simdata[0][5]; // Ntg
//		(*xaq)[4] = 0.0;                                              // N excluding Ntg
//		(*xaq)[5] = 0.0;                                              // N(-3), i.e. dissolved NH3 and NH4+
		(*xaq)[6] = simdata[0][12] * (*mass_w)/nAir0 / simdata[0][5]; // Mg
		(*xaq)[7] = simdata[0][8]  * (*mass_w)/nAir0 / simdata[0][5]; // Ca
		(*xaq)[8] = simdata[0][10] * (*mass_w)/nAir0 / simdata[0][5]; // Fe
		(*xaq)[9] = simdata[0][17] * (*mass_w)/nAir0 / simdata[0][5]; // Si
		(*xaq)[10]= simdata[0][14] * (*mass_w)/nAir0 / simdata[0][5]; // Na
		(*xaq)[11]= simdata[0][16] * (*mass_w)/nAir0 / simdata[0][5]; // S
		(*xaq)[12]= simdata[0][14] * (*mass_w)/nAir0 / simdata[0][5]; // Cl
		// Alternatively, track moles of carbonate precipitated
//		(*xaq)[215] = simdata[0][215];          // Aragonite CaCO3
	}

	for (i=0;i<kinsteps;i++) free (simdata[i]);
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
		for (i=0;i<nvar;i++) {                           // Rest of parameters
			GetSelectedOutputValue(instance,line,i,&v);
			if (fabs(v.dVal) < 1e-50) (*data)[0][i] = 0.0;
			else (*data)[0][i] = v.dVal;
		}
	}
	else if (nvar == nvarKin) { // Kinetic simulation
		for (i=0;i<nvar;i++) {
			for (j=0;j<line;j++) {
				GetSelectedOutputValue(instance,j,i,&v);
				if (fabs(v.dVal) < 1e-50) (*data)[j][i] = 0.0;
				else (*data)[j][i] = v.dVal;
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
		double *xgas, double *xaq, double *xriver, int forcedPP, double kintime, int kinsteps, char **tempinput, int readRockVol) {

	int i = 0;
	int solution1 = 0; // Switch: is the SOLUTION block being read?
	int solution2 = 0; // Switch: is the SOLUTION 2 block being read? (for mixing calculation)
	int eqphases = 0; // Switch: is the EQUILIBRIUM_PHASES block being read?
	int gasphase = 0; // Switch: is the GAS_PHASE block being read?
	int mix = 0;	  // Switch: is the MIX block being read?
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
	char molMassCrust[10];
	char rockVol[10];

	char **gas_str1 = (char**)malloc(nAtmSpecies*sizeof(char*));
	for (i=0;i<nAtmSpecies;i++) gas_str1[i] = (char*)malloc(1024);

	char **gas_str2 = (char**)malloc(nAtmSpecies*sizeof(char*));
	for (i=0;i<nAtmSpecies;i++) gas_str2[i] = (char*)malloc(1024);

	char **gas_str3 = (char**)malloc(nAtmSpecies*sizeof(char*));
	for (i=0;i<nAtmSpecies;i++) gas_str3[i] = (char*)malloc(1024);

	char **aq_str = (char**) malloc(nAqSpecies*sizeof(char*));
	for (i=0;i<nAqSpecies;i++) aq_str[i] = (char*)malloc(1024);

	char (**river_str) = (char**) malloc(nvarKin*sizeof(char*));
	for (i=0;i<nvarKin;i++) river_str[i] = (char*)malloc(1024);

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
	for (i=0;i<nAqSpecies;i++) {
		aq_str[i][0] = '\0';
	}
	for(i=0;i<nvarKin;i++) {
		river_str[i][0] = '\0';
	}

	int line_length = 300;
	char line[line_length]; // Individual line

	// Convert to PHREEQC input units
	pressure = pressure/atm2bar; // Convert from bar to atm
	gasvol = gasvol*1000.0;      // Convert from m3 to L

	// Kinetic simulation applies only to continental weathering, for which the pH of rain is used, not that of the ocean
	if (kintime) pH = 7.0;

	sprintf(itime_str, "%d", itime);
	sprintf(temp_str, "%g", temp);
	sprintf(pressure_str, "%g", pressure);
	sprintf(vol_str, "%g", gasvol);
	sprintf(pH_str, "%g", pH);
	sprintf(pe_str, "%g", pe);
	sprintf(mass_w_str, "%g", mass_w);

	// Used if gases are in EQUILIBRIUM_PHASES block, which is not the case here but could be for SeafWeather
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

	// Used in OceanStart, ContWeather, and MixRiverOcean where all or some PP are forced
	for (i=0;i<nAtmSpecies;i++) {
		if (xgas[i] > 0.0)
			sprintf(gas_str2[i], "%g", log(xgas[i]*pressure)/log(10.0)); // Convert from mixing ratio to -log partial pressure
		else
			sprintf(gas_str2[i], "%g", 0.0);
	}

	// Use in OceanDiss for GAS_PHASE block
	for (i=0;i<nAtmSpecies;i++) sprintf(gas_str3[i], "%g", xgas[i]*pressure); // Partial pressure

//	if (forcedPP) {
//		for (i=0;i<nAqSpecies;i++) sprintf(aq_str[i],"1.0"); // Doesn't matter what xaq is since PHREEQC will set molalities from partial pressures
//	}
//	else {
		for (i=0;i<nAqSpecies;i++) {
			sprintf(aq_str[i], "%g mol/kgw", xaq[i]);
		}
		for (i=0;i<nvarKin;i++) {
			sprintf(river_str[i], "%g mol/kgw", xriver[i]);
		}
//	}

	strcpy(*tempinput,TemplateFile);
//	strcat(*tempinput,itime_str);
	strcat(*tempinput,"Exec.txt"); // File title complete

	sprintf(steps_str1,"%g in ", kintime);    // Duration of kinetic simulation
	sprintf(steps_str2,"%d steps", kinsteps); // Number of time steps of kinetic simulation
	strcat(steps_str1, steps_str2);           // "kintime in kinsteps steps"

	fin = fopen (TemplateFile,"r"); 	// Open input file
	if (fin == NULL) printf("WritePHREEQCInput: Missing input file. Path: %s\n", TemplateFile);
	fout = fopen (*tempinput,"w");      // Open output file
	if (fout == NULL) printf("WritePHREEQCInput: Missing output file. Path: %s\n", *tempinput);

	while (fgets(line, line_length, fin)) {

		// Block switches
		if (kintime) { // Continental weathering simulation only
			if (line[0] == 'T' && line[1] == 'I' && line[2] == 'T' && line[3] == 'L') { // TITLE line with crust molar mass value
				for (i=0;i<10;i++) molMassCrust[i] = line[i+64];
				mass_w *= strtod((const char*)molMassCrust, NULL)/1000.0; // Multiply W:R by crust molar mass to get mass of water in kg
				sprintf(mass_w_str, "%g", mass_w);
			}
		}
		if (readRockVol) { // Read mass and density from TITLE line
			if (line[0] == 'T' && line[1] == 'I' && line[2] == 'T' && line[3] == 'L') { // TITLE line with crust molar mass value
				for (i=0;i<10;i++) rockVol[i] = line[i+40];
				mass_w *= strtod((const char*)rockVol, NULL); // Multiply W:R by input rock volume to get mass of water in kg
				sprintf(mass_w_str, "%g", mass_w);
			}
		}
		if (line[0] == 'S' && line[1] == 'O' && line[2] == 'L' && line[9] == '1') solution1 = 1; // SOLUTION block
		if (line[0] == 'S' && line[1] == 'O' && line[2] == 'L' && line[9] == '2') { // SOLUTION 2 block
			solution2 = 1;
			solution1 = 0;
		}
		if (line[0] == 'E' && line[1] == 'Q' && line[2] == 'U' && line[3] == 'I') { // EQUILIBRIUM_PHASES block
			eqphases = 1;
			solution1 = 0;
		}
		if (line[0] == 'G' && line[1] == 'A' && line[2] == 'S' && line[3] == '_') gasphase = 1; // GAS_PHASE block
		if (line[0] == 'M' && line[1] == 'I' && line[2] == 'X') mix = 1;                        // MIX block
		if (line[0] == 'S' && line[1] == 'E' && line[2] == 'L' && line[3] == 'E') sel = 1;      // SELECTED_OUTPUT block

		// SOLUTION
		if (line[1] == 'p' && line[2] == 'H') {
			if      (solution1) fprintf(fout, "%s charge\n", ConCat("\tpH \t \t",pH_str));
			else if (solution2) fprintf(fout, "%s charge\n", ConCat("\tpH \t \t",river_str[3]));
		}
		else if (line[1] == 't' && line[2] == 'e' && line[3] == 'm' && line[4] == 'p') {
			if      (solution1) fprintf(fout, "%s\n", ConCat("\ttemp \t \t",temp_str));
			else if (solution2) fprintf(fout, "%s\n", ConCat("\ttemp \t \t",river_str[5]));
		}
		else if (line[1] == 'p' && line[2] == 'r' && line[3] == 'e' && line[4] == 's')         fprintf(fout, "%s\n", ConCat("\tpressure \t",pressure_str));
		else if (line[1] == 'p' && line[2] == 'e') {
			if      (solution1) fprintf(fout, "%s\n", ConCat("\tpe \t \t",pe_str));
			else if (solution2) fprintf(fout, "%s\n", ConCat("\tpe \t \t",river_str[4]));
		}
		// Common water mass even for mixing calculation because proportions of solutions mixed are set in MIX block. These proportions could be set here instead, the outcome is the same.
		else if (!sel && line[1] == '-' && line[2] == 'w' && line[3] == 'a' && line[4] == 't') fprintf(fout, "%s\n", ConCat("\t-water \t \t",mass_w_str));
		else if (!eqphases && line[1] == 'C' && line[2] == '(' && line[3] == '4' && line[4] == ')') {
			if (forcedPP && xgas[0] > 0.0) {
				strcat(aq_str[0], "\tCO2(g) \t");
				strcat(aq_str[0], gas_str2[0]);
				fprintf(fout, "%s\n", ConCat("\tC(4) \t\t", aq_str[0]));
			}
			else {
				if      (solution1) fprintf(fout, "%s\n", ConCat("\tC(4) \t\t",aq_str[0]));
				else if (solution2) fprintf(fout, "%s\n", ConCat("\tC(4) \t\t",river_str[9])); // TODO this is C, get C(4) instead
			}
		}
		else if (!eqphases && line[1] == 'C' && line[2] == '(' && line[3] == '-' && line[4] == '4') {
			if (forcedPP && xgas[1] > 0.0) {
				strcat(aq_str[1], "\tCH4(g) \t");
				strcat(aq_str[1], gas_str2[1]);
				fprintf(fout, "%s\n", ConCat("\tC(-4) \t\t", aq_str[1]));
			}
			else {
				if (solution1) fprintf(fout, "%s\n", ConCat("\tC(-4) \t\t",aq_str[1])); // TODO also set C(-4) from river for river-ocean mixing calculation
			}
		}
		else if (!eqphases && line[1] == 'O' && line[2] == '('  && line[3] == '0' && line[4] == ')') {
			if ((forcedPP && xgas[2] > 0.0) || solution2) {
				strcat(aq_str[2], "\tO2(g) \t");
				strcat(aq_str[2], gas_str2[2]); // Doesn't matter if it's river_str (river-ocean mixing calculation) or aq_str (ocean start, continental weathering) because the value should be overwritten by gas_str at PHREEQC execution
				fprintf(fout, "%s\n", ConCat("\tO(0) \t\t", aq_str[2]));
			}
			else fprintf(fout, "%s\n", ConCat("\tO(0) \t\t", aq_str[2]));
		}
		else if (!eqphases && line[1] == 'N' && line[2] == 't'  && line[3] == 'g' && line[4] == '\t') {
			if ((forcedPP && xgas[3] > 0.0) || solution2) {
				strcat(aq_str[3], "\tNtg(g) \t");
				strcat(aq_str[3], gas_str2[3]); // Doesn't matter if it's river_str (river-ocean mixing calculation) or aq_str (ocean start, continental weathering) because the value should be overwritten by gas_str at PHREEQC execution
				fprintf(fout, "%s\n", ConCat("\tNtg \t\t", aq_str[3]));
			}
			else fprintf(fout, "%s\n", ConCat("\tNtg \t\t", aq_str[3]));
		}

		// Next set of rows for ocean dissolution, mixing, and seafloor weathering, i.e., !forcedPP
		else if (!forcedPP && line[1] == 'M' && line[2] == 'g') {
			if      (solution1) fprintf(fout, "%s\n", ConCat("\tMg \t\t", aq_str[6]));
			else if (solution2) fprintf(fout, "%s\n", ConCat("\tMg \t\t", river_str[11]));
		}
		else if (!forcedPP && line[1] == 'C' && line[2] == 'a') {
			if      (solution1) fprintf(fout, "%s\n", ConCat("\tCa \t\t", aq_str[7]));
			else if (solution2) fprintf(fout, "%s\n", ConCat("\tCa \t\t", river_str[13]));
		}
		else if (!forcedPP && line[1] == 'F' && line[2] == 'e') {
			if      (solution1) fprintf(fout, "%s\n", ConCat("\tFe \t\t", aq_str[8]));
			else if (solution2) fprintf(fout, "%s\n", ConCat("\tFe \t\t", river_str[14]));
		}
		else if (!forcedPP && line[1] == 'S' && line[2] == 'i') {
			if      (solution1) fprintf(fout, "%s\n", ConCat("\tSi \t\t", aq_str[9]));
			else if (solution2) fprintf(fout, "%s\n", ConCat("\tSi \t\t", river_str[12]));
		}
		else if (!forcedPP && line[1] == 'N' && line[2] == 'a') {
			if      (solution1) fprintf(fout, "%s\n", ConCat("\tNa \t\t", aq_str[10]));
			else if (solution2) fprintf(fout, "%s\n", ConCat("\tNa \t\t", river_str[10]));
		}
		else if (!forcedPP && line[1] == 'S' && line[2] == '(') {
			if      (solution1) fprintf(fout, "%s\n", ConCat("\tS(6) \t\t", aq_str[11]));
			else if (solution2) fprintf(fout, "%s\n", ConCat("\tS(6) \t\t", river_str[10]));
		}
		else if (!forcedPP && line[1] == 'C' && line[2] == 'l') {
			if      (solution1) fprintf(fout, "%s\n", ConCat("\tCl \t\t", aq_str[12]));
			else if (solution2) fprintf(fout, "%s\n", ConCat("\tCl \t\t", river_str[16]));
		}

		// MIX -- could remove this line and set proportions of solutions to be mixed with specific water masses instead.
		else if (mix && line[1] == '1') fprintf(fout, "\t1 \t %s\n", mass_w_str); // Need to reprint the "1" otherwise it is overwritten

		// KINETICS
		else if (line[0] == '-' && line[1] == 's' && line[2] == 't' && line[3] == 'e')         fprintf(fout, "%s\n", ConCat("-steps\t",steps_str1));

		// GAS_PHASE
		else if (gasphase && line[1] == 'v' && line[2] == 'o' && line[3] == 'l') fprintf(fout, "%s\n", ConCat("\tvolume \t\t",vol_str));
		else if (gasphase && line[1] == 'C' && line[2] == 'O' && line[3] == '2' && line[4] == '(') fprintf(fout, "%s\n", ConCat("\tCO2(g) \t\t",gas_str3[0]));
		else if (gasphase && line[1] == 'C' && line[2] == 'H' && line[3] == '4' && line[4] == '(') fprintf(fout, "%s\n", ConCat("\tCH4(g) \t\t",gas_str3[1]));
		else if (gasphase && line[1] == 'O' && line[2] == '2' && line[3] == '(' && line[4] == 'g') fprintf(fout, "%s\n", ConCat("\tO2(g) \t\t",gas_str3[2]));
		else if (gasphase && line[1] == 'N' && line[2] == 't' && line[3] == 'g' && line[4] == '(') fprintf(fout, "%s\n", ConCat("\tNtg(g) \t\t",gas_str3[3]));
		else if (gasphase && line[1] == 'H' && line[2] == '2' && line[3] == 'O' && line[4] == '(') fprintf(fout, "%s\n", ConCat("\tH2O(g) \t\t",gas_str3[4]));

		// EQUILIBRIUM_PHASES
		else if (eqphases && line[2] == 'C' && line[3] == 'O' && line[4] == '2' && line[5] == '(' && line[6] == 'g') fprintf(fout, "%s\n", ConCat("\tCO2(g) \t\t",gas_str1[0]));
		else if (eqphases && line[2] == 'C' && line[3] == 'H' && line[4] == '4' && line[5] == '(' && line[6] == 'g') fprintf(fout, "%s\n", ConCat("\tCH4(g) \t\t",gas_str1[1]));
		else if (eqphases && line[2] == 'O' && line[3] == '2' && line[4] == '(' && line[5] == 'g' && line[6] == ')') fprintf(fout, "%s\n", ConCat("\tO2(g) \t\t",gas_str1[2]));
		else if (eqphases && line[2] == 'N' && line[3] == 't' && line[4] == 'g' && line[5] == '(' && line[6] == 'g') fprintf(fout, "%s\n", ConCat("\tNtg(g) \t\t",gas_str1[3]));
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
	for (i=0;i<nAqSpecies;i++) {
		free(aq_str[i]);
	}
	for (i=0;i<nvarKin;i++) {
		free(river_str[i]);
	}
	free(gas_str1);
	free(gas_str2);
	free(gas_str3);
	free(aq_str);
	free(river_str);

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
	if (cmdline == 1) strcat(cmd,"Debug/sel*");
	else strcat(cmd,"sel*");
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
			 + xgas[3]*2.0*M_N         // Ntg (i.e. inert N2)
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
//	strcat(aMELTStmp,"alphaMELTS-1.9/");
	strcat(aMELTStmp,"alphaMELTS-1.9/ > /dev/null"); // Prevent alphaMELTS from printing out to terminal
	strcat(aMELTSsys, aMELTStmp);

	// --- Run alphaMELTS ---
	printf("Running alphaMELTS\n");
	system(aMELTSsys);
	printf("alphaMELTS ran successfully\n");
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
