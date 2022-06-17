/*
 * ExoCcycleGeo.h
 *
 *  Created on: Jun 13, 2019
 *      Author: Marc Neveu (marc.f.neveu@nasa.gov)
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
#include <Var.h>                        // To use the external PHREEQC geochemical code
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
#define MPa2Pa 1.0e6					// Convert MPa to Pa
#define atm2bar 1.01325                 // Convert atm to bar, same conversion factor as for PHREEQC
#define Kelvin 273.15                   // Convert between Celsius and Kelvin, same conversion factor as for PHREEQC
#define Gyr2sec 3.15576e16              // Convert 1 billion years to s
#define Myr2sec 3.15576e13              // Convert 1 million years to s
#define Yr2sec 3.15576e7                // Convert 1 year to s

// Earth scalings
//#define TsurfEarth 288.0                // Mean equilibrium surface temperature on Earth (K)
//#define RmeltEarth 3.8e-19              // Rate of melt generation on Earth (s-1) Kite et al. 2009 Fig. 15; http://dx.doi.org/10.1088/0004-637X/700/2/1732
//#define deltaCvolcEarth 2.2e5           // Surface C flux from subaerial+submarine volcanic outgassing (mol C s-1) Donnadieu et al. 2006; http://dx.doi.org/10.1029/2006GC001278
//#define deltaCcontwEarth 8.4543e-10     // Surface C flux from continental weathering on Earth (mol C m-2 s-1)
//#define MORlength 60000.0*km2m          // Length of mid-ocean ridges, unconstrained parameter, default 60000 km (Earth today), likely did not vary monotonically in the past

// Accretion parameters
#define chi 2.0                         // Ratio of planetesimal velocity to escape velocity
#define h_frac 0.1                      // Fraction of impactor energy deposited at depth

// Geodynamics parameters
#define k 4.18                          // Mantle thermal conductivity (W m-1 K-1) (Turcotte and Schubert 2002)
#define alpha 3.0e-5                    // Mantle thermal expansivity (K-1) (Turcotte and Schubert 2002)
#define Cp 1295.0                       // Mantle specific heat capacity (J kg-1 K-1) (Akaogi & Ito 1993)
#define beta 1.0/3.0                    // Exponent for scaling Nusselt number to Rayleigh number, 1/4 to 1/3 (Schubert et al. 2001)
#define Ra_c 1707.762                   // Critical Rayleigh number for convection, http://home.iitk.ac.in/~sghorai/NOTES/benard/node15.html, see also Koschmieder EL (1993) Benard cells and Taylor vortices, Cambridge U Press, p. 20.
#define TminMELTS 750.0                 // Don't run MELTS below 750ºC to avoid it crashing (ºC)

// Atmosphere parameters
//#define xCO2g0 355.0e-6               // Reference atmospheric CO2 mixing ratio (ppmv)
#define runoff_Earth 7.75e-9            // Reference runoff (m s-1) = 0.67e-3 m day-1 (Edson et al. 2012, http://dx.doi.org/10.1089/ast.2011.0762)

// Geochem parameters
#define nAtmSpecies 5                   // Number of atmospheric species whose abundances are passed between the physical and chemical models
#define nAqSpecies 13                   // Number of aqueous species whose concentrations are passed between the physical and chemical models
#define nvarEq 1036                     // Number of geochemical variables stored in each PHREEQC equilibrium simulation
#define nvarKin 140                     // Number of geochemical variables stored in each PHREEQC kinetic simulation
#define naq 258                         // Number of aqueous species (+ physical parameters)
#define nmingas 389                     // Number of minerals and gases
#define nelts 31                        // 30 elements + 1 extra column in WaterRock/Molar_masses.txt
#define rhoH2O 1000.0                   // Density of water (kg m-3)
#define tcirc (1.0e7*Yr2sec)			// Timescale of ocean cycling through seafloor hydrothermal systems (Mottl 1983; Kadko et al. 1995) (s)

#ifndef EXOCCYCLEGEO_H_
#define EXOCCYCLEGEO_H_

double *exoCinput (double *input, char path[1024]);

//-------------------------------------------------------------------
//                  Read ExoCcycleGeo input file
//-------------------------------------------------------------------

double *exoCinput (double *input, char path[1024]) {

	FILE *f;
	int i = 0;
	int scan = 0;

	int line_length = 300;
	char line[line_length]; // Individual line
	int line_no = 0;        // Line number
	int tab = 51;           // Column number of inputs
//	fpos_t pos;

	char *title = (char*)malloc(1024);
	title[0] = '\0';
	if (cmdline == 1) strncat(title,path,strlen(path)-20);
	else strncat(title,path,strlen(path)-18);
	strcat(title,"Inputs/ExoCcycleGeoInput.txt");

	i = 0;
	f = fopen (title,"r");
	if (f == NULL) {
		printf("ExoCcycleGeo: Cannot find ExoCcycleGeo.txt file.\n");
		printf("Was ExoCcycleGeo launched from the right folder?\n");
		printf("The following option is active: command line %d\n", cmdline);
		exit(0);
	}
	else {
		while (fgets(line, line_length, f)) {
			line_no++;
			// Grid
			if (line_no == 5) {
				fgets(line, tab, f);  // Simulation time step (years)
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 6) {
				fgets(line, tab, f);  // Simulation start time (years)
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 7) {
				fgets(line, tab, f);  // Simulation end time (years)
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			// Interior
			else if (line_no == 11) {
				fgets(line, tab, f);  // Planet mass (Earth masses)
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 12) {
				fgets(line, tab, f);  // Core mass fraction
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 13) {
				fgets(line, tab, f);  // Inner layer material code
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 14) {
				fgets(line, tab, f);  // Mid layer material code
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 15) {
				fgets(line, tab, f);  // Outer layer material code
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 16) {
				fgets(line, tab, f);  // Radionuclide inventory
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 17) {
				fgets(line, tab, f);  // Rheology
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 18) {
				fgets(line, tab, f);  // Upper mantle redox
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 19) {
				fgets(line, tab, f);  // Mass fraction of carbon in the mantle
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 20) {
				fgets(line, tab, f);  // Mole fraction of C outgassed as CH4, relative to CH4+CO2
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			// Surface
			else if (line_no == 24) {
				fgets(line, tab, f);  // Mass of surface ocean (kg)
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 25) {
				fgets(line, tab, f);  // Areal land fraction
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 26) {
				fgets(line, tab, f);  // Initial temperature (K)
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 27) {
				fgets(line, tab, f);  // Initial pressure (bar)
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 28) {
				fgets(line, tab, f);  // Runoff rate (m/day)
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 29) {
				fgets(line, tab, f);  // Water residence time on continents (years)
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			// Atmosphere
			else if (line_no == 33) {
				fgets(line, tab, f);  // Atmospheric CO2 mixing ratio
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 34) {
				fgets(line, tab, f);  // Atmospheric CH4 mixing ratio
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 35) {
				fgets(line, tab, f);  // Atmospheric O2 mixing ratio
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 36) {
				fgets(line, tab, f);  // Atmospheric N2 mixing ratio
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
			else if (line_no == 37) {
				fgets(line, tab, f);  // Atmospheric H2O mixing ratio
				scan = fscanf(f, "%lg", &input[i]), i++;
				if (scan != 1) printf("Error scanning ExoCcycleGeo input file at entry i = %d\n",i);
			}
		}
	}
	fclose(f);

	free (title);
	return input;
}

#endif /* EXOCCYCLEGEO_H_ */
