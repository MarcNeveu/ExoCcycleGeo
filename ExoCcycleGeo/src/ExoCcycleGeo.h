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
#define TsurfEarth 288.0                // Mean equilibrium surface temperature on Earth (K)
#define RmeltEarth 3.8e-19              // Rate of melt generation on Earth (s-1) Kite et al. 2009 Fig. 15; http://dx.doi.org/10.1088/0004-637X/700/2/1732
#define deltaCvolcEarth 2.2e5           // Surface C flux from subaerial+submarine volcanic outgassing (mol C s-1) Donnadieu et al. 2006; http://dx.doi.org/10.1029/2006GC001278
#define deltaCcontwEarth 8.4543e-10     // Surface C flux from continental weathering on Earth (mol C m-2 s-1)
#define MORlength 60000.0*km2m          // Length of mid-ocean ridges, unconstrained parameter, default 60000 km (Earth today), likely did not vary monotonically in the past

// Mantle parameters
#define k 4.18                          // Mantle thermal conductivity (W m-1 K-1)
#define alpha 3.0e-5                    // Mantle thermal expansivity (K-1)
#define Cp 914.0                        // Mantle specific heat capacity (J kg-1 K-1)

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

#ifndef EXOCCYCLEGEO_H_
#define EXOCCYCLEGEO_H_



#endif /* EXOCCYCLEGEO_H_ */
