/*
 * Geodynamics.h
 *
 *  Created on: Jun 14, 2019
 *      Author: Marc Neveu (marc.f.neveu@nasa.gov)
 */

#include "ExoCcycleGeo.h"

#ifndef GEODYNAMICS_H_
#define GEODYNAMICS_H_

//-------------------------------------------------------------------
// SUBROUTINE DECLARATIONS
//-------------------------------------------------------------------

double viscosity (double T, double P, double flowLaw[6], double grainSize, double C_OH, double dtime);
double combVisc (double T, double P, double flowLawDiff[6], double flowLawDisl[6], double grainSize, double C_OH, double dtime);
double Psolidus (double T);
double dPsolidusdT (double T);

/*--------------------------------------------------------------------
 *
 * Subroutine viscosity
 *
 * Calculates non-Newtonian viscosity in Pa s.
 *
 *--------------------------------------------------------------------*/

double viscosity (double T, double P, double flowLaw[6], double grainSize, double C_OH, double dtime) {

	double visc = 0.0;

	visc = MPa2Pa*pow(pow(10.0,-flowLaw[2]) * pow(grainSize,flowLaw[4]) * pow(C_OH,-flowLaw[5]) * exp((flowLaw[0] + P*flowLaw[1])/(R_G*T)), 1.0/flowLaw[3]) * pow(1.0/dtime, 1.0/flowLaw[3]-1.0);

	return visc;
}

/*--------------------------------------------------------------------
 *
 * Subroutine combVisc
 *
 * Calculates non-Newtonian viscosity in Pa s combined from parallel
 * diffusion and dislocation flow laws (e.g. Cizkova et al. 2012 eq. 1).
 *
 *--------------------------------------------------------------------*/

double combVisc (double T, double P, double flowLawDiff[6], double flowLawDisl[6], double grainSize, double C_OH, double dtime) {

	double visc = 0.0;
	double viscDryDiff = 0.0;
	double viscDryDisl = 0.0;

	viscDryDiff = viscosity(T, P, flowLawDiff, grainSize, C_OH, dtime); // Korenaga & Karato (2008), dry diffusion
	viscDryDisl = viscosity(T, P, flowLawDisl, grainSize, C_OH, dtime); // Korenaga & Karato (2008), dry dislocation
	visc = 1.0/(1.0/viscDryDiff + 1.0/viscDryDisl); // Parallel combination

	return visc;
}

/*--------------------------------------------------------------------
 *
 * Subroutine Psolidus
 *
 * Calculates solidus pressure as a function of temperature using
 * expression of McKenzie & Bickle (1988)
 *
 *--------------------------------------------------------------------*/

double Psolidus (double T) {

	double Psol = 0.0;

	Psol = (T-Kelvin - 1100.0)/136.0 + 4.968e-4*exp(1.2e-2*(T-Kelvin-1100.0)); // in GPa
	Psol = Psol*1.0e9;

	return Psol;
}

/*--------------------------------------------------------------------
 *
 * Subroutine dPsolidusdT
 *
 * Calculates derivative of solidus pressure as a function of
 * temperature using expression of McKenzie & Bickle (1988)
 *
 *--------------------------------------------------------------------*/

double dPsolidusdT (double T) {

	double dPsoldT = 0.0;

	dPsoldT = (T-Kelvin)/136.0 + 4.968e-4*1.2e-2*(T-Kelvin)*exp(1.2e-2*(T-Kelvin-1100.0));
	dPsoldT = dPsoldT*1.0e9;

	return dPsoldT;
}

#endif /* GEODYNAMICS_H_ */
