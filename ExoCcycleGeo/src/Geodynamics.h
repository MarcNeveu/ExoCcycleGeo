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

double brittleDuctile (double T, double rhoLith, double zLith, double gsurf, double Tmantle, double Tsurf, double flowLawDiff[5], double flowLawDisl[5], double grainSize,
		double dtime);
double brittleDuctile_prime (double T, double rhoLith, double zLith, double gsurf, double Tmantle, double Tsurf, double flowLawDiff[5], double flowLawDisl[5],
		double grainSize, double dtime);
double viscosity (double T, double P, double flowLaw[5], double grainSize, double dtime);
double combVisc (double T, double P, double flowLawDiff[5], double flowLawDisl[5], double grainSize, double dtime);
double dviscdT (double T, double dPdT, double flowLaw[5], double grainSize, double dtime, double Tsurf);
double Psolidus (double T);
double dPsolidusdT (double T);

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

double brittleDuctile (double T, double rhoLith, double zLith, double gsurf, double Tmantle, double Tsurf, double flowLawDiff[5], double flowLawDisl[5], double grainSize,
		double dtime) {

	double f = 0.0;

	double P = rhoLith*gsurf*zLith*(T-Tsurf)/(Tmantle-Tsurf); // Pressure (Pa)
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

double brittleDuctile_prime (double T, double rhoLith, double zLith, double gsurf, double Tmantle, double Tsurf, double flowLawDiff[5], double flowLawDisl[5],
		double grainSize, double dtime) {

	double f_prime = 0.0;

	double P = rhoLith*gsurf*zLith*(T-Tsurf)/(Tmantle-Tsurf); // Pressure (Pa)
	double dPdT = rhoLith*gsurf*zLith/(Tmantle-Tsurf);        // Geotherm (Pa K-1), independent of T
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
