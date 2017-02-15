/*
 ============================================================================
 Name        : ExoCcycleGeo.c
 Author      : Marc Neveu
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI_greek 3.14159265358979323846
#define G 6.67e-11                      // Gravitational constant (SI)
#define rEarth 6378000.0                // Earth radius (m)
#define mEarth 6.0e24                   // Earth mass (kg)
#define km2m 1000.0                     // Convert from km to m

#define TsurfEarth 288                  // Mean equilibrium surface temperature on Earth (K)
#define RmeltEarth 3.8e-19              // Rate of melt generation on Earth (s-1 kg-1) Kite et al. 2009 Fig. 15; http://dx.doi.org/10.1088/0004-637X/700/2/1732
#define deltaCvolcEarth 2.2e5           // Surface C flux from subaerial+submarine volcanic outgassing (mol C s-1) Donnadieu et al. 2006; http://dx.doi.org/10.1029/2006GC001278
#define deltaCcontwEarth 8.4543e-10*4.0*PI_greek*rEarth*rEarth // Surface C flux from continental weathering on Earth (mol C s-1)

#define Ra_c 1707.762                   // Critical Rayleigh number for convection, http://home.iitk.ac.in/~sghorai/NOTES/benard/node15.html, see also Koschmieder EL (1993) Benard cells and Taylor vortices, Cambridge U Press, p. 20.
#define Cp 914.0                        // Specific heat capacity of mantle rock (J kg-1 K-1)
#define k 4.18                          // Thermal conductivity of crust (W m-1 K-1)
#define A0 70000.0                      // Activation temperature for thermal model (K)

#define atmCO2_0 355                    // Reference atmospheric CO2 mixing ratio (ppmv)
#define runoff_0 0.665                  // Reference runoff (mm/day) Edson et al. 2012, http://dx.doi.org/10.1089/ast.2011.0762

int main(int argc, char *argv[]) {

	int cmdline = 0;

	double r_p = 0.0;               // Planet radius (m)

	double deltaCvolc = 0.0;        // Surface C flux from volcanic outgassing (subaerial+submarine)
	double deltaCcontw = 0.0;       // Surface C flux from continental weathering
	double deltaCseafw = 0.0;       // Surface C flux from seafloor weathering
	double deltaC = 0.0;            // Net surface C flux from all geological processes

	double Rmelt = 0.0;             // Rate of melt generation
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

	double m_p = 2.0*mEarth;        // Planet mass (kg)
	double L = 0.3;             // Fraction of planet surface covered by land

	// Atmospheric inputs
	double Tsurf = 288.0;       // Surface temperature (K)
	double atmCO2 = 400;        // Atmospheric CO2 mixing ratio (ppmv)
	double runoff = 0.7;        // Atmospheric runoff (mm/day)

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
	printf("Atmospheric CO2 mixing ratio \t \t %g ppmv\n", atmCO2);
	printf("Surface runoff \t \t \t \t %g mm d-1\n", runoff);
	printf("-------------------------------------------------------------------\n");

	printf("\n");
	printf("Computing geo C fluxes...\n");

	r_p = pow(m_p/mEarth,0.27)*rEarth;  // Planet radius (m)
	gsurf = G*m_p/r_p/r_p;
	Pf = rhoCrust*gsurf*zCrust;
	TbaseLid = Tmant - 2.23*Tmant*Tmant/A0;
	Nu = pow(Ra/Ra_c,0.25);

	//-------------------------------------------------------------------
	// Calculate surface C flux from volcanic outgassing (Kite et al. 2009)
	//-------------------------------------------------------------------

	// Assumes that all melt generated reaches the surface
	vConv = 2.0*(Nu-1.0) * (k/(rhoMagma*Cp*(r_p-r_c))) * (Tmant-Tsurf) / (Tmant-TbaseLid);
	Rmelt = 4.0*PI_greek*r_p*r_p*vConv*rhoMagma/m_p * (rhoCrust*gsurf*zCrust/(Pf-P0));
	deltaCvolc = deltaCvolcEarth * Rmelt/RmeltEarth;

	//-------------------------------------------------------------------
	// Calculate surface C flux from continental weathering (Edson et al. 2012)
	//-------------------------------------------------------------------

	deltaCcontw = -L * 0.5*deltaCcontwEarth*r_p*r_p/rEarth/rEarth * pow(atmCO2/atmCO2_0,0.3) * runoff/runoff_0 * exp((Tsurf-TsurfEarth)/17.7);

	//-------------------------------------------------------------------
	// Calculate surface C flux from ocean dissolution
	//-------------------------------------------------------------------

	//-------------------------------------------------------------------
	// Calculate surface C flux from seafloor weathering
	//-------------------------------------------------------------------

	deltaCreac = 0.0; // TODO call PHREEQC to get deltaCreac, net mol C leached/precipitated per kg of rock
	tcirc = 1.0; // TODO compute tcirc based on Nu(Ra(basal heat flux))
	deltaCseafw = -(1.0-L) * 4.0/3.0*PI_greek*(pow(r_p,3)-pow(r_p-zCrack,3))/tcirc*deltaCreac*rhoCrust;

	//-------------------------------------------------------------------
	// Compute net geo C flux
	//-------------------------------------------------------------------

	deltaC = deltaCvolc + deltaCcontw + deltaCseafw;

	printf("\n");
	printf("Net surface C flux from... \t mol C s-1\n");
	printf("-------------------------------------------\n");
	printf("volcanic outgassing \t \t %g\n", deltaCvolc);
	printf("continental weathering \t \t %g\n", deltaCcontw);
	printf("seafloor weathering \t \t %g\n", deltaCseafw);
	printf("-------------------------------------------\n");
	printf("all geo processes \t \t %g\n", deltaC);

	printf("\nExiting ExoCcycleGeo...\n");
	return EXIT_SUCCESS;
}
