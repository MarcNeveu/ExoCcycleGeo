# ExoCcycleGeo
Calculates net C fluxes from a terrestrial planet's interior to its atmosphere with a coupled geophysical-geochemical, time-evolution box model of outgassing, weathering, and ocean-atmosphere equilibria. Input parameters include planet size (0.5 to 2 Earth masses), detailed chemical coomposition and interior structure, and land coverage and patchiness. The weathering and ocean-atmosphere equilibrium model computes kinetically limited and equilibrium water-rock reactions using the [USGS PHREEQC v3](https://www.usgs.gov/software/phreeqc-version-3) code. The outgassing model scales outgassing from melting, which is computed based on planet size and age (i.e., temperature-pressure profile) and composition using the [alphaMELTS](https://magmasource.caltech.edu/alphamelts/) implementation of the MELTS code. Currently only the outgassing model is fully benchmarked, so this is the only capability accessible in this repository.

# Installation

The installation steps outlined below are valid for Mac OS 10.9+. *ExoCcycleGeo* could also plausibly run on Windows and Linux, but compilation instructions are not set up and external I/O handling needs to be modified in the source code. 

## Install *R*
*R* is needed to run the geochemistry package *CHNOSZ*.
Go to http://www.r-project.org and follow instructions.

## Install *CHNOSZ*
Open *R* using either the installed application icon or in a terminal by typing

	R
	
In *R*, type the command

	install.packages("CHNOSZ")

## Install *Rcpp* and *RInside*
*Rcpp* and *RInside* are libraries that allow *R* applications to be embedded in C or C++ codes. From a Terminal window, open *R* and install the *Rcpp* and *RInside* packages:

	install.packages("Rcpp")
	install.packages("RInside")

## Install *IPHREEQC*
The *IPHREEQC* library is a module that allows the *PHREEQC* application to be embedded in C or C++ codes. Go to http://wwwbrr.cr.usgs.gov/projects/GWC_coupled/phreeqc to download *IPHREEQC* (for Linux), unzip it, and follow the default installation instructions (you need admin credentials on your machine):

	./configure
	make
	make install

In v3.8.6, the *./configure* script omits copying a header file *PHRQ_exports.h* to the include folder. You can manually remedy this by copying, from the unzipped folder, *src/phreeqcpp/common/PHRQ_exports.h* to the */usr/local/include/* folder.

## Install *ExoCcycleGeo*
Go to https://github.com/MarcNeveu/ExoCcycleGeo. Click the green *Clone or download* button to the right of the page, then either:
- *Download ZIP* on the bottom right. Unzip ExoCcycleGeo-master.zip. Move the *ExoCcycleGeo* folder inside to any folder you would like, we will call it *Path_to_GitFolder* here.
- if you are familiar with GitHub, you can clone the directory with your favorite tool (I use Git within the Eclipse developing environment).

The folders within */Path_to_GitFolder/ExoCcycleGeo/* contain:
- *src*: All source code files
- *Debug*: The executable program *ExoCcycleGeo*
- *alphaMELTS-1.9*: The *alphaMELTS* code executed by the outgassing subroutine within *ExoCcycleGeo*
- *PHREEQC-3.1.2*: The *PHREEQC* code executed by the weathering and ocean-atmosphere equilibrium subroutines within *ExoCcycleGeo*
- *Data*: Lookup tables. Currently there is only one used for calculated the planet's self-compression.

At this time, and I apologize for the clunkiness, you also need to manually edit three files to change the */Path_to_GitFolder/* from mine (*/Users/mneveu/eclipse-workspace/ExoCcycleGeo/ExoCcycleGeo/*) to yours (*/Path_to_GitFolder/ExoCcycleGeo/*):
- In *alphaMELTS-1.9/run_alphameltsExoC.command* at lines 10, 432 (3 instances), and 457 (2 instances)
- In *alphaMELTS-1.9/ExoC/ExoC_env.txt* at l. 6
- In *alphaMELTS-1.9/ExoCbatch.txt* at l. 2.
    
# Running the code

## Start the code
To execute the code, I usually open a Terminal window, get to the folder */Path_to_GitFolder/ExoCcycleGeo/Debug* which contains the executable by typing

    cd /Path_to_GitFolder/ExoCcycleGeo/Debug
    
and start the execution by typing 

    ./ExoCcycleGeo &
    
## Changing inputs

### Planet composition

**Outgassing**: the composition that controls carbon outgassing in the model is that of the upper mantle, which is set in *alphaMELTS-1.9/ExoC/ExoCcycleGeo.melts*. To change it, you can either edit the "Initial Composition" lines of this file or create a duplicate, *Duplicate.melts*, and tell *ExoCbatch.txt* in the same folder to run *Duplicate.melts* rather than *ExoCcycleGeo.melts*. The default composition is that of a mid-ocean ridge basalt (MORB).

**Ocean-atmosphere equilibrium**: the composition that controls ocean-atmosphere equilibrium in the model is set in *PHREEQC-3.1.2/io/OceanDiss.txt*. At the first time step, it is initiated instead with *PHREEQC-3.1.2/io/OceanStart.txt*. The default composition is that of Earth's ocean.

**Continental weathering**: the composition that controls continental weathering in the model is set in *PHREEQC-3.1.2/io/ContWeather.txt*.

### Other inputs

Other inputs that can be modified in *Inputs/ExoCcycleGeoInput.txt* include:
- **Simulation time step** (default 10 million years, which achieves numerical convergence for the outgassing part of the code);
- **Simulation start time** (default 0.6 billion years (Gyr) after formation, i.e. 4.57-0.6=3.97 Gyr ago);
- **Simulation end time** (default 5 Gyr after formation, i.e. 0.4 Gyr into the future);
- **Planet mass** (default 1 Earth mass, min 0.5 Earth masses, max 2 Earth masses);
- **Core mass fraction** (default 0.325), which sets the mantle thickness, vigor of convection, and ultimately the near-surface pressure-temperature profile;
- **Rough core and mantle compositions**: three 3-digit codes used to identify three materials used to calculate core and mantle radii based on mass from compression equations of state (default codes set to approximate Earth's radius, solid metal inner core, liquid metal outer core, and silicate mantle). Materials and material parameters are provided in *Data/Compression_planmat.txt*;
- **Radionuclide abundances**: 5 choices of radionuclide abundances: enter `0` for case 3 of Lyubetskaya & Korenaga (2007), `1` for the highest heating case from Turcotte & Schubert (2002) as reported in Table 1 of Kite et al. (2009), `2` for an intermediate heating case from Ringwood (1991) as reported in Table 1 of Kite et al. (2009), `3` for the lowest heating case from Table 2 of Lyubetskaya & Korenaga (2007), other entries default to abundances of McDonough & Sun (1995) as reported in Table 2 of Lyubetskaya & Korenaga (2007) that result in intermediate heating similar to Ringwood (1991);
- **Rheology** (of the mantle rock): enter `0` for the dry olivine rheology of Korenage & Karato (2008), or `1` for their wet olivine rheology. Other entries are not accepted;
- **Redox state** (of the ocean crust): enter `1` for the present-day Earth surface (log fO2 = -0.7 here fO2 is the oxygen fugacity), `2` to set fO2 at the hematite-magnetite (HM) redox buffer, `3` to set fO2 at the fayalite-magnetite-quartz (FMQ) buffer, `4` to set fO2 at the iron-wüstite (IW) buffer;
- **Mass fraction of carbon in the mantle** (default 0.2%);
- **Mass of the ocean** (default 1.4e21 kg, the mass of Earth's ocean);
- **Areal land fraction** (default 29%);
- **Initial surface temperature** (default 288 K);
- **Initial surface pressure** (default 1 bar);
- **Runoff rate**, i.e. rate at which water is delivered from rivers into the ocean (default 0.7 mm/day, i.e., 70 cm/year precipitation with 35% of that rainwater reaching the ocean and the rest evaporated back to the atmosphere before reaching the ocean);
- **Residence time of water on continents** (default 10 years), i.e. continental patchiness;
- **Initial atmospheric composition** (default 280 ppm CO2, 0.0001 ppm CH4, 20% O2, 79% N2, 1% H2O).

## Reading outputs

- *CompoAtmosph.txt* provides the composition of the atmosphere at each time step.
- *CompoOcean.txt* provides the size and composition of the ocean at each time step.
- *ContweatherFluxes.txt* provides the bicarbonate trapping capacity of cations dissolved in rivers at each time step.
- *Geotherm.txt* provides the pressure, temperature, melt fraction, and density profiles at each time step, in radial grid zones from the upper mantle up to the surface.
- *Outgassing.txt* provides geodynamic parameters of the planet at each time step.
- *ReservoirsFluxes.txt* provides reservoirs and fluxes of carbon at each time step.

These files are all tab-separated. Parameters and units are provided in the first 1-2 rows of each file.
    
# Modifying the source code

If you wish to modify the code, set up your compiler and linker so that all the relevant flags are added. My compiling and linking commands (executed in a Terminal window from within the *ExoCcycleGeo/Debug* folder) look like this:
 
    gcc -I/usr/local/include -I/Library/Frameworks/R.framework/Versions/Current/Resources/include -I/Library/Frameworks/R.framework/Versions/Current/Resources/library/RInside/include -I/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include -O0 -g3 -Wall -c -fmessage-length=0 -o src/ExoCcycleGeo.o ../src/ExoCcycleGeo.c
    gcc -L/usr/lib -L/usr/local/lib -L/Library/Frameworks/R.framework/Versions/Current/Resources/lib /usr/local/lib/libiphreeqc-3.8.6.dylib /usr/local/lib/libiphreeqc.dylib /usr/local/lib/libiphreeqc.a -o ExoCcycleGeo src/ExoCcycleGeo.o -lR

Email me if you have any issues.

# References

Antoshechkina, P. M., and P. D. Asimow (2010), alphaMELTS 3.0 and the MAGMA website: educational and research tools for studying the petrology and geochemistry of plate margins, Abstract ED41B-0644 presented at 2010 Fall Meeting, AGU, San Francisco, Calif., 13-17 Dec.

Antoshechkina, P. M., P. D. Asimow, E. H. Hauri, and P. I. Luffi (2010), Effect of water on mantle melting and magma differentiation, as modeled using alphaMELTS 3.0, Abstract V53C-2264 presented at 2010 Fall Meeting, AGU, San Francisco, Calif., 13-17 Dec.

Asimow, P.D., and M.S. Ghiorso (1998), Algorithmic modifications extending MELTS to calculate subsolidus phase relations, American Mineralogist, 83 (9-10), 1127-1132.

Ghiorso, M.S., and R.O. Sack (1995), Chemical Mass-Transfer in Magmatic Processes IV. A Revised and Internally Consistent Thermodynamic Model for the Interpolation and Extrapolation of Liquid-Solid Equilibria in Magmatic Systems at Elevated-Temperatures and Pressures, Contributions to Mineralogy and Petrology, 119 (2-3), 197-212.

Ghiorso, M.S., M.M. Hirschmann, P.W. Reiners, and V.C. Kress (2002), The pMELTS: A revision of MELTS for improved calculation of phase relations and major element partitioning related to partial melting of the mantle to 3 GPa, Geochemistry Geophysics Geosystems, 3(5), art. no. 1030, doi:10.1029/2001GC000217.

Neveu, M., Bartlow, T.L.C., Felton, R., Domagal-Goldman, S.D., and Desch, S. (2023) Geophysical and geochemical controls on abiotic carbon cycling on Earth-like planets, ESS Open Archive preprint (submitted to JGR: Planets), https://doi.org/10.22541/essoar.167751627.79222663/v1.

Parkhurst, D.L. and Appelo, C.A.J., 2013. Description of input and examples for PHREEQC version 3: a computer program for speciation, batch-reaction, one-dimensional transport, and inverse geochemical calculations (No. 6-A43). US Geological Survey.

Smith, P. M., and P. D. Asimow (2005), Adiabat_1ph: A new public front-end to the MELTS, pMELTS, and pHMELTS models, Geochem. Geophys. Geosyst., 6, art. no. Q02004, doi:10.1029/2004GC000816.

Thompson, R. N., A. J. V. Riches, P. M. Antoshechkina, D. G. Pearson, G. M. Nowell, C. J. Ottley, A. P. Dickin, V. L. Hards, A.-K. Nguno and V. Niku-Paavola (2007), Origin of CFB magmatism: Multi-tiered intracrustal picrite-rhyolite magmatic plumbing at Spitzkoppe, western Namibia, during early-Cretaceous Etendeka magmatism, J. Petrol., 48, 1119-1154.
