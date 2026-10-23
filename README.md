# Eyjafjallajokull_VolcanicAsh_Modeling
CESM model code changes to incorporate volcanic ash into CAM6 used to produce the results in Wu et al. "Modeling Volcanic Ash Impacts of the 2010 Eyjafjallajökull Eruption"
Version CESM2.1.0 for CAM6

## First Author Information
Name: Lily Wu\
ORCID: https://orcid.org/0009-0001-4019-508X \
Institution: Cornell University\
Address: 112 Hollister Dr.\
Email: lw784@cornell.edu

## Principal Investigator Information
Name: Natalie Mahowald\
ORCID: https://orcid.org/0000-0002-2873-997X \
Institution: Cornell University\
Address: 112 Hollister Dr.\
Email: mahowald@cornell.edu

# DATA & FILE OVERVIEW
A list of SourceMods changes for the src.cam module in CESM CAM6. Here, we highlight some important modifications for files that require changes!!!

***PLEASE READ THROUGH THIS SECTION IF YOU ARE INTERESTED IN MODIFYING THE CODE!!!!!***

* **modal_aer_opt.F90** - 

* **modal_aer_opt.F90.orig.Lily** - original modal_aer_opt file that included more iron and ash variables that were decomposed into the tracer components. 
