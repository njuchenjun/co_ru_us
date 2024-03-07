# co_ru_us

Umbrella sampling along the desorption path
CO @ Ru(0001), using a 6D or 42D PES.

## input files

the current INPUT is for 6D US-free energy calculations,
i.e., all surface atoms are fixed.

the INPUT.12ru is for 42D calculations,
with 3 layers of surface Ru atoms (4x4 per layer) relaxed.

## pes codes

the pes is called via "pes\_interface.f90", from "MD\_VV.f90".

the 6D pes contains 2 files,
"nnfit\_3.14.txt" and "pes\_ruFIXco.f90"

the HD pes contains 2 files and 1 folder,
"pes\_hd.f90", "pes\_hd.inc1.f90", and "pes8.dat"


