# myqcs

myqcs is a quantum circuit simulator written in Julia.

The package is developed to work in the jupyter notebook. 

The circuit can be built using commands or by reading in an OpenQASM file. the user can also specify whether they want the simulation done in density matrices or as a wave-function. 

Noise can be added to the density matrices. The noise can be added uniformly or can be specify for each gate.

Example for the use of the package can be found in the nb folder.

The package uses the libraries ITensors, Random, DelimitedFiles
