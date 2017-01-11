# husky.3

Husky computes the electronic transport properties of molecular junctions based on an extended Huckel electronic structure and using a non-equilibrium Green function formalism. 

# Installation

Once you've dowloaded the core, create a OBJ directory that will store the .o files. 
You have then to generate the Makefile using the script createMakefile.sh. 
Then type make and in all logic it shouldn't work. You may have to change a few #includes to locate the different library the code uses.

# Example

One example is provided. The conductance of a benzene di-thiol is computed using husky. 

![alt text](./junction.png)
