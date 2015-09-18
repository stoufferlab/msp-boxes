## msp-boxes

This package implements the functions neeeded to order a matrix and identify it's "communities" using a combination of simulated annealing and BIC (Sales-Pardo et al. 2007) in C++.

Marta Sales-Pardo, Roger Guimerà, André A. Moreira, and Luís A. Nunes Amaral (2007) "Extracting the hierarchical organization of complex systems."
*PNAS* 104(39): 15224-15229 doi:[10.1073/pnas.0703740104][doi].

## Build instructions
After cloning the repository, you can build the executables by:

        sudo apt-get install autoconf libtool
	autoreconf --install
	./configure
	make

This produces:

* `src/ordering/ordering-kernal` what does this do?
* `src/ordering/ordering-exhaustive` what does this do?
* `src/boxes/main_boxes_greedy_nbox_ic_wdiag_kernel` what does this do?

You can view the output of any of these scripts with the `--help` flag.

[doi]: http://dx.doi.org/10.1073/pnas.0703740104
