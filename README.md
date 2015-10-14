## msp-boxes

This package implements the functions needed to order a matrix and identify
it's "communities" using a combination of simulated annealing and BIC
(Sales-Pardo et al. 2007) in C++. As input, it takes a non-negative adjacency
or similarity matrix. Full details of the algorithm are available at:

Marta Sales-Pardo, Roger Guimerà, André A. Moreira, and Luís A. Nunes Amaral (2007) "Extracting the hierarchical organization of complex systems."
*PNAS* 104(39): 15224-15229 doi:[10.1073/pnas.0703740104][doi].

## Build instructions
After cloning the repository, you can build the executables by running:

    sudo apt-get install autoconf libtool
    autoreconf --install
    ./configure
    make

This produces:

* `src/ordering/ordering-kernal`.This is the main script you should use. It
  uses simulated annealing to identify the ordering of rows/columns that
  optimizes the BIC of the box model.

* `src/ordering/ordering-exhaustive`. This exhaustively searches the ordering
  of rows/columns that optimizes the BIC of the box model.

* `src/boxes/main_boxes_greedy_nbox_ic_wdiag_kernel`. Starting from an initial
  row ordering, this optimizes the row ordering using a steepest decent
  approach. This is effectively a zero temperature quenching of the simulating
  annealing procedure.

You can view the options of any of these scripts with the `--help` flag.

[doi]: http://dx.doi.org/10.1073/pnas.0703740104
