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

* `src/ordering/ordering-kernel`.This is the script you should use first. It groups identical rows/columns into kernels and then
  uses simulated annealing to identify the ordering of those kernels that
  optimizes the BIC of the box model.

* `src/ordering/ordering-exhaustive`. This is the script you should use second. It exhaustively checks if reordering of adjacent kernels can further improve the the BIC of the box model.

* `src/boxes/main_boxes_greedy_nbox_ic_wdiag_kernel`. Starting from an initial
  ordering, this identifies break points between "communities" in the matrix and accepts them based on steepest decent and whether or not the BIC of the clustered model is significantly improved.

You can view the options of any of these scripts with the `--help` flag.

## Usage

The primary way to use this package is to run the `ordering-kernel` executable
with something like this:

```sh
./src/ordering/ordering-kernel -f path/to/matrix.dat
```

where `path/to/matrix.dat` is a text file that stores the information in the
matrix. By default, `ordering-kernel` assumes that the matrix will be in a
format similar to [data/test-matrix.dat]. Using the `-d` flag you can
optionally specify that the data is stored in a sparse format.

After this program is run, it outputs several files in the working directory,
the most important of which are:

* `coclas-final.dat` which prints the final reordering of the matrix.

* `transtable-final.dat` which lists how each input index is mapped to the
  output index.

This can then be followed up by the additional optimization step:

```sh
./src/ordering/ordering-exhaustive -f path/to/coclas-final.dat
```

which will generate two additional output files of particular interest:

* `coclas-exhaustive.dat` which prints the exhaustive reordering of the matrix.

* `transtable-exhaustive.dat` which lists how each input index is mapped to the
  output index. Note that, following the typical workflow, this translates from `coclas-final.dat` to `coclas-exhaustive.dat` and not from the *original* matrix (i.e. `matrix.dat` in this example).

Given either of these reordered output matrices, it is then possible to identify the most parsimonious grouping of rows/columns into "communities" based on BIC using:

```sh
./src/boxes/main_boxes_greedy_nbox_ic_wdiag_kernel N path/to/coclas-exhaustive.dat 1 0
```

where `N` is the size of the matrix and the final two arguments are poorly commented at present (though these default values should work in the vast majority of cases).


[doi]: http://dx.doi.org/10.1073/pnas.0703740104
