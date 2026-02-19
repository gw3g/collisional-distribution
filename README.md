# collisional-distribution

Calculation of energy loss/gain distribution for fast parton within the HTL approximation.

This code uses features of the [GNU Scientific Library](http://www.gnu.org/software/gsl/), 
which is free software under the GNU GPL.  Some of the examples also use the 
[OpenMP](http:://openmp.org) API for parallelizing data generation (although this is not 
needed for the main algorithm). 


Usage
-----

For the most general use, you should compile `src/base.cpp` and link it with your program, 
and include the corresponding header file `inc/base.hpp`. This will provide access to useful subroutines 
described below. Several examples are provided and described in the next section (see the provided 
makefile for an example of how to compile with `g++`). 

The necessary ingredients for the single scattering spectrum is:

```c++
double S_full(double epsilon, double T, double mu, double nf);
double R_hard(double epsilon, double T, double mu, double nf);
```
where `epsilon`, `T` and `mu` are all give in units of the Debye mass.

For calculating the probability distribution itself, interpolation grids need to be precomputed first for the functions G(eta) and H(eta). 

```c++
void tabulate_G_and_H(int N_eta, double T, double mu, double nf);
```
where `N_eta` is the number of points in the grid. 

Finally, to actually carry out the eta-integration, the `itp` namespace shoud be used and initialized, e.g.

```c++
using namespace itp;
init(double T, double mu, double nf);
```

The probability distribution f is computed through `lev_int(double xbar,double Delta)` ...


Examples
--------

### single scattering spectrum

May be compiled with `make spectrum`, and run with `./bin/spectrum`.

### probability distribution

May be compiled with `make distribution`, and run with `./bin/distribution`.

