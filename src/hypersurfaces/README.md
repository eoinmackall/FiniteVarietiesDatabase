# Directory structure

There are two subdirectories.

## 1. mprocs

This directory contains an abandoned project that aims to make the hypersurface enumeration algorithms distributed, in order to rely on GAP functionality for orbit calculations.

The main difficulty in completing this project is serializing data to be transferred to distributed processes.

## 2. mthreads

This directory contains three files. All of the functions in this directory can be ran with one or multiple threads. Higher thread counts does not necessarily provide better performance.

### 2.1) projective_equivalence_classes.jl

This file contains code for enumerating representatives for equivalence classes of hypersurfaces over finite field under projective equivalence. The main function uses a union-find algorithm.

### 2.2) orbit_size.jl

This file contains code for calculating the size of the nonzero orbits of GL(n+1,q) on the space of homogeneous polynomials of degree d in n+1 variables over GF(q). If gcd(q-1,d)=1, then this is the same as the size of the set of orbit representatives for the linear action on projective hypersurfaces.

### 2.3) equiv_classes_filtration_method.jl

This file contains code for calculating orbit representatives for the GL(n+1,q) action on the space of homogeneous polynomials of degree d in n+1 variables over GF(q). The algorithm used filters the vector space of homogeneous polynomials by GL(n+1,q) stable subspaces, and lifts representatives along successive quotients.

___

# Examples

Here's an example computing a list of representatives for the set of 1,732,563 equivalence classes of quartic surfaces in P^3 over GF(2).

```julia-repl
julia> using FiniteVarietiesDB
julia> using Oscar
julia> F=GF(2);
...
julia> @time projective_hypersurface_equivalence_classes_from_filtration(F,3,4,verbose=true);
Starting chain collection...
Found the following chains:
35 34 30 24 4 0 
35 34 30 24 20 0 
35 34 30 10 4 0 
35 34 14 10 4 0 
Found chain at position #4 with maximal relative dimension = 20
Beginning orbit collection
Starting stage #1 -- finding orbits in V/V_1
Starting stage #2 -- lifting orbits along chain
Lifting step: #1
Lifting step: #2
Lifting step: #3
Starting stage #3 -- finding representatives in V
111.210500 seconds (2.73 G allocations: 177.072 GiB, 43.63% gc time, 2 lock conflicts, 1.28% compilation time)
```

