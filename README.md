# FiniteVarietiesDB.jl

This project aims to provide a free and accessible database for varities over finite fields.
This project makes critical use of Oscar.jl (the Open Source Computer Algebra Research system) as a 
foundation for mathematical algorithms and infrastructure and of DuckDB as a foundation for its database
functionality.

## Functionality

- Currently, this project provides functionality for computing a set of class representatives for equivalence classes of hypersurfaces in projective space over a finite field under projective equivalence.

## TODO:

The following is a list of future plans for this project.

1. Rewrite the data collection scripts to collect into individual (script local) database files that are then added to a larger database pool after data collection.
2. Write functionality for transfering the polynomial elements of the database into Oscar(+) objects.
3. Add computations of invariants for elements of the database (e.g. smoothness, zeta function, number of lines, etc.).
4. Rewrite the hypersurface collection algorithms in C using FLINT.
5. Add functionality for computing lists of isomorphism classes of curves of genus g.

