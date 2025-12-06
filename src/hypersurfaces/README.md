Here's an example.

```julia-repl
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
