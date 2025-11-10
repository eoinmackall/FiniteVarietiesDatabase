module FiniteVarietiesDB

using Oscar
using Base.Threads
import Base: ==, hash

include("./hypersurfaces/equiv_classes_filtration_method.jl")
include("./hypersurfaces/orbit_size.jl")
include("./hypersurfaces/projective_equivalence_classes.jl")


export DChainNode, add_subnode!, delete_subnode!, add_subspace_descending

export homogeneous_monomial_basis
export HomogeneousExponents
export normal_forms
export orbit_size
export projective_hypersurface_equivalence_classes
export ProjectiveCoefficients
export symmetric_representation

end
