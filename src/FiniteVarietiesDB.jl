module FiniteVarietiesDB

using Oscar
using Base.Threads
import Base: ==, hash

include("./hypersurfaces/equiv_classes_filtration_method.jl")
include("./hypersurfaces/orbit_size.jl")
include("./hypersurfaces/projective_equivalence_classes.jl")


export DChainNode
export add_subnode!
export delete_subnode!
export is_in
export collect_chains

export homogeneous_monomial_basis
export HomogeneousExponents
export normal_forms
export orbit_size
export PGL
export projective_hypersurface_equivalence_classes
export projective_hypersurface_equivalence_classes_from_filtration
export ProjectiveCoefficients
export symmetric_representation

end
