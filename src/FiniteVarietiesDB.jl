module FiniteVarietiesDB

using Oscar
using Base.Threads
using Distributed
import Base: ==, hash

if nprocs() > 1
    @everywhere using Oscar
    @everywhere using Distributed
    @everywhere using Base.Threads
end

include("./misc/chain_maker.jl")
include("./hypersurfaces/mthreads/equiv_classes_filtration_method.jl")
include("./hypersurfaces/mthreads/orbit_size.jl")
include("./hypersurfaces/mthreads/projective_equivalence_classes.jl")
#include("./hypersurfaces/mprocs/equiv_classes_filtration_method2.jl")

if nprocs() > 1
    include("./hypersurfaces/mprocs/equiv_classes_filtration_method2.jl")
end

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
