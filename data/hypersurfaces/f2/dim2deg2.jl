using Pkg
Pkg.activate(joinpath(@__DIR__, "../../.."))
Pkg.instantiate()

using DataFrames
using DuckDB
using DBInterface

using FiniteVarietiesDB
using Oscar


function main()

    F = GF(2)
    (f, S) = projective_hypersurface_equivalence_classes_from_filtration(F, 3, 2, verbose=true)

    function poly_to_string(f,s)
        return string(forget_grading(f(s)))
    end

    n_poly = length(S)
    classes = Vector{String}(undef, n_poly)
    for i in 1:n_poly
        s = pop!(S)
        classes[i] = poly_to_string(f,s)
    end

    df = DataFrame(
        field = 2,
        polynomial = classes,
        dimension = 2,
        degree = 2
    )

    output_filename = joinpath(@__DIR__, "hypersurfaces_f2_dim2_deg2.parquet")
    con = DBInterface.connect(DuckDB.DB)

    try
        DuckDB.register_data_frame(con, df, "df_temp")
        DBInterface.execute(con, "COPY df_temp TO '$output_filename' (FORMAT PARQUET)")
    finally
        DBInterface.close!(con)
        println("Wrote $(nrow(df)) rows to $output_filename")
    end
end

main()
