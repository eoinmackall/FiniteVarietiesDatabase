# Completed

using Pkg
Pkg.activate(joinpath(@__DIR__, "../../.."))

using DataFrames
using DuckDB
using DBInterface

using FiniteVarietiesDB
using Oscar

F=GF(2)
classes=[string(f) for f in projective_hypersurface_equivalence_classes_from_filtration(F, 2, 6, verbose=true)]
df=DataFrame(
    field = 2,
    polynomial = classes,
    dimension = 1,
    degree = 6
)

hypersurfaces_path = joinpath(@__DIR__, "..", "hypersurfaces.db")

con = DBInterface.connect(DuckDB.DB, hypersurfaces_path)
DuckDB.register_data_frame(con, df, "df_temp")

create_db = """
CREATE TABLE IF NOT EXISTS hypersurfaces (
    field INTEGER,
    polynomial VARCHAR,
    dimension INTEGER,
    degree INTEGER
);
"""

DBInterface.execute(con, create_db)
DBInterface.execute(con, "INSERT INTO hypersurfaces SELECT * FROM df_temp")
DBInterface.close!(con)

println("Wrote $(nrow(df)) rows to $hypersurfaces_path")
