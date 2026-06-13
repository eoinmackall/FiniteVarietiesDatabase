using Pkg
Pkg.activate(joinpath(@__DIR__, "../../.."))
Pkg.instantiate()

using DuckDB
using DBInterface

function compile_hypersurfaces()
    hypersurfaces_path = joinpath(@__DIR__, "..", "hypersurfaces.db")
    parquet_pattern = joinpath(@__DIR__, "*.parquet")
    
    con = DBInterface.connect(DuckDB.DB, hypersurfaces_path)

    try
        # Create the table if it doesn't exist
        create_db = """
        CREATE TABLE IF NOT EXISTS hypersurfaces (
            field INTEGER,
            polynomial VARCHAR,
            dimension INTEGER,
            degree INTEGER
        );
        """
        DBInterface.execute(con, create_db)

        # Bulk insert all parquet files matching the pattern
        DBInterface.execute(con, """
            INSERT INTO hypersurfaces 
            SELECT * FROM read_parquet('$parquet_pattern')
        """)
        
    finally
        DBInterface.close!(con)
        println("Data compiled successfully into $hypersurfaces_path")
    end
end

compile_hypersurfaces()
