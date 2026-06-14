using Oscar
using DuckDB
using DataFrames

function parse_poly_builder(poly_str::String, R)
    builder = MPolyBuildCtx(R)
    vars = gens(R)
    var_idx = Dict(string(v) => i for (i, v) in enumerate(vars))
    num_vars = length(vars)
    
    clean_str = replace(poly_str, " " => "")
    terms = split(replace(clean_str, "-" => "+-"), "+", keepempty=false)
    
    for term in terms
        exp_vec = zeros(Int, num_vars)
        coeff = 1
        factors = split(term, "*")
        
        for factor in factors
            if occursin("^", factor)
                v_str, e_str = split(factor, "^")
                exp_vec[var_idx[v_str]] += parse(Int, e_str)
            elseif haskey(var_idx, factor)
                exp_vec[var_idx[factor]] += 1
            elseif factor == "-"
                coeff *= -1
            else
                coeff *= parse(Int, factor)
            end
        end
        push_term!(builder, base_ring(R)(coeff), exp_vec)
    end
    return finish(builder)
end

# Database Configuration
db_path = "data/hypersurfaces/hypersurfaces.db"
con = DBInterface.connect(DuckDB.DB, db_path)

# 1. Read entire table into memory (moved outside the try block)
df = DBInterface.execute(con, "SELECT * FROM hypersurfaces") |> DataFrame
n_rows = nrow(df)

# 2. Pre-generate rings serially
unique_configs = unique(zip(df.field, df.dimension))
rings = Dict{Tuple{Int, Int}, Any}()

for (q, n) in unique_configs
    F = GF(q)
    var_names = ["x$i" for i in 0:n+1] 
    R, _ = graded_polynomial_ring(F, var_names)
    rings[(q, n)] = R
end

# 3. Compute smoothness in parallel
is_smooth_col = Vector{Bool}(undef, n_rows)

for i in 1:n_rows
    q = df.field[i]
    n = df.dimension[i]
    poly_str = df.polynomial[i]
    
    R = rings[(q, n)]
    f = parse_poly_builder(poly_str, R)
    
    V = variety(f, check=false)
    is_smooth_col[i] = is_smooth(V)
end

# 4. Attach new column to DataFrame
df.is_smooth = is_smooth_col

# 5. Safe Bulk Database Update via Transaction
try
    # Transaction now only wraps the actual update
    DBInterface.execute(con, "BEGIN TRANSACTION")
    
    DuckDB.register_data_frame(con, df, "updated_hypersurfaces")
    DBInterface.execute(con, "CREATE OR REPLACE TABLE hypersurfaces AS SELECT * FROM updated_hypersurfaces")
    
    DBInterface.execute(con, "COMMIT")
catch e
    DBInterface.execute(con, "ROLLBACK")
    rethrow(e)
finally
    close(con)
end
