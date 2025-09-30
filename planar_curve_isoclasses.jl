using Nemo

if length(ARGS) != 2
    error("Please provide a prime power p^n and a degree d.")
end

x = parse(BigInt, ARGS[1])
y = parse(BigInt, ARGS[2])


function reduce(x::BigInt)

    if x < typemax(Int128)
        x = Int128(x)
    end

    if x < typemax(Int64)
        x = Int64(x)
    end

    if x < typemax(Int)
        x = Int(x)
    end

    return x
end

x = reduce(x)
y = reduce(y)


"""
Produces a ?? containing a unique isomorphism class for each curve of degree d
over a finite field of p^n elements
"""
function plane_curve_isomorphism_classes(x, y)









end

