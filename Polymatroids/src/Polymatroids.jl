module Polymatroids

using Oscar
using Combinatorics: powerset

export PolymatroidCone, h, Δ, dim, facets, rays, is_matroid

struct PolymatroidCone
    n::Int
    vars::Vector

    function PolymatroidCone(n::Int)
        vars = collect(powerset(1:n)) # TODO: Order fixed?
        new(n, vars)
    end
end

function getvar(H::PolymatroidCone, A::Vector)::Int
    A = sort(unique(A))
    findfirst(B -> A == B, H.vars)
end

function h(H::PolymatroidCone, A::Vector)::Vector
    v = [0 for x in H.vars]
    v[getvar(H, A)] = 1
    v
end

function Δ(H::PolymatroidCone, i::Int, j::Int, K::Vector)::Vector
    h(H, [i, K...]) + h(H, [j, K...]) - h(H, [i, j, K...]) - h(H, K)
end

function cone(H::PolymatroidCone)
    N = 1:H.n
    A = Vector()
    # Normalized
    push!(A, +h(H, []))
    push!(A, -h(H, []))
    # Submodular
    for i in N
        push!(A, -Δ(H, i, i, setdiff(N, i)))
    end
    for ij in powerset(N, 2)
        (i, j) = ij
        for K in powerset(setdiff(N, ij))
            push!(A, -Δ(H, i, j, K))
        end
    end
    Polyhedron(transpose(hcat(A...)), [0 for x in A])
end

# Make sure all rays have coprime integer coordinates.
normalization_factor(r)    = lcm([denominator(x) for x in r]...) // gcd([numerator(x) for x in r])
rays(H::PolymatroidCone)   = [normalization_factor(r) .* Vector(r) for r in Oscar.rays(cone(H))]
facets(H::PolymatroidCone) = [h.a for h in Oscar.facets(cone(H))]
dim(H::PolymatroidCone)    = Oscar.dim(cone(H))

# Given coprime integer coordinates, a polymatroid is a matroid
# if and only if it is bounded by the cardinality function.
is_matroid(H::PolymatroidCone, v::Vector) = all(A -> v[getvar(H, A)] <= length(A), H.vars)

end

# vim: set expandtab ts=4 sts=-1 sw=4 tw=0:
