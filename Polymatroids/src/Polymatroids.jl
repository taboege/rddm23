module Polymatroids

using Oscar
using Combinatorics: powerset

export PolymatroidCone, h, Δ, dim, facets, rays, is_matroid

"""
    PolymatroidCone(n::Int)

Creates an object representing the set of polymatroids on ground set `1:n`.
This object stores the ground set and establishes a linear ordering on its
powerset which is used for indexing.
"""
struct PolymatroidCone
    n::Int
    vars::Vector

    function PolymatroidCone(n::Int)
        vars = collect(powerset(1:n)) # TODO: Order fixed?
        new(n, vars)
    end
end

# Return the index of a set `A` in `H`.
function getvar(H::PolymatroidCone, A::Vector)::Int
    A = sort(unique(A))
    findfirst(B -> A == B, H.vars)
end

"""
    h(H::PolymatroidCone, A::Vector)::Vector

Return the unit vector in the ambient space of `H` where the coordinate
corresponding to the set `A` is `1` and all other coordinates are zero.
"""
function h(H::PolymatroidCone, A::Vector)::Vector
    v = [0 for x in H.vars]
    v[getvar(H, A)] = 1
    v
end

"""
    Δ(H::PolymatroidCone, i::Int, j::Int, K::Vector)::Vector

Return the submodular difference expression `h({i} ∪ K) + h({j} ∪ K) - h({i,j} ∪ K) - h(K)`.
These functionals being non-negative for all pairs `i, j` and disjoint
sets `K` implies submodularity and can therefore be used to define
the polyhedral cone of all polymatroids.
"""
function Δ(H::PolymatroidCone, i::Int, j::Int, K::Vector)::Vector
    h(H, [i, K...]) + h(H, [j, K...]) - h(H, [i, j, K...]) - h(H, K)
end

"""
    cone(H::PolymatroidCone)::Polyhedron

Return an OSCAR `Polyhedron` representation of `H`.
"""
function cone(H::PolymatroidCone)::Polyhedron
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

"""
    rays(H::PolymatroidCone)

Return the extreme rays of `H`, each with coprime integer coordinates.
"""
normalization_factor(r) = lcm([denominator(x) for x in r]...) // gcd([numerator(x) for x in r])
rays(H::PolymatroidCone) = [normalization_factor(r) .* Vector(r) for r in Oscar.rays(cone(H))]

"""
    facets(H::PolymatroidCone)

Return the facet-defining functionals of `H`.
"""
facets(H::PolymatroidCone) = [h.a for h in Oscar.facets(cone(H))]

"""
    dim(H::PolymatroidCone)

Return the dimension of `H` as a polyhedron.
"""
dim(H::PolymatroidCone) = Oscar.dim(cone(H))

"""
    is_matroid(H::PolymatroidCone, v::Vector)

Checks whether a given integer polymatroid `v` is a matroid.
This is the case if and only if `v` is bounded by the cardinality
function.

NOTE: It is not checked whether `v` is an integer-valued polymatroid.
"""
is_matroid(H::PolymatroidCone, v::Vector) = all(A -> v[getvar(H, A)] <= length(A), H.vars)

end

# vim: set expandtab ts=4 sts=-1 sw=4 tw=0:
