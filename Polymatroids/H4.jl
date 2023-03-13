push!(LOAD_PATH, "src")
using Polymatroids

H = PolymatroidCone(4)
println("H4 is a cone of dimension ", dim(H), ".")
println(length(facets(H)), " facets:")
for h in sort(facets(H))
    println("[", join(h, ", "), "]")
end
print("\n")

println(length(rays(H)), " extreme rays:")
for r in sort(rays(H))
    print("[", join(r, ", "), "]")
    if is_matroid(H, r)
    	print(" (matroid)")
    end
    print("\n")
end

# vim: set expandtab ts=4 sts=-1 sw=4 tw=0:
