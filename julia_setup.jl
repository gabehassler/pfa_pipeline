using Pkg

Pkg.add(["CSV", "DataFrames", "LinearAlgebra", "Statistics", "UnPack", "RCall"])
Pkg.add(Pkg.PackageSpec(url="https://github.com/gabehassler/BeastUtils.jl.git"))
Pkg.add(Pkg.PackageSpec(url="https://github.com/Julia-BEAST-utils/PosteriorSummary.jl.git"))
