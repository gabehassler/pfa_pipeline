import Random
Random.seed!(666)



################################################################################
## meta-parameters
################################################################################

const SIM_DIRECTORY = joinpath(@__DIR__, "simulations")


push!(LOAD_PATH, @__DIR__)


function sim_name(name::String, n::Int, k::Int, p::Int, rep::Int)
    return name * "_n$n" * "_k$k" * "_p$p" * "s$rep"
end

beast_path = ""
beast_sle = 1000
final_chain_length = 1000
final_file_freq = 10

repeats = 2

sparsity = 0.1
selection_burnin = 0.5

keep_threshold = 0.9
plot_burnin = 0.9

julia_seed = 666
beast_seed = 666

constrain_loadings = false


################################################################################
## simulate data
################################################################################

using Distributions, LinearAlgebra
using BeastUtils.Simulation, BeastUtils.TreeUtils, BeastUtils.DataStorage

ks = [2, 4]
ps = [10, 20]
ns = [100, 200]
n_sims = 2

λ_dist = Gamma(1, 1)
mults_dist = Gamma(2.0, 2.0)

sim_name = "sim"


for n in ns
    taxa = ["taxon$i" for i = 1:n]
    for k in ks
        for p in ps
            for i = 1:n_sims
                Vt = svd(randn(k, p)).Vt
                λ = rand(λ_dist, p)
                mults = rand(mults_dist, k)
                precs = [prod(mults[1:i]) for i = 1:k]
                sq_norms = rand(Chisq(p), k)
                norms = [sqrt(sq_norms[i] / precs[i]) for i = 1:k]

                L = Diagonal(norms) * Vt
                tree = rtree(taxa, ultrametric = true)
                standardize_height!(tree)

                lfm = LatentFactorModel(L, Diagonal(λ))

                tsm = TraitSimulationModel(taxa, tree, lfm)

                data = simulate(tsm)

                name = sim_name(sim_name, n, k, p, n_sims)

                data_path = joinpath(SIM_DIRECTORY, name, "$name.csv")
                store_data(data_path, taxa, data)
            end
        end
    end
end






# TODO


################################################################################
## Model selection under different selection statistics
################################################################################

# TODO

################################################################################
## Evaluate performance of selection statistics
################################################################################

# TODO: loadings
# TODO: residual precision