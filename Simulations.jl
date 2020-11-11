module Simulations

export simulate_data, run_sim_pipelines, process_simulations

using Distributions, LinearAlgebra, PhyloNetworks, DataFrames, CSV
using BeastUtils.Simulation, BeastUtils.TreeUtils, BeastUtils.DataStorage
using PipelineFunctions, ImportVariables

const INSTRUCTIONS = joinpath(@__DIR__, "sim_selection.csv")
const SIM_DIRECTORY = joinpath(@__DIR__, "simulations")
const DATA_NAME = "data.csv"
const NEWICK_NAME = "newick.txt"
const METADATA_NAME = "metadata.csv"
const STATS = ["MSE", "LPD"]
const BATCH = false


beast_path = joinpath(@__DIR__, "beast.jar")
beast_sle = 100
final_chain_length = 1000
final_file_freq = 10

repeats = 1

sparsity = 0.1
selection_burnin = 0.5

keep_threshold = 0.9
plot_burnin = 0.9

julia_seed = 666
beast_seed = 666

constrain_loadings = false



function parse_initial(file::AbstractString, base_name::String)
    reg_string = "$(base_name)_n(\\d+)_k(\\d+)_p(\\d+)_s.*"
    regex = Regex(reg_string)
    m = match(regex, file)

    return parse(Int, m.captures[1]),
           parse(Int, m.captures[2]),
           parse(Int, m.captures[3])
end


function make_sim_name(name::String, n::Int, k::Int, p::Int, rep::Int)
    return name * "_n$n" * "_k$k" * "_p$p" * "_s$rep"
end


function simulate_data(sim_name::String,
                        ns::Vector{Int}, ks::Vector{Int}, ps::Vector{Int},
                        n_sims::Int;
                        shrinkage_shape::Float64 = 2.0,
                        shrinkage_scale::Float64 = 2.0,
                        res_shape::Float64 = 2.0,
                        res_scale::Float64 = 1.0,
                        overwrite::Bool = false)

    λ_dist = Gamma(res_shape, res_scale)
    mults_dist = Gamma(shrinkage_shape, shrinkage_scale)

    this_dir = joinpath(SIM_DIRECTORY, sim_name)

    if !overwrite && isdir(this_dir) && !isempty(this_dir)
        error("Directory $this_dir is not empty. Set OVERWRITE=true to overwrite.")
    end

    mkpath(this_dir)

    for file in readdir(this_dir)
        rm(joinpath(this_dir, file), recursive=true)
    end


    for n in ns
        taxa = ["taxon$i" for i = 1:n]
        for k in ks
            for p in ps
                for i = 1:n_sims
                    Vt = svd(randn(k, p)).Vt
                    λ = 1.0 ./ rand(λ_dist, p)
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

                    name = make_sim_name(sim_name, n, k, p, i)
                    name_dir = joinpath(this_dir, name)
                    mkpath(name_dir)

                    data_path = joinpath(name_dir, DATA_NAME)
                    DataStorage.store_data(data_path, taxa, data)

                    newick_path = joinpath(name_dir, NEWICK_NAME)
                    write(newick_path, writeTopology(tree))

                    metadata_path = joinpath(name_dir, METADATA_NAME)
                    df = DataFrame()
                    traits = ["trait_$i" for i = 1:p]
                    df.label = traits
                    df.pretty = traits
                    df.cat = fill("NA", p)
                    CSV.write(metadata_path, df)
                end
            end
        end
    end
    return this_dir
end



function run_sim_pipelines(this_dir::String)
    old_dir = pwd()
    vars_dict = Dict{String, Dict{String, PipelineVariables}}()
    for dir in readdir(this_dir)
        data_path = joinpath(this_dir, dir, DATA_NAME)
        newick_path = joinpath(this_dir, dir, NEWICK_NAME)
        metadata_path = joinpath(this_dir, dir, METADATA_NAME)
        cd(joinpath(this_dir, dir))

        vars_dict[dir] = Dict{String, PipelineVariables}()

        for stat in STATS

            vars = PipelineVariables(stat,
                                    data_path,
                                    newick_path,
                                    INSTRUCTIONS,
                                    metadata_path,
                                    true,
                                    true,
                                    true,
                                    true,
                                    false,
                                    true,
                                    beast_path,
                                    beast_sle,
                                    final_chain_length,
                                    final_file_freq,
                                    repeats,
                                    sparsity,
                                    selection_burnin,
                                    stat,
                                    keep_threshold,
                                    plot_burnin,
                                    julia_seed,
                                    beast_seed,
                                    constrain_loadings,
                                    false
                                    )
            vars_dict[dir][stat] = vars
            if BATCH
                # stat_dir = joinpath(this_dir, dir, stat)
                # mkdir(stat_dir)
                # jld_path = joinpath(stat_dir, JLD_NAME)
                # save(jld_path, "vars", vars)



                # run(`julia `)

                error("not implemented")

            else
                run_pipeline(vars)
            end
        end
    end
    cd(old_dir)
    return vars_dict
end

function process_simulations(this_dir::String,
                    vars_dict::Dict{String, Dict{String, PipelineVariables}},
                    sim_name::String)
    n_sims = length(readdir(this_dir))
    df = DataFrame([String, Int, Int, Int, Int], ["stat", "n", "k", "p", "k_inferred"], n_sims * length(STATS))


    ind = 1

    for dir in readdir(this_dir)
        (n, k, p) = parse_initial(dir, sim_name)
        for stat in STATS
            cd(joinpath(this_dir, dir))

            vars = vars_dict[dir][stat]
            vars.make_selection_xml = false
            vars.run_selection_xml = false
            vars.make_final_xml = false
            vars.run_final_xml = false

            k_effective = run_pipeline(vars)

            df.stat[ind] = stat
            df.n[ind] = n
            df.k[ind] = k
            df.p[ind] = p
            df.k_inferred[ind] = k_effective
            ind += 1
        end
    end

    CSV.write(joinpath(this_dir, "results.csv"), df)
end

end