module Simulations

export simulate_data, run_sim_pipelines, process_simulations,
       check_status

using Distributions, LinearAlgebra, PhyloNetworks, DataFrames, CSV
using BeastUtils.Simulation, BeastUtils.TreeUtils, BeastUtils.DataStorage
using PipelineFunctions, ImportVariables, BatchSetup

const INSTRUCTIONS = Dict{Bool, String}(
                        true => joinpath(@__DIR__, "sim_shrink.csv"),
                        false => joinpath(@__DIR__, "sim_noShrink.csv")
                                        )

const SIM_DIRECTORY = joinpath(@__DIR__, "simulations")
const DATA_NAME = "data.csv"
const NEWICK_NAME = "newick.txt"
const METADATA_NAME = "metadata.csv"
const STATS = ["MLPD", "CLPD"]
const SHRINKS = [true, false]
const BATCH = false
const BATCH_DIR = joinpath(@__DIR__, "batch")

const ALL = 0
const MAKE_XML = 1
const PROCESS_LOGS = 2
const MAKE_DF = 3
const DONE = 4

const BASE_NAME = "X"
const RUN_NAMES = ["X_shrink$s" for s in SHRINKS]

beast_path = joinpath(@__DIR__, "beast.jar")
beast_sle = 100
final_chain_length = 20000
final_file_freq = 20

repeats = 1

sparsity = 0.1
selection_burnin = 0.5

keep_threshold = 0.9
plot_burnin = 0.9

julia_seed = 666
beast_seed = 666

function make_batch(this_dir::String)
    bis = BatchInfo[]

    for dir in readdir(this_dir, join = true)
        if isdir(dir)
            dirs = ["\$HOME", "pfa_pipeline", "simulations",
                    splitpath(this_dir)[end],
                    splitpath(dir)[end],
                    RUN_NAME,
                    "selection"]

            xml_dir = join([dirs; "xml"], "/")
            logs_dir = join([dirs; "logs"], "/")
            for xml in readdir(joinpath(this_dir, dir, RUN_NAME, "selection", "xml"))
                nm = string(split(xml, '.')[1])
                bi = BatchInfo(nm)
                bi.name = splitpath(dir)[end] * "_" * bi.name
                bi.run_time = "4:00:00"
                bi.source_dir = xml_dir
                bi.dest_dir = logs_dir
                bi.email = false
                bi.tmp_dir = true

                push!(bis, bi)
            end
        end
    end

    setup_sh(this_dir, bis)
end

function check_status(this_dir::String; batch::Bool = false)

    status = DONE

    if !isdir(this_dir) || length(readdir(this_dir)) == 0
        return return_status(MAKE_XML, batch)
    end

    if isfile(joinpath(this_dir, "results.csv"))
        return DONE
    end

    for dir in readdir(this_dir, join = true)
        if isdir(dir)
            dir_name = splitpath(dir)[end]

            for i = 1:length(RUN_NAMES)
                run_dir = joinpath(dir, RUN_NAMES[i])
                for file in readdir(run_dir)
                    if endswith(file, ".log")
                        return MAKE_DF
                    end
                end

                selection_dir = joinpath(run_dir, "selection")
                xml_dir = joinpath(selection_dir, "xml")
                logs_dir = joinpath(selection_dir, "logs")

                if !isdir(xml_dir) || isempty(xml_dir)
                    return return_status(MAKE_XML, batch)
                end

                if isdir(logs_dir)
                    if length(readdir(logs_dir)) != 0
                        return return_status(PROCESS_LOGS, batch)
                    else
                        return return_status(MAKE_XML, batch)
                    end
                end
            end
        end
    end

    error("unknown condition")

end

function return_status(status::Int, batch::Bool)
    if !batch
        if !(status == ALL || status == DONE || status == MAKE_DF)
            status = ALL
        end
    end

    return status
end


function parse_initial(file::AbstractString, base_name::String)
    reg_string = "$(base_name)_n(\\d+)_k(\\d+)_p(\\d+)_s.*"
    regex = Regex(reg_string)
    m = match(regex, file)

    @show file

    return parse(Int, m.captures[1]),
           parse(Int, m.captures[2]),
           parse(Int, m.captures[3])
end


function make_sim_name(name::String, n::Int, k::Int, p::Int, rep::Int)
    return name * "_n$n" * "_k$k" * "_p$p" * "_s$rep"
end


function random_orthonormal(k::Int, p::Int;
                         dense_dim::Int = k,
                         sparse_count::Int = dense_dim + 1)

    sparse_dim = k - dense_dim
    L = zeros(k, p)
    L[1:dense_dim, :] .= svd(randn(dense_dim, p)).Vt

    total_sparse = sparse_count * sparse_dim
    if total_sparse > sparse_dim * p
        error("cannot create sparse orthogonal matrix with " *
              "sparse_count=$sparse_count, dense_dim=$dense_dim, and p=$p.")
    end

    sparse_inds = sample(1:p, total_sparse, replace=false)
    sparse_inds = reshape(sparse_inds, sparse_count, sparse_dim)

    for i = 1:sparse_dim
        cols = @view sparse_inds[:, i]
        n = nullspace(L[1:dense_dim, cols])

        L[i + dense_dim, cols] .= vec(n)
    end

    return L
end



function simulate_data(sim_name::String,
                        ns::Vector{Int}, ks::Vector{Int}, ps::Vector{Int},
                        n_sims::Int;
                        shrinkage_shape::Float64 = 2.0,
                        shrinkage_scale::Float64 = 1.0,
                        res_shape::Float64 = 2.0,
                        res_scale::Float64 = 1.0,
                        overwrite::Bool = false,
                        status::Int = ALL,
                        sparse_loadings::Bool = false)

    λ_dist = Gamma(res_shape, res_scale)
    mults_dist = Gamma(shrinkage_shape, shrinkage_scale)

    this_dir = joinpath(SIM_DIRECTORY, sim_name)
    mkpath(this_dir)


    if !(status == ALL || status == DONE || status == MAKE_XML)
        return this_dir
    end

    if !overwrite && isdir(this_dir) && !(length(readdir(this_dir)) == 0)
        error("Directory $this_dir is not empty. Set OVERWRITE=true to overwrite.")
    end


    for file in readdir(this_dir)
        rm(joinpath(this_dir, file), recursive=true)
    end


    for n in ns
        taxa = ["taxon$i" for i = 1:n]
        for k in ks
            for p in ps
                for i = 1:n_sims
                    dense_dim = sparse_loadings ? div(k, 2) : k
                    Vt = random_orthonormal(k, p, dense_dim=dense_dim)
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



function run_sim_pipelines(this_dir::String;
                           status::Int = NONE,
                           compare_shrinkage::Bool = true)
    old_dir = pwd()
    vars_dict = Dict{String, Vector{PipelineVariables}}()

    make_selection = true
    run_selection = true
    process_logs = true
    do_final = true

    if status == MAKE_XML
        make_selection = true
        run_selection = false
        process_logs = false
        do_final = false
    elseif status == PROCESS_LOGS
        make_selection = false
        run_selection = false
        process_logs = true
        do_final = true
    end

    for dir in readdir(this_dir, join = true)
        if isdir(dir)
            data_path = joinpath(dir, DATA_NAME)
            newick_path = joinpath(dir, NEWICK_NAME)
            metadata_path = joinpath(dir, METADATA_NAME)
            cd(dir)

            n = length(SHRINKS)

            vars_dict[dir] = Vector{PipelineFunctions.PipelineVariables}(undef, n)


            for i = 1:n
                shrink = SHRINKS[i]
                instructions = INSTRUCTIONS[shrink]

                vars = PipelineVariables(RUN_NAMES[i],
                                        data_path,
                                        newick_path,
                                        instructions,
                                        metadata_path,
                                        make_selection,
                                        run_selection,
                                        process_logs,
                                        do_final,
                                        do_final,
                                        false,
                                        false,
                                        true,
                                        beast_path,
                                        beast_sle,
                                        final_chain_length,
                                        final_file_freq,
                                        0,
                                        repeats,
                                        sparsity,
                                        selection_burnin,
                                        STATS,
                                        keep_threshold,
                                        plot_burnin,
                                        julia_seed,
                                        beast_seed,
                                        false,
                                        shrink,
                                        true,
                                        false
                                        )
                vars_dict[dir][i] = vars

                if status != MAKE_DF
                    run_pipeline(vars)
                end
            end
        end
    end

    if status == MAKE_XML
        make_batch(this_dir)
    end

    cd(old_dir)
    return vars_dict
end

function process_simulations(this_dir::String,
                    vars_dict::Dict{String, Vector{PipelineVariables}},
                    sim_name::String)
    n_sims = count(isdir(dir) for dir in readdir(this_dir, join = true))
    df = DataFrame([String, Bool, Int, Int, Int, Int],
                   ["stat", "shrink", "n", "k", "p", "k_inferred"],
                   n_sims * length(STATS) * length(SHRINKS))


    ind = 1

    for dir in readdir(this_dir, join = true)
        if isdir(dir)
            (n, k, p) = parse_initial(dir, sim_name)
            cd(dir)

            all_vars = vars_dict[dir]
            for vars in all_vars
                vars.make_selection_xml = false
                vars.run_selection_xml = false
                vars.make_final_xml = false
                vars.run_final_xml = false
                vars.compute_k = true

                eff_ks = run_pipeline(vars)

                n_stats = length(STATS)
                rng = ind:(ind + n_stats - 1)

                for i = 1:length(rng)
                    df.stat[rng[i]] = STATS[i]
                    df.k_inferred[rng[i]] = eff_ks[i]
                end

                df.n[rng] .= n
                df.k[rng] .= k
                df.p[rng] .= p
                df.shrink[rng] .= vars.shrink
                ind += n_stats
            end
        end
    end

    CSV.write(joinpath(this_dir, "results.csv"), df)
    return df
end

end