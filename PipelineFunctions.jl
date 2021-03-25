module PipelineFunctions

export run_pipeline, safe_mkdir

using BEASTXMLConstructor, BEASTXMLConstructor.BeastNames,
      BeastUtils.DataStorage, BeastUtils.MatrixUtils,
      BeastUtils.RunBeast, BeastUtils.Logs, PosteriorSummary

push!(LOAD_PATH, @__DIR__)
using SpecialSVDProcessing, ImportVariables

using CSV, DataFrames, LinearAlgebra, Statistics, UnPack, RCall
import Random


const FIX_GLOBAL = false
const FIX_FIRST = true
const INITIALZE_SHRINKAGE = true
const BASE_SHAPE = 1.0
const BASE_SCALE = 1.0
const SHAPE = "shape"
const SCALE = "scale"
const BOTH = "both"
const SCALE_BY = SHAPE
const SCALE_BASES = true
const K_FOLD = false

# const LPD_STAT = "LPD"
const MSE_STAT = "MSE"
const LPD_COND = "CLPD"
const LPD_MARG = "MLPD"
const NO_STAT = ""
const JOINT = "joint"
const REMOVED = "removed"
const partition_dict = Dict(LPD_MARG => true, LPD_COND => false, MSE_STAT => false)
const trait_dict = Dict(LPD_MARG => REMOVED, LPD_COND => JOINT, MSE_STAT => JOINT)
const label_dict = Dict(LPD_MARG => "$(trait_dict[LPD_MARG]).", LPD_COND => "$(trait_dict[LPD_COND]).", MSE_STAT => "traitValidation.TotalSum")
const mult_dict = Dict(LPD_MARG => 1, LPD_COND => 1, MSE_STAT => -1)

const INIT_SHRINK_WITH_SVD = true
const SAMPLED_HMC = false



mutable struct XMLRun
    newick::String
    taxa::Vector{String}
    data::Matrix{Float64}
    missing_data::Union{Nothing, Matrix{Float64}}
    joint_data::Matrix{Float64}
    k::Int
    shrink::Bool
    shape_multiplier::Float64
    scale_multiplier::Float64
    rep::Int
    base_name::String
    filename::String
    chain_length::Int
    file_freq::Int
    L_init::Matrix{Float64}
    selection_stats::Vector{String}
end

function XMLRun(name::String, newick::String, taxa::Vector{String}, data::Matrix{Float64})
    return XMLRun(name, newick, taxa, data, 0)
end

function XMLRun(name::String, newick::String, taxa::Vector{String}, data::Matrix{Float64},
                k::Int)
    n, p = size(data)
    if length(taxa) != n
        error("The length of the taxa vector must equal the number of rows " *
                "in the data matrix.")
    end
    return XMLRun(newick, taxa, data, nothing, data, k, false, NaN, NaN,
                  0, name, "", 100, 10, default_loadings(k, p), [MSE_STAT])
end

function default_loadings(k::Int, p::Int)
    x = zeros(k, p)
    for i = 1:min(k, p)
        x[i, i] = 1.0
    end
    return x
end

import Base.copy
function Base.copy(run::XMLRun)
    return XMLRun(run.newick, run.taxa, run.data, run.missing_data,
                    run.joint_data,
                    run.k,
                    run.shrink, run.shape_multiplier, run.scale_multiplier,
                    run.rep, run.base_name, run.filename, run.chain_length,
                    run.file_freq, run.L_init, run.selection_stats)
end

function name_run!(run::XMLRun)
    @unpack shrink, rep, base_name = run
    nm = [base_name]
    if shrink
        @unpack scale_multiplier, shape_multiplier = run
        scale_string = float2string(log10(scale_multiplier))
        shape_string = float2string(log10(shape_multiplier))
        push!(nm, "scale$(scale_string)_shape$shape_string")
    # else

    end
    @unpack k = run
    push!(nm, "k$k")
    push!(nm, "r$rep")
    run.filename = join(nm, '_')
    return run.filename
end


function float2string(x::Float64; keep_decimals::Int = 1)
    s = string(x)
    s1, s2 = split(s, '.')
    if keep_decimals > 0
        return "$(s1)p$(s2[1:keep_decimals])"
    else
        return s1
    end
end


function make_xml(run::XMLRun, vars::PipelineVariables, dir::String;
                  standardize::Bool = false, log_factors::Bool = false)

    taxa = run.taxa
    data = run.data
    removed_data = run.missing_data
    newick = run.newick
    shrink = run.shrink
    chain_length = run.chain_length
    fle = run.file_freq
    k = run.k
    L_init = run.L_init
    mult_shape = run.shape_multiplier
    mult_scale = run.scale_multiplier
    selection_stats = run.selection_stats

    shapes = fill(mult_shape, k)
    shapes[1] = BASE_SHAPE
    scales = fill(mult_scale, k)
    scales[1] = BASE_SCALE

    if SCALE_BASES
        trV = sum(tr(MatrixUtils.missing_cov(data)[1]))
        if SCALE_BY == SHAPE
            shapes[1] = 1.0 / trV
        elseif SCALE_BY == SCALE
            scales[1] = 1.0 / trV
        elseif SCALE_BY == BOTH
            v = sqrt(1.0 / trV)
            shapes[1] = v
            scales[1] = v
        end
    end

    if shrink
        if SAMPLED_HMC
            bx = BEASTXMLConstructor.make_sampled_pfa_xml(data, taxa, newick, k,
                                            chain_length = chain_length,
                                            log_factors = log_factors,
                                            timing = true,
                                            force_ordered = true,
                                            hmc_scale=false,
                                            scale_together=false,
                                            always_draw_factors=false)
            set_screen_logEvery(bx, vars.beast_sle)
            set_file_logEvery(bx, fle)
        else

            bx = BEASTXMLConstructor.make_orthogonal_pfa_xml(data, taxa, newick, k,
                                            chain_length = chain_length,
                                            sle = vars.beast_sle,
                                            fle = fle,
                                            log_factors = log_factors,
                                            timing = true,
                                            force_ordered = true)
        end
    else
        bx = BEASTXMLConstructor.make_pfa_xml(data, taxa, newick, k,
                                        chain_length = chain_length,
                                        sle = vars.beast_sle,
                                        fle = fle,
                                        log_factors = log_factors,
                                        useHMC = false,
                                        timing = true)
    end

    if vars.full_eval != -1
        BEASTXMLConstructor.set_full_eval!(bx, vars.full_eval)
    end


    mbd = BEASTXMLConstructor.get_mbd(bx)
    facs = BEASTXMLConstructor.get_integratedFactorModel(bx)
    facs.standardize_traits = standardize

    BEASTXMLConstructor.set_loadings!(facs, L_init)
    if shrink
        set_scale = true
        if INITIALZE_SHRINKAGE
            set_scale = false
        end

        BEASTXMLConstructor.set_shrinkage_mults!(facs,
                                            shapes = shapes,
                                            scales = scales,
                                            set_scale = set_scale,
                                            set_mults = true)
    end
    like = BEASTXMLConstructor.get_traitLikelihood(bx)

    ind_dict = Dict{String, Int}()

    if !isnothing(removed_data)
        nms = String[]
        stat_ind = 2
        for i = 1:length(selection_stats)
            stat = selection_stats[i]
            trait_name = trait_dict[stat]
            if !(trait_name in nms)
                if trait_name == JOINT
                    BEASTXMLConstructor.add_trait!(bx, run.joint_data, trait_name)
                elseif trait_name == REMOVED
                    BEASTXMLConstructor.add_trait!(bx, removed_data, trait_name)
                else
                    error("unknown trait name")
                end

                ind_dict[trait_name] = stat_ind
                stat_ind += 1
            end
        end
    end

    ops = BEASTXMLConstructor.get_operators(bx)

    # ops[2].weight = 5.0

    lgo = BEASTXMLConstructor.get_loadings_op(bx, component = "matrix")
    lgo.weight = 3.0

    if vars.constrain_loadings
        if shrink
            p = size(data, 2)
            mask = zeros(p, k)
            mask[2:end, 2:end] .= 1.0
            lgo.mask = vec(mask)

            lgo2 = BEASTXMLConstructor.LoadingsGibbsOperatorXMLElement(facs, like)
            lgo2.sparsity = "firstRow"
            push!(ops, lgo2)
        else
            lgo.sparsity = "hybrid"
        end
    end

    if FIX_GLOBAL && shrink
        error("not implemented")
    end
    if shrink
        if FIX_GLOBAL
            mgo = findfirst(x ->
                    typeof(x) == BEASTXMLConstructor.ShrinkageScaleOperators,
                    ops)
            mgo.fix_globals = FIX_GLOBAL
        else
            indices = collect(1:k)
            if FIX_FIRST
                indices = collect(2:k)
            end
            BEASTXMLConstructor.set_muliplicative_gamma_indices(bx, indices)
        end
    end
    # if FIX_GLOBAL && shrink
    #     msop = ops[findfirst(x -> isa(x, BEASTXMLConstructor.ShrinkageScaleOperators), ops)]
    #     msop.fix_globals = true
    # end
    # if FIX_FIRST && shrink
    #     msop = ops[findfirst(x -> isa(x, BEASTXMLConstructor.ShrinkageScaleOperators), ops)]
    #     msop.fix_first = true
    # end


    if !isnothing(removed_data) && selection_stats != [NO_STAT]

        for stat in selection_stats

            if stat == LPD_COND || stat == LPD_MARG
                flpd = BEASTXMLConstructor.FactorLogPredictiveDensity(facs, like, trait_ind = ind_dict[trait_dict[stat]])
                BEASTXMLConstructor.add_loggable(bx, flpd, already_made = false)
            elseif stat == MSE_STAT
                tree_model = BEASTXMLConstructor.get_treeModel(bx)

                trait_validation =
                    BEASTXMLConstructor.TraitValidationXMLElement(tree_model, like)

                cross_validation =
                    BEASTXMLConstructor.CrossValidationXMLElement(trait_validation)

                BEASTXMLConstructor.set_validation_type!(cross_validation,
                                                    BeastNames.SQUARED_ERROR)

                BEASTXMLConstructor.set_log_sum!(cross_validation, true)

                BEASTXMLConstructor.add_loggable(bx, cross_validation,
                                            already_made=false)
            else
                error("unknown statistic")
            end
        end
    end
    path = joinpath(dir, run.filename * ".xml")
    BEASTXMLConstructor.save_xml(path, bx)
    return path
end

function remove_observations(data::Matrix{Float64}, inds::Vector{Int};
                             partition::Bool = PARTITION_MISSING)
    obs_data = copy(data)
    mis_data = copy(data)
    obs_data[inds] .= NaN
    if partition
        obs_inds = setdiff(1:length(data), inds)
        mis_data[obs_inds] .= NaN
    end

    return obs_data, mis_data
end


function safe_mkdir(dir::String)
    paths = splitpath(dir)
    n = length(paths)
    path = ""
    for i = 1:n
        path = joinpath(path, paths[i])
        if !isdir(path)
            mkdir(path)
        end
    end
end

function safe_csvwrite(path::String, df::DataFrame; overwrite::Bool = true)
    if isfile(path) && !overwrite
        error("Cannot replace file '$path' with `overwrite` set to `false`.")
    end
    CSV.write(path, df)
end


# cd(@__DIR__)
# data_path = joinpath(@__DIR__, "data", data_filename)
# newick_path = joinpath(@__DIR__, "data", newick_filename)
# instructions_path = joinpath(@__DIR__, instructions_filename)

struct TreeData
    newick::String
    taxa::Vector{String}
    data::Matrix{Float64}
    N::Int
    P::Int

    function TreeData(newick::String, taxa::Vector{String}, data::Matrix{Float64})
        n, p = size(data)
        @assert length(taxa) == n
        return new(newick, taxa, data, n, p)
    end
end





function run_pipeline(vars::PipelineVariables)
    t1 = time()
    old_dir = pwd()
    if (vars.julia_seed != -1)
        Random.seed!(vars.julia_seed)
    end

    safe_mkdir(vars.name)

    cd(vars.name)

    newick = read(vars.newick_path, String)
    taxa, data = DataStorage.csv_to_data(vars.data_path)

    td = TreeData(newick, taxa, data)

    ## model selection
    best_models = model_selection(vars, td)

    eff_ks = nothing

    if vars.make_final_xml || vars.run_final_xml || vars.plot_loadings || vars.compute_k

        ## final run
        n_models = length(best_models)
        svd_paths = Vector{String}(undef, n_models)

        if n_models == 1
            nms = [vars.name]
        else
            nms = [vars.name * stat for stat in vars.selection_statistics]
        end
        for i = 1:n_models
            svd_paths[i] = run_final_xml(vars, best_models[i], data, name = nms[i])
        end

        if (vars.plot_loadings || vars.compute_k)


            ## plot results
            eff_ks = fill(-1, n_models)
            for i = 1:n_models
                eff_ks[i] = plot_loadings(vars, best_models[i], svd_paths[i])
            end

            if length(eff_ks) == 1
                eff_ks = fill(eff_ks[1], length(vars.selection_statistics))
            end
            @show eff_ks
        end
    end
    t2 = time()
    write("time.txt", "$(t2 - t1) seconds")

    cd(old_dir)
    return eff_ks
end

function model_selection(vars::PipelineVariables, tree_data::TreeData)
    mod_select_name = "selection"
    select_dir = "selection"
    safe_mkdir(select_dir)
    select_xml_dir = joinpath(select_dir, "xml")
    select_log_dir = joinpath(select_dir, "logs")
    safe_mkdir(select_xml_dir)
    safe_mkdir(select_log_dir)


    ## Read instructions and setup xml files

    df = DataFrame(CSV.File(vars.instructions_path))
    n_opts = size(df, 1)
    selection_data = copy(tree_data.data)
    if vars.standardize
        standardize_data!(selection_data)
    end

    xmlruns = [XMLRun(vars.name * "_" * mod_select_name, tree_data.newick,
                    tree_data.taxa, selection_data)
                    for i = 1:n_opts, j = 1:vars.repeats]

    L_init = default_loadings(maximum(df.k), tree_data.P)
    if vars.shrink && INITIALZE_SHRINKAGE && vars.make_selection_xml
        init_run = XMLRun(vars.name, tree_data.newick,
                          tree_data.taxa, selection_data)

        init_run.k = maximum(df.k)
        init_run.L_init = default_loadings(init_run.k, tree_data.P)
        init_run.chain_length = 1000
        init_run.filename = vars.name * "_init"
        xml_path = make_xml(init_run, vars, select_dir)

        RunBeast.run_beast(xml_path, seed = vars.beast_seed,
                           overwrite = vars.overwrite,
                           beast_jar = vars.beast_path)
        log_file = basename(xml_path)[1:(end - 4)] * ".log"
        cols, L_all = get_log_match(log_file, "L")
        L_final = L_all[end, :]
        L_final = reshape(L_final, tree_data.P, init_run.k)'
        L_svd = svd(L_final)
        L_init .= Diagonal(L_svd.S) * L_svd.Vt
    end


    for i = 1:n_opts
        k_max = df.k[i]
        selection_exp = df.shape_exp[i]
        shape = BASE_SHAPE
        scale = BASE_SCALE

        if SCALE_BY == SHAPE
            shape = 10.0^selection_exp
        elseif SCALE_BY == SCALE
            scale = 10.0^selection_exp
        elseif SCALE_BY == BOTH
            scale = 10.0^(0.5 * selection_exp)
            shape = scale
        end

        chain_length = df.chain_length[i]
        fle = df.file_freq[i]
        for j = 1:vars.repeats
            run = xmlruns[i, j]
            run.k = k_max
            run.L_init = L_init[1:k_max, :]
            run.shrink = vars.shrink
            run.shape_multiplier = shape
            run.scale_multiplier = scale
            run.rep = j
            run.chain_length = chain_length
            run.file_freq = fle
            run.selection_stats = vars.selection_statistics

            name_run!(run)
        end
    end

    if vars.make_selection_xml
        partitions = Int[]

        n_points = length(selection_data)
        if K_FOLD
            reps = vars.repeats
            sparsity = vars.sparsity
            fold = Int(round(1.0 / sparsity))
            # if abs(sparsity * reps - 1) > 1e-10
            #     @warn "Sparsity is set to $sparsity but repeats are set " *
            #         "to $reps. Ignoring sparsity and doing $reps-fold cross " *
            #         "validation."
            # end

            fold_repeats, r = divrem(reps, fold)
            @assert r == 0


            partitions = rand(1:fold, n_points, fold_repeats)
        end

        for j = 1:vars.repeats
            if K_FOLD
                fold_i = div(j - 1, fold) + 1
                fold_val = j - fold * (fold_i - 1)

                inds = findall(isequal(j), partitions[:, fold_i])
            else
                inds = findall(x -> x < vars.sparsity, rand(n_points))
            end
            observed_data, missing_data =
                    remove_observations(selection_data,
                            inds,
                            partition = true)

            for i = 1:n_opts
                run = xmlruns[i, j]
                run.data = observed_data
                run.missing_data = missing_data
                run.joint_data = selection_data
            end

        end

        for run in xmlruns
            make_xml(run, vars, select_xml_dir)
        end
    end

    ## Run selection xml and store statistics

    statistic_paths = [joinpath(select_dir, stat * ".csv")
                            for stat in vars.selection_statistics]

    nms = [Symbol("run$i") for i = 1:vars.repeats]
    statistic_dfs = Dict(stat => DataFrame(zeros(n_opts, vars.repeats), nms)
                                for stat in vars.selection_statistics)
    # safe_csvwrite(statistic_path, statistic_df, overwrite = vars.overwrite)
    for j = 1:vars.repeats
        for i = 1:n_opts
            run = xmlruns[i, j]
            xml_path = joinpath(select_xml_dir, "$(run.filename).xml")
            log_name = "$(run.filename).log"
            log_path = joinpath(select_log_dir, "$(run.filename).log")

            if vars.run_selection_xml


                RunBeast.run_beast(xml_path, seed = vars.beast_seed,
                                overwrite = vars.overwrite,
                                beast_jar = vars.beast_path)
                mv(log_name, log_path, force = vars.overwrite)
            end

            if vars.process_selection_logs
                for stat_ind = 1:length(vars.selection_statistics)
                    stat = vars.selection_statistics[stat_ind]
                    col, log_data = Logs.get_log_match(log_path,
                                            label_dict[stat],
                                            burnin = vars.selection_burnin)

                    if stat == LPD_COND
                        col, like_data = Logs.get_log_match(log_path,
                                            "likelihood",
                                            burnin = vars.selection_burnin)
                        log_data .-= like_data
                    end

                    statistic_dfs[stat][i, j] = mean(log_data)
                    safe_csvwrite(statistic_paths[stat_ind], statistic_dfs[stat],
                                  overwrite = true)
                end
            end


        end
    end


    ## Figure out which model is best

    n_stats = length(vars.selection_statistics)
    @show statistic_paths

    for path in statistic_paths
        if !isfile(path)
            return nothing
        end
    end


    stat_means = [mean(Matrix{Float64}(DataFrame(CSV.File(path))), dims=2) for path in statistic_paths]
    best_inds = [findmax(vec(mult_dict[vars.selection_statistics[i]] *
                                        stat_means[i]))[2] for i = 1:n_stats]


    if length(Set(best_inds)) == 1
        best_inds = [best_inds[1]]
    end

    final_runs = [copy(xmlruns[ind, 1]) for ind in best_inds]
    return final_runs
end

function run_final_xml(vars::PipelineVariables, final_run::XMLRun, data::Matrix{Float64}; name::String = vars.name)
    final_run.rep = 0
    final_run.data = data
    final_run.missing_data = nothing
    final_run.filename = name
    final_run.chain_length = vars.final_chain_length
    final_run.file_freq = vars.final_file_freq

    final_filename = final_run.filename * ".xml"

    if vars.make_final_xml
        xml = make_xml(final_run, vars, pwd(); standardize = vars.standardize, log_factors = true)
    end

    if vars.run_final_xml
        RunBeast.run_beast(final_run.filename * ".xml", seed = vars.beast_seed,
                           overwrite = vars.overwrite,
                           beast_jar = vars.beast_path)
    end

    fn = final_run.filename
    svd_path = "$(fn)_processed.log"

    start_ind = vars.constrain_loadings ? 2 : 1

    log_path = "$fn.log"
    if isfile(log_path)

        SpecialSVDProcessing.svd_logs("$fn.log", svd_path, final_run.k, size(data, 2),
                                    rotate_factors = true,
                                    do_svd = !vars.shrink,
                                    relevant_rows = collect(start_ind:final_run.k),
                                    relevant_cols = collect(start_ind:size(data, 2)))
        if vars.standardize
            n, p = size(data)
            partition_variance(log_path, svd_path, n, final_run.k, p)
        end
    end

    return svd_path
end

function plot_loadings(vars::PipelineVariables, final_run::XMLRun, svd_path::String)
    metadata_path = vars.labels_path

    final_log = svd_path
    k_effective = 0


    # cd(@__DIR__)
    # nms = string.(names(CSV.read(data_path)))

    labels_df = DataFrame(CSV.File(metadata_path))
    cat_levs = unique(labels_df.cat)

    nm = final_run.filename
    # final_log = nm * ".log"
    csv_path = "$nm.csv"
    k_effective = process_log(vars, final_log, csv_path, vars.data_path,
                                metadata_path)

    if vars.plot_loadings

        # display([nms pretty_names categories])
        @rput csv_path
        plot_path = "$nm.pdf"
        @rput plot_path
        # @rput pretty_names
        @rput cat_levs
        fact = 1:k_effective
        @rput fact

        # below needed to avoid issues with '°' character for temperatures
        tmp_path = "tmp.csv"
        @rput tmp_path
        CSV.write(tmp_path, DataFrame(levs = labels_df.pretty))

        import_r_functions()

        R"""
        pretty_names <- read.csv(tmp_path, encoding="UTF-8")$levs
        plot_loadings(csv_path, plot_path, pretty_names, cat_levs)
        """

        rm(tmp_path)
    end

    taxa, F = process_for_factors(final_log, final_run.filename * "_factors.txt")
    return k_effective
end






################################################################################
## Step 1 - Determine the appropriate shrinkage
################################################################################


## Read instructions and stup XML files

# sv_threshold = 1e-1




## Run XML files





## Find "best" shrinkage



################################################################################
## Actual run
################################################################################



################################################################################
## Plot results
################################################################################



function process_log(vars::PipelineVariables, log_path::String,
                     csv_path::String, data_path::String,
                     labels_path::String; scale_by_factors::Bool = true)

    cols, data = get_log(log_path, burnin = vars.plot_burnin)
    L_header = "L"
    sv_header = "sv"
    fac_header = "factors."

    L_inds = findall(x -> startswith(x, L_header), cols)
    sv_inds = findall(x -> startswith(x, sv_header), cols)
    L_cols = cols[L_inds]
    L_data = @view data[:, L_inds]

    k = length(sv_inds)
    p, r = divrem(length(L_cols), k)
    @assert r == 0

    n = size(data, 1)
    upper_threshold = Int(floor(vars.keep_threshold * n))
    lower_threshold = Int(ceil((1.0 - vars.keep_threshold) * n))

    df = DataFrame([String, String, Int, Int, Float64, Float64, Float64, Float64],
                   ["trait", "cat", "col", "row", "L", "perc", "hpdu", "hpdl"],
                   k * p)

    labels_df = DataFrame(CSV.File(labels_path))
    labels = labels_df.label
    new_names = labels_df.pretty
    trait_types = labels_df.cat

    original_labels = string.(CSV.File(data_path).names)[2:end] # TODO: this is super inefficient, just get labels

    @assert length(labels) == length(original_labels)
    perm = indexin(original_labels, labels)

    @assert length(findall(isnothing, perm)) == 0

    row_counts = zeros(Int, k)

    stdev_adjustments = ones(n, k)

    if scale_by_factors
        fac_inds = findall(startswith(fac_header), cols)
        n_taxa, r = divrem(length(fac_inds), k)
        @assert r == 0
        F_cols = @view cols[fac_inds]
        F_data = @view data[:, fac_inds]

        for i = 1:n
            F = reshape(F_data[i, :], k, n_taxa)
            stds = std(F, dims=2, corrected=false)
            stdev_adjustments[i, :] .= vec(stds)
        end
    end



    for i = 1:k
        for j in 1:p
            ind = (i - 1) * p + j
            @assert L_cols[ind] == "$L_header$i$(j)"

            df.trait[ind] = new_names[perm[j]]
            df.cat[ind] = trait_types[perm[j]]
            df.col[ind] = perm[j]
            df.row[ind] = i

            vals = @view(L_data[:, ind]) .* @view(stdev_adjustments[:, i])
            μ = mean(vals)
            df.L[ind] = μ

            hpd = hpd_interval(vals)
            df.hpdu[ind] = hpd.upper
            df.hpdl[ind] = hpd.lower

            perc_same = count(x -> sign(μ) * x > 0.0, vals) / n
            df.perc[ind] = perc_same
            if perc_same > vars.keep_threshold
                row_counts[i] += 1
            end
        end
    end

    keep_rows = findall(x -> x > 1, row_counts)
    k_effective = length(keep_rows)

    safe_csvwrite(csv_path, df, overwrite = vars.overwrite)

    return k_effective
    # return levs[perm]
end

function process_for_factors(svd_path::String, out_path::String)
    cols, data = get_log(svd_path)
    k = length(findall(x -> startswith(x, "sv"), cols))
    fac_inds = findall(x -> startswith(x, "factors."), cols)
    fac_cols = cols[fac_inds]
    fac_means = vec(mean(data[:, fac_inds], dims = 1))
    nk = length(fac_cols)
    n, r = divrem(nk, k)
    @assert r == 0

    F = zeros(n, k)
    taxa = Vector{String}(undef, n)

    ind = 1
    for i = 1:n
        taxon = split(fac_cols[ind], '.')[2]
        taxa[i] = taxon
        for j = 1:k
            s = split(fac_cols[ind], '.')
            @assert s[2] == taxon
            @assert s[3] == "$j"
            F[i, j] = fac_means[ind]
            ind += 1
        end
    end

    df = DataFrame(taxon = taxa)

    origins = [x[1:2] for x in taxa]
    df.origin = origins
    for i = 1:k
        df[!, Symbol("f$i")] = @view F[:, i]
    end

    CSV.write(out_path, df, delim = '\t')
    return taxa, F
end

function import_r_functions()
    R"""
    library(ggplot2)
    library(wesanderson)
    plot_loadings <- function(csv_path, plot_name, trait_levs, cat_levs){
        df  <- read.csv(csv_path, header=TRUE, encoding="UTF-8")

        df$trait <- factor(df$trait, levels=trait_levs)
        df$cat <- factor(df$cat, levels=cat_levs)
        df$L <- sapply(df$L, as.numeric)
        df$sign_perc <- df$perc
        for (i in 1:length(df$perc)) {
          if (df$L[i] < 0.0) {
            df$sign_perc[i] <- 1 - df$sign_perc[i]
          }
          if (df$hpdl[i] == 0 && df$hpdu[i] == 0) {
            df$sign_perc[i] <- 0.5
          }
        }

        # df$sign <- sapply(df$L, sign)
        # n <- dim(df)[[1]]
        # for (i in 1:n) {
        #     if (df$hpdl[i] <= 0 && df$hpdu[i] >= 0) {
        #     df$sign[i] = 0
        #     }
        # }
        # df$sign <- factor(df$sign, levels=c(1, 0, -1))
        ymin = min(df$hpdl)
        ymax = max(df$hpdu)

        k <- max(df$row)
        facet_labels <- character(k)
        facet_names <- integer(k)
        for (i in 1:k) {
            facet_labels[i] = paste("factor", i)
            facet_names[i] = i
        }
        names(facet_labels) <- facet_names

        # ps <- c()
        # for (i in 1:max(df$row)) {
        #   ps <- c(ps, plot_single_row(df, 1, ymin, ymax))
        # }

        wes_pal <- wes_palette("GrandBudapest2")
        pal <- c(wes_pal[1], "grey", wes_pal[4])

        p <- ggplot(df) +
            geom_hline(yintercept=0, linetype="dashed") +
            geom_point(aes(x=trait, y=L, color=sign_perc), size=1.5) +
            geom_errorbar(aes(x=trait, ymin=hpdl, ymax=hpdu, color=sign_perc), width=0.0, size=1) +
            scale_color_gradient2(low="blue", mid="grey", high="red", midpoint=0.5, limits=c(0, 1), name="probability > 0") +
            # scale_color_gradient2(low="orange", mid="white", high="purple", limits=c(-1, 1), name="L") +
            #facet_grid(~ cat, scales="free_x", space="free_x") +
            #geom_tile() +
            #scale_fill_gradient2(low="orange", mid="white", high="purple", midpoint=0) +
            #scale_x_discrete(position = "top") +
            # scale_y_discrete() +
            labs(y="loadings", x="") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle=90, hjust=1, size=4),
                panel.border = element_rect(colour = "black", fill=NA),
                panel.grid.major.y = element_blank(),
                panel.grid.minor.y = element_blank()#,
                # legend.position = "none"
                ) +
            ylim(ymin, ymax) +
            facet_grid(cols=vars(cat),
                    rows=vars(row),
                    labeller = labeller(row=facet_labels),
                    scales="free_x",
                    space="free_x",
                    switch = "x")

        # axis.title.y = element_text())

        ggsave(plot_name, width=10, height=7, units="in")
        gc()
    }

    plot_factors <- function(csv_path, plot_name){
        df  <- read.csv(csv_path, header=TRUE)

        ggplot(df, aes(x=f1, y=f2, color=class)) +
                geom_point() +
                labs(y="f2", x="f1") +
        theme_minimal()

        ggsave(plot_name, width=6, height=4, units="in")

    gc()
    }
    """
end

function partition_variance(F::AbstractArray{Float64, 2},
                            L::AbstractArray{Float64, 2},
                            λ::AbstractArray{Float64, 1})

    k, p = size(L)
    n, kf = size(F)
    @assert k == kf
    FtF = F' * F
    LtL = L * L'
    vs_single = zeros(k)
    vs_pairwise = zeros(k, k)
    for k1 = 1:k
        vs_single[k1] = FtF[k1, k1] * LtL[k1, k1]
        for k2 = (k1 + 1):k
            vs_pairwise[k1, k2] = FtF[k1, k2] * LtL[k1, k2]
            vs_pairwise[k2, k1] = vs_pairwise[k1, k2]
        end
    end

    v_single = sum(vs_single)
    v_pairwise = sum(vs_pairwise)

    v_fac = v_single + v_pairwise
    v_total = v_fac + sum(1.0 ./ λ) * n
    perc_fac = v_fac / v_total
    percs_fac_rel = [x / v_single for x in vs_single]
    percs_fac_abs = [x / v_total for x in vs_single]

    return (factor_percent = perc_fac,
            individual_percents_abs = percs_fac_abs,
            individual_percent = v_single / v_fac,
            interaction_percent = v_pairwise / v_fac,
            individual_percents_rel = percs_fac_rel,
            )
end

function partition_variance(log_path::String, svd_path::String,
                            n::Int, k::Int, p::Int;
                            fac_start::String = "factors.",
                            load_start::String = "L",
                            prec_start::String = "factorPrecision",
                            overwrite::Bool = true)

    cols, data = get_log(svd_path)
    fac_inds = findall(x -> startswith(x, fac_start), cols)
    load_inds = findall(x -> startswith(x, load_start), cols)

    prec_cols, prec_data = get_log_match(log_path, prec_start)

    @assert length(fac_inds) == k * n
    @assert length(load_inds) == k * p
    @assert length(prec_cols) == p

    states = size(data, 1)

    new_cols = ["variancePercent.factors";
                ["variancePercent.factor$i" for i = 1:k]
                "relativeFactorContribution.marginal";
                "relativeFactorContribution.covariance";
                ["relativePercent.factor$i" for i = 1:k]
                ]

    new_data = zeros(states, length(new_cols))

    for state = 1:states
        F_state = @view data[state, fac_inds]
        L_state = @view data[state, load_inds]
        F = reshape(F_state, k, n)'
        L = reshape(L_state, p, k)'
        λ = @view prec_data[state, :]
        partition = partition_variance(F, L, λ)

        new_data[state, :] .= [partition.factor_percent;
                              partition.individual_percents_abs;
                              partition.individual_percent;
                              partition.interaction_percent;
                              partition.individual_percents_rel]
    end

    cols = overwrite ? [cols; new_cols] : new_cols
    data = overwrite ? [data new_data] : new_data

    new_path = overwrite ? svd_path :
                            log_path[1:(end - 4)] * "varianceExplained.log"
    make_log(new_path, data, cols, includes_states = true)
end



end