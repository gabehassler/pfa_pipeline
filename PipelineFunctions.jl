module PipelineFunctions

export run_pipeline, safe_mkdir

using BeastUtils.XMLConstructor, BeastUtils.DataStorage, BeastUtils.MatrixUtils,
        BeastUtils.RunBeast, BeastUtils.Logs, BeastUtils.BeastNames

push!(LOAD_PATH, @__DIR__)
using SpecialSVDProcessing, ImportVariables

using CSV, DataFrames, LinearAlgebra, Statistics, UnPack, RCall
import Random

cd(@__DIR__)





const FIX_GLOBAL = false
const FIX_FIRST = false
const BASE_SHAPE = 2.0

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
        @unpack shape_multiplier = run
        shape_string = float2string(log10(shape_multiplier))
        push!(nm, "shape$shape_string")
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

    if shrink

        bx = XMLConstructor.make_orthogonal_pfa_xml(data, taxa, newick, k,
                                        chain_length = chain_length,
                                        sle = vars.beast_sle,
                                        fle = fle,
                                        log_factors = log_factors)
    else
        bx = XMLConstructor.make_pfa_xml(data, taxa, newick, k,
                                        chain_length = chain_length,
                                        sle = vars.beast_sle,
                                        fle = fle,
                                        log_factors = log_factors,
                                        useHMC = false)
    end

    if vars.full_eval != -1
        XMLConstructor.set_full_eval(bx, vars.full_eval)
    end


    mbd = XMLConstructor.get_mbd(bx)
    facs = XMLConstructor.get_integratedFactorModel(bx)
    facs.standardize_traits = standardize
    XMLConstructor.set_loadings!(facs, L_init)
    if shrink
        XMLConstructor.set_shrinkage_mults!(facs,
                                            shapes = shapes,
                                            scales = scales)
    end
    like = XMLConstructor.get_traitLikelihood(bx)

    ind_dict = Dict{String, Int}()

    if !isnothing(removed_data)
        nms = String[]
        stat_ind = 2
        for i = 1:length(selection_stats)
            stat = selection_stats[i]
            trait_name = trait_dict[stat]
            if !(trait_name in nms)
                if trait_name == JOINT
                    XMLConstructor.add_trait!(bx, run.joint_data, trait_name)
                elseif trait_name == REMOVED
                    XMLConstructor.add_trait!(bx, removed_data, trait_name)
                else
                    error("unknown trait name")
                end

                ind_dict[trait_name] = stat_ind
                stat_ind += 1
            end
        end
    end

    ops = XMLConstructor.get_operators(bx)

    lgo = XMLConstructor.get_loadings_op(bx)
    lgo.weight = 3.0

    if vars.constrain_loadings
        if shrink
            p = size(data, 2)
            mask = zeros(p, k)
            mask[2:end, 2:end] .= 1.0
            lgo.mask = vec(mask)

            lgo2 = XMLConstructor.LoadingsGibbsOperatorXMLElement(facs, like)
            lgo2.sparsity = "firstRow"
            push!(ops, lgo2)
        else
            lgo.sparsity = "hybrid"
        end
    end

    if FIX_GLOBAL || FIX_FIRST
        error("not implemented")
    end
    # if FIX_GLOBAL && shrink
    #     msop = ops[findfirst(x -> isa(x, XMLConstructor.ShrinkageScaleOperators), ops)]
    #     msop.fix_globals = true
    # end
    # if FIX_FIRST && shrink
    #     msop = ops[findfirst(x -> isa(x, XMLConstructor.ShrinkageScaleOperators), ops)]
    #     msop.fix_first = true
    # end


    if !isnothing(removed_data) && selection_stats != [NO_STAT]

        for stat in selection_stats

            if stat == LPD_COND || stat == LPD_MARG
                flpd = XMLConstructor.FactorLogPredictiveDensity(facs, like, trait_ind = ind_dict[trait_dict[stat]])
                XMLConstructor.add_loggable(bx, flpd, already_made = false)
            elseif stat == MSE_STAT
                tree_model = XMLConstructor.get_treeModel(bx)

                trait_validation =
                    XMLConstructor.TraitValidationXMLElement(tree_model, like)

                cross_validation =
                    XMLConstructor.CrossValidationXMLElement(trait_validation)

                XMLConstructor.set_validation_type!(cross_validation,
                                                    BeastNames.SQUARED_ERROR)

                XMLConstructor.set_log_sum!(cross_validation, true)

                XMLConstructor.add_loggable(bx, cross_validation,
                                            already_made=false)
            else
                error("unknown statistic")
            end
        end
    end
    path = joinpath(dir, run.filename * ".xml")
    # @show path
    XMLConstructor.save_xml(path, bx)
    return path
end

function remove_observations(data::Matrix{Float64}, sparsity::Float64;
                             partition::Bool = PARTITION_MISSING)
    if !(0.0 <= sparsity <= 1.0)
        error("Must specify a sparsity between 0.0 and 1.0")
    end

    obs_data = copy(data)
    mis_data = copy(data)

    for i = 1:length(data)
        if rand() < sparsity
            obs_data[i] = NaN
        elseif partition
            mis_data[i] = NaN
        end
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

    df = CSV.read(vars.instructions_path)
    n_opts = size(df, 1)
    selection_data = copy(tree_data.data)
    if vars.standardize
        standardize_data!(selection_data)
    end

    xmlruns = [XMLRun(vars.name * "_" * mod_select_name, tree_data.newick,
                    tree_data.taxa, selection_data)
                    for i = 1:n_opts, j = 1:vars.repeats]

    for i = 1:n_opts
        k_max = df.k[i]
        shape_exp = df.shape_exp[i]
        shape = 10.0^shape_exp
        # scale_exp = df.scale_exp[i]
        # scale = 10.0^scale_exp
        scale = 1.0
        chain_length = df.chain_length[i]
        fle = df.file_freq[i]
        for j = 1:vars.repeats
            run = xmlruns[i, j]
            run.k = k_max
            run.L_init = default_loadings(k_max, tree_data.P)
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
        for j = 1:vars.repeats
            observed_data, missing_data =
                    remove_observations(selection_data,
                            vars.sparsity,
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

    for path in statistic_paths
        if !isfile(path)
            return nothing
        end
    end


    stat_means = [mean(Matrix{Float64}(CSV.read(path)), dims=2) for path in statistic_paths]
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

    SpecialSVDProcessing.svd_logs("$fn.log", svd_path, final_run.k, size(data, 2),
                                  rotate_factors = true,
                                  do_svd = !vars.shrink,
                                  relevant_rows = collect(start_ind:final_run.k),
                                  relevant_cols = collect(start_ind:size(data, 2)))

    return svd_path
end

function plot_loadings(vars::PipelineVariables, final_run::XMLRun, svd_path::String)
    metadata_path = vars.labels_path

    final_log = svd_path
    k_effective = 0


    # cd(@__DIR__)
    # nms = string.(names(CSV.read(data_path)))

    labels_df = CSV.read(metadata_path)
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
        pretty_names <- read.csv(tmp_path)$levs
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
                     labels_path::String)

    cols, data = get_log(log_path, burnin = vars.plot_burnin)
    L_header = "L"
    sv_header = "sv"

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

    L = Matrix{Union{Missing, Float64}}(undef, k, p)
    fill!(L, missing)

    labels_df = CSV.read(labels_path)
    labels = labels_df.label
    new_names = labels_df.pretty
    trait_types = labels_df.cat

    original_labels = string.(names(CSV.read(data_path))[2:end]) # TODO: this is super inefficient, just get labels

    @assert length(labels) == length(original_labels)
    perm = indexin(original_labels, labels)

    @assert length(findall(isnothing, perm)) == 0


    for i = 1:k
        for j in 1:p
            col = (i - 1) * p + j
            @assert L_cols[col] == "$L_header$i$(j)"
            vals = @view(L_data[:, col])
            n_pos = count(x -> x > 0.0, vals)
            n_neg = count(x -> x < 0.0, vals)
            if n_pos > upper_threshold || n_neg > upper_threshold
                L[i, perm[j]] = mean(@view L_data[:, col])
            end
        end
    end

    m_count = count(x -> ismissing(x), L)
    @show m_count / (k * p)

    row_counts = [count(x -> !ismissing(x), @view L[i, :]) for i = 1:k]
    col_counts = [count(x -> !ismissing(x), @view L[:, j]) for j = 1:p]

    keep_rows = findall(x -> x > 1, row_counts)
    k_effective = length(keep_rows)

    df = DataFrame()
    df.L = vec(L[keep_rows, :])
    df.col = repeat(1:p, inner=k_effective)
    df.row = repeat(keep_rows, outer=p)

    trait_names = string.(names(CSV.read(data_path))[2:end])

    @assert length(trait_names) == p

    df.trait = repeat(new_names, inner=k_effective)
    df.cat = repeat(trait_types, inner=k_effective)
    safe_csvwrite(csv_path, df, overwrite = vars.overwrite)

    levs = df.trait
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
    plot_loadings <- function(csv_path, plot_name, trait_levs, cat_levs){
    df  <- read.csv(csv_path, header=TRUE)

    df$trait <- factor(df$trait, levels=trait_levs)
    df$cat <- factor(df$cat, levels=cat_levs)
    df$L <- sapply(df$L, as.numeric)
    ggplot(df, aes(x=trait, y=row, fill=L)) +
            facet_grid(~ cat, scales="free_x", space="free_x") +
            geom_tile() +
            scale_fill_gradient2(low="orange", mid="white", high="purple", midpoint=0) +
            scale_x_discrete(position = "top") +
            scale_y_discrete(limits=fact) +
      labs(y="Loadings Row", x="Trait") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle=90, hjust=0),
            strip.text.x = element_blank()
            )
            # axis.title.y = element_text())

    ggsave(plot_name, width=15, height=10, units="in")
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




end