using BeastUtils.XMLConstructor, BeastUtils.DataStorage, BeastUtils.MatrixUtils,
        BeastUtils.RunBeast, BeastUtils.Logs, BeastUtils.BeastNames

push!(LOAD_PATH, @__DIR__)
using SpecialSVDProcessing

using CSV, DataFrames, LinearAlgebra, Statistics, UnPack

cd(@__DIR__)
include("instructions.jl")

import Random
if (JULIA_SEED != -1)
    Random.seed!(JULIA_SEED)
end

const FIX_GLOBAL = false
const FIX_FIRST = false
const BASE_SHAPE = 2.0

const LPD_STAT = "LPD"
const MSE_STAT = "MSE"
const NO_STAT = ""
partition_dict = Dict(LPD_STAT => true, MSE_STAT => false)
label_dict = Dict(LPD_STAT => "removed.", MSE_STAT => "traitValidation.TotalSum")
mult_dict = Dict(LPD_STAT => 1, MSE_STAT => -1)
const PARTITION_MISSING = partition_dict[SELECTION_STATISTIC]
const STAT_LABEL = label_dict[SELECTION_STATISTIC]
const STAT_MULT = mult_dict[SELECTION_STATISTIC]


mutable struct XMLRun
    newick::String
    taxa::Vector{String}
    data::Matrix{Float64}
    missing_data::Union{Nothing, Matrix{Float64}}
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
    selection_stat::String
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
    return XMLRun(newick, taxa, data, nothing, k, false, NaN, NaN,
                  0, name, "", 100, 10, default_loadings(k, p), LPD_STAT)
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
    return XMLRun(run.newick, run.taxa, run.data, run.missing_data, run.k,
                    run.shrink, run.shape_multiplier, run.scale_multiplier,
                    run.rep, run.base_name, run.filename, run.chain_length,
                    run.file_freq, run.L_init, run.selection_stat)
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


function make_xml(run::XMLRun, dir::String; standardize::Bool = false, log_factors::Bool = false)

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
    selection_stat = run.selection_stat

    shapes = fill(mult_shape, k)
    shapes[1] = BASE_SHAPE
    scales = fill(mult_scale, k)

    bx = XMLConstructor.make_orthogonal_pfa_xml(data, taxa, newick, k,
                                     chain_length = chain_length,
                                     sle = SLE,
                                     fle = fle,
                                     log_factors = log_factors)


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

    if !isnothing(removed_data)
        XMLConstructor.add_trait!(bx, removed_data, "removed")
    end

    ops = XMLConstructor.get_operators(bx)

    lgo = ops[findfirst(x -> isa(x, XMLConstructor.HMCOperatorXMLElement), ops)]
    lgo.weight = Float64(3)
    if CONSTRAIN_LOADINGS
        error("not implemented")
        sparsity_constraint = CONSTRAIN_LOADINGS ? "hybrid" : "none"
        XMLConstructor.sparsity_constraint!(lgo, sparsity_constraint)
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


    if !isnothing(removed_data) && selection_stat != NO_STAT
        if selection_stat == LPD_STAT
            flpd = XMLConstructor.FactorLogPredictiveDensity(facs, like)
            XMLConstructor.add_loggable(bx, flpd, already_made = false)
        elseif selection_stat == MSE_STAT
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
    if !isdir(dir)
        mkdir(dir)
    end
end

function safe_csvwrite(path::String, df::DataFrame; overwrite::Bool = true)
    if isfile(path) && !overwrite
        error("Cannot replace file '$path' with `overwrite` set to `false`.")
    end
    CSV.write(path, df)
end


cd(@__DIR__)
data_path = joinpath(@__DIR__, "data", data_filename)
newick_path = joinpath(@__DIR__, "data", newick_filename)
instructions_path = joinpath(@__DIR__, instructions_filename)

safe_mkdir(name)

cd(name)



newick = read(newick_path, String)
taxa, data = DataStorage.csv_to_data(data_path)
const N = length(taxa)
const P = size(data, 2)

################################################################################
## Step 1 - Determine the appropriate shrinkage
################################################################################


## Read instructions and stup XML files
mod_select_name = "selection"
select_dir = "selection"
safe_mkdir(select_dir)
select_xml_dir = joinpath(select_dir, "xml")
select_log_dir = joinpath(select_dir, "logs")
safe_mkdir(select_xml_dir)
safe_mkdir(select_log_dir)



df = CSV.read(instructions_path)
n_opts = size(df, 1)
selection_data = copy(data)
standardize_data!(selection_data)

xmlruns = [XMLRun(name * "_" * mod_select_name, newick, taxa, selection_data)
            for i = 1:n_opts, j = 1:REPEATS]

# sv_threshold = 1e-1

for i = 1:n_opts
    k_max = df.k[i]
    shape_exp = df.shape_exp[i]
    shape = 10.0^shape_exp
    # scale_exp = df.scale_exp[i]
    # scale = 10.0^scale_exp
    scale = 1.0
    chain_length = df.chain_length[i]
    fle = df.file_freq[i]
    for j = 1:REPEATS
        run = xmlruns[i, j]
        run.k = k_max
        run.L_init = default_loadings(k_max, P)
        run.shrink = true
        run.shape_multiplier = shape
        run.scale_multiplier = scale
        run.rep = j
        run.chain_length = chain_length
        run.file_freq = fle
        run.selection_stat = SELECTION_STATISTIC

        name_run!(run)
    end
end

if MAKE_SELECTION_XML

for j = 1:REPEATS
    observed_data, missing_data = remove_observations(selection_data, SPARSITY)

    for i = 1:n_opts
        run = xmlruns[i, j]
        run.data = observed_data
        run.missing_data = missing_data
    end

end

for run in xmlruns
    make_xml(run, select_xml_dir)
end

end

## Run XML files

lpd_path = joinpath(select_dir, "lpd.csv")



if RUN_SELECTION_XML
    nms = [Symbol("run$i") for i = 1:REPEATS]
    lpd_df = DataFrame(zeros(n_opts, REPEATS), nms)
    safe_csvwrite(lpd_path, lpd_df, overwrite = OVERWRITE)
    for j = 1:REPEATS
        for i = 1:n_opts
            run = xmlruns[i, j]
            xml_path = joinpath(select_xml_dir, "$(run.filename).xml")
            RunBeast.run_beast(xml_path, seed = BEAST_SEED,
                               overwrite = OVERWRITE,
                               beast_jar = joinpath(BEAST_HOME, "beast.jar"))
            log_name = "$(run.filename).log"

            col, log_data = Logs.get_log_match(log_name, STAT_LABEL, burnin = SELECTION_BURNIN)
            @assert length(col) == 1
            lpd_df[i, j] = mean(log_data)
            safe_csvwrite(lpd_path, lpd_df, overwrite = true)

            mv(log_name,
                joinpath(select_log_dir, "$(run.filename).log"),
                force = OVERWRITE)
        end
    end
end


## Find "best" shrinkage

lpds = Matrix{Float64}(CSV.read(lpd_path))
best_ind = findmax(vec(STAT_MULT * mean(lpds, dims = 2)))[2]
display(lpds)
display(best_ind)

################################################################################
## Actual run
################################################################################

final_run = copy(xmlruns[best_ind, 1])
final_run.rep = 0
final_run.data = data
final_run.missing_data = nothing
final_run.filename = name
final_run.chain_length = FINAL_CHAIN_LENGTH
final_run.file_freq = FINAL_FILE_FREQUENCY

final_filename = final_run.filename * ".xml"

if MAKE_FINAL_XML
    xml = make_xml(final_run, pwd(); standardize = true, log_factors = true)
end
if RUN_FINAL_XML
    RunBeast.run_beast(final_run.filename * ".xml", seed = BEAST_SEED,
                       overwrite = OVERWRITE,
                       beast_jar = joinpath(BEAST_HOME, "beast.jar"))
end

fn = final_run.filename
svd_path = "$(fn)_processed.log"

start_ind = CONSTRAIN_LOADINGS ? 2 : 1

SpecialSVDProcessing.svd_logs("$fn.log", svd_path, final_run.k, size(data, 2),
                              rotate_factors = true,
                              do_svd = false,
                              relevant_rows = collect(start_ind:final_run.k),
                              relevant_cols = collect(start_ind:size(data, 2)))



################################################################################
## Plot results
################################################################################

using BeastUtils.Logs

using RCall, Statistics, DataFrames

metadata_path = joinpath(@__DIR__, "data", labels_filename)


function process_log(log_path::String, csv_path::String, data_path::String,
                     labels_path::String; burnin = 0.1)

    cols, data = get_log(log_path, burnin= PLOT_BURNIN)
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
    upper_threshold = Int(floor(KEEP_THRESHOLD * n))
    lower_threshold = Int(ceil((1.0 - KEEP_THRESHOLD) * n))

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
            if n_pos > upper_threshold || n_pos < lower_threshold
                L[i, perm[j]] = mean(@view L_data[:, col])
            end
        end
    end

    m_count = count(x -> ismissing(x), L)
    @show m_count / (k * p)

    row_counts = [count(x -> !ismissing(x), @view L[i, :]) for i = 1:k]
    col_counts = [count(x -> !ismissing(x), @view L[:, j]) for j = 1:p]

    last_row = 1
    for i = 2:k
        if row_counts[i] > 1
            last_row = i
        end
    end

    @show row_counts
    @show last_row

    keep_rows = 1:last_row

    df = DataFrame()
    df.L = vec(L[keep_rows, :])
    df.col = repeat(1:p, inner=last(keep_rows))
    df.row = repeat(keep_rows, outer=p)

    trait_names = string.(names(CSV.read(data_path))[2:end])

    @assert length(trait_names) == p

    df.trait = repeat(new_names, inner=last(keep_rows))
    df.cat = repeat(trait_types, inner=last(keep_rows))
    safe_csvwrite(csv_path, df, overwrite = OVERWRITE)

    levs = df.trait
    return last_row
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

if PLOT_LOADINGS

# cd(@__DIR__)
# nms = string.(names(CSV.read(data_path)))

labels_df = CSV.read(metadata_path)
cat_levs = unique(labels_df.cat)

nm = final_run.filename
# final_log = nm * ".log"
final_log = svd_path
csv_path = "$nm.csv"
k_effective = process_log(final_log, csv_path, data_path, metadata_path, burnin = PLOT_BURNIN)
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

R"""
pretty_names <- read.csv(tmp_path)$levs
plot_loadings(csv_path, plot_path, pretty_names, cat_levs)
"""

rm(tmp_path)
end

taxa, F = process_for_factors(final_log, name * "_factors.txt");
# @show maximum(abs.(F[:, 1]))
# @show maximum(abs.(F[:, 2]))
cd(@__DIR__)