module ImportVariables

export PipelineVariables, import_variables

mutable struct PipelineVariables
    name::String
    data_path::String
    newick_path::String
    instructions_path::String
    labels_path::String

    make_selection_xml::Bool
    run_selection_xml::Bool
    process_selection_logs::Bool
    make_final_xml::Bool
    run_final_xml::Bool
    compute_k::Bool
    plot_loadings::Bool

    overwrite::Bool

    beast_path::String
    beast_sle::Int
    final_chain_length::Int
    final_file_freq::Int
    full_eval::Int

    repeats::Int
    sparsity::Float64
    selection_burnin::Float64

    selection_statistics::Vector{String}

    keep_threshold::Float64
    plot_burnin::Float64

    julia_seed::Int
    beast_seed::Int

    constrain_loadings::Bool
    shrink::Bool
    standardize::Bool

    batch::Bool
end


function string_or_vector(v::String)
    return [v]
end

function string_or_vector(v::AbstractArray{T}) where T <: AbstractString
    return v
end

function import_variables(path::String)
    include(path)
    selections_stats = SELECTION_STATISTIC
    if !(typeof(selections_stats) <: Vector{String})
        selections_stats = [SELECTION_STATISTIC]
    end

    names = string_or_vector(name)
    data = string_or_vector(data_filename)
    newick = string_or_vector(newick_filename)
    instructions = string_or_vector(instructions_filename)
    labels = string_or_vector(labels_filename)
    files = Dict("names" => names, "data files" => data, "newick files" => newick,
             "instruction files" => instructions, "label files" => labels)
    ks = keys(files)
    ns = [length(files[k]) for k in ks]
    max_n = maximum(ns)
    if max_n > 1
        for k in ks
            if length(files[k]) == 1
                if k == "names"
                    files[k] = [files[k][1] * "$i" for i = 1:max_n]
                else
                    files[k] = repeat(files[k], max_n)
                end
            elseif length(files[k]) != max_n
                max_inds = findall(isequal(max_n), ns)
                error("You have supplied $(length(files[k])) $k but " *
                      "$max_n " *
                      join(ks[max_inds[1:(end - 1)]], ", ") *
                      ", and $(ks[max_inds[end]]).")
            end
        end
    end


    pvs = Vector{PipelineVariables}(undef, max_n)
    DATA_DIR = joinpath(@__DIR__, "data")

    for i = 1:max_n
        pvs[i] = PipelineVariables(files["names"][i],
                                joinpath(DATA_DIR, files["data files"][i]),
                                joinpath(DATA_DIR, files["newick files"][i]),
                                joinpath(@__DIR__, files["instruction files"][i]),
                                joinpath(DATA_DIR, files["label files"][i]),
                                MAKE_SELECTION_XML,
                                RUN_SELECTION_XML,
                                MAKE_FINAL_XML,
                                MAKE_FINAL_XML,
                                RUN_FINAL_XML,
                                PLOT_LOADINGS,
                                PLOT_LOADINGS,
                                OVERWRITE,
                                BEAST_HOME,
                                SLE,
                                FINAL_CHAIN_LENGTH,
                                FINAL_FILE_FREQUENCY,
                                LIKELIHOOD_CHECK_COUNT,
                                REPEATS,
                                SPARSITY,
                                SELECTION_BURNIN,
                                selections_stats,
                                KEEP_THRESHOLD,
                                PLOT_BURNIN,
                                JULIA_SEED,
                                BEAST_SEED,
                                CONSTRAIN_LOADINGS,
                                SHRINK_LOADINGS,
                                STANDARDIZE_DATA,
                                false
                                )




    end
    return pvs
end


function import_variables()
    return import_variables(joinpath(@__DIR__, "instructions.jl"))
end


end