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

    batch::Bool
end


function import_variables(path::String)
    include(path)
    DATA_DIR = joinpath(@__DIR__, "data")
    return PipelineVariables(name,
                            joinpath(DATA_DIR, data_filename),
                            joinpath(DATA_DIR, newick_filename),
                            joinpath(@__DIR__, instructions_filename),
                            joinpath(DATA_DIR, labels_filename),
                            MAKE_SELECTION_XML,
                            RUN_SELECTION_XML,
                            RUN_SELECTION_XML,
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
                            [SELECTION_STATISTIC],
                            KEEP_THRESHOLD,
                            PLOT_BURNIN,
                            JULIA_SEED,
                            BEAST_SEED,
                            CONSTRAIN_LOADINGS,
                            SHRINK_LOADINGS,
                            false
                            )
end


function import_variables()
    return import_variables(joinpath(@__DIR__, "instructions.jl"))
end


end