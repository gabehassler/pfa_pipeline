push!(LOAD_PATH, @__DIR__)

using PipelineFunctions, ImportVariables

try
    vars = import_variables()
    run_pipeline(vars)
catch e
    cd(@__DIR__)
    throw(e)
end

cd(@__DIR__)