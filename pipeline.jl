push!(LOAD_PATH, @__DIR__)

using PipelineFunctions, ImportVariables

instructions_paths = String[]
if length(ARGS) > 0
    for file in ARGS
        path1 = joinpath(@__DIR__, file)
        if isfile(path1)
            push!(instructions_paths, path1)
        elseif isfile(file)
            push!(instructions_paths, file)
        else
            error("could not find instructions file $file")
        end
    end
else
    push!(instructions_paths, joinpath(@__DIR__, "instructions.jl"))
end


try
    for path in instructions_paths
        vars = import_variables(path)
        run_pipeline(vars)
    end
catch e
    cd(@__DIR__)
    throw(e)
end

cd(@__DIR__)