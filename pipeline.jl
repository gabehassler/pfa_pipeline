push!(LOAD_PATH, @__DIR__)


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

using PipelineFunctions, ImportVariables

try
    for path in instructions_paths
        all_vars = import_variables(path)
        for vars in all_vars
            run_pipeline(vars)
        end
    end
catch e
    @error "Something went wrong" exception=(e, catch_backtrace())
    cd(@__DIR__)
end

cd(@__DIR__)