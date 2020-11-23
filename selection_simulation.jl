import Random
Random.seed!(666)

# using JLD, BatchSetup


# function submit_pipeline(starting_dir::String, jld_path:String)
#     cd(starting_dir)



################################################################################
## meta-parameters
################################################################################


const OVERWRITE = true
const BATCH = false
const SPARSE_LOADINGS = true


const JLD_NAME = "vars.jld"

cd(@__DIR__)
push!(LOAD_PATH, @__DIR__)

using Simulations






################################################################################
## simulate data
################################################################################



ks = [2, 4]
ps = [80]
ns = [150]
n_sims = 2
sim_name = "diff3"
status = check_status(joinpath(Simulations.SIM_DIRECTORY, sim_name), batch = BATCH)

this_dir = simulate_data(sim_name, ns, ks, ps, n_sims, overwrite = OVERWRITE,
                         status = status, sparse_loadings = SPARSE_LOADINGS)



################################################################################
## Model selection under different selection statistics
################################################################################

vars_dict = run_sim_pipelines(this_dir; status = status)



################################################################################
## Evaluate performance of selection statistics
################################################################################

if status != Simulations.MAKE_XML

    df = process_simulations(this_dir, vars_dict, sim_name)
end