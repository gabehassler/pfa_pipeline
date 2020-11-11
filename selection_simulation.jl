import Random
Random.seed!(666)

# using JLD, BatchSetup

using ImportVariables

# function submit_pipeline(starting_dir::String, jld_path:String)
#     cd(starting_dir)



################################################################################
## meta-parameters
################################################################################


const OVERWRITE = true
const BATCH = false

const JLD_NAME = "vars.jld"

cd(@__DIR__)
push!(LOAD_PATH, @__DIR__)

using Simulations






################################################################################
## simulate data
################################################################################



ks = [2]
ps = [10]
ns = [100]
n_sims = 1
sim_name = "sim2"

this_dir = simulate_data(sim_name, ns, ks, ps, n_sims, overwrite = OVERWRITE)


################################################################################
## Model selection under different selection statistics
################################################################################

vars_dict = run_sim_pipelines(this_dir)



################################################################################
## Evaluate performance of selection statistics
################################################################################

df = process_simulations(this_dir, vars_dict, sim_name)