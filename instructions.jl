# Instructions for factor analysis pipeline

################################################################################
## variables you MUST specify
################################################################################

name = "sarah_mse" # this will be the folder name that the results are stored
                    # in as well as the name of the final xml and log files

data_filename = "Chameleons_nop_data.csv" # csv file where the data are stored (must be in ./data directory)
newick_filename = "tree_nop.txt" # file where the newick tree is stored (must be in ./data directory)
instructions_filename = "chameleon_selection_nop.csv" # csv file with relevant instructions for model selection
labels_filename = "Chameleons_nop_labels.csv" # stores labeling information (and order) for loadings plot


################################################################################
## variables you can optionally specify (defaults are ok for most purposes)
################################################################################

## Which parts of the pipeline get run
const MAKE_SELECTION_XML = true # create selection xml files
const RUN_SELECTION_XML = true # run model selection xml files
const MAKE_FINAL_XML = true # make final xml file
const RUN_FINAL_XML = true # run final xml file
const PLOT_LOADINGS = true # make plot summarizing the loadings matrix
const OVERWRITE = true

## BEAST-specific instructions
const BEAST_HOME = @__DIR__ # set this to the directory where your beast.jar file is located
const SLE = 1000 # frequency at which BEAST logs to screen
const FINAL_CHAIN_LENGTH = 10000 # chain length for final xml run
const FINAL_FILE_FREQUENCY = 10 # log-to-file frequency for final xml run


## Relevant files specific to this run


## Model selection variables
const REPEATS = 1 # how many cross-validation sets do you want to run
const SPARSITY = 0.1 # what proportion of data to withold for cross-validation
const SELECTION_BURNIN = 0.5 # burnin for model selection (i.e. the proportion of
                             # states that will be ignored when calculating the
                             # posterior mean predictive likelihood)

const SELECTION_STATISTIC = "MSE" # specific statistic you're trying to maximize (or minimize)
                                  # Options are:
                                  #     1) "LPD" - log predictive density
                                  #     2) "MSE" - mean squared error

## Plotting files and variables

const KEEP_THRESHOLD = 0.90 # proportion of posterior samples that must be on
                            # the same side of 0.0 to be inlcuded in plot
const PLOT_BURNIN = 0.5 # burnin for loadings plot


## Random number seeds
const JULIA_SEED = 435094328702954 # random number seed for Julia (set to -1 for a random seed)
const BEAST_SEED = 666 # the number of the BEAST (set to -1 for a random seed)

## General modeling choices

const CONSTRAIN_LOADINGS = true # set to `true` to enforce the constraint that the
                          # first trait only loads onto the first factor