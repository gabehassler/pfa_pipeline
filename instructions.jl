# Instructions for factor analysis pipeline

################################################################################
## variables you MUST specify
################################################################################

name = "cham_trees" # this will be the folder name that the results are stored
                    # in as well as the name of the final xml and log files

data_filename = "Chamaeleonidae_data_181.csv" # csv file where the data are stored (must be in ./data directory)
newick_filename = ["Chamaeleonidae_tree_1.nwk",
                    "Chamaeleonidae_tree_2.nwk",
                    "Chamaeleonidae_tree_3.nwk",
                    "Chamaeleonidae_tree_4.nwk",
                    "Chamaeleonidae_tree_5.nwk",
                    "Chamaeleonidae_tree_6.nwk",
                    "Chamaeleonidae_tree_7.nwk",
                    "Chamaeleonidae_tree_8.nwk",
                    "Chamaeleonidae_tree_9.nwk",
                    "Chamaeleonidae_tree_10.nwk",
                    "Chamaeleonidae_tree_11.nwk",
                    "Chamaeleonidae_tree_12.nwk",
                    "Chamaeleonidae_tree_13.nwk",
                    "Chamaeleonidae_tree_14.nwk",
                    "Chamaeleonidae_tree_15.nwk"] # file where the newick tree is stored (must be in ./data directory)
instructions_filename = "dummy_selection.csv" # csv file with relevant instructions for model selection
labels_filename = "Chamaeleonidae_labels_181.csv" # stores labeling information (and order) for loadings plot


################################################################################
## variables you can optionally specify (defaults are ok for most purposes)
################################################################################

## Which parts of the pipeline get run
MAKE_SELECTION_XML = true # create selection xml files
RUN_SELECTION_XML = true # run model selection xml files
MAKE_FINAL_XML = true # make final xml file
RUN_FINAL_XML = true # run final xml file
PLOT_LOADINGS = true # make plot summarizing the loadings matrix
OVERWRITE = true


## BEAST-specific instructions
BEAST_HOME = joinpath(@__DIR__, "beast.jar") # set this to the directory where your beast.jar file is located
SLE = 100 # frequency at which BEAST logs to screen
FINAL_CHAIN_LENGTH = 10000 # chain length for final xml run
FINAL_FILE_FREQUENCY = 10 # log-to-file frequency for final xml run
LIKELIHOOD_CHECK_COUNT = -1 # how many steps in the MCMC chain will BEAST
                            # double check the likelihood.
                            # Set to -1 for default BEAST settings.


## Model selection variables
REPEATS = 1 # how many cross-validation sets do you want to run
SPARSITY = 0.1 # what proportion of data to withold for cross-validation
SELECTION_BURNIN = 0.5 # burnin for model selection (i.e. the proportion of
                             # states that will be ignored when calculating the
                             # posterior mean predictive likelihood)

SELECTION_STATISTIC = ["MSE"] # specific statistic you're trying to maximize (or minimize)
                                # Options are:
                                #     1) "CLPD" - log predictive density conditional on observed data
                                #     2) "MLPD" - marginal log predictive density
                                #     3) "MSE" - mean squared error


## Plotting files and variables
KEEP_THRESHOLD = 0.90 # proportion of posterior samples that must be on
                            # the same side of 0.0 to be inlcuded in plot
PLOT_BURNIN = 0.5 # burnin for loadings plot


## Random number seeds
JULIA_SEED = -1 # random number seed for Julia (set to -1 for a random seed)
BEAST_SEED = -1 # the number of the BEAST (set to -1 for a random seed)


## General modeling choices
CONSTRAIN_LOADINGS = true # set to `true` to enforce the constraint that the
                          # first trait only loads onto the first factor
SHRINK_LOADINGS = false # set to `true` to use a multiplicative gamma shrinkage
                        # prior on the loadings
STANDARDIZE_DATA = true # set to `true` to standardize each trait to have
                        # mean 0 and variance 1
