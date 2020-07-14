# pfa_pipeline

## Summary

The pipeline (as coded in the 'pipeline.jl' script) has the following steps:

1.	Read data from a csv file and a newick tree from a text file
2.	Create several xml (you determine exactly how many) for running model selection to determine the appropriate number of factors
3.	Run model selection xml
4.	Process log files from model selection xml and determine the "best" model
5.	Make and run final xml
6.	Process log file from final xml
7.	Plot results from final xml run.

## Computing Environment
You'll need to do a few things to set up your computing environment for this all to work (a lot of this is probably already set up, but I just want to make sure I'm explicit about everything).

1.	Make sure Java is installed (it almost certainly is if you can run BEAST).
    a.	From the command line, type `java -version`.
2.	Checkout the 'repeated_measures' branch of BEAST and build the 'beast.jar' file using Ant (in case you don't usually build BEAST this way, I've included a 'beast.jar' file in the repo that should work).
3.	Create a `BEAST_HOME` environment variable that points to the directory where you have BEAST installed (or the location of this repo if you're using the included 'beast.jar' file).
4.	Make sure you have R installed with the ggplot2 package.
5.	Install Julia 1.4 or later (https://julialang.org/downloads/).
6.	Unpack the zip file I sent you (or copy its contents into another folder).
7.	Install relevant Julia packages by running the 'julia_setup.jl' script. To do this on the command line, navigate to the unpacked zip folder and type `julia julia_setup.jl`.
    a.	One of the Julia packages lets Julia run R code. This can be tricky to install depending on how R is installed, but usually works fine. If it doesn't work, let me know and I can send you a separate R script for plotting.
8.	Make sure Julia can run BEAST by entering `julia -e "using BeastUtils.RunBeast; RunBeast.check_beast()"` into the command line. If everything is setup correctly, you should see the BEAST intro text followed by "Java and BEAST checks suceeded." If it doesn't work, please send me the error message you get.
9.	Finally, to make sure everything works, run the 'pipeline.jl' in the this repo. To do this from the command line, navigate to the folder containing the 'pipeline.jl' file and enter `julia pipeline.jl` into the command line. This should run a small example to make sure everything is working. If it works, you should see a 'mammals_example' folder appear with a 'mammals_example.pdf' containing a plot of the loadings. This should take a few minutes.

## Instructions
To run your own examples, you'll need to add the following to the 'data' folder in the unpacked zip file.

1.	A csv file containing your data. The first column in the csv file should be labeled "taxon" and have the taxa labels. It doesn't need to be standardized as the script will automatically handle that.
2.	A text file containing a newick tree
3.	A third file that maps the column names in your data file (excluding the first 'taxon' column) to the names you want to be displayed in the plot. See the 'mammals_labels.csv' file in the 'data' folder for an example. This file needs to have 3 columns labeled 'label' (the labels column labels from the original data csv file), 'pretty' (the labels you want displayed on the plot), and 'cat' (just set everything to 'NA'). Note that the order of the labels doesn't need to be the same as the ordering of the columns in the data csv file.

Because the specific parameters for the model selection vary from analysis to analysis, you'll also need to create another csv file that supplies instructions for the model selection part of the pipeline. You can look at the 'mammals_selection.csv' file as an example. It needs to have 4 columns with labels:

1.	'k' - the maximum number of factors. For smaller examples with few traits, you can just set this to the number of traits
2.	'shape_exp' - the exponent on the shrinkage parameter, which will end up determining the number of factors (hopefully). Basically we put increasingly strong priors driving each successive row of the loadings matrix closer to zero, and a bigger number here means stronger shrinkage. In theory, you can set this to be any non-negative number, but I usually have 5 rows with 'shape_exp' set to 1, 2, 3, 4, and 5.
3.	'chain_length' - How many iterations you want to run the MCMC chain for.
4.	'file_freq' - How often you want BEAST to write to the log file.

Finally, you need to set up the 'instructions.jl' file. Open this file in any text editor and replace the 'name' variable with whatever you want to. If you're running the amphibians dataset for the first time, then you could change it to "amphibians1", for example. You also have to rename the 'data_filename', 'newick_filename', 'instructions_filename', and 'labels_filename' variables to the appropriate filenames.


Once all that is done, just navigate to this repo in the command line and type `julia pipeline.jl`. For your first run, I'd recommend setting the `chain_length` in the model selection instructions to something low (like 100) just to make sure everything works without wasting a lot of time. Once you do that, you can delete the folder it just created (the default behavior is to not overwrite files you've already created), set the chain length longer, and rerun.


The log file ending in 'svd.log' is what you should use to determine if the chain ran for long enough. That won't get created until the end though. If you want to get a sense of intermediate process, you can look at any log file, just ignore any entries corresponding to the loadings since you can't get a good estimate of those until post-processing.
