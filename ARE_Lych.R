#######################################################
# This is a script for the R package "BioGeoBEARS" based in the script created by Nick Matzke. See the original scrip link below:
# http://phylo.wikidot.com/biogeobears#script
# 
# All scripts are copyright Nicholas J. Matzke, please cite if you use. License: GPL-3
# http://cran.r-project.org/web/licenses/GPL-3
# 
# Please, submit your questions on the BioGeoBEARS google group.
# https://groups.google.com/forum/#!forum/biogeobears
#
# The package is designed for ML and Bayesian inference of 
# 
# (a) ancestral geographic ranges, and 
# 
# (b) perhaps more importantly, models for the evolution of geographic range across a phylogeny.
#
# The script below implements and compares:
# 
# (1) The standard 2-parameter DEC model implemented in the program LAGRANGE (Ree & Smith 2008); users will notice that the ML
#     parameter inference and log-likelihoods are identical
#
# (2) A DEC+J model implemented in BioGeoBEARS, wherein a third parameter, j, is added, representing the relative per-event
#     weight of founder-event / jump speciation events at cladogenesis events.  The higher j is, the more probability these
#     events have, and the less probability the standard LAGRANGE cladogenesis events have.
#
# (3) Some standard model-testing (LRT and AIC) is implemented at the end so that users may compare models
#
# (4) The script does similar tests of a DIVA-like model (Ronquist 1997)
#     and a BAYAREA-like model (Landis, Matzke, Moore, & Huelsenbeck, 2013)
# 
#######################################################

#######################################################
# The script below is freely available and does not offer any guarantee. The use of it is entirely the responsibility of the
# user, as well as any possible eventuality resulting from its misuse.
# This code was tested using R version 3.6.3 and "BioGeoBEARS" 1.1.2, on Windows 10.

# To facilitate its use, the script has been divided into blocks and models as described below. You can search by term and 
# find the desired section easily. For the statistical calculations section look for the term "STATISTICS". You can find the plot
# section searching for the term "PLOT" and the number of the model, e.g., "PLOT 12". Note that several segments of the script 
# are enclosed in {} to speed up its execution. The script was built so that the user can make a few manual adjustments as 
# possible, however, feel free to improve it and share it with the community. I take this opportunity to leave my thanks to 
# Nicholas J. Matzke for the script provided in PhyloWiki webpage, which I used as base.

# Script author: Fábio Alves.

# SUMMARY

# BLOCK 1: DEC, DEC+X, DEC+J AND DEC+J+X ANALYSIS
# MODEL 1: DEC
# MODEL 2: DEC+X
# MODEL 3: DEC+J
# MODEL 4: DEC+J+X

# BLOCK 2: DIVALIKE, DIVALIKE+X, DIVALIKE+J AND DIVALIKE+J+X ANALYSIS
# MODEL 5: DIVALIKE
# MODEL 6: DIVALIKE+X
# MODEL 7: DIVALIKE+J
# MODEL 8: DIVALIKE+J+X

# BLOCK 3: BAYAREALIKE, BAYAREALIKE+X, BAYAREALIKE+J AND BAYAREALIKE+J+X ANALYSIS
# MODEL 9: BAYAREALIKE
# MODEL 10: BAYAREALIKE+X
# MODEL 11: BAYAREALIKE+J
# MODEL 12: BAYAREALIKE+J+X

#######################################################

#######################################################
# Installing BioGeoBEARS and required packages # More details see https://github.com/nmatzke/BioGeoBEARS
#######################################################
#install.packages(c("GenSA", "FD", "parallel", "rexpokit", "cladoRcpp", "devtools"))
#devtools::install_github(repo="nmatzke/BioGeoBEARS")
# If the BioGeoBEARS instalation keeps trying to reinstall rexpokit and cladoRcpp from scratch (including demanding a compiler
# like gcc or gfortran) even though they are already installed, try adding "dependencies=FALSE":
# devtools::install_github(repo="nmatzke/BioGeoBEARS", dependencies=FALSE)
#
#######################################################

#######################################################
# SETUP: Increase memory and load packages
#######################################################
{memory.size(2*memory.limit()) #Increased size of available memory
memory.limit(2*memory.limit()) #Increased limit of available memory
Packages <- c("GenSA", "FD", "parallel", "rexpokit", "cladoRcpp", "BioGeoBEARS")
lapply(Packages, library, character.only = T)} #Load the packages
# "GenSA" is better than "optimx" (although somewhat slower)
# use "snow" if you want to use multicore functionality. Some systems/R versions prefer "parallel", try either)

#######################################################
# SETUP: YOUR WORKING DIRECTORY
#######################################################
# You will need to set your working directory to match your local system

# Note these very handy functions!
# Command "setwd(x)" sets your working directory
# Command "getwd()" gets your working directory and tells you what it is.
# Command "list.files()" lists the files in your working directory
# To get help on any command, use "?".  E.g., "?list.files"

# Set your working directory for output files. Default here is your home directory ("~"). Change this as you like
{wd = np("./BioGeoBEARS")
setwd(wd)
getwd() # Double-check your working directory

#######################################################
# SETUP: Extension data directory
#######################################################
# When R packages contain extra files, they are stored in the "extdata" directory inside the installed package.
# BioGeoBEARS contains various example files and scripts in its extdata directory.
# Each computer operating system might install BioGeoBEARS in a different place, depending on your OS and settings. 
# However, you can find the extdata directory like this:
extdata_dir = np(system.file("extdata", package="BioGeoBEARS"))
extdata_dir
list.files(extdata_dir)

# "system.file" looks in the directory of a specified package (in this case BioGeoBEARS)
# The function "np" is just a shortcut for normalizePath(), which converts the path to the format appropriate for your system
# (e.g., Mac/Linux use "/", but Windows uses "\\", if memory serves).

# Even when using your own data files, you should KEEP these commands in your script, since the plot_BioGeoBEARS_results function
# needs a script from the extdata directory to calculate the positions of "corners" on the plot. This cannot be made into a
# straight up BioGeoBEARS function because it uses C routines from the package APE which do not pass R CMD check for some reason.

#######################################################
# SETUP: YOUR TREE FILE AND GEOGRAPHY FILE
#######################################################
# Phylogeny file
# Notes: 
# 1. Must be binary/bifurcating: no polytomies
# 2. No negative branchlengths (e.g. BEAST MCC consensus trees sometimes have negative branchlengths)
# 3. Be careful of very short branches, as BioGeoBEARS will interpret ultrashort branches as direct ancestors
# 4. You can use non-ultrametric trees, but BioGeoBEARS will interpret any tips significantly below the top of the tree as
#    fossils! This is only a good idea if you actually do have fossils in your tree, as in e.g. Wood, Matzke et al. (2013),
#    Systematic Biology.
# 5. The default settings of BioGeoBEARS make sense for trees where the branchlengths are in units of millions of years, and the
#    tree is 1-1000 units tall. If you have a tree with a total height of e.g. 0.00001, you will need to adjust e.g. the max
#    values of d and e, or (simpler) multiply all your branchlengths to get them into reasonable units.
# 6. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#######################################################

# Load the tree in Newick format. "trfn" = "tree file name"
trfn = "tree.tre"

# Look at the raw Newick file:
moref(trfn)

# Look at your phylogeny (plots to a PDF, which avoids issues with multiple graphics in same window):
pdffn = "tree.pdf"
pdf(file=pdffn, width=9, height=12)

tr = ladderize(read.tree(trfn))
tr
plot(tr)
title("Lychnophorinae phylogeny")
axisPhylo() # plots timescale

dev.off()
cmdstr = paste0("open ", pdffn)
#system(cmdstr)

#######################################################
# Geography file
# Notes:
# 1. This is a PHYLIP-formatted file. This means that in the first line, 
#    - the 1st number equals the number of rows (species)
#    - the 2nd number equals the number of columns (number of areas)
#    - after a tab, put the areas in parentheses, with spaces: (A B C D)
#
# 1.5. Example first line:
#    10    4    (A B C D)
# 
# 2. The second line, and subsequent lines:
#    speciesA    0110
#    speciesB    0111
#    speciesC    0001
#         ...
# 
# 2.5a. This means a TAB between the species name and the area 0/1s
# 2.5b. This also means NO SPACE AND NO TAB between the area 0/1s.
# 
# 3. See example files at:
#    http://phylo.wikidot.com/biogeobears#files
# 
# 4. Make you understand what a PLAIN-TEXT EDITOR is:
#    http://phylo.wikidot.com/biogeobears#texteditors
#
# 3. The PHYLIP format is the same format used for C++ LAGRANGE geography files.
#
# 4. All names in the geography file must match names in the phylogeny file.
#
# 5. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#
# 6. Operational taxonomic units (OTUs) should ideally be phylogenetic lineages, 
#    i.e. genetically isolated populations.  These may or may not be identical 
#    with species.  You would NOT want to just use specimens, as each specimen 
#    automatically can only live in 1 area, which will typically favor DEC+J 
#    models.  This is fine if the species/lineages really do live in single areas,
#    but you wouldn't want to assume this without thinking about it at least. 
#    In summary, you should collapse multiple specimens into species/lineages if 
#    data indicates they are the same genetic population.
######################################################

geogfn = "lych_geog.txt" #Load the geography file

# Look at the raw geography text file:
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))

# Set the maximum number of areas any species may occupy. This cannot be larger than the number of areas you set up, but it can
# be smaller.
max_range_size = 8
}

####################################################
# KEY HINT: The number of states (= number of different possible geographic ranges) depends on (a) the number of areas and (b)
# max_range_size. If you have more than about 500-600 states, the calculations will get REALLY slow, since the program has to
# exponentiate a matrix of e.g. 600x600. Often the computer will just sit there and crunch, and never get through the calculation
# of the first likelihood.
# 
# (this is also what is usually happening when LAGRANGE hangs: you have too many states!)
#
# To check the number of states for a given number of ranges, try:
#numstates_from_numareas(numareas=4, maxareas=4, include_null_range=TRUE)
#numstates_from_numareas(numareas=4, maxareas=4, include_null_range=FALSE)
#numstates_from_numareas(numareas=4, maxareas=3, include_null_range=TRUE)
#numstates_from_numareas(numareas=4, maxareas=2, include_null_range=TRUE)

# Large numbers of areas have problems:
#numstates_from_numareas(numareas=10, maxareas=10, include_null_range=TRUE)

# ...unless you limit the max_range_size:
#numstates_from_numareas(numareas=10, maxareas=2, include_null_range=TRUE)
####################################################

run_ARE <- function(..., DEC=F, DECX=F, DECJ=F, DECJX=F, DIVALIKE=F, DIVALIKEX=F, DIVALIKEJ=F, DIVALIKEJX=F, BAYAREALIKE=F, 
                    BAYAREALIKEX=F, BAYAREALIKEJ=F, BAYAREALIKEJX=F) {

#######################################################
# BLOCK 1: DEC, DEC+X, DEC+J AND DEC+J+X ANALYSIS
#######################################################
# NOTE: The BioGeoBEARS "DEC" model is identical with the Lagrange DEC model, and should return identical ML estimates of
# parameters, and the same log-likelihoods, for the same datasets.
#
# Ancestral state probabilities at nodes will be slightly different, since BioGeoBEARS is reporting the ancestral state
# probabilities under the global ML model, and Lagrange is reporting ancestral state probabilities after re-optimizing the
# likelihood after fixing the state at each node. These will be similar, but not identical.
# See Matzke (2014), Systematic Biology, for discussion.
#
# Also see Matzke (2014) for presentation of the founder-event speciation model (+J). https://www.jstor.org/stable/43700617
# Also see Van Dam & Matzke (2016) for presentation of the distance-dependent dispersal model (+X).
# http://dx.doi.org/10.1111/jbi.12727
#######################################################

#######################################################
# MODEL 1: DEC
#######################################################
# Intitialize a default model (DEC model)
if (DEC) {
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, Jeremy M.; Matzke, Nicholas J.; O'Meara, Brian C. 
# (2015). Non-null Effects of the Null Range in Biogeographic Models: Exploring Parameter Estimation in the DEC Model. 
# bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
# 1. Here, un-comment ONLY the files you want to use.
# 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
# 3. For example files see (a) extdata_dir, or (b) http://phylo.wikidot.com/biogeobears#files and BioGeoBEARS Google Group posts 
# for further hints)
#
# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = detectCores() - 4
# (use more cores to speed it up; this requires "library(parallel)" and/or "library(snow)". The package "parallel" is now default
# on Macs in R 3.0+, but apparently still has to be typed on some Windows machines. Note: apparently "parallel" works on Mac 
# command-line R, but not R.app. BioGeoBEARS checks for this and resets to 1 core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+). I have experimented with sparse matrix 
# exponentiation in EXPOKIT/rexpokit, but the results are imprecise and so I haven't explored it further. In a Bayesian analysis,
# it might work OK, but the ML point estimates are not identical. Also, I have not implemented all functions to work with 
# force_sparse=TRUE. Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE # get ancestral states from "optim" run

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

  resDEC = bears_optim_run(BioGeoBEARS_run_object)
  save(resDEC, file="Lychnophorinae_DEC.Rdata")
} else {
  # Loads saved analysis results
  load("Lychnophorinae_DEC.Rdata", envir = .GlobalEnv)
}

#######################################################
# MODEL 2: DEC+X
#######################################################
if (DECX) {
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = detectCores() - 4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE # get ancestral states from optim run

# Set up DEC+x model
# Add distance-dependent dispersal (+x)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0.0

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resDECx = bears_optim_run(BioGeoBEARS_run_object)
save(resDECx, file="Lychnophorinae_DEC+X.Rdata")
} else {
  # Loads saved analysis results
  load("Lychnophorinae_DEC+X.Rdata", envir = .GlobalEnv)
}

#######################################################
# MODEL 3: DEC+J
#######################################################
if (DECJ) {
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = detectCores() - 4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resDECj = bears_optim_run(BioGeoBEARS_run_object)
save(resDECj, file="Lychnophorinae_DEC+J.Rdata")
} else {
  # Loads saved analysis results
  load("Lychnophorinae_DEC+J.Rdata", envir = .GlobalEnv)
}

#######################################################
# MODEL 4: DEC+J+X
#######################################################
if (DECJX) {
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = FALSE         # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = detectCores() - 4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC+J+X model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Add distance-dependent dispersal (+x)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0.0

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resDECjx = bears_optim_run(BioGeoBEARS_run_object)
save(resDECjx, file="Lychnophorinae_DEC+J+X.Rdata")
} else {
  # Loads saved analysis results
  load("Lychnophorinae_DEC+J+X.Rdata", envir = .GlobalEnv)
}

#######################################################
# BLOCK 2: DIVALIKE, DIVALIKE+X, DIVALIKE+J AND DIVALIKE+J+X ANALYSIS
#######################################################
# NOTE: The BioGeoBEARS "DIVALIKE" model is not identical with Ronquist (1997)'s parsimony DIVA. It is a likelihood
# interpretation of DIVA, constructed by modelling DIVA's processes the way DEC does, but only allowing the processes DIVA 
# allows (widespread vicariance: yes; subset sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
#
# DIVALIKE is a likelihood interpretation of parsimony DIVA, and it is "like DIVA" -- similar to, but not identical to, 
# parsimony DIVA.
#
# I thus now call the model "DIVALIKE", and you should also. ;-)
#######################################################

#######################################################
# MODEL 5: DIVALIKE
#######################################################
if (DIVALIKE) {
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = detectCores() -4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resDIVALIKE = bears_optim_run(BioGeoBEARS_run_object)
save(resDIVALIKE, file="Lychnophorinae_DIVALIKE.Rdata")
} else {
  # Loads saved analysis results
  load("Lychnophorinae_DIVALIKE.Rdata", envir = .GlobalEnv)
}

#######################################################
# MODEL 6: DIVALIKE+X
#######################################################
if (DIVALIKEX) {
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = detectCores() -4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE+X model
# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add distance-dependent dispersal (+x)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resDIVALIKEx = bears_optim_run(BioGeoBEARS_run_object)
save(resDIVALIKEx, file="Lychnophorinae_DIVALIKE+X.Rdata")
} else {
  # Loads saved analysis results
  load("Lychnophorinae_DIVALIKE+X.Rdata", envir = .GlobalEnv)
}

#######################################################
# MODEL 7: DIVALIKE+J
#######################################################
if (DIVALIKEJ) {
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = detectCores() -4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resDIVALIKEj = bears_optim_run(BioGeoBEARS_run_object)
save(resDIVALIKEj, file="Lychnophorinae_DIVALIKE+J.Rdata")
} else {
  # Loads saved analysis results
  load("Lychnophorinae_DIVALIKE+J.Rdata", envir = .GlobalEnv)
}

#######################################################
# MODEL 8: DIVALIKE+J+X
#######################################################
if (DIVALIKEJX) {
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = FALSE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = detectCores() -4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DIVALIKE+J+X model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

# Add distance-dependent dispersal (+x)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0.0

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resDIVALIKEjx = bears_optim_run(BioGeoBEARS_run_object)
save(resDIVALIKEjx, file="Lychnophorinae_DIVALIKE+J+X.Rdata")
} else {
  # Loads saved analysis results
  load("Lychnophorinae_DIVALIKE+J+X.Rdata", envir = .GlobalEnv)
}

#######################################################
# BLOCK 3:BAYAREALIKE, BAYAREALIKE+X, BAYAREALIKE+J AND BAYAREALIKE+J+X ANALYSIS
#######################################################
# NOTE: As with DIVA, the BioGeoBEARS BayArea-like model is not identical with the full Bayesian model implemented in the 
# "BayArea" program of Landis et al. (2013). 
#
# Instead, this is a simplified likelihood interpretation of the model. Basically, in BayArea and BioGeoBEARS-BAYAREALIKE, 
# "d" and "e" work like they do in the DEC model of Lagrange (and BioGeoBEARS), and then BayArea's cladogenesis assumption
# (which is that nothing in particular happens at cladogenesis) is replicated by BioGeoBEARS.
#
# This leaves out 3 important things that are in BayArea:
# 1. Distance dependence (you can add this with a distances matrix + the "x" parameter in BioGeoBEARS, however)
# 2. A correction for disallowing "e" events that drive a species extinct (a null geographic range)
# 3. The neat Bayesian sampling of histories, which allows analyses on large numbers of areas.
#
# The main purpose of having a "BAYAREALIKE" model is to test the importance of the cladogenesis model on particular datasets. 
# Does it help or hurt the data likelihood if there is no special cladogenesis process?
# 
# BAYAREALIKE is a likelihood interpretation of BayArea, and it is "like BayArea" -- similar to, but not identical to, Bayesian 
# BayArea. I thus now call the model "BAYAREALIKE", and you should also. ;-)
#######################################################

#######################################################
# MODEL 9: BAYAREALIKE
#######################################################
if (BAYAREALIKE) {
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = detectCores() -4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check the inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resBAYAREALIKE = bears_optim_run(BioGeoBEARS_run_object)
save(resBAYAREALIKE, file="Lychnophorinae_BAYAREALIKE.Rdata")
} else {
  # Loads saved analysis results
  load("Lychnophorinae_BAYAREALIKE.Rdata", envir = .GlobalEnv)
}

#######################################################
# MODEL 10: BAYAREALIKE+X
#######################################################
if (BAYAREALIKEX) {
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = detectCores() -4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE+X model
# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# No jump dispersal/founder-event speciation
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = 0.01
# BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = 0.01

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Add distance-dependent dispersal (+x)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0.0

# Check the inputs
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resBAYAREALIKEx = bears_optim_run(BioGeoBEARS_run_object)
save(resBAYAREALIKEx, file="Lychnophorinae_BAYAREALIKE+X.Rdata")
} else {
  # Loads saved analysis results
  load("Lychnophorinae_BAYAREALIKE+X.Rdata", envir = .GlobalEnv)
}

#######################################################
# MODEL 11: BAYAREALIKE+J
#######################################################
if (BAYAREALIKEJ) {
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"
BioGeoBEARS_run_object$num_cores_to_use = detectCores() -4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows machines. I can't replicate this on my 
# Mac machines, but it is almost certainly just some precision under-run issue, when optim/optimx tries some parameter value 
# just below zero.  The "min" and "max" options on each parameter are supposed to prevent this, but apparently optim/optimx 
# sometimes go slightly beyond these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" slightly for each 
# parameter:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resBAYAREALIKEj = bears_optim_run(BioGeoBEARS_run_object)
save(resBAYAREALIKEj, file="Lychnophorinae_BAYAREALIKE+J.Rdata")
} else {
  # Loads saved analysis results
  load("Lychnophorinae_BAYAREALIKE+J.Rdata", envir = .GlobalEnv)
}

#######################################################
# MODEL 12: BAYAREALIKE+J+X
#######################################################
if (BAYAREALIKEJX) {
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = FALSE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"
BioGeoBEARS_run_object$num_cores_to_use = detectCores() -4
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE+J+X model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No subset sympatry
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Add distance-dependent dispersal (+x)
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["x","max"] = 0.0

# NOTE (NJM, 2014-04): BAYAREALIKE+J seems to crash on some computers, usually Windows machines. I can't replicate this on my 
# Mac machines, but it is almost certainly just some precision under-run issue, when optim/optimx tries some parameter value 
# just below zero.  The "min" and "max" options on each parameter are supposed to prevent this, but apparently optim/optimx 
# sometimes go slightly beyond these limits.  Anyway, if you get a crash, try raising "min" and lowering "max" slightly for each 
# parameter:
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resBAYAREALIKEjx = bears_optim_run(BioGeoBEARS_run_object)
save(resBAYAREALIKEjx, file="Lychnophorinae_BAYAREALIKE+J+X.Rdata")
} else {
  # Loads saved analysis results
  load("Lychnophorinae_BAYAREALIKE+J+X.Rdata", envir = .GlobalEnv)
}}

plot_ARE <- function(..., plotDEC=T, plotDECX=T, plotDECJ=T, plotDECJX=T, plotDIVALIKE=T, plotDIVALIKEX=T, plotDIVALIKEJ=T, 
                     plotDIVALIKEJX=T, plotBAYAREALIKE=T, plotBAYAREALIKEX=T, plotBAYAREALIKEJ=T, plotBAYAREALIKEJX=T) {

#######################################################
# PLOT 1: ancestral states - DEC
#######################################################
if (plotDEC) {
  analysis_titletxt ="BioGeoBEARS Lychnophorinae DEC Model"
  
  # Setup
  results_object = resDEC
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

  tiff("1.DEC_Square.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.05, 
                           tipcex=0.6, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
  
  # Pie chart
  tiff("13.DEC_Pie.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.05, 
                           tipcex=0.6, statecex=0.25, splitcex=0.25, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
}
#######################################################
# PLOT 2: ancestral states - DEC+X
#######################################################
if (plotDECX) {
  analysis_titletxt ="BioGeoBEARS Lychnophorinae DEC+X Model"
  
  # Setup
  results_object = resDECx
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

  tiff("2.DEC+X_Square.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.05, 
                           tipcex=0.6, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
  
  # Pie chart
  tiff("14.DEC+X_Pie.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.05, 
                           tipcex=0.6, statecex=0.25, splitcex=0.25, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
}
#######################################################
# PLOT 3: ancestral states - DEC+J
#######################################################
if (plotDECJ) {
  analysis_titletxt ="BioGeoBEARS Lychnophorinae DEC+J Model"
  
  # Setup
  results_object = resDECj
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

  tiff("3.DEC+J_Square.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.05, 
                           tipcex=0.6, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
  
  # Pie chart
  tiff("15.DEC+J_Pie.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.05, 
                           tipcex=0.6, statecex=0.25, splitcex=0.25, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
}
#######################################################
# PLOT 4: ancestral states - DEC+J+X
#######################################################
if (plotDECJX) {
  analysis_titletxt ="BioGeoBEARS Lychnophorinae DEC+J+X Model"
  
  # Setup
  results_object = resDECjx
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

  tiff("4.DEC+J+X_Square.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.05, 
                           tipcex=0.6, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
  
  # Pie chart
  tiff("16.DEC+J+X_Pie.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.05, 
                           tipcex=0.6, statecex=0.25, splitcex=0.25, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
}

#######################################################
# PLOT 5: ancestral states - DIVALIKE
#######################################################
if (plotDIVALIKE) {
  analysis_titletxt ="BioGeoBEARS Lychnophorinae DIVALIKE Model"
  
  # Setup
  results_object = resDIVALIKE
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

  tiff("5.DIVALIKE_Square.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.05, 
                           tipcex=0.6, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
  
  # Pie chart
  tiff("17.DIVALIKE_Pie.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.05, 
                           tipcex=0.6, statecex=0.25, splitcex=0.25, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
}
#######################################################
# PLOT 6: ancestral states - DIVALIKE+X
#######################################################
if (plotDIVALIKEX) {
  analysis_titletxt ="BioGeoBEARS Lychnophorinae DIVALIKE+X Model"
  
  # Setup
  results_object = resDIVALIKEx
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

  tiff("6.DIVALIKE+X_Square.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.05, 
                           tipcex=0.6, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
  
  # Pie chart
  tiff("18.DIVALIKE+X_Pie.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.05, 
                           tipcex=0.6, statecex=0.25, splitcex=0.25, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
}
#######################################################
# PLOT 7: ancestral states - DIVALIKE+J
#######################################################
if (plotDIVALIKEJ) {
  analysis_titletxt ="BioGeoBEARS Lychnophorinae DIVALIKE+J Model"
  
  # Setup
  results_object = resDIVALIKEj
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

  tiff("7.DIVALIKE+J_Square.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.05, 
                           tipcex=0.6, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
  
  # Pie chart
  tiff("19.DIVALIKE+J_Pie.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.05, 
                           tipcex=0.6, statecex=0.25, splitcex=0.25, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
}
#######################################################
# PLOT 8: ancestral states - DIVALIKE+J+X
#######################################################
if (plotDIVALIKEJX) {
  analysis_titletxt ="BioGeoBEARS Lychnophorinae DIVALIKE+J+X Model"
  
  # Setup
  results_object = resDIVALIKEjx
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  tiff("8.DIVALIKE+J+X_Square.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.05, 
                           tipcex=0.6, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
  
  # Pie chart
  tiff("20.DIVALIKE+J+X_Pie.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.05, 
                           tipcex=0.6, statecex=0.25, splitcex=0.25, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
}

#######################################################
# PLOT 9: ancestral states - BAYAREALIKE
#######################################################
if (plotBAYAREALIKE) {
  analysis_titletxt ="BioGeoBEARS Lychnophorinae BAYAREALIKE Model"
  
  # Setup
  results_object = resBAYAREALIKE
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

  tiff("9.BAYAREALIKE_Square.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.05, 
                           tipcex=0.6, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
  
  # Pie chart
  tiff("21.BAYAREALIKE_Pie.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.05, 
                           tipcex=0.6, statecex=0.25, splitcex=0.25, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
}
#######################################################
# PLOT 10: ancestral states - BAYAREALIKE+X
#######################################################
if (plotBAYAREALIKEX) {
  analysis_titletxt ="BioGeoBEARS Lychnophorinae BAYAREALIKE+X Model"
  
  # Setup
  results_object = resBAYAREALIKEx
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

  tiff("10.BAYAREALIKE+X_Square.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.05, 
                           tipcex=0.6, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
  
  # Pie chart
  tiff("22.BAYAREALIKE+X_Pie.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.05, 
                           tipcex=0.6, statecex=0.25, splitcex=0.25, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
}
#######################################################
# PLOT 11: ancestral states - BAYAREALIKE+J
#######################################################
if (plotBAYAREALIKEJ) {
  analysis_titletxt ="BioGeoBEARS Lychnophorinae BAYAREALIKE+J Model"
  
  # Setup
  results_object = resBAYAREALIKEj
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

  tiff("11.BAYAREALIKE+J_Square.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.05, 
                           tipcex=0.6, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
  
  # Pie chart
  tiff("23.BAYAREALIKE+J_Pie.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.05, 
                           tipcex=0.6, statecex=0.25, splitcex=0.25, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
}
#######################################################
# PLOT 12: ancestral states - BAYAREALIKE+J+X
#######################################################
if (plotBAYAREALIKEJX) {
  analysis_titletxt ="BioGeoBEARS Lychnophorinae BAYAREALIKE+J+X Model"
  
  # Setup
  results_object = resBAYAREALIKEjx
  scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
  
  tiff("12.BAYAREALIKE+J+X_Square.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="text", label.offset=0.05, 
                           tipcex=0.6, statecex=0.8, splitcex=0.8, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
  
  # Pie chart
  tiff("24.BAYAREALIKE+J+X_Pie.tiff", units = "px", width = 6000, height = 3375, res = 300) #Activate to export
  plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list(c("j","x")), plotwhat="pie", label.offset=0.05, 
                           tipcex=0.6, statecex=0.25, splitcex=0.25, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                           include_null_range=TRUE, tr=tr, tipranges=tipranges, xlims = 6.5, ylims = 75)
  dev.off()
}}

#########################################################################
# 
# CALCULATE SUMMARY STATISTICS TO COMPARE
# DEC, DEC+X, DEC+J, DEC+J+X, DIVALIKE, DIVALIKE+X, DIVALIKE+J, DIVALIKE+J+X, BAYAREALIKE, BAYAREALIKE+X, BAYAREALIKE+J,
# BAYAREALIKE+J+X
# 
#########################################################################

#########################################################################
# REQUIRED READING:
#
# Practical advice / notes / basic principles on statistical model comparison in general, and in BioGeoBEARS:
# http://phylo.wikidot.com/advice-on-statistical-model-comparison-in-biogeobears
#########################################################################

# Set up empty tables to hold the statistical results
{restable = NULL
teststable = NULL

#######################################################
# Statistics -- DEC vs. DEC+X
#######################################################
# We have to extract the log-likelihood differently, depending on the version of optim/optimx
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECx)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

# DEC, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(resDEC, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
# DEC+X, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(resDECx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

# The null hypothesis for a Likelihood Ratio Test (LRT) is that two models confer the same likelihood on the data. 
# See: Brian O'Meara's webpage: http://www.brianomeara.info/tutorials/aic ...for an intro to LRT, AIC, and AICc

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DEC vs. DEC+J
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

# DEC+J, alternative model for Likelihood Ratio Test (LRT)
res = extract_params_from_BioGeoBEARS_results_object(resDECj, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DEC vs. DEC+J+X
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDEC)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECjx)

numparams1 = 4
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

# DEC+J+X, alternative model for Likelihood Ratio Test (LRT)
res = extract_params_from_BioGeoBEARS_results_object(resDECjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DEC+X vs. DEC+J
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDECx)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECj)

numparams1 = 3
numparams2 = 3
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

tmp_tests = conditional_format_table(stats)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DEC+X vs. DEC+J+X
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDECx)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECjx)

numparams1 = 4
numparams2 = 3
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

tmp_tests = conditional_format_table(stats)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DEC+J vs. DEC+J+X
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDECj)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDECjx)

numparams1 = 4
numparams2 = 3
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

tmp_tests = conditional_format_table(stats)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DIVALIKE vs. DIVALIKE+X
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEx)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

# DIVALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(resDIVALIKE, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
# DIVALIKE+X, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(resDIVALIKEx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DIVALIKE vs. DIVALIKE+J
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

# DIVALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res = extract_params_from_BioGeoBEARS_results_object(resDIVALIKEj, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DIVALIKE vs. DIVALIKE+J+X
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEjx)

numparams1 = 4
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

# DIVALIKE+J+X, alternative model for Likelihood Ratio Test (LRT)
res = extract_params_from_BioGeoBEARS_results_object(resDIVALIKEjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DIVALIKE+X vs. DIVALIKE+J
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEx)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)

numparams1 = 3
numparams2 = 3
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

tmp_tests = conditional_format_table(stats)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DIVALIKE+X vs. DIVALIKE+J+X
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEx)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEjx)

numparams1 = 4
numparams2 = 3
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

tmp_tests = conditional_format_table(stats)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- DIVALIKE+J vs. DIVALIKE+J+X
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEj)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resDIVALIKEjx)

numparams1 = 4
numparams2 = 3
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

tmp_tests = conditional_format_table(stats)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- BAYAREALIKE vs. BAYAREALIKE+X
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEx)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

# BAYAREALIKE, null model for Likelihood Ratio Test (LRT)
res2 = extract_params_from_BioGeoBEARS_results_object(resBAYAREALIKE, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
# BAYAREALIKE+X, alternative model for Likelihood Ratio Test (LRT)
res1 = extract_params_from_BioGeoBEARS_results_object(resBAYAREALIKEx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

rbind(res2, res1)
tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res2, res1)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- BAYAREALIKE vs. BAYAREALIKE+J
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

numparams1 = 3
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

# BAYAREALIKE+J, alternative model for Likelihood Ratio Test (LRT)
res = extract_params_from_BioGeoBEARS_results_object(resBAYAREALIKEj, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- BAYAREALIKE vs. BAYAREALIKE+J+X
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKE)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEjx)

numparams1 = 4
numparams2 = 2
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

# BAYAREALIKE+J+X, alternative model for Likelihood Ratio Test (LRT)
res = extract_params_from_BioGeoBEARS_results_object(resBAYAREALIKEjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

tmp_tests = conditional_format_table(stats)

restable = rbind(restable, res)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- BAYAREALIKE+X vs. BAYAREALIKE+J
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEx)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)

numparams1 = 3
numparams2 = 3
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

tmp_tests = conditional_format_table(stats)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- BAYAREALIKE+X vs. BAYAREALIKE+J+X
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEx)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEjx)

numparams1 = 4
numparams2 = 3
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

tmp_tests = conditional_format_table(stats)
teststable = rbind(teststable, tmp_tests)

#######################################################
# Statistics -- BAYAREALIKE+J vs. BAYAREALIKE+J+X
#######################################################
LnL_2 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEj)
LnL_1 = get_LnL_from_BioGeoBEARS_results_object(resBAYAREALIKEjx)

numparams1 = 4
numparams2 = 3
stats = AICstats_2models(LnL_1, LnL_2, numparams1, numparams2)

tmp_tests = conditional_format_table(stats)
teststable = rbind(teststable, tmp_tests)

#########################################################################
# ASSEMBLE RESULTS TABLES: DEC, DEC+J, DIVALIKE, DIVALIKE+J, BAYAREALIKE, BAYAREALIKE+J
#########################################################################
teststable$alt = c("DEC+X", "DEC+J", "DEC+J+X", "DEC+J", "DEC+J+X", "DEC+J+X", 
                   "DIVALIKE+X", "DIVALIKE+J", "DIVALIKE+J+X", "DIVALIKE+J", "DIVALIKE+J+X", "DIVALIKE+J+X", 
                   "BAYAREALIKE+X", "BAYAREALIKE+J", "BAYAREALIKE+J+X", "BAYAREALIKE+J", "BAYAREALIKE+J+X", "BAYAREALIKE+J+X")
teststable$null = c("DEC", "DEC", "DEC", "DEC+X", "DEC+X", "DEC+J", 
                    "DIVALIKE", "DIVALIKE", "DIVALIKE", "DIVALIKE+X", "DIVALIKE+X", "DIVALIKE+J", 
                    "BAYAREALIKE", "BAYAREALIKE", "BAYAREALIKE", "BAYAREALIKE+X", "BAYAREALIKE+X", "BAYAREALIKE+J")
row.names(restable) = c("DEC", "DEC+X", "DEC+J", "DEC+J+X", "DIVALIKE", "DIVALIKE+X", "DIVALIKE+J", "DIVALIKE+J+X",
                        "BAYAREALIKE", "BAYAREALIKE+X", "BAYAREALIKE+J", "BAYAREALIKE+J+X")
restable = put_jcol_after_ecol(restable)

# Look at the results!!
restable
teststable

#######################################################
# Save the results tables for later -- check for, e.g., convergence issues
#######################################################

# Loads to "restable"
save(restable, file="restable.Rdata")
#load(file="restable.Rdata")

# Loads to "teststable"
save(teststable, file="teststable.Rdata")
#load(file="teststable.Rdata")

# Also save to text files
write.table(restable, file="restable.txt", quote=FALSE, sep="\t")
write.table(unlist_df(teststable), file="teststable.txt", quote=FALSE, sep="\t")

#######################################################
# Model weights of all models
#######################################################
restable2 = restable

# With AICs:
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC = put_jcol_after_ecol(restable_AIC)
restable_AIC

# With AICcs -- factors in sample size
samplesize = length(tr$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc = put_jcol_after_ecol(restable_AICc)
restable_AICc

# Also save to text files
write.table(restable_AIC, file="restable_AIC.txt", quote=FALSE, sep="\t")
write.table(restable_AICc, file="restable_AICc.txt", quote=FALSE, sep="\t")

# Save with nice conditional formatting
write.table(conditional_format_table(restable_AIC), file="restable_AIC_formatted.txt", quote=FALSE, sep="\t")
write.table(conditional_format_table(restable_AICc), file="restable_AICc_formatted.txt", quote=FALSE, sep="\t")
}

#######################################################
# BSM: Biogeographical Stochastic Mapping
#######################################################
{res = resBAYAREALIKEjx
clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

#######################################################
# Get the inputs for Biogeographical Stochastic Mapping
# Note: this can be slow for large state spaces and trees, since the independent likelihoods for each branch are being 
# pre-calculated, e.g., for 10 areas, this requires calculation of a 1024x1024 matrix for each branch.  On a tree with ~800 tips 
# and thus ~1600 branches, this was about 1.6 gigs for storage of "BSM_inputs_file.Rdata".
# Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis, the same settings will be used for 
# get_inputs_for_stochastic_mapping().
#######################################################
runInputsSlow = F
if (runInputsSlow)
{
  stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
  save(stochastic_mapping_inputs_list, file="BSM_inputs_file.Rdata")
} else {
  # Loads to "stochastic_mapping_inputs_list"
  load("BSM_inputs_file.Rdata")
} # END if (runInputsSlow)

# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))

runBSMslow = F
if (runBSMslow)
{
  # Saves to: RES_clado_events_tables.Rdata
  # Saves to: RES_ana_events_tables.Rdata
  BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=5, nummaps_goal=5, 
                      maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
  
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {
  # Load previously saved...
  
  # Loads to: RES_clado_events_tables
  load(file="RES_clado_events_tables.Rdata")
  # Loads to: RES_ana_events_tables
  load(file="RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
} # END if (runBSMslow)

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
head(clado_events_tables[[1]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)

include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
max_range_size = 8

# Note: If you did something to change the states_list from the default given the number of areas, you would
# have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)

colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, 
                                                           max_range_size=max_range_size, plot_null_range=TRUE)

#######################################################
# Summarize stochastic map tables
#######################################################
length(clado_events_tables)
length(ana_events_tables)

head(clado_events_tables[[1]][,-20])
tail(clado_events_tables[[1]][,-20])

head(ana_events_tables[[1]])
tail(ana_events_tables[[1]])

areanames = names(tipranges@df)
actual_names = areanames

# Get the dmat and times (if any)
dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables

# Simulate the source areas
BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables = BSMs_w_sourceAreas$ana_events_tables

# Count all anagenetic and cladogenetic events
counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)

summary_counts_BSMs = counts_list$summary_counts_BSMs
print(conditional_format_table(summary_counts_BSMs))

# Histogram of event counts
hist_event_counts(counts_list, pdffn="BAYAREALIKE+J+X_histograms_of_event_counts.pdf")

#######################################################
# Print counts to files
#######################################################
tmpnames = names(counts_list)
cat("\n\nWriting tables* of counts to tab-delimited text files:\n(* = Tables have dimension=2 (rows and columns). 
    Cubes (dimension 3) and lists (dimension 1) will not be printed to text files.) \n\n")
for (i in 1:length(tmpnames))
{
  cmdtxt = paste0("item = counts_list$", tmpnames[i])
  eval(parse(text=cmdtxt))
  
  # Skip cubes
  if (length(dim(item)) != 2)
  {
    next()
  }
  
  outfn = paste0(tmpnames[i], ".txt")
  if (length(item) == 0)
  {
    cat(outfn, " -- NOT written, *NO* events recorded of this type", sep="")
    cat("\n")
  } else {
    cat(outfn)
    cat("\n")
    write.table(conditional_format_table(item), file=outfn, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
  } # END if (length(item) == 0)
} # END for (i in 1:length(tmpnames))
cat("...done.\n")

#######################################################
# Check that ML ancestral state/range probabilities and
# the mean of the BSMs approximately line up
#######################################################
library(MultinomialCI)    # For 95% CIs on BSM counts
check_ML_vs_BSM(res, clado_events_tables, "BAYAREALIKE+J+X", tr=NULL, plot_each_node=FALSE, linreg_plot=TRUE, MultinomialCI=TRUE)
}
