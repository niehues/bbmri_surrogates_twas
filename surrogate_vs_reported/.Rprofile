# R libraries are assumed to be stored in r_dir
r_dir = "~/data/volume_2/R"
.libPaths(c(r_dir, .libPaths()))
# set location of tar binaries
Sys.setenv(TAR = "/bin/tar")

## when using renv
# set renv cache to writable location
renv_dir <- file.path("~/data/volume_2/renv/")
if (!dir.exists(renv_dir)) {dir.create(renv_dir)}
Sys.setenv(RENV_PATHS_ROOT = renv_dir)
# use packages installed in system library of workspace instead of installing into project-local environment
renv::settings$external.libraries(c("/etc/miniconda/lib/R/library", r_dir))
# do not snapshot BBMRIomics package
renv::settings$ignored.packages("BBMRIomics")

source("renv/activate.R") # run renv::init() # if activate.R does not exist
