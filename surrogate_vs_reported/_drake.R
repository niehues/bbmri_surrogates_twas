################################################################################
## Title:         _drake.R
## Description:   This file is required when by executing workflow using
##                `drake::r_make()`
## Author:        Anna Niehues
## Date created:  2020-11-20
## Email:         anna.niehues@radboudumc.nl 
################################################################################
## Notes:
## 
################################################################################
source("R/functions.R")
source("input.R")
source("R/plan.R")

# options(clustermq.scheduler = "multicore") # for parallel computing, requires libzmq / zeromq

# configure drake with plan
config <- drake::drake_config(
  plan,
  # parallelism = "clustermq",
  # jobs = 3,
  memory_strategy = "autoclean", 
  garbage_collection = TRUE
)



