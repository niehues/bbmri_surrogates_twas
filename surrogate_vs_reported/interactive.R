################################################################################
## Title:         interactive.R
## Description:   Commands for running workflow in interactive session
## Author:        Anna Niehues
## Date created:  2020-11-20
## Email:         anna.niehues@radboudumc.nl
################################################################################
## Notes:
## These commands are demonstrated at
## https://github.com/wlandau/drake-example/blob/main/main/interactive.R
################################################################################
#renv::activate()
renv::status()

library(drake)
r_outdated() # list outdated analysis steps
drake::r_make() # run drake plan
r_outdated()
r_vis_drake_graph()
# clean(purge = TRUE)
# clean(destroy = T, purge = T, garbage_collection = T)

# build single target
source("_drake.R"); drake::drake_build(report, plan)
# load a target into the environment
loadd(model_parameters)

