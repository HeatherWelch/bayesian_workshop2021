### dataframe for workshop
library(devtools)
devtools::install_github('cjcarlson/embarcadero')

library(embarcadero)
library(blockCV)
library(modEvA)

source("/Users/heatherwelch/Dropbox/OLE/github/OLE_Projects/utilities/load_libraries.R")
date="03_22_21"
dat=read.csv(glue("/Users/heatherwelch/Dropbox/OLE/data/species_data_w_bkgd_envtData/sp_dat_bkgd_envt_{date}.csv"))
test=dat %>% group_by(id) %>% summarise(n=n())
df=dat %>% filter(id=="black-footedAlbatross_TOPP")
