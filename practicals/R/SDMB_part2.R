# R scripts for workshop "Take your species distribution models to the next level with Bayesian nonparametric regression"
# A. Marcia Barbosa (https://modtools.wordpress.com), March 2021


# NOTE: most commands on this script can take a LONG TIME to run - plan ahead and prepare to be patient!!!


# FIRST, DO THIS!!!
# RStudio menu "Session - Set Working Directory - To Source File Location"


# LOAD THE NECESSARY PACKAGES IN THE CURRENT R SESSION ####

library(embarcadero)


# if you have (re)started a clean R session, import the necessary data to the current R workspace:

# dat <- read.csv("../data/mammal_data.csv")
# var_stack <- stack(list.files("../data/wclim_spain", full.names = TRUE))
# mod_BART <- readRDS("mod_BART.rds")


# SELECT RELEVANT PREDICTOR VARIABLES WITH BART ####

# check out your data:

names(dat)
# spc_col <- 8  # a species' presence/absence column IN THIS DATASET (change as appropriate!)
# var_cols <- 66:85  # variables are in these columns IN THIS DATASET (change as appropriate!)
names(dat)[spc_col]
names(dat)[var_cols]  # check that everything looks OK


varimp(mod_BART, plots = TRUE)


# get a diagnostic plot of variable importance for our data and variables:
# (you can skip this step if it takes way too long, as results are just visual and not needed further in the script)

?varimp.diag

varimp_BART <- varimp.diag(x.data = dat[ , var_cols], y.data = dat[ , spc_col], iter = 10)  # the recommended default is 50 iterations; I dropped it to 10 to save some time here, but WHEN YOU DO THIS FOR REAL ANALYSIS, REMOVE 'iter = 10' !!!!!

varimp_BART


# select minimal subset of relevant variables with BART:
# (we selected variables before in script 02, but using GLM)

?variable.step

set.seed(789)
varselect <- variable.step(x.data = dat[ , var_cols], y.data = dat[ , spc_col], iter = 10)  # again, I reduced the number of iterations to 10 to save time here, but the sensible default is 50, so REMOVE 'iter = 10' WHEN YOU DO MORE SERIOUS WORK !!!!!

varselect
# for the meanings of these variables, see https://www.worldclim.org/data/bioclim.html


# BUILD A FINAL MODEL WITH THE BART-SELECTED VARIABLES ####

set.seed(4)
mod_BART_varselect <- dbarts::bart(y.train = dat[ , spc_col], x.train = dat[ , varselect], keeptrees = TRUE)

summary(mod_BART_varselect)


# if you want to use this BART model e.g. for prediction or for plotting response curves IN FUTURE R SESSIONS, you need to explicitly ask for the full information to be included when you next save the model object:

invisible(mod_BART_varselect$fit$state)


# GET PREDICTIONS ON THE DATA FRAME ####
# (including confidence intervals)
mod_BART_varselect=mod_BART
head(dat)

source("https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df")  # I edited the BART predict function to work on data frames rather than raster layers

bart_pred <- predict_bart_df(mod_BART_varselect, dat, quantiles = c(0.025, 0.975))
head(bart_pred)
bart_pred$uncertainty <- bart_pred[ , 3] - bart_pred[ , 2]
head(bart_pred)

# plot predicted values with confidence intervals:
plot(sort(bart_pred$pred), pch = ".")
lines(sort(bart_pred$q0025), col = "darkgrey")
lines(sort(bart_pred$q0975), col = "darkgrey")

# plot predicted value against prediction uncertainty:
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
plot(x = bart_pred$pred, y = bart_pred$uncertainty)


# map the predictions

# bart_pred <- data.frame(bart_pred, dat[ , c("longitude", "latitude")])
bart_pred <- data.frame(bart_pred, dat[ , c("lon", "lat")])
# coordinates(bart_pred) <- ~ longitude + latitude
coordinates(bart_pred) <- ~ lon + lat
spplot(bart_pred, cex = 0.1)


# GET PREDICTIONS ON RASTER MAPS ####
# (if you have also raster maps of your variables)

var_stack <- stack(list.files("/Users/heatherwelch/Dropbox/OLE/data/environmental_rasters1.25/2012-01-01", full.names = TRUE,pattern=".grd"))
names_var=list.files("/Users/heatherwelch/Dropbox/OLE/data/environmental_rasters1.25/2012-09-01",pattern=".grd") %>% 
  gsub(".grd","",.)
names(var_stack)=names_var
e <- as(extent(to360(120), to360(-100), 10, 62), 'SpatialPolygons')
var_stack=crop(var_stack,e)

BART_varselect_P <- predict2.bart(mod_BART_varselect, x.layers = var_stack, quantiles = c(0.025, 0.975))  # may take very long for large or high-resolution rasters!
BART_varselect_P  # RasterStack with 3 layers
names(BART_varselect_P)
names(BART_varselect_P) <- c("predicted", "q0025", "q0975")

BART_varselect_uncert <- BART_varselect_P[["q0975"]] - BART_varselect_P[["q0025"]]

par(mfrow = c(2, 2), mar = c(3, 2, 2, 2))  # set 2x2 plots per window
plot(BART_varselect_P[["q0025"]], main = 'Lower 95% CI bound', zlim = c(0, 1))
plot(BART_varselect_P[["q0975"]], main = 'Upper 95% CI bound', zlim = c(0, 1))
plot(BART_varselect_P[["predicted"]], main = 'Posterior mean', zlim = c(0, 1))
plot(BART_varselect_uncert, main = "Credible interval width", zlim = c(0, 1))

df_maphigh=rasterToPoints(BART_varselect_uncert)%>% as.data.frame()
colnames(df_maphigh)=c("rows","cols","value")

plot_sp=ggplot()+
  geom_tile(data=df_maphigh,aes(x = rows, y = cols, fill=value))+
  # geom_point(data=sp_dat,aes(x=lon,y=lat),color="red",size=.3)+
  scale_fill_gradientn(colours = pals::parula(100),na.value="black")+
  geom_polygon(data = fortify(maps::map("world2",plot=FALSE,fill=TRUE)), aes(x=long, y = lat, group=group),color="black",fill="grey")+
  # geom_polygon(data = fortify(boxx,plot=FALSE,fill=TRUE), aes(x=long, y = lat, group=group),color="red",fill=NA)+
  # geom_polygon(data=fortify(boxx,plot=F,fill=F), aes(x=long, y = lat, group=group),color="black",fill=NA,size=1)+
  geom_sf(data=bbox,color="grey",fill=NA)+
  theme_classic()+xlab(NULL)+ylab(NULL)+
  coord_sf(xlim = c(120, 260), ylim = c(10,62),expand=F)+
   ggtitle(glue("Black-footed Albs BART uncertainty"))

plot_sp

df_maphigh=rasterToPoints(BART_varselect_P[["predicted"]])%>% as.data.frame()
colnames(df_maphigh)=c("rows","cols","value")

plot_sp=ggplot()+
  geom_tile(data=df_maphigh,aes(x = rows, y = cols, fill=value))+
  # geom_point(data=sp_dat,aes(x=lon,y=lat),color="red",size=.3)+
  scale_fill_gradientn(colours = pals::parula(100),na.value="black")+
  geom_polygon(data = fortify(maps::map("world2",plot=FALSE,fill=TRUE)), aes(x=long, y = lat, group=group),color="black",fill="grey")+
  # geom_polygon(data = fortify(boxx,plot=FALSE,fill=TRUE), aes(x=long, y = lat, group=group),color="red",fill=NA)+
  # geom_polygon(data=fortify(boxx,plot=F,fill=F), aes(x=long, y = lat, group=group),color="black",fill=NA,size=1)+
  geom_sf(data=bbox,color="grey",fill=NA)+
  theme_classic()+xlab(NULL)+ylab(NULL)+
  coord_sf(xlim = c(120, 260), ylim = c(10,62),expand=F)+
  ggtitle(glue("Black-footed Albs BART prediction"))

plot_sp

# plot predicted values with confidence intervals:
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
plot(sort(getValues(BART_varselect_P[["predicted"]])), pch = ".")
lines(sort(getValues(BART_varselect_P[["q0025"]])), col = "darkgrey")
lines(sort(getValues(BART_varselect_P[["q0975"]])), col = "darkgrey")

# plot predicted value against prediction uncertainty:
plot(x = getValues(BART_varselect_P[["predicted"]]), y = getValues(BART_varselect_uncert))


# PARTIAL DEPENDENCE PLOTS ####

# get partial dependence plots (with Bayesian confidence intervals) for different variables, e.g. the first 2 in 'varselect':

?embarcadero::partial

varselect

partial_BART_1 <- partial(mod_BART_varselect, x.vars = varselect[1])  # takes time!
partial_BART_1

partial_BART_2 <- partial(mod_BART_varselect, x.vars = varselect[2], smooth = 5, trace = FALSE)
partial_BART_2

# you can get the partial dependence plots for all variables in the model at once, but notice this may take a very long time!
#par(mfrow = c(4, 2))
#partials <- partial(mod_BART_varselect)
#partials


# get a two-dimensional dependence plot, e.g. using the first 2 selected variables:

?pd2bart
partial_BART_1and2 <- pd2bart(mod_BART_varselect, xind = c(varselect[1], varselect[2]))
par(mfrow = c(1, 1), mar = c(5, 5, 2, 1))
plot(partial_BART_1and2)

# you can also choose variables by name (if they are in your model):
partial_BART_bio4_bio5 <- pd2bart(mod_BART_varselect, xind = c("bio4", "bio5"))


# OBTAIN SPATIAL PARTIAL DEPENDENCE PLOTS ####

# these are a cartographic version of partial dependence plots, showing where a particular variable most favours species presence

?spartial
spartial_BART_1 <- spartial(mod_BART_varselect, envs = var_stack, x.vars = varselect[1])
plot(spartial_BART_1, main = paste("partial prediction:", varselect[1]))

spartial_BART_2 <- spartial(mod_BART_varselect, envs = var_stack, x.vars = varselect[2])
plot(spartial_BART_2, main = paste("partial prediction:", varselect[2]))


# you can generate plots with different variables and combinations, to see how each variable or interaction affects your species
# try interpreting the results in light of your species' ecology
# read each function's help file and possibly try alternative options


# plot the MCMC draws (BART learning process):
source("https://raw.githubusercontent.com/cjcarlson/embarcadero/master/R/plot_mcmc.R")  # (while the function name doesn't stabilize in the package)
plot.mcmc(mod_BART_varselect, var_stack, wait = 0.5)


# save results for future use:
saveRDS(mod_BART_varselect, "mod_BART_varselect.rds")
write.csv(bart_pred, "BART_preds_part2.csv", row.names = FALSE)


# for further info, see "mee313389-sup-0001-Supinfo.pdf" under "Supporting Information" at the bottom of https://doi.org/10.1111/2041-210X.13389
