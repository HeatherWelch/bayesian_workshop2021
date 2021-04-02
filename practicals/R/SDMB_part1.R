# R scripts for workshop "Take your species distribution models to the next level with Bayesian nonparametric regression"
# A. Marcia Barbosa (https://modtools.wordpress.com), March 2021


# FIRST, DO THIS!!!
# RStudio menu "Session - Set Working Directory - To Source File Location"


# LOAD THE NECESSARY PACKAGES IN THE CURRENT R SESSION ####

library(embarcadero)  # for building BART models
library(modEvA)  # for evaluating model performance
library(blockCV)  # for block cross-validation


# DEFINE THE MODELLING COLUMNS ####

# we need a data frame with presence/(pseudo)absence of records of a species and a set of environmental variables:

dat <- read.csv("/Users/heatherwelch/Dropbox/bayesian_workshop2021/IBSwkshp/practicals/data/mammal_data.csv")
dat=df
# presence/absence of records of non-flying terrestrial mammals on UTM 10x10 km2 grid cells of mainland Spain
# source: Villares J.M. (2018) Inventario Espanol de Especies Terrestres (MAGRAMA). Version 1.4. Ministry of Agriculture, Food and Environment.
# accessed from R via 'rgbif::occ_data' on 2021-03-23
# you can find "species_names.csv" in the "data" folder
library(zoo)
dat=dat %>% mutate(day=yday(dt))
head(dat)
names(dat)
#gbm.x=c("bathy_sd","eke","l.chl","mld","oxy200m","PPupper200m","sla","sst_sd","random","sst","bathy","day")
spc_col <- 8  # a species' presence/absence column IN THIS DATASET (change as appropriate!)
var_cols <- 66:85  # variables are in these columns IN THIS DATASET (change as appropriate!)
spc_col <- 11  # a species' presence/absence column IN THIS DATASET (change as appropriate!)
var_cols <- c(30,19,20,21,22,24,25,32,26,15,35)
names(dat)[spc_col]  # check that it's OK
names(dat)[var_cols]  # check that it's OK

myspecies <- names(dat)[spc_col]

# if you have spatial coordinates in the data frame, map some species and variables to check everything's in place:
dat_spatial <- dat
coordinates(dat_spatial) <- dat_spatial[ , c("longitude", "latitude")]
coordinates(dat_spatial) <- dat_spatial[ , c("lon", "lat")]
crs(dat_spatial) <- "+proj=longlat"
spplot(dat_spatial, zcol = myspecies, cex = 0.5)
spplot(dat_spatial, zcol = "bio1", cex = 0.5)  # temperature
spplot(dat_spatial, zcol = "bio12", cex = 0.5)  # precipitation
spplot(dat_spatial, zcol = "adt_sd", cex = 0.5)  # temperature
spplot(dat_spatial, zcol = "bathy_sd", cex = 0.5)  # precipitation
# for the meanings of the 'bio' variables, see https://www.worldclim.org/data/bioclim.html


# COMPUTE BAYESIAN ADDITIVE REGRESSION TREES (BART) ####

?bart
set.seed(123)  # set a seed of random numbers so next command yields the same result in different runs of the script
mod_BART <- bart(y.train = dat[ , myspecies], x.train = dat[ , var_cols], keeptrees = TRUE)
summary(mod_BART)

# if you want to use this BART model e.g. for prediction or for plotting response curves IN FUTURE R SESSIONS, you need to run the following command to explicitly ask for the full information to be included when you save the model object:
invisible(mod_BART$fit$state)


# GET PREDICTIONS ON THE DATA TABLE ####

source("https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df")  # I edited the 'embarcadero' predict function to work on data frames rather than raster layers

dat$BART_P <- predict_bart_df(mod_BART, dat)  # can take time for large datasets!
head(dat)  # predictions now on the table

# map the predictions:
dat_spatial$BART_P <- dat$BART_P
spplot(dat_spatial, zcol = "BART_P", cex = 0.5)


# GET PREDICTIONS ON A RASTER STACK ####
# (if you have also raster maps of your variables)

var_stack <- stack(list.files("/Users/heatherwelch/Dropbox/bayesian_workshop2021/IBSwkshp/practicals/data/wclim_spain", full.names = TRUE))
var_stack <- stack(list.files("/Users/heatherwelch/Dropbox/OLE/data/environmental_rasters1.25/2012-01-01", full.names = TRUE,pattern=".grd"))
names_var=list.files("/Users/heatherwelch/Dropbox/OLE/data/environmental_rasters1.25/2012-09-01",pattern=".grd") %>% 
  gsub(".grd","",.)
names(var_stack)=names_var
template=raster("/Users/heatherwelch/Dropbox/OLE/spatial_data/template.grd")
library(sf)
bbox=st_read("/Users/heatherwelch/Dropbox/OLE/bounding_boxes/species_03_22_21/black-footedAlbatross_TOPP.shp")
var_stack2=var_stack %>% mask(.,bbox)
var_stack=var_stack2
a=crop(var_stack2,bbox)
plot(a[[1:16]])
plot(a[[16:22]])
# citation for these data:
# Fick S.E. & Hijmans R.J. (2017) Worldclim 2: New 1-km spatial resolution climate surfaces for global land areas. International Journal of Climatology 37:4302-4315. Accessed from R via 'raster::getData' on 2021-03-23. Cut to study area using 'raster::crop' and 'raster::mask'.
var_stack
plot(var_stack)
BART_P <- predict2.bart(mod_BART, x.layers = var_stack)  # can take time for large rasters!
test=dismo::predict(mod_BART,var_stack)
plot(BART_P)
## predicting with missing data: https://onlinelibrary.wiley.com/doi/abs/10.1002/cjs.11248

df_maphigh=rasterToPoints(BART_P)%>% as.data.frame()
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
  coord_sf(xlim = c(180, 260), ylim = c(10,62),expand=F)+
  # ggtitle(glue("{names_sp[ii]}_{datelist[i]} continuous"))

# save the created objects to disk for future use:
saveRDS(mod_BART, "mod_BART.rds")
writeRaster(BART_P, "BART_pred.tif")  # you can open this in GIS software


# EVALUATE MODEL PREDICTIONS ####

# area under the ROC Curve:
?AUC
par(mfrow = c(1, 1), mar = c(5, 4, 2, 1))
AUC(obs = dat[ , spc_col], pred = dat[ , "BART_P"])  # see the 'plots' pane

# area under the precision-recall Curve:
AUC(obs = dat[ , spc_col], pred = dat[ , "BART_P"], curve = "PR")


# threshold-based classification metrics:

# first, let's convert our predictions of presence probability into binary presence or absence
# you don't normally need to do this, but be aware that when using threshold-based measures, the binarized prediction is what you are evaluating:
dat$BART_P_01 <- ifelse(dat$BART_P >= prevalence(dat[ , myspecies]), 1, 0)
head(dat)
dat_spatial$BART_P_01 <- dat$BART_P_01
spplot(dat_spatial, zcol = "BART_P_01")
# let's do the same for the raster predictions:
BART_P_01 <- BART_P
BART_P_01[BART_P >= prevalence(dat[ , myspecies])] <- 1
BART_P_01[BART_P < prevalence(dat[ , myspecies])] <- 0
par(mfrow = c(2, 1))
plot(BART_P)
plot(BART_P_01)

?threshMeasures
classif_metrics <- c("CCR", "Sensitivity", "Specificity", "Precision", "Recall", "TSS", "kappa")
threshMeasures(obs = dat[ , spc_col], pred = dat[ , "BART_P"], thresh = "preval", measures = classif_metrics, ylim = c(0, 1))


# now let's look at some measures of calibration, which assess the accuracy of the continuous predictions:

# Miller calibration line:
?MillerCalib
MillerCalib(obs = dat[ , spc_col], pred = dat[ , "BART_P"])

# Hosmer-Lemeshow goodness-of-fit:
?HLfit
HLfit(obs = dat[ , spc_col], pred = dat[ , "BART_P"], bin.method = "round.prob")
# beware: results may strongly depend on the arbitrary choice of 'bin.method'!

# BUT THIS JUST EVALUATES HOW THE MODELS FIT THE SAME DATA ON WHICH THEY WERE TRAINED
# YOU CAN SET ASIDE SOME DATA TO LEAVE OUT OF MODEL TRAINING AND USE FOR TESTING OUT-OF-SAMPLE PREDICTIVE CAPACITY
# BLOCK CROSS-VALIDATION (below) IS CURRENTLY THE MOST APPROPRIATE METHOD


# DIVIDE STUDY AREA INTO SPATIAL BLOCKS ####
# for model cross-validation

?spatialBlock
names(dat)

# if you have your variables as maps in a raster stack, you can calculate their range of spatial autocorrelation
# but note that this can be too stringent, i.e. make blocks too large for the model training sets to contain sufficient information
?spatialAutoRange
#
sarange <- spatialAutoRange(var_stack)  # see the 'plots' pane
sarange$range

# get spatial blocks of a given size (e.g. 150 km2):
set.seed(321)  # set a seed of random numbers so next command yields the same result in different runs of the script
blocks <- spatialBlock(speciesData = dat_spatial, rasterLayer = var_stack[[1]], theRange = 150000, species = myspecies, k = 5)  # argument 'species' is optional and makes the process slower - see ?spatialBlock

blocks$folds
blocks$foldID

dat$foldID <- blocks$foldID
head(dat)

# map each fold:
dat_spatial$foldID <- dat$foldID
folds <- sort(unique(dat$foldID))
par(mfrow = c(3, 2), mar = c(1, 1, 1, 1))
for (f in folds) {
  plot(dat_spatial, col = "grey")
  points(subset(dat_spatial, foldID == f), col = "blue")
  title(paste("Fold", f))
}


# COMPUTE MODELS AND GET PREDICTIONS LEAVING OUT EACH FOLD IN TURN ####

names(dat)

for (f in folds) {
  print(f)
  dat_train <- subset(dat, foldID != f)
  mod_BART_fold <- bart(y.train = dat_train[ , spc_col], x.train = dat_train[ , var_cols], keeptrees = TRUE, verbose = FALSE)
  dat[ , paste0("BART_fold", f, "_P")] <- predict_bart_df(mod_BART_fold, dat)
}  # end for s for f

# see the new predictions added to the the data frame:
head(dat)


# EVALUATE EACH MODEL ON ITS VALIDATION FOLD ####

# identify columns containing the predictions of different folds:
fold_cols <- grep("_fold", names(dat))
fold_cols
names(dat)[fold_cols]

# choose the measures to calculate:
measures <- c("AUC", "AUPR", "TSS", "MCS")

# create an empty table to receive the cross-validation results:
crossval <- as.data.frame(matrix(nrow = length(folds), ncol = length(measures)))
colnames(crossval) <- measures
crossval  # for now it's only filled with NAs

par(mfrow = c(2, 2), mar = c(4, 3, 1, 1), oma = c(0, 0, 2, 0))
for (f in folds) {
  fold_col <- names(dat)[grep(paste0("BART_fold", f), names(dat))]
  fold_dat <- subset(dat, foldID == f)
  crossval[f, "AUC"] <- AUC(obs = fold_dat[ , spc_col], pred = fold_dat[ , fold_col], simplif = TRUE, plot = TRUE, main = "AUC")
  crossval[f, "AUPR"] <- AUC(obs = fold_dat[ , spc_col], pred = fold_dat[ , fold_col], curve = "PR", simplif = TRUE, plot = TRUE, main = "AUPR")
  crossval[f, "TSS"] <- threshMeasures(obs = fold_dat[ , spc_col], pred = fold_dat[ , fold_col], thresh = "preval", measures = "TSS", simplif = TRUE, standardize = FALSE, main = "TSS")
  crossval[f, "MCS"] <- MillerCalib(obs = fold_dat[ , spc_col], pred = fold_dat[ , fold_col], main = "Miller line")$slope
  HLfit(obs = fold_dat[ , spc_col], pred = fold_dat[ , fold_col], bin.method = "quantiles", main = "H-L fit")
  mtext(paste("Fold", f), outer = TRUE, cex = 1.5)
}
# press the back arrow on the top left of your plotting pane to see the different fold evaluations

crossval

# boxplots of cross-validation metrics:
par(mfrow = c(1, 1), mar = c(7, 3, 2, 1))
boxplot(crossval, las = 2)
abline(h = 1, col = "darkgrey", lty = 2)  # remember Miller calibration slope (MCS) should ideally be close to 1 (not bigger = better)

# save this session's predictions to a CSV file:
head(dat)
names(dat)
write.csv(dat[ , c("cells", "BART_P", "foldID", "BART_fold1_P", "BART_fold2_P", "BART_fold3_P", "BART_fold4_P", "BART_fold5_P")], "BART_preds_part1.csv", row.names = FALSE)
