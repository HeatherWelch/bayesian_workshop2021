### dataframe for workshop
library(devtools)
devtools::install_github('cjcarlson/embarcadero')

library(embarcadero)
library(blockCV)
library(modEvA)
library(dbarts)

# library(embarcadero)  # for building BART models
# library(modEvA)  # for evaluating model performance
# library(blockCV)  # for block cross-validation

source("/Users/heatherwelch/Dropbox/OLE/github/OLE_Projects/utilities/load_libraries.R")
date="03_22_21"
dat=read.csv(glue("/Users/heatherwelch/Dropbox/OLE/data/species_data_w_bkgd_envtData/sp_dat_bkgd_envt_{date}.csv"))
test=dat %>% group_by(id) %>% summarise(n=n())
df=dat %>% filter(id=="black-footedAlbatross_TOPP")

dat=df
dat=dat %>% mutate(day=yday(dt))
spc_col <- 11  # a species' presence/absence column IN THIS DATASET (change as appropriate!)
#gbm.x=c("bathy_sd","eke","l.chl","mld","oxy200m","PPupper200m","sla","sst_sd","random","sst","bathy","day")
var_cols <- c(30,19,20,21,22,24,25,32,26,15,33)
names(dat)[spc_col]  # check that it's OK
names(dat)[var_cols]  # check that it's OK
myspecies <- names(dat)[spc_col]

dat_spatial <- dat
coordinates(dat_spatial) <- dat_spatial[ , c("lon", "lat")]
crs(dat_spatial) <- "+proj=longlat"

?bart
set.seed(123)  # set a seed of random numbers so next command yields the same result in different runs of the script
mod_BART <- bart(y.train = dat[ , myspecies], x.train = dat[ , var_cols], keeptrees = TRUE)
summary(mod_BART)

invisible(mod_BART$fit$state)

source("https://raw.githubusercontent.com/AMBarbosa/unpackaged/master/predict_bart_df")  # I edited the 'embarcadero' predict function to work on data frames rather than raster layers
dat$BART_P <- predict_bart_df(mod_BART, dat)  # can take time for large datasets!
head(dat)  # predictions now on the table

# map the predictions:
dat_spatial$BART_P <- dat$BART_P
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

BART_P <- predict2.bart(mod_BART, x.layers = var_stack)  # can take time for large rasters!
saveRDS(mod_BART, "mod_BART_albatross.rds")

### day 2 #####

main=dat %>% 
  dplyr::select(c(id,date,day,species,lat,lon,presAbs,adt_sd,adt,bathy_sd,bathy,dist_shore,eke,l.chl,mld,oxy200m,oxycline,PPupper200m,sla_sd,sla,sst_sd,sst)) 
# daylight=daylength(main$lat,main$lon,JD(main$date),-9) %>% as.data.frame()
# main=main %>%  mutate(daylight=daylight$daylen)

main=main %>% mutate(random=sample(1:nrow(main)))
main=main %>% mutate(presAbs=as.integer(round(presAbs))) %>% .[complete.cases(.),]

#full model no reception
gbm.x=c("bathy_sd","eke","l.chl","mld","oxy200m","PPupper200m","sla","sst_sd","sst","bathy","day")
family="bernoulli"
if(main %>% filter(presAbs==1) %>% nrow()<500){lr=0.001}
if(main %>% filter(presAbs==1) %>% nrow()>=500){lr=0.01}
tc=3
bf=0.6
tolerance = .00001
# type=glue("{sp[i]}")
counter=1

    gap_brt_poiss_step = gbm.fixed(main,gbm.x=gbm.x,gbm.y="presAbs",family=family,learning.rate = lr, tree.complexity =tc, bag.fraction = bf,n.trees = 2000)
    if(gap_brt_poiss_step$n.trees<2000){
      print("N trees too low, refitting with smaller LR")
      lr=lr/10*5
      gap_brt_poiss_step = gbm.fixed(main,gbm.x=gbm.x,gbm.y="sum_gaps",family=family,learning.rate = lr, tree.complexity =tc, bag.fraction = bf,n.trees = 2000)
    }
    
    if(gap_brt_poiss_step$n.trees<2000){
      print("N trees too low, refitting with smaller LR")
      lr=lr/10*5
      gap_brt_poiss_step = gbm.fixed(main,gbm.x=gbm.x,gbm.y="sum_gaps",family=family,learning.rate = lr, tree.complexity =tc, bag.fraction = bf,n.trees = 2000)
    }
    
    name=glue("/Users/heatherwelch/Dropbox/bayesian_workshop2021/models/albatross_lr{lr}_tc{tc}_bf{bf}_{counter}.rds")
    write_rds(gap_brt_poiss_step,name)
    counter=counter+1

    ## brt predictions
    mod=readRDS("/Users/heatherwelch/Dropbox/bayesian_workshop2021/models/albatross_lr0.01_tc3_bf0.6_2.rds")
    a=predict(var_stack,mod,n.trees=mod$gbm.call$best.trees,type="response",na.rm=F)
    mod=readRDS("/Users/heatherwelch/Dropbox/bayesian_workshop2021/models/albatross_lr0.01_tc3_bf0.6_3.rds")
    b=predict(var_stack,mod,n.trees=mod$gbm.call$best.trees,type="response",na.rm=F)
    mod=readRDS("/Users/heatherwelch/Dropbox/bayesian_workshop2021/models/albatross_lr0.01_tc3_bf0.6_4.rds")
    c=predict(var_stack,mod,n.trees=mod$gbm.call$best.trees,type="response",na.rm=F)
    mod=readRDS("/Users/heatherwelch/Dropbox/bayesian_workshop2021/models/albatross_lr0.01_tc3_bf0.6_5.rds")
    d=predict(var_stack,mod,n.trees=mod$gbm.call$best.trees,type="response",na.rm=F)
    mod=readRDS("/Users/heatherwelch/Dropbox/bayesian_workshop2021/models/albatross_lr0.01_tc3_bf0.6_6.rds")
    e=predict(var_stack,mod,n.trees=mod$gbm.call$best.trees,type="response",na.rm=F)

    all=stack(a,b,c,d,e)
    meann=mean(all)
    sdd=calc(all,sd)

    df_maphigh=rasterToPoints(meann)%>% as.data.frame()
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
      ggtitle(glue("Black-footed Albs BRT prediction"))
    
    plot_sp
    
    df_maphigh=rasterToPoints(sdd)%>% as.data.frame()
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
      ggtitle(glue("Black-footed Albs BRT uncertainty (SD)"))
    
    plot_sp
    
    