#load packages
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(hypervolume)

### 1   *******************************************************************
###            ***  loading data ***
###    -------------------------------------------------------------------
### loading data
data<-read.csv('data.csv')%>%
  #z-scored data
  dplyr::mutate(d15N_scaled=scale(d15N),d13C_scaled=scale(d13C),
                C20.4n6_scaled=scale(C20.4n6),C20.5n3_scaled=scale(C20.5n3),C22.6n3_scaled=scale(C22.6n3))
head(data)

## ======================================================
##     =******************************************=
##      1.     calculate holes based on isotopes
##  *****************************************************

ind_ID<-unique(data$individual_ID)

hole_volume_ratio<-{}
for (i in 1:length(ind_ID)) {
  temp0<-data[data$individual_ID==ind_ID[i],]%>%dplyr::select(d15N_scaled,d13C_scaled)
  #get overall hypercolume using default parameters
  hv_1<-hypervolume_gaussian(temp0)
  hv_1@Name <- "hv_1_whole"
  # compute convex expectation
    # first thin the hypervolume
  hv_1_thinned = hypervolume_thin(hv_1, num.points=500)
  hv_1_ec <- expectation_convex(hv_1_thinned, check.memory=FALSE, use.random=TRUE)
  hv_1_ec@Name <- "Convex expectation"
  # find holes
  hv_1_holes <- hypervolume_holes(hv_1, hv_1_ec, set.check.memory=FALSE)
  hv_1_holes@Name <- "Holes"
  # extract volume statistics
  volumes <- get_volume(hypervolume_join(hv_1, hv_1_ec, hv_1_holes))
  # calculate fraction of volume that is holey
  hole_volume_ratio_0 <- data.frame(ind_ID=ind_ID[i],hole_ratio=volumes["Holes"] / volumes["Convex expectation"])
  #print(hole_volume_ratio)
  hole_volume_ratio<-rbind(hole_volume_ratio,hole_volume_ratio_0)
}
hole_volume_ratio
#plot
hole_volume_ratio%>%
  ggplot(aes(hole_ratio))+
  geom_histogram(aes(y=after_stat(c(count/sum(count))*100)),
                 size=0.1,binwidth=0.002,fill='blue',col='white')+
  geom_vline(xintercept = mean(hole_volume_ratio$hole_ratio),color='red',size=1.5,lty=2)+
  scale_x_continuous(name = 'Ratio of volumes between\ndetected holes and convex expectation',
                     limits = c(0.019,0.037),breaks = seq(0,1,0.004),expand = c(0,0))+
  scale_y_continuous(name = 'Frequency (%)',limits = c(0,37),breaks = seq(0,50,10),expand = c(0,0))+
  theme_bw()+
  labs(subtitle = 'Detection based on stable isotopes data')->is_hole_volume_plot
is_hole_volume_plot


## =================================================================================
##     =******************************************=
##      2.     calculate holes based on essential fatty acids
##  *******************************************************************************

###     ***************************************************
####     (1)       ### estimate MDS
###  ------------------------------------------------------
library(ade4)
head(data)

fa_data<-data_energy_is_fa%>%dplyr::select(individual_ID,C20.4n6_scaled,C20.5n3_scaled,C22.6n3_scaled)%>%
  ##remove rows with NA value in any column
  na.omit()%>%
  dplyr::group_by(index)%>%
  dplyr::mutate(ind_no=paste0(index,'_', 1:n()))%>%as.data.frame()

# use the dist.ktab function in the ade4 package to calculate trait distance matrix
fa_dist <- ktab.list.df(list(fa_data[,c(2:4)]))
fa_dist <- dist.ktab(fa_dist, type = "Q")
# calculate PCOA
fa_cmd <- cmdscale(fa_dist, k = 2, add = TRUE, eig = TRUE)
# extract taxa coordinates and format them
fa_pcoa <- fa_cmd$points %>%
  as.data.frame() %>%
  dplyr::rename(dim_1 = V1, dim_2 = V2) %>%
  tibble::rownames_to_column("ind_ID")
# evaluate the goodness of fit
fa_cmd$GOF #0.8084489 0.8084489

# prepare lab names
x_lab <- paste("DIM 1 (", round(fa_cmd$eig[1] / sum(fa_cmd$eig) * 100, 2), "%)", sep = "")
y_lab <- paste("DIM 2 (", round(fa_cmd$eig[2] / sum(fa_cmd$eig) * 100, 2), "%)", sep = "")

## merge MDS dimension data with fatty acids transformed data
library(tibble)
fa_data_standardized <- fa_pcoa%>%dplyr::left_join(fa_data%>%
                                                   tibble::rownames_to_column("ind_ID"),by=c('ind_ID'))
### (2) detected holes
ind_ID<-unique(fa_data_standardized$individual_ID)
fa_hole_volume_ratio<-{}
for (i in 1:length(ind_ID)) {
  temp0<-fa_data_standardized[fa_data_standardized$individual_ID==ind_ID[i],]%>%dplyr::select(dim_1,dim_2)
  #get overall hypercolume using default parameters
  hv_1<-hypervolume_gaussian(temp0)
  hv_1@Name <- "hv_1_whole"
  # compute convex expectation
  # first thin the hypervolume
  hv_1_thinned = hypervolume_thin(hv_1, num.points=500)
  hv_1_ec <- expectation_convex(hv_1_thinned, check.memory=FALSE, use.random=TRUE)
  hv_1_ec@Name <- "Convex expectation"
  # find holes
  hv_1_holes <- hypervolume_holes(hv_1, hv_1_ec, set.check.memory=FALSE)
  hv_1_holes@Name <- "Holes"
  # extract volume statistics
  volumes <- get_volume(hypervolume_join(hv_1, hv_1_ec, hv_1_holes))
  # calculate fraction of volume that is holey
  hole_volume_ratio_0 <- data.frame(ind_ID=ind_ID[i],hole_ratio=volumes["Holes"] / volumes["Convex expectation"])
  #print(fa_hole_volume_ratio)
  fa_hole_volume_ratio<-rbind(fa_hole_volume_ratio,hole_volume_ratio_0)
}
fa_hole_volume_ratio

#plot
fa_hole_volume_ratio%>%
  ggplot(aes(hole_ratio))+
  geom_histogram(aes(y=after_stat(c(count/sum(count))*100)),
                 size=0.1,binwidth=0.002,fill='orange',col='white')+#
  geom_vline(xintercept = mean(fa_hole_volume_ratio$hole_ratio),color='red',size=1.5,lty=2)+
  scale_x_continuous(name = 'Ratio of volumes between\ndetected holes and convex expectation',
                     limits = c(0.019,0.037),breaks = seq(0,1,0.004),expand = c(0,0))+
  scale_y_continuous(name = 'Frequency (%)',limits = c(0,37),breaks = seq(0,50,10),expand = c(0,0))+
  theme_bw()+
  labs(subtitle = 'Detection based on essential fatty acids data')->fa_hole_volume_plot
fa_hole_volume_plot

library(ggpubr)
hole_volume_plot<-ggarrange(is_hole_volume_plot,fa_hole_volume_plot,
                                     vjust = 2.5,
                                     labels = c("(a)",'(b)'), 
                                     nrow = 1,ncol = 2, 
                                     #widths = c(0.55,0.45),
                                     font.label = list(size = 8, color = "black", face = "bold", family = "serif"))
hole_volume_plot

#ggsave





