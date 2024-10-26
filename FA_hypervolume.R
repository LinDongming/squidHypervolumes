# require library
library(dplyr)
library(ggplot2)
library(ggpubr) 
library(plyr) 
library(tidyr)

#####  FA   #######

### 1   *******************************************************************
###            ***  loading data ***
###    -------------------------------------------------------------------
### loading data
data<-read.csv('data.csv')%>%
  #z-scored data
  dplyr::mutate(d15N_scaled=scale(d15N),d13C_scaled=scale(d13C),
                C20.4n6_scaled=scale(C20.4n6),C20.5n3_scaled=scale(C20.5n3),C22.6n3_scaled=scale(C22.6n3))
head(data)

### 2 ***************************************************************************
##       ** obtain testing results and make mean±sd table for plotting **
###  ---------------------------------------------------------------------------
 #long format dataset for testing between tissues
data_longformat<-data%>%
  dplyr::select(tissue,C20.4n6,C20.5n3,C22.6n3)%>%na.omit()%>%
  tidyr::pivot_longer(cols = c(C20.4n6,C20.5n3,C22.6n3),names_to = 'variable',values_to = 'value')

testFr<-c(unique(data_longformat$variable))
test_tissues<-{}
#test_CE_sex_forPlot<-{}
for (i in 1:length(testFr)) {
  temp0<-data_longformat%>%dplyr::filter(variable==testFr[i])
  # normal distribution testing
  kt<-with(data=temp0,ks.test(value,'pnorm',mean=mean(value,na.rm=T),sd=sd(value,na.rm=T))) 
  if(kt$p.value>0.05){ ## satisfy normality
    # anova test
    lm_is<-lm(value~tissue,data=temp0)
    fa_ano<-anova(lm_is)  # anova test
    fa_hsd.lm.fr<- data.frame(test_variable=paste0(testFr[i]),test='ANOVA',
                              statisticV=fa_ano$`F value`[1],P=fa_ano$`Pr(>F)`[1])
    
    test_tissues<-rbind(test_tissues,fa_hsd.lm.fr)
  }else{ # reject normality
    kw.t<-with(data=temp0,kruskal.test(x=value,g=tissue)) #K-W test (Kruskal-Wallis test)
    fa_hsd.lm.fr<- data.frame(test_variable=paste0(testFr[i]),
                              test='Kruskal-Wallis test',statisticV=kw.t$statistic[[1]],P=kw.t$p.value)
    
    test_tissues<-rbind(test_tissues,fa_hsd.lm.fr)
  }
}
test_tissues

### 2   *******************************************************************
###            *** estimate individual hypervolume using gausian ***
###    -------------------------------------------------------------------

##(1) estimate MDS
###  ------------------------------------------------------
library(ade4)
head(data)
#Append a sequential count of a column based on tissue sampling
fa_data<-data%>%dplyr::select(individual_ID,C20.4n6_scaled,C20.5n3_scaled,C22.6n3_scaled)%>%
  ##remove rows with NA value in any column
  na.omit()%>%
  dplyr::group_by(individual_ID)%>%
  dplyr::mutate(ind_no=paste0(individual_ID,'_', 1:n()))%>%as.data.frame()

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
#fa_mds_data$index<-as.character(fa_mds_data$index)
head(fa_data)
library(tibble)
fa_data_standardized <- fa_pcoa%>%
  dplyr::left_join(fa_data%>%
                     tibble::rownames_to_column("ind_ID"),by=c('ind_ID'))
head(fa_data_standardized)

## (2) esitimate individual hypervolume
#Calculate bandwidth from all individuals for use in the Gaussian hypervolumes
library(hypervolume)

ia_ind<-unique(fa_data_standardized$individual_ID)
#ia_overall_bandwidth<-estimate_bandwidth(IS_illex.longFormat[,10:11][which(IS_illex.longFormat$index%in%ia_ind),])
#Using gaussian method, varying bandwidth and quantile:
fa_gausian_0.95<-{}
fa_hv_ind<-list()
for (i in 1:length(ia_ind)) {
  print(i)
  individual_i<-ia_ind[i]
  temp0<-fa_data_standardized[fa_data_standardized$individual_ID==ia_ind[i],]%>%dplyr::select(dim_1,dim_2)
  #get overall hypercolume using default parameters
  fa_hv_ind[[i]]<-hypervolume_gaussian(temp0)
  val_0.95<-get_volume(fa_hv_ind[[i]])%>%as.data.frame()
  temp1<-data.frame(individual_ID=individual_i,volume_0.95=val_0.95[1,])
  fa_gausian_0.95<-rbind(fa_gausian_0.95,temp1)
}

### (3) estimate the union of hypervolume to calculate specialization index
fa_hv_ind_join<-do.call(hypervolume_join,fa_hv_ind)#convert the list to a HypervolumeList

hv1 <- fa_hv_ind_join@HVList[[1]]
for(j in 2:length(fa_hv_ind_join@HVList)) {
  #j <- 2
  hv2 <- fa_hv_ind_join@HVList[[j]]
  hvSet <- hypervolume_set(hv1, hv2, check.memory=FALSE, verbose=FALSE)
  hv1 <- hvSet@HVList$Union
}
# estimate individual specialization index as individual niche divided by the population hypervolume
fa_union_volume<-get_volume(hv1)%>%as.data.frame()
fa_niche<-fa_gausian_0.95%>%dplyr::mutate(fa_union_volume=fa_union_volume[1,],
                                          specialization_index=1-volume_0.95/fa_union_volume)
fa_niche
summary(fa_niche)

##(3) estimate the pairwise overlap
fa_Jaccard_similarity<-{} #The Jaccard Similarity will be 0 if the two sets don't share any values and 1
#if the two sets are identical. The set may contain either numerical values or strings.
for (i in 1:length(fa_hv_ind)) {
  hv1=fa_hv_ind[[i]]
  ia_index_i=ia_ind[i]
  z=i+1
  if(z<52){
    for (j in z:length(fa_hv_ind)) {
      hv2=fa_hv_ind[[j]]
      ia_index_j=ia_ind[j]
      #Computes the intersection, union, and unique components of two hypervolumes
      hv_set <- hypervolume_set(hv1, hv2, check.memory=FALSE)
      jac_sim<-hv_set@HVList$Intersection@Volume/hv_set@HVList$Union@Volume#niche overlap; 
      #same to hypervolume_overlap_statistics(hv_set)
      jac_nested<-hv_set@HVList$Intersection@Volume/min(hv1@Volume,hv2@Volume)#nestedness index
      #hypervolume_overlap_statistics(hv_set)
      temp0<-data.frame(comp_index=paste0(i,'_',j),
                        comp_index_ind=paste0(ia_index_i,'_',ia_index_j),
                        jac_similarity=jac_sim,jac_nestedIndex=jac_nested)
      fa_Jaccard_similarity<-rbind(fa_Jaccard_similarity,temp0)
    }}
}

### 3   *******************************************************************
###            *** plot individual specialization index ***
###    -------------------------------------------------------------------
#(1) plot fatty acids variation for individuals
head(data)
cc <- scales::seq_gradient_pal("darkgrey", "lightgrey", "Lab")(seq(0,1,length.out=100))
#(a) plot C20.4n6 density distribution using kernel method for each individual
data_energy_is_fa%>%dplyr::select(individual_ID,C20.4n6)%>%na.omit()%>%
  ggplot(aes(x = C20.4n6)) + 
  geom_density(aes(color=as.factor(individual_ID)),alpha = 0.4,show.legend = F,linewidth=0.2) + 
  scale_colour_manual(values=cc)+
  geom_density(aes(x = mean(C20.4n6)),color='red',linewidth=0.5,show.legend = F) + 
  scale_x_continuous(name = expression(paste('20:4n6')),expand=c(0, 0),limits = c(-0.01,2.8),breaks = seq(0,20,0.5))+
  scale_y_continuous(name = expression(paste('Kernel density')),expand=c(0, 0),limits = c(-0.01,10),breaks = seq(0,20,2))+
  theme(axis.line = element_line(color = "black",linewidth=0.3),
        panel.grid = element_line(linewidth = 0.1))->C20.4n6_density.plot
C20.4n6_density.plot

#(b) plot C20.5n3 density distribution using kernel method for each individual
data_energy_is_fa%>%dplyr::select(individual_ID,C20.5n3)%>%na.omit()%>%
  ggplot(aes(x = C20.5n3)) + 
  geom_density(aes(color=as.factor(individual_ID)),alpha = 0.4,show.legend = F,linewidth=0.2) + 
  scale_colour_manual(values=cc)+
  geom_density(aes(x = mean(C20.5n3)),color='red',linewidth=0.5,show.legend = F) + 
  scale_x_continuous(name = expression(paste('20:5n3')),expand=c(0, 0),limits = c(5.5,18.5),breaks = seq(0,20,3))+
  scale_y_continuous(name = expression(paste('Kernel density')),expand=c(0, 0),limits = c(-0.01,3.5),breaks = seq(0,20,1))+
  theme(axis.line = element_line(color = "black",linewidth=0.3),
        panel.grid = element_line(linewidth = 0.1))->C20.5n3_density.plot
C20.5n3_density.plot

#(c) plot C22.6n3 density distribution using kernel method for each individual
data_energy_is_fa%>%dplyr::select(individual_ID,C22.6n3)%>%na.omit()%>%
  ggplot(aes(x = C22.6n3)) + 
  geom_density(aes(color=as.factor(individual_ID)),alpha = 0.4,show.legend = F,linewidth=0.2) + 
  scale_colour_manual(values=cc)+
  geom_density(aes(x = mean(C22.6n3)),color='red',linewidth=0.5,show.legend = F) + 
  scale_x_continuous(name = expression(paste('22:6n3')),expand=c(0, 0),limits = c(12,38),breaks = seq(0,50,5))+
  scale_y_continuous(name = expression(paste('Kernel density')),expand=c(0, 0),limits = c(-0.001,0.305),breaks = seq(0,20,0.1))+
  theme(axis.line = element_line(color = "black",linewidth=0.3),
        panel.grid = element_line(linewidth = 0.1))->C22.6n3_density.plot
C22.6n3_density.plot

#(d) combine plots
fa_density.plot<-ggarrange(C20.4n6_density.plot,C20.5n3_density.plot,C22.6n3_density.plot,
                           labels=c("(a)",'(b)','(c)'),
                           hjust = -3.5,#More negative values move the label further to the right on the plot canvas
                           vjust=3.5, #More positive values move the label further down on the plot canvas
                           font.label = list(size = 10, color = "black", face = "plain", family = "serif"),
                           nrow = 1,ncol = 3)
fa_density.plot

#(2) plot individual specialization index distribution
#fa_specialization<-read.csv('01_fa_individualSpecialization2024-01-20.csv')
head(fa_niche)
#create dataframe for polar plotting
fa_niche_gp<-fa_niche%>%dplyr::select(individual_ID,specialization_index)%>%
  dplyr::mutate(specialization_index1=specialization_index-0.4)%>%
  dplyr::mutate(ID=as.numeric(row.names(.)),
                angle0=90-360*(as.numeric(row.names(fa_niche))-0.5)/(nrow(fa_niche)+6))%>%
  dplyr::mutate(hjust=ifelse(as.numeric(angle0)<0,ifelse(abs(angle0)>90,1,0),0),
                angle=ifelse(as.numeric(angle0)<0,ifelse(abs(angle0)>90,angle0+180,angle0),angle0))%>%
  dplyr::mutate(index_gp=case_when(specialization_index1<0.15~'gp_1',
                                   specialization_index1<0.3~'gp_2',
                                   specialization_index1<0.45~'gp_3',
                                   specialization_index1<0.6~'gp_4'))

fa_niche_gp%>%ggplot()+
  geom_linerange(aes(xmin=0,xmax=56.5, y = 0),lty="solid", color = "black")+
  geom_linerange(aes(xmin=0,xmax=56.5, y = 0.15),lty="solid", size=0.2, color = "grey80")+
  geom_linerange(aes(xmin=0,xmax=56.5, y = 0.3),lty="solid", size=0.5, color = "grey80")+
  geom_linerange(aes(xmin=0,xmax=56.5, y = 0.45),lty="solid",size=0.2,  color = "grey80")+
  geom_linerange(aes(xmin=0,xmax=56.5, y = 0.6),lty="solid", size=0.2, color = "grey80")+
  #手动添加坐标及标题
  geom_text(x=-1.1,y=-0.02,label="0.40",color="black",size=2)+
  geom_text(x=-1,y=0.15,label="0.55",color="black",size=2)+
  geom_text(x=-1,y=0.3,label="0.70",color="black",size=2)+
  geom_text(x=-1,y=0.45,label="0.85",color="black",size=2)+
  geom_text(x=-1,y=0.6,label="1.0",color="black",size=2)+
  geom_col(aes(ID, specialization_index1, fill=index_gp), alpha = 0.6)+
  #极坐标转换
  coord_polar(direction=1)+
  #x,y轴范围确定
  scale_y_continuous(limits = c(-0.25,0.6),breaks = seq(0,1,0.2))+#,0.5
  scale_x_continuous(limits = c(0,58))+
  #主题
  labs(title='Individual specialization index',subtitle = "(based on fatty acid data)")+
  theme_void()+p8.8+
  theme(legend.position = 'none',
        plot.title = element_text(color="black", size=9, face="bold"),#bold.italic
        plot.subtitle = element_text(color="black", size=7.5, face="plain"),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y.left = element_blank())->fa_specialization_plot
fa_specialization_plot

###(3) plot pairwise overlap distribution
head(fa_Jaccard_similarity)

fa_Jaccard_similarity%>%ggplot(aes(x=jac_similarity))+
  geom_histogram(aes(y=after_stat(c(count/sum(count))*100)),
                 size=0.3,binwidth = 0.05,color='lightgrey',fill='lightblue',
                 center = TRUE,show.legend = F)+
  xlab("Niche overlap for pairwise individuals")+ylab("Frequency (%)")+
  scale_x_continuous(expand = c(0, 0),limits = c(-0.01,0.9),breaks = seq(0,1,0.2))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,15),breaks = seq(0,50,5))+
  theme_bw()+
  theme(panel.grid.major = element_line(linewidth=0.1))->fa_Jaccard_similarity.plot
fa_Jaccard_similarity.plot

## save plots
fa_individual_specialization.plot<-ggarrange(fa_density.plot,
                                             ggarrange(fa_specialization_plot,fa_Jaccard_similarity.plot,
                                                       labels=c("(d)",' (e)'),
                                                       hjust = -2.5,#More negative values move the label further to the right on the plot canvas
                                                       vjust=3.5, #More positive values move the label further down on the plot canvas
                                                       font.label = list(size = 10, color = "black", face = "plain", family = "serif"),
                                                       widths = c(0.55,0.45),
                                                       nrow = 1,ncol = 2),
                                             heights = c(0.4,0.6),
                                             nrow = 2)
                                             
fa_individual_specialization.plot

#ggsave

### 4 ***************************************************************************
###      *** correlated specialization index between IS and FA ***
###  ---------------------------------------------------------------------------
head(is_specialization_niche)#esimated from "FA_hypervolume.R" code file
head(fa_niche)
 # merge data
is_fa_specialization<-is_specialization_niche[,c('individual_ID','specialization_index')]%>%
  dplyr::rename(is_specialization_index=specialization_index)%>%
  dplyr::left_join(fa_niche[,c('individual_ID','specialization_index')]%>%
                     dplyr::rename(fa_specialization_index=specialization_index),by='individual_ID')%>%
  na.omit()

summary(lm(is_specialization_index~fa_specialization_index,data = is_fa_specialization))#2.666   0.0103 *

is_fa_specialization%>%ggplot(aes(fa_specialization_index,is_specialization_index))+
  geom_point(size=0.8)+
  geom_smooth(method = 'lm')+
  annotate("text", x = 0.65, y = 0.90, label = substitute(paste(italic(y) == a,' + ',b, italic(x), ';  ',italic(r^2) == c), 
                                                      list(a ='0.54',b = "0.26",c='0.41')),hjust=0,
           size=2,family = "serif", fontface = "plain")+
  scale_x_continuous('Inidividual specialization index\nbased on fatty acids data',limits = c(0.54,0.98),breaks = seq(0,1,0.1))+
  scale_y_continuous('Inidividual specialization index\nbased on isotopes data',limits = c(0.55,0.92),breaks = seq(0,1,0.1))+
  theme_bw()+
  theme(text = element_text(size = 10))->specialization_is_fa_plot
specialization_is_fa_plot

#ggsave




