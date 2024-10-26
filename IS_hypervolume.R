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

### ======================================================================
###    *******************************************************************
###         1.  estimate individual hypervolume using gausian 
###    -------------------------------------------------------------------

#(1) esitimate individual hypervolume

ia_ind<-unique(data$individual_ID)
is_gausian_0.95<-{}
is_hv_ind<-list()
for (i in 1:length(ia_ind)) {
  print(i)
  is_ind_i<-data[which(data$individual_ID%in%ia_ind[i]),]%>%dplyr::select(d15N_scaled,d13C_scaled)
  is_hv_ind[[i]] = hypervolume_gaussian(na.omit(is_ind_i),name = 'gaussian')
  val_0.95<-get_volume(is_hv_ind[[i]])%>%as.data.frame()
  temp0<-data.frame(index=ia_ind[i],volume_0.95=val_0.95[1,])
  is_gausian_0.95<-rbind(is_gausian_0.95,temp0)
}

##(2) estimate the union of hypervolume to calculate specialization index
is_hv_ind_join<-do.call(hypervolume_join,is_hv_ind)#convert the list to a HypervolumeList

hv1 <- is_hv_ind_join@HVList[[1]]
for(j in 2:length(is_hv_ind_join@HVList)) {
  hv2 <- is_hv_ind_join@HVList[[j]]
  hvSet <- hypervolume_set(hv1, hv2, check.memory=FALSE, verbose=FALSE)
  hv1 <- hvSet@HVList$Union
}
# estimate individual specialization index as individual niche divided by the population hypervolume
is_union_volume0<-get_volume(hv1)%>%as.data.frame()
is_specialization_niche<-is_gausian_0.95%>%dplyr::mutate(is_union_volume=is_union_volume0[1,],
                                          specialization_index=1-volume_0.95/is_union_volume)#used for polar plot

### (3)  estimate the individual hypervolume overlap 
is_hv_ind ## the list of hypervolume for each individual
is_Jaccard_similarity<-{} #The Jaccard Similarity will be 0 if the two sets don't share any values and 1
                           #if the two sets are identical. 

for (i in 1:length(is_hv_ind)) {
  hv1=is_hv_ind[[i]]
  ia_index_i=ia_ind[i]
  z=i+1
  if(z<66){
    for (j in z:length(is_hv_ind)) {
      hv2=is_hv_ind[[j]]
      ia_index_j=ia_ind[j]
      #Computes the intersection, union, and unique components of two hypervolumes
      hv_set <- hypervolume_set(hv1, hv2, check.memory=FALSE)
      jac_sim<-hv_set@HVList$Intersection@Volume/hv_set@HVList$Union@Volume#niche overlap; 
      #same to hypervolume_overlap_statistics(hv_set)
      jac_nested<-hv_set@HVList$Intersection@Volume/min(hv1@Volume,hv2@Volume)#nestedness index
      #hypervolume_overlap_statistics(hv_set)
      temp0<-data.frame(comp_index=paste0(i,'_',j),
                        #set variable to store the pair index
                        comp_index_ind=paste0(ia_index_i,'_',ia_index_j), 
                        jac_similarity=jac_sim,
                        jac_nestedIndex=jac_nested)
      is_Jaccard_similarity<-rbind(is_Jaccard_similarity,temp0)
    }}}

### 3   *******************************************************************
###            *** plot individual specialization index ***
###    -------------------------------------------------------------------

#(1) polar plot individual specialization index distribution

head(is_specialization_niche)
#minus index by 0.3 to plot the polar chart
is_specialization_niche_03<-is_specialization_niche%>%
  dplyr::mutate(specialization_index1=specialization_index-0.3)
summary(is_specialization_niche_03)
str(is_specialization_niche_03)
is_specialization_sel<-is_specialization_niche_03%>%dplyr::select(index,specialization_index1)%>%
  dplyr::mutate(ID=as.numeric(row.names(.)),
                angle0=90-360*(as.numeric(row.names(is_specialization_niche_03))-0.5)/(nrow(is_specialization_niche_03)+6))%>%
  dplyr::mutate(hjust=ifelse(as.numeric(angle0)<0,ifelse(abs(angle0)>90,1,0),0),
                angle=ifelse(as.numeric(angle0)<0,ifelse(abs(angle0)>90,angle0+180,angle0),angle0))%>%
  dplyr::mutate(index_gp=case_when(specialization_index1<.15~'gp_45', #0.45
                                   specialization_index1<.3~'gp_60', #0.65
                                   specialization_index1<0.45~'gp_75',#0.8
                                   specialization_index1<0.6~'gp_90',#0.95
                                   TRUE~'gp_10'))

is_specialization_sel$index_gp<-factor(is_specialization_sel$index_gp,levels = c('gp_45','gp_60','gp_75','gp_90','gp_10'))
is_specialization_sel%>%ggplot()+
  geom_linerange(aes(xmin=0,xmax=66, y = 0.0),lty="solid", color = "black")+
  geom_linerange(aes(xmin=0,xmax=66, y = 0.15),lty="solid", color = "grey80",size=0.3)+
  geom_linerange(aes(xmin=0,xmax=66, y = 0.3),lty="solid", color = "grey80",size=0.3)+
  geom_linerange(aes(xmin=0,xmax=66, y = 0.4),lty=3, color = "blue")+
  geom_linerange(aes(xmin=0,xmax=66, y = 0.45),lty="solid", color = "grey80",size=0.3)+
  geom_linerange(aes(xmin=0,xmax=66, y = 0.6),lty="solid", color = "grey80",size=0.5)+
  #手动添加坐标及标题
  geom_text(x=-2,y=-0.02,label="0.30",color="black",size=2)+
  geom_text(x=-1.5,y=0.15,label="0.45",color="black",size=2)+
  #geom_text(x=-0.4,y=0.3,label="0.60",color="black",size=1.5)+
  geom_text(x=-1,y=0.3,label="0.60",color="black",size=2)+
  geom_text(x=-1,y=0.39,label="0.70",color="black",size=2)+
  geom_text(x=-1,y=0.45,label="0.75",color="black",size=2)+
  geom_text(x=-1,y=0.6,label="0.90",color="black",size=2)+
  #柱状图
  geom_col(aes(ID, specialization_index1, fill=index_gp), alpha = 0.9)+
  #极坐标转换
  coord_polar(direction=1)+
  #x,y轴范围确定
  scale_y_continuous(limits = c(-0.25,0.7),breaks = seq(0,1,0.3))+
  scale_x_continuous('Individual specialization index',limits = c(0,72))+
  #主题
  labs(title='Individual specialization index',subtitle = "(based on stable isotopes data)")+
  #主题
  theme_void()+
  theme(legend.position = 'none',
        plot.title = element_text(color="black", size=9, face="bold"),#bold.italic
        plot.subtitle = element_text(color="black", size=7.5, face="plain"),
        panel.grid.major = element_blank(),
        axis.line = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y.left = element_blank())->is_spe_index.plot
is_spe_index.plot

###(2) plot specialization index distribution by maturity stage, figure 3b
### combined niche data with maturity stage
head(data)
head(is_specialization_niche)

data$individual_ID<-as.integer(data$individual_ID)

#create dataframe
is_specialization_mat<-is_specialization_niche%>%
  dplyr::left_join(distinct(data[,c('individual_ID','mat')]),by=c('index'='individual_ID'))%>%
  dplyr::mutate(mat_num=case_when(mat=='III'~1,mat=='IV'~2,mat=='V'~3,mat=='VI'~4,mat=='VII'~5))

is_specialization_mat%>%ggplot(aes(x=mat_num, y=specialization_index)) + 
  stat_boxplot(geom ='errorbar',aes(group=mat),width=0.1)+ 
  geom_boxplot(aes(fill=mat,group=mat),width=0.4,outlier.shape = NA,show.legend = F)+
  #add mean value in the boxplot
  stat_summary(fun.y = "mean",geom = "point",shape=21,size=1.5)+#,fill="grey"
  annotate("text", x = 4, y = 0.95, label = substitute(paste(italic(F) == a,',  ',italic(p)==b),  #italic(r^2) == a,";  ", 
                                                       list(a ='3.37',b = " 0.014")),hjust=0,
           size=2,family = "serif", fontface = "plain")+
  geom_smooth(size=0.8,color='blue')+
  scale_x_continuous(name='Maturity stage',limits=c(0.5,5.5),breaks=seq(1,5,1),
                   labels=c('III','IV','V','VI','VII'))+
  scale_y_continuous(name='Individual specialization index',limits = c(0.4,1),breaks = seq(0,1000,.2))+
  theme_bw()->is_specialization_mat.plot
is_specialization_mat.plot

###(3) plot pairwise overlap distribution 
head(is_Jaccard_similarity)
is_Jaccard_similarity_mat%>%dplyr::group_by(mat_comb)%>%
  summarise(n=n())%>%as.data.frame()
sort(unique(Jaccard_similarity_mat$mat_comb))

is_Jaccard_similarity_mat%>%ggplot(aes(x=jac_similarity))+
  geom_histogram(aes(y=after_stat(c(count/sum(count))*100)),#count[group==1]/sum(count[group==1]))*100)),#
                 #count[group==2]/sum(count[group==2]))*100)),
                 linewidth=0.3,color='lightgrey',fill='lightblue',center = TRUE,breaks=seq(0,1,0.05),show.legend = F)+#breaks=seq(0.15,1,0.015), ,position=position_dodge2(width = 0.1),binwidth = 0.05,
  #scale_fill_manual(name="性别",values=c("black","white"))+
  #geom_vline(aes(xintercept=0.5),lty=2,size=1.5,alpha=0.6,color='red')+
  xlab("Parewise niche overlap for individuals")+ylab("Frequency (%)")+
  scale_x_continuous(expand = c(0, 0),limits = c(0.12,1),breaks = seq(-0.1,1,0.2))+
  scale_y_continuous(expand = c(0, 0),limits = c(0,18),breaks = seq(0,50,5))+
  theme_bw()+p8.8->is_Jaccard_similarity.plot

is_Jaccard_similarity.plot


#### (4) plot pairwise overlap for each sequential pair stage, figure 3d
###combined jaccard similarity with maturity
head(data_energy_is_fa)
str(data_energy_is_fa)
data_energy_is_fa$index<-as.character(data_energy_is_fa$index)
unique(is_Jaccard_similarity$comp_index_ind)
mat_ch<-c('III_III','IV_IV','V_V','VI_VI','VII_VII')

#combined 
is_Jaccard_similarity_mat<-is_Jaccard_similarity%>%
  dplyr::mutate(comp_index_ind1=comp_index_ind)%>%
  tidyr::separate(comp_index_ind1,c('index1','index2'),'_')%>%
  dplyr::left_join(data_energy_is_fa[,c('index','mat')],by=c('index1'='index'))%>%
  dplyr::left_join(data_energy_is_fa[,c('index','mat')],by=c('index2'='index'))%>%
  dplyr::mutate(mat1_comb=paste0(mat.x,'_',mat.y))%>%
  dplyr::filter(!mat1_comb%in%c("III_III","IV_IV","V_V","VI_VI","VII_VII",'IV_III','VII_III','VII_IV','VII_V','VII_VI','VI_III','VI_IV','VI_V','V_III','V_IV'))%>%
  dplyr::mutate(mat1_comb_num=case_when(mat1_comb=="III_IV"~1,mat1_comb=="III_V"~2,mat1_comb=="III_VI"~3,mat1_comb=="III_VII"~4,
                                        mat1_comb=="IV_V"~5,mat1_comb=="IV_VI"~6,mat1_comb=="IV_VII"~7,
                                        mat1_comb=="V_VI"~8,mat1_comb=="V_VII"~9,
                                        mat1_comb=="VI_VII"~10))
#calculate mean value
is_Jaccard_similarity_mat_mean<-is_Jaccard_similarity_mat%>%dplyr::group_by(mat1_comb_num)%>%#,mat1_comb_num mat1_comb
  dplyr::summarise(jac_similarity_mean=mean(jac_similarity))

is_Jaccard_similarity_mat%>%ggplot(aes(x=mat1_comb_num, y=jac_similarity)) + 
  #add errorbar for the box, and width argument to set the width of the whisker end caps
  stat_boxplot(geom ='errorbar',aes(group=mat_comb),width=0.1)+ 
  geom_boxplot(aes(fill=mat.x,group=mat_comb),width=0.4,outlier.shape = NA,show.legend = F)+
  #add mean value in the boxplot
  stat_summary(fun.y = "mean",geom = "point",shape=21,size=1.5)+
  geom_smooth(aes(x=mat1_comb_num, y=jac_similarity_mean),formula = y~x,
              data = is_Jaccard_similarity_mat_mean,size=0.8,color='blue')+
  scale_x_continuous(name='Pairwise saturity stages',limits = c(0.5,10.5),breaks = seq(1,10,1),
                     labels = c("III~IV","III~V","III~VI","III~VII","IV~V","IV~VI","IV~VII","V~VI","V~VII","VI~VII"))+
  scale_y_continuous(name='Parewise niche overlap',limits = c(0.15,1),breaks = seq(0,1000,.2))+
  theme_bw()->is_Jaccard_similarity_mat.plot

is_Jaccard_similarity_mat.plot

## combine plots
is_specialization_index.plot<-ggarrange(ggarrange(is_spe_index.plot,is_specialization_mat.plot,
                                                          labels=c("(a)",' (b)'),
                                                          hjust = -2.5,#More negative values move the label further to the right on the plot canvas
                                                          vjust=3.5, #More positive values move the label further down on the plot canvas
                                                          font.label = list(size = 10, color = "black", face = "plain", family = "serif"),
                                                          widths = c(0.6,0.4),
                                                          nrow = 1,ncol = 2),
                                                ggarrange(is_Jaccard_similarity.plot,is_Jaccard_similarity_mat.plot,
                                                          labels=c("(c)",'(d)'),
                                                          hjust = -3.5,#More negative values move the label further to the right on the plot canvas
                                                          vjust=3.5, #More positive values move the label further down on the plot canvas
                                                          font.label = list(size = 10, color = "black", face = "plain", family = "serif"),
                                                          #widths = c(0.45,0.4),
                                                          nrow = 1,ncol = 2),
                                                nrow = 2,ncol = 1,
                                                heights = c(0.6,0.4)
                                                )
is_specialization_index.plot

#ggsave

### 5   LMM  ML\gonad energy\sst\chla ~specialization

#  ##********************************************************************
##                *** LMM of gonadEnergy\sst\chla on specialization ****
###----------------------------------------------------------------------

library(lme4)
library(lmerTest)
library(effects)
library(MuMIn)

#===========================================================
# read environmental data
Env<-read.csv("env.csv")%>%
  dplyr::filter(individual_ID%in%c(unique(data$individual_ID)))%>%
  dplyr::select(individual_ID,caught_month,long.dec,lat.dec,mat,
                MLD,bottomT,ssh,surface_sal,surface_temperature,CHL)
head(Env)

# combined specialization index with reproductive energy and enviromental data
head(data)
head(is_specialization_niche)
data$individual_ID<-as.integer(data$individual_ID)

is_energy_specialization_env<-is_specialization_niche%>%
  dplyr::left_join(distinct(data_energy_is_fa[,c('individual_ID','mat','gonadEnergy_percentage')]),by='individual_ID')%>%
  dplyr::left_join(IEnv[,c(1,6:10)],by='individual_ID')%>%
  dplyr::mutate(gonadEnergy_percentage_scaled=scale(gonadEnergy_percentage),
                MLD_scaled=scale(MLD),bottomT_scaled=scale(bottomT),ssh_scaled=scale(ssh), 
                surface_sal_scaled=scale(surface_sal),surface_temperature_scaled=scale(surface_temperature),
                chla_scaled=scale(CHL))%>%
  dplyr::filter(!is.na(MLD_scaled))

head(is_energy_specialization_env)

####
library(buildmer)
library(lme4)
# explanatory variables selection 
f<-specialization_index~gonadEnergy_percentage_scaled+
  MLD_scaled+ 
  bottomT_scaled+ 
  surface_sal_scaled+
  surface_temperature_scaled+
  chla_scaled+
  (1|mat1)
m <- buildmer(f,data=is_energy_specialization_env,
              buildmerControl=buildmerControl(direction='order',
                                              args=list(control=lmerControl(optimizer='bobyqa'))))
m@model
f1 <- formula(m@model)

m1 <- buildmer(f1,data=is_energy_specialization_env,
               buildmerControl=list(direction='backward',
                                    args=list(control=lmerControl(optimizer='bobyqa'))))
summary(m1)

m1@summary
m1@p
###  lmm
is_energy_env.lme_all<-lmerTest::lmer(specialization_index~gonadEnergy_percentage_scaled+
                                        bottomT_scaled+ 
                                        surface_temperature_scaled+(1|mat),
                                      data=is_energy_specialization_env,control = lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))) 
summary(is_energy_env.lme_all)
MuMIn::r.squaredGLMM(is_energy_env.lme_all)

###predictions
## predicted effect
is_energy_specialization_env$pred_all_rand <-predict(is_energy_env.lme_all, type='response')  

specialization_gonadEnergy_plot<-is_energy_specialization_env%>%ggplot()+
  geom_point(aes(gonadEnergy_percentage,specialization_index),size=0.5)+
  geom_smooth(aes(gonadEnergy_percentage,pred_all_rand),method = 'lm',color='blue',linewidth=1)+
  annotate("text", x =14.5, y = 0.95, label = substitute(paste(italic(β) == a,',  ',italic(t)==b,',  ',italic(p)==c),  #italic(r^2) == a,";  ", 
                                                      list(a ='0.023',b = " 2.67",c=" 0.0009")),hjust=0,
           size=2,family = "serif", fontface = "plain")+
  scale_x_continuous(name='Reproductive energy (%)',limits = c(0,45),breaks = seq(0,100,10))+
  scale_y_continuous(name='Individual specialization index',limits = c(0.4,1),breaks = seq(0,1000,.2))+
  theme_bw()
specialization_gonadEnergy_plot

specialization_surface_temperature_plot<-is_energy_specialization_env%>%ggplot()+
  geom_point(aes(surface_temperature,specialization_index),size=0.5)+
  geom_smooth(aes(surface_temperature,pred_all_rand),method = 'lm',color='blue',linewidth=1)+
  annotate("text", x =12.2, y = 0.95, label = substitute(paste(italic(β) == a,',  ',italic(t)==b,',  ',italic(p)==c),  #italic(r^2) == a,";  ", 
                                                       list(a ='0.038',b = " 3.04",c=" 0.0035")),hjust=0,
           size=2,family = "serif", fontface = "plain")+
  scale_x_continuous(name='Sea surface temperature (℃)',limits = c(11,15),breaks = seq(0,100,1))+
  scale_y_continuous(name='Individual specialization index',limits = c(0.4,1),breaks = seq(0,1000,.2))+
  theme_bw()
specialization_surface_temperature_plot

specialization_bt_plot<-is_energy_specialization_env%>%ggplot()+
  geom_point(aes(bottomT,specialization_index),size=0.5)+
  geom_smooth(aes(bottomT,pred_all_rand),method = 'lm',color='blue',linewidth=1)+
  annotate("text", x =6.5, y = 0.95, label = substitute(paste(italic(β) == a,',  ',italic(t)==b,',  ',italic(p)==c),  #italic(r^2) == a,";  ", 
                                                       list(a ='0.029',b = " 2.45",c=" 0.017")),hjust=0,
           size=2,family = "serif", fontface = "plain")+
  scale_x_continuous(name = 'Bottom temerature (℃)',limits = c(5.5,9),breaks = seq(0,500,1))+
  scale_y_continuous(name='Individual specialization index',limits = c(0.4,1),breaks = seq(0,1000,.2))+
  theme_bw()
specialization_bt_plot

## save plots
Specialization_lmer.plot<-ggarrange(specialization_gonadEnergy_plot,specialization_surface_temperature_plot,specialization_bt_plot,
                                                labels=c("(a)",'(b)','(c)'),
                                                hjust = -3.5,#More negative values move the label further to the right on the plot canvas
                                                vjust=4.5, #More positive values move the label further down on the plot canvas
                                                font.label = list(size = 10, color = "black", face = "plain", family = "serif"),
                                                #widths = c(0.45,0.4),
                                                nrow = 1,ncol = 3)
Specialization_lmer.plot
#ggsave



