# load packages
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggpubr)


### 1   *******************************************************************
###            ***  loading data ***
###    -------------------------------------------------------------------
### loading data
data<-read.csv('data.csv')
head(data)
dim(data)

  
#### 2 ************************************************************
##           ***  test differences  ***
#### ==============================================================
# (1)testing between tissues
head(data)[1:2,]
data_test<-data%>%dplyr::select(tissue,d15N,d13C)%>% 
  #formate data
  tidyr::pivot_longer(cols = c(d15N,d13C),names_to='isotope_factor',values_to='value')%>%as.data.frame()
head(data_test)

isotopeFr<-c(unique(data_test$isotope_factor))
test_IS_tissue<-{}
#test_CE_sex_forPlot<-{}
for (i in 1:length(isotopeFr)) {
  temp0<-data_test%>%dplyr::filter(isotope_factor==isotopeFr[i])
  kt<-with(data=temp0,ks.test(value,'pnorm',mean=mean(value,na.rm=T),sd=sd(value,na.rm=T))) # normal distribution test
  if(kt$p.value>0.05){ ## satisfy normality
    # anova test
    lm_is<-lm(value~tissue,data=temp0)
    is.ano<-anova(lm_is)  # anova test
    #is.hsd.lm<-agricolae::HSD.test(lm_is, 'tissue') #Tukey HSD
    is.hsd.lm.fr<- data.frame(isotope_factor=paste0(isotopeFr[i]),test='ANOVA',statisticV=is.ano$`F value`[1],P=is.ano$`Pr(>F)`[1])
    
    test_IS_tissue<-rbind(test_IS_tissue,is.hsd.lm.fr)
    }else{ # reject normality
    kw.t<-with(data=temp0,kruskal.test(x=value,g=tissue)) #K-W test (Kruskal-Wallis test)
    is.hsd.lm.fr<- data.frame(isotope_factor=paste0(isotopeFr[i]),test='Kruskal-Wallis test',statisticV=kw.t$statistic[[1]],P=kw.t$p.value)
    
    test_IS_tissue<-rbind(test_IS_tissue,is.hsd.lm.fr)
  }
}
test_IS_tissue

### 3   *******************************************************************
###            ***  plot individual IS variations ***
###    -------------------------------------------------------------------
library(patchwork)

#(1) transform data frame to estimate mean and sd values
is_indiviual_mn_sd<-data%>%
  dplyr::group_by(individual_ID)%>%
  dplyr::summarise(N15_range=max(d15N)-min(d15N),N15_mn=mean(d15N),N15_sd=sd(d15N),
                   C13_range=max(d13C)-min(d13C),C13_mn=mean(d13C),C13_sd=sd(d13C))%>%as.data.frame
  
is_indiviual_mn_sd

# (2) plot mean and sd for each individuals
is_mnSd_plot<-ggplot(data = is_indiviual_mn_sd,aes(x=C13_mn,y=N15_mn))+
  geom_errorbar(aes(xmin=C13_mn-C13_sd, xmax=C13_mn+C13_sd),color='red', width=.05,size=0.1,#alpha=0.5,
                position=position_dodge(.9))+
  geom_errorbar(aes(ymin=N15_mn-N15_sd, ymax=N15_mn+N15_sd), width=.05,color='blue',size=0.1,#alpha=0.5,
                position=position_dodge(.9))+
  geom_point(stat="identity",size=1.5,color="grey50")+#position=position_dodge(0.9)
  xlab(expression(paste(delta^{13},'C')))+
  ylab(expression(paste(delta^{15},'N')))
is_mnSd_plot

#(3) add mean and sd by tissues to individual IS plot
is_tissues_NC<-data%>%dplyr::group_by(tissue)%>%
  dplyr::summarise(n=n(),
                   d15N_mn=sprintf(mean(d15N),fmt = '%#.2f'),
                   d15N_sd=sprintf(sd(d15N),fmt = '%#.2f'),
                   d13C_mn=sprintf(mean(d13C),fmt = '%#.2f'),
                   d13C_sd=sprintf(sd(d13C),fmt = '%#.2f'),
                   d15N_range=paste0(paste(sprintf(min(d15N),fmt = '%#.2f'),sprintf(max(d15N),fmt = '%#.2f'),sep='~')),
                   d15N_span=sprintf(max(d15N)-min(d15N),fmt = '%#.2f'),
                   d15N_mn.sd=paste(sprintf(mean(d15N),fmt = '%#.2f'),sprintf(sd(d15N),fmt = '%#.2f'),sep='±'),
                   d13C_range=paste0(paste(sprintf(min(d13C),fmt = '%#.2f'),sprintf(max(d13C),fmt = '%#.2f'),sep='~')),
                   d13C_span=sprintf(max(d13C)-min(d13C),fmt = '%#.2f'),
                   d13C_mn.sd=paste(sprintf(mean(d13C),fmt = '%#.2f'),sprintf(sd(d13C),fmt = '%#.2f'),sep='±')
  )%>%
  dplyr::mutate_at(vars(n,d15N_mn,d15N_sd,d13C_mn,d13C_sd,d15N_span,d13C_span),as.numeric)
is_tissues_NC

#plot
is_mnsd_plot_tissue<-is_mnSd_plot+
  geom_errorbar(aes(x=d13C_mn,y=d15N_mn,xmin=d13C_mn-d13C_sd, xmax=d13C_mn+d13C_sd),data=is_tissues_NC,color='black', width=.05,size=0.5,
                position=position_dodge(.9))+
  geom_errorbar(aes(x=d13C_mn,y=d15N_mn,ymin=d15N_mn-d15N_sd, ymax=d15N_mn+d15N_sd),data=is_tissues_NC, width=.05,color='black',size=0.5,
                position=position_dodge(.9))+
  geom_point(aes(x=d13C_mn,y=d15N_mn,shape=tissue),data=is_tissues_NC,stat="identity",color='black',size=2.5)+
  scale_shape_manual('Tissue',values = c(19,17,15),labels = c("Digestive gland", "Ovary", "Mantle muscle"))+
  theme(legend.position = c(0.85,0.15))

is_mnsd_plot_tissue  

#(4) plot C13 density distribution using kernel method for each individual
#summary(data)
cc <- scales::seq_gradient_pal("darkgrey", "lightgrey", "Lab")(seq(0,1,length.out=100))

is_C13_density_plot<-ggplot(data = data, aes(x = d13C)) + 
  geom_density(aes(color=factor(individual_ID)),alpha = 0.5,show.legend = F,linewidth=0.2) + 
  scale_colour_manual(values=cc)+
  geom_density(data = is_indiviual_mn_sd, aes(x = C13_mn),color='blue',linewidth=1,show.legend = F)+
  theme_void()
  
is_C13_density_plot

#(5) plot N15 density distribution using kernel method for each individual
is_N15_density_plot<-ggplot(data = data, aes(x = d15N)) + 
  geom_density(aes(color=factor(individual_ID)),alpha = 0.4,show.legend = F,linewidth=0.2) + 
  scale_colour_manual(values=cc)+
  geom_density(data = is_indiviual_mn_sd, aes(x = N15_mn),color='red',linewidth=1,show.legend = F) + 
  theme_void()+coord_flip()
  
is_N15_density_plot

#(6) combine mean+sd plot with IS density distribution plots 
is_plot_density<-is_C13_density_plot + 
  patchwork::plot_spacer() + #The function plot_spacer() adds an empty plot to the top right corner.
  (is_mnsd_plot_tissue+theme(plot.margin = unit(c(0,0,5,0), "pt"))) + (is_N15_density_plot +theme(plot.margin = unit(c(0,5,0,0), "pt")))+
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1.5), heights = c(1.5, 4))
is_plot_density


##  =================================================================
##    *************************************************************
##            plotting for the isotopes by maturity stages
##  ----------------------------------------------------------------
head(data)
anova(lm(d15N~mat,data = subset(data,tissue%in%'digestive')))#5.6451 0.0006224 ***

data%>%dplyr::filter(tissue=='digestive')%>%
  ggplot(aes(mat,d15N))+
  stat_boxplot(geom ='errorbar',aes(group=mat),width=0.1)+ 
  geom_boxplot(aes(fill=mat,group=mat),width=0.4,outlier.shape = NA,show.legend = F)+
  #add mean value in the boxplot
  stat_summary(fun.y = "mean",geom = "point",shape=21,size=1.5)+
  annotate("text", x = 1.5, y = 13.5, 
           label = substitute(paste(italic(F) == a,',  ',italic(p)==b),  #italic(r^2) == a,";  ", 
                              list(a ='5.65',b = " 0.0006")),hjust=0,
           size=3,family = "serif", fontface = "plain")+
  scale_y_continuous(name=expression(paste(delta^{15},'N')),limits = c(10,14),breaks = seq(0,1000,1),expand = c(0,0))+
  labs(x='Maturity stage')+
  theme_bw()+
  ggtitle("Digestive gland")+
  theme(plot.title = element_text(color="black", size=10, face="bold"))->dg_d15N_plot#"bold.italic"
dg_d15N_plot

### plot figure S2b
# test difference between maturity for ovary d15N
anova(lm(d15N~mat,data = subset(data_energy_is_fa,tissue%in%'ovary'))) #7.4656 5.81e-05 ***

data%>%dplyr::filter(tissue=='ovary')%>%
  ggplot(aes(mat,d15N))+
  stat_boxplot(geom ='errorbar',aes(group=mat),width=0.1)+ 
  geom_boxplot(aes(fill=mat,group=mat),width=0.4,outlier.shape = NA,show.legend = F)+
  #add mean value in the boxplot
  stat_summary(fun.y = "mean",geom = "point",shape=21,size=1.5)+#,fill="grey"
  #geom_smooth(size=0.8,color='blue')+
  annotate("text", x = 1.5, y = 16.2, label = substitute(paste(italic(F) == a,',  ',italic(p)==b),  #italic(r^2) == a,";  ", 
                                                         list(a ='7.46',b = " 5.81e-05")),hjust=0,
           size=3,family = "serif", fontface = "plain")+
  scale_y_continuous(name=expression(paste(delta^{15},'N')),limits = c(12.3,16.5),breaks = seq(0,1000,1),expand = c(0,0))+
  labs(x='Maturity stage')+
  theme_bw()+
  ggtitle("Ovary")+
  theme(plot.title = element_text(color="black", size=10, face="bold"))->ovary_d15N_plot
ovary_d15N_plot

####################################################################
#####   *********************************************************
#(7) combine mean+sd plot with IS distribution plots by maturity (dg+ovary)
IS_combine.plot<-ggarrange(is_plot_density,
                           ggarrange(dg_d15N_plot,ovary_d15N_plot,
                                     labels=c("(b)",'(c)'),
                                     hjust = -1.5,#More negative values move the label further to the right on the plot canvas
                                     vjust=2, #More positive values move the label further down on the plot canvas
                                     font.label = list(size = 10, color = "black", face = "plain", family = "serif"),
                                     nrow = 2,ncol = 1),
                           labels=c('(a)'),
                           hjust = -3.5,#More negative values move the label further to the right on the plot canvas
                           vjust=4.5, #More positive values move the label further down on the plot canvas
                           font.label = list(size = 10, color = "black", face = "plain", family = "serif"),
                           widths = c(0.65,0.35),
                           nrow = 1,ncol = 2)
IS_combine.plot

#save







