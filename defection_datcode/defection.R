###examining ONT data for the SARS integration manuscript

###prep env
library(tidyverse)
library(gridExtra)
library(data.table)
library(ggridges)
library(ggirl)
theme_set(theme_classic())

###read data
IndsPath="./"
QhistList=list.files(IndsPath,pattern="*.qhist.txt")
for( x in QhistList) {
  name=gsub(".qhist.txt","",x)
  df=read.delim(file=paste0(IndsPath,x)) %>% 
    mutate(Read=name,bps=n(),tq=sum(Read1_linear),mq=tq/bps)
  assign(paste0(name,"_iqh"),df) 
}
NPChim_Ind_Qhist=as.data.frame(rbindlist(mget(ls(pattern = "_iqh")),fill=T))
#
NPChim_Ind_Qhist %>% 
  ggplot()+
  geom_ridgeline(aes(x=X.BaseNum,y=fct_reorder(Read,bps),height=log2(Read1_linear),color=mq),fill=NA)+
  scale_color_viridis_c("Mean read quality", option="cividis")+
  theme(axis.text.y = element_blank(),axis.ticks.y=element_blank())+
  ggtitle("Base quality for each of 61 chimeric ONT reads")+xlab("Position")+ylab("Quality (log2) for each read")

Strings=NPChim_Ind_Qhist %>% 
  mutate(rand=rnorm(n(),mean=1,sd=0.1),newq=Read1_linear*rand) %>% 
  mutate(newq_spline=predict(smooth.spline(newq,nknots=100000))$y) %>% 
  ggplot()+
  geom_ridgeline(aes(x=log10(X.BaseNum),y=fct_reorder(Read,bps),height=log2(newq_spline),color=mq),fill=NA,size=1)+
  scale_color_viridis_c(option="cividis")+
  theme_void()+theme(legend.position="none")

#ggsave("Defection.png",Strings,units="cm",width=30,height=50,bg="transparent")
#ggsave("Defection.svg",Strings,units="cm",width=30,height=50,bg="transparent")
#ggsave("Defection.pdf",Strings,units="cm",width=30,height=50,bg="transparent")
