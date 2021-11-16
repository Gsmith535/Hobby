#!/usr/bin/Rscript
#Packages
library(tidyverse)
library(data.table)
library(ggridges)
library(gridExtra)
theme_set(theme_void())

#All size combo - multiple cutoff lengths for soft-shredding put back together and compared using k=4
INPath="./"
FilList=list.files(INPath,pattern="*self.csv")

### Read in MAG data
Check=fread(file=paste0(INPath,"check.txt"))
colnames(Check)=c("Access","Lineage","UID","Numb_Genomes","Numb_markers","Numb_markersets","0","1","2","3","4","5ormore","Complete","Contam","Hetero")
Check=Check %>% 
  mutate(Genome=gsub("\\.[12]","",Access))

GTDB=fread(file=paste0(INPath,"gtdb.tsv"))
GTDB=GTDB %>% 
  mutate(Genome=gsub("\\.[12]","",user_genome))

Genomes=merge(Check,GTDB,by="Genome",all=T) %>% 
  separate(classification,c("GTDB_Kin","GTDB_Phy","GTDB_Cla","GTDB_Ord","GTDB_Fam","GTDB_Gen","GTDB_Spe"),sep=";",remove=F)

### Read in sourmash data
for( x in FilList) {
  name=gsub("_[12]_.*","",x)
  df=fread(file=paste0(INPath,x),sep=",")
  rownames(df)=colnames(df)
  df2=df %>% 
    rownames_to_column(var="Seq1") %>% 
    mutate(Genome=name) %>% 
    select(Genome,Seq1,everything()) %>% 
    gather(Seq2,Sim,-Genome,-Seq1) %>% 
    filter(Seq1!=Seq2) %>% 
    mutate(Seq1_tmp=gsub(".*frag","",Seq1),Seq1_tmp=gsub("_.*pos_","&",Seq1_tmp),Seq1_tmp=gsub("_frg_","&",Seq1_tmp),Seq1_tmp=gsub("-","&",Seq1_tmp),Seq1_tmp=ifelse(!grepl("&",Seq1_tmp),NA,Seq1_tmp)) %>% 
    mutate(Seq2_tmp=gsub(".*frag","",Seq2),Seq2_tmp=gsub("_.*pos_","&",Seq2_tmp),Seq2_tmp=gsub("_frg_","&",Seq2_tmp),Seq2_tmp=gsub("-","&",Seq2_tmp),Seq2_tmp=ifelse(!grepl("&",Seq2_tmp),NA,Seq2_tmp)) %>% 
    separate(Seq1_tmp,c("Seq1_Frag_size","Seq1_Frag_position_start","Seq1_Frag_position_end","Seq1_Frag_numb","Seq1_Frag_total"),sep="&") %>% 
    separate(Seq2_tmp,c("Seq2_Frag_size","Seq2_Frag_position_start","Seq2_Frag_position_end","Seq2_Frag_numb","Seq2_Frag_total"),sep="&") %>% 
    select(-matches("Frag_size")) %>% 
    mutate_at(vars(matches("_Frag_")),as.character) %>% 
    mutate_at(vars(matches("_Frag_")),as.integer) %>% 
    select(everything(),Sim) 
  df4=merge(Genomes,df2,by="Genome",all.x = F) 
  assign(paste0(name,"_dat"),df4) 
  df5=df4 %>% 
    select(Genome,Seq1_Frag_position_end,Sim) %>% 
    group_by(Genome) %>% 
    summarise(Genome_size=max(Seq1_Frag_position_end),Sim_min=min(Sim),Sim_mean=mean(Sim),Sim_med=median(Sim),Sim_sd=sd(Sim),Sim_3xsd=3*sd(Sim),Sim_cut3xsd=mean(Sim)-(3*sd(Sim)),Rem_cut3xsd=sum(Sim<(mean(Sim)-(3*sd(Sim)))),Kep_cut3xsd=sum(Sim>mean(Sim)-(3*sd(Sim))),Sim_se=(sd(Sim)/sqrt(n())),Sim_comps=n()) %>% 
    select(Genome,Genome_size,everything())
  assign(paste0(name,"_sumry"),df5)
  rm(df,df2,df4,df5,name,x)
}

### Combine some objects

Summaries_temp=as.data.frame(rbindlist(mget(ls(pattern = "_sumry")),fill=T))
Soursumrys=merge(Genomes,Summaries_temp,by="Genome")

Soursumrys_l=Soursumrys %>%
  mutate(Rem_cut3xsd=Rem_cut3xsd/2,Kep_cut3xsd=Kep_cut3xsd/2,Sim_comps=Sim_comps/2) %>% 
  mutate(Kep_frac_cut3xsd=Kep_cut3xsd/(Rem_cut3xsd+Kep_cut3xsd)) %>% 
  select(-Sim_se) %>% 
  gather("Stat","val",Sim_min:Sim_cut3xsd,-Sim_comps) %>% 
  mutate(Stat=factor(Stat,levels=unique(Stat))) 

Sourdat=as.data.frame(rbindlist(mget(ls(pattern = "_dat"))))

Sourdat_l=Sourdat %>% 
  mutate(Sim_oppadd=(Sim*-1)+2) %>% 
  gather("Sim_cat","Sim_m",c(Sim,Sim_oppadd)) %>% 
  mutate(Sim_s=Sim_m^2) %>% 
  mutate(Nudge=sqrt(log(Seq1_Frag_total))/10) %>% 
  mutate(Rand=rnorm(Sim_m,mean=1,sd=0.5))

rm(list=c(ls(pattern = "_dat"),ls(pattern = "_sumry"),"Summaries_temp","Sourdat"))

### Look at variation within genomes
#p_ridg_noise = Sourdat_l %>% 
#  filter(Sim_cat=="Sim") %>%
#  ggplot()+
#  geom_density_ridges(aes(x=Sim_m,y=interaction(Genome,classification),fill=Seq1_Frag_total),bandwidth=0.005,rel_min_height=0.0001,scale=5,color="#0571b0",size=0.1)+
#  theme(legend.position = "none")
#p_ridg_noisier = Sourdat_l %>% 
#  sample_frac(0.01) %>% 
#  ggplot()+
#  geom_density_ridges(aes(x=Sim_m*Rand,y=interaction(Genome,classification),fill=Seq1_Frag_total),bandwidth=0.01,rel_min_height=0.0001,scale=5,color="#0571b0",size=0.1)+
#  theme(legend.position = "none")+
#  scale_fill_distiller(palette="BrBG")
#p_ridg_mir =  Sourdat_l %>% 
#  ggplot()+
#  geom_density_ridges(aes(x=Sim_m,y=interaction(Genome,classification),fill=Seq1_Frag_total),bandwidth=0.03,rel_min_height=0.001,scale=5,color="#0571b0",size=0.1)+
#  theme(legend.position = "none")+
#  scale_fill_distiller(palette="BrBG")
#p_ridg_nudge = Sourdat_l %>% 
#  ggplot()+
#  geom_density_ridges(aes(x=Sim_m+Nudge,y=interaction(Genome,classification),fill=Seq1_Frag_total),bandwidth=0.03,rel_min_height=0.001,scale=3,color="#0571b0",size=0.1)+
#  theme(legend.position = "none")+
#  scale_fill_distiller(palette="BrBG")
#
#p_ridg = grid.arrange(p_ridg_mir,p_ridg_nudge,p_ridg_noise,p_ridg_noisier,nrow=2,ncol=2)
#
#p_ridg_noise_nudge = Sourdat_l %>% 
#  ggplot()+
#  geom_density_ridges(aes(x=(Sim_m*Rand)+(1/(1-Nudge))^4,y=interaction(Genome,classification),fill=Seq1_Frag_total),bandwidth=0.005,rel_min_height=0.0001,scale=5,color="#0571b0",size=0.1)+
#  theme(legend.position = "none")+
#  scale_fill_distiller(palette="BrBG")
#p_ridg_noisier_nudge = Sourdat_l %>% 
#  sample_frac(0.01) %>% 
#  ggplot()+
#  geom_density_ridges(aes(x=(Sim_m*Rand)+(1/(1-Nudge))^4,y=interaction(Genome,classification),fill=Seq1_Frag_total),bandwidth=0.01,rel_min_height=0.0001,scale=5,color="#0571b0",size=0.1)+
#  theme(legend.position = "none")+
#  scale_fill_distiller(palette="BrBG")
#
#p_ridg2 = grid.arrange(p_ridg_noise,p_ridg_noisier,p_ridg_noise_nudge,p_ridg_noisier_nudge,nrow=2,ncol=2)

#Best is p_ridg_noise_nudge
p_ridg_noise_nudge_keep_A = Sourdat_l %>% 
  ggplot()+
  geom_density_ridges(aes(x=(Sim_m*Rand)+(1/(1-Nudge))^4,y=interaction(Genome,classification),fill=Seq1_Frag_total),bandwidth=0.005,rel_min_height=0.0001,scale=5,color="gray67",size=0.1)+
  theme(legend.position = "none")+
  scale_fill_distiller(palette="BrBG")
p_ridg_noise_nudge_keep_A_ind = Sourdat_l %>% 
  ggplot()+
  geom_density_ridges(aes(x=(Sim_m*Rand)+(1/(1-Nudge))^4,y=interaction(Genome,classification),fill=Seq1_Frag_total),bandwidth=0.005,rel_min_height=0.0001,scale=5,color="gray67",size=0.05)+
  theme(legend.position = "none")+
  scale_fill_distiller(palette="BrBG")
#ggsave("Dissoc.pdf",p_ridg_noise_nudge_keep_A_ind,units="cm",width=10,height=10,bg="transparent")
#ggsave("Dissoc.png",p_ridg_noise_nudge_keep_A_ind,units="cm",width=10,height=10,bg="transparent")
#ggsave("Dissoc.svg",p_ridg_noise_nudge_keep_A_ind,units="cm",width=10,height=10,bg="transparent")

p_ridg_noise_nudge_keep_D = Sourdat_l %>% 
  ggplot()+
  geom_density_ridges(aes(x=(Sim_m*Rand)+(1/(1-Nudge))^4,y=fct_rev(interaction(Genome,classification)),fill=Seq1_Frag_total),bandwidth=0.005,rel_min_height=0.0001,scale=5,color="gray67",size=0.05)+
  theme(legend.position = "none")+
  scale_x_reverse()+
  scale_fill_distiller(palette="BrBG",direction=1)

p_ridg_noise_nudge_keep_B = Sourdat_l %>% 
  ggplot()+
  geom_density_ridges(aes(x=(Sim_m*Rand)+(1/(1-Nudge))^4,y=interaction(Genome,classification),fill=as.numeric(interaction(Genome,classification)),color=as.numeric(interaction(Genome,classification))),bandwidth=0.005,rel_min_height=0.0001,scale=3,size=0.3)+
  theme(legend.position = "none", plot.background = element_rect(color = "#238b45", fill = "black"))+
  scale_fill_distiller(palette="YlOrBr",direction=1)+scale_color_distiller(palette="BuPu",direction=1)
p_ridg_noise_nudge_keep_B_ind = Sourdat_l %>% 
  ggplot()+
  geom_density_ridges(aes(x=(Sim_m*Rand)+(1/(1-Nudge))^4,y=interaction(Genome,classification),fill=as.numeric(interaction(Genome,classification)),color=as.numeric(interaction(Genome,classification))),bandwidth=0.005,rel_min_height=0.0001,scale=3,size=0.1)+
  theme(legend.position = "none", plot.background = element_rect(color = NA, fill = "black"))+
  scale_fill_distiller(palette="YlOrBr",direction=1)+scale_color_distiller(palette="BuPu",direction=1)
#ggsave("DunDissoc.png",p_ridg_noise_nudge_keep_B_ind,units="cm",width=10,height=10,bg="transparent")
#ggsave("DunDissoc.svg",p_ridg_noise_nudge_keep_B_ind,units="cm",width=10,height=10,bg="transparent")

p_ridg_noise_nudge_keep_C = Sourdat_l %>% 
  ggplot()+
  geom_density_ridges(aes(x=(Sim_m*Rand)+(1/(1-Nudge))^4,y=interaction(Genome,classification)),bandwidth=0.005,rel_min_height=0.0001,scale=3,fill="black",color="white",size=0.3)+
  theme(legend.position = "none", plot.background = element_rect(color = "black", fill = "black"))
p_ridg_noise_nudge_keep_C_ind = Sourdat_l %>% 
  ggplot()+
  geom_density_ridges(aes(x=(Sim_m*Rand)+(1/(1-Nudge))^4,y=interaction(Genome,classification)),bandwidth=0.005,rel_min_height=0.0001,scale=3,fill="black",color="white",size=0.1)+
  theme(legend.position = "none", plot.background = element_rect(color = "black", fill = "black"))
#ggsave("JoyDissoc.png",p_ridg_noise_nudge_keep_C_ind,units="cm",width=10,height=10,bg="transparent")
#ggsave("JoyDissoc.svg",p_ridg_noise_nudge_keep_C_ind,units="cm",width=10,height=10,bg="transparent")

p_grid = grid.arrange(p_ridg_noise_nudge_keep_A,p_ridg_noise_nudge_keep_B,p_ridg_noise_nudge_keep_C,p_ridg_noise_nudge_keep_D,nrow=2,ncol=2)

#ggsave("MultDissoc.png",p_grid,units="cm",width=20,height=20,bg="transparent")
#ggsave("MultDissoc.svg",p_grid,units="cm",width=20,height=20,bg="transparent")
#ggsave("MultDissoc.pdf",p_grid,units="cm",width=20,height=20,bg="transparent")
