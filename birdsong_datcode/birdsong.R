###Packages
library(rtweet)
library(tidyverse)
library(zoo)
library(ggirl)

load("birdsong.Rdata")

#u1_fol=get_followers("#######",n=Inf,retryonratelimit = TRUE) ### removed to protect identities
#u2_fol=get_followers("#######",n=Inf,retryonratelimit = TRUE) ### removed to protect identities
#save.image("birdsong.Rdata")

u1_foldat=lookup_users(u1_fol$user_id)
u2_foldat=lookup_users(u2_fol$user_id)


u1_foldat %>% 
  select(account_created_at) %>% 
  separate(account_created_at," ",into=c("YMD","Time")) %>% 
  select(-Time) %>% 
  separate(YMD,"-",into=c("Year","Month","Day")) %>% 
  group_by(Year) %>% 
  tally() %>% 
  mutate(pr=n/sum(n))

u1_foldat %>% 
  select(account_created_at) %>% 
  separate(account_created_at," ",into=c("YMD","Time")) %>% 
  select(-Time) %>% 
  group_by(YMD) %>% 
  tally() %>% 
  separate(YMD,"-",into=c("Year","Month","Day"),remove=F) %>% 
  unite("MD",Month:Day,remove=F) %>% 
  #filter(Year>2012) %>% 
  ggplot()+
  geom_rect(aes(ymax=max(n)),xmin="03_01",xmax="03_31",ymin=0,fill="gray95",color=NA)+
  geom_line(aes(x=MD,y=n,color=as.numeric(as.character(Year)),group=Year),size=1)+
  scale_color_distiller(palette = "BuGn",direction=1)+
  theme_void()+theme(legend.position="none")



u1_rand=u1_foldat %>% 
  select(account_created_at) %>% 
  separate(account_created_at," ",into=c("YMD","Time")) %>% 
  select(-Time) %>% 
  group_by(YMD) %>% 
  tally() %>% 
  separate(YMD,"-",into=c("Year","Month","Day"),remove=F) %>% 
  unite("MD",Month:Day,remove=F) %>% 
  filter(Year>2016) %>% 
  mutate(Year=factor(Year,levels=c(2021,2020,2019,2018,2017,2016))) %>% 
  mutate(user="u1") %>% 
  mutate(rand=rnorm(n,mean=1,sd=0.1),nrand=n*rand,nrandroll=rollmean(nrand, 7, na.pad=T))
u2_rand=u2_foldat %>% 
  select(account_created_at) %>% 
  separate(account_created_at," ",into=c("YMD","Time")) %>% 
  select(-Time) %>% 
  group_by(YMD) %>% 
  tally() %>% 
  separate(YMD,"-",into=c("Year","Month","Day"),remove=F) %>% 
  unite("MD",Month:Day,remove=F) %>% 
  filter(Year>2016) %>% 
  mutate(Year=factor(Year,levels=c(2021,2020,2019,2018,2017,2016))) %>% 
  mutate(user="u2") %>% 
  mutate(rand=rnorm(n,mean=1,sd=0.1),nrand=sqrt(n)*rand,nrandroll=rollmean(nrand, 7, na.pad=T))

randdat=rbind(u1_rand,u2_rand) %>% 
  filter(!is.na(nrandroll)) %>% 
  filter(MD!="01_01",
         MD!="01_02",
         MD!="01_03") %>% 
  mutate(newn=predict(smooth.spline(nrandroll,nknots=1000))$y)

randdat %>%
  ggplot()+
  geom_ribbon(aes(x=MD,ymax=log2(nrandroll),fill=interaction(Year,user),group=interaction(Year,user)),ymin=0,alpha=0.9)+
  geom_line(aes(x=MD,y=log2(nrandroll),color=interaction(Year,user),group=interaction(Year,user)),alpha=0.3,stroke=1)+
  #scale_fill_manual(values=c("#006d2c","#2ca25f","#66c2a4","#b2e2e2","#bdc9e1","#74a9cf","#2b8cbe","#045a8d"))+
  scale_fill_manual(values=c("#006d2c","#2ca25f","#66c2a4","#99d8c9","#ccece6",
                             "#d0d1e6","#a6bddb","#74a9cf","#2b8cbe","#045a8d"))+
  scale_color_manual(values=c("#ccece6","#99d8c9","#66c2a4","#2ca25f","#006d2c",
                              "#045a8d","#2b8cbe","#74a9cf","#a6bddb","#d0d1e6"))+
  ##scale_fill_brewer(palette="BrBG")+
  theme_void()+theme(legend.position="none")
Scene=randdat %>%
    ggplot()+
    geom_ribbon(aes(x=MD,ymax=log2(newn),fill=interaction(Year,user),group=interaction(Year,user)),ymin=0,alpha=0.9)+
    geom_line(aes(x=MD,y=log2(newn),color=interaction(Year,user),group=interaction(Year,user)),alpha=0.5,lwd=0.3)+
    scale_fill_manual(values=c("#006d2c","#2ca25f","#66c2a4","#99d8c9","#ccece6",
                             "#d0d1e6","#a6bddb","#74a9cf","#2b8cbe","#045a8d"))+
    scale_color_manual(values=c("#ccece6","#99d8c9","#66c2a4","#2ca25f","#006d2c",
                                "#045a8d","#2b8cbe","#74a9cf","#a6bddb","#d0d1e6"))+
    theme_void()+theme(legend.position="none")

#ggsave("Birdsong.png",Scene,units="cm",width=40,height=20,bg="transparent")
#ggsave("Birdsong.svg",Scene,units="cm",width=40,height=20,bg="transparent")
#ggsave("Birdsong.pdf",Scene,units="cm",width=40,height=20,bg="transparent")
