###Packages
library(rtweet)
library(tidyverse)
library(zoo)
library(ggirl)

load("birdsong.Rdata") ### in lieu of re-acquiring these data, which can change, it is stored in this object
#u1_fol=get_followers("#######",n=Inf,retryonratelimit = TRUE) ### user names removed here to protect identities
#u2_fol=get_followers("#######",n=Inf,retryonratelimit = TRUE) 
#save.image("birdsong.Rdata") ### a more complete dataset would have been used, but too large for Github

u1_foldat=lookup_users(u1_fol$user_id) ### get user data
u2_foldat=lookup_users(u2_fol$user_id) 

u1_rand=u1_foldat %>% 
  select(account_created_at) %>% ### this is the data I am interested in
  separate(account_created_at," ",into=c("YMD","Time")) %>% ### break up time and date
  select(-Time) %>% ### remove unnecessary information
  group_by(YMD) %>% ### perform following functions per YMD entry, ie specific day
  tally() %>% ### count number of followers whose accounts were created that day
  separate(YMD,"-",into=c("Year","Month","Day"),remove=F) %>% ### break down date into its components, keeping original
  unite("MD",Month:Day,remove=F) %>% ### make new variable that I will use as an axis so year can be a layer
  filter(Year>2016) %>% ### remove older accounts
  mutate(Year=factor(Year,levels=c(2021,2020,2019,2018,2017,2016))) %>% ### reorder the years, affecting plotting characteristics
  mutate(user="u1") %>% ### add a tag specific for this user to combine later
  mutate(nroll=rollmean(n, 7, na.pad=T), ### calculate 7-day rolling average of the tally
         rand=rnorm(n,mean=1,sd=0.1), ### create random noise
         nrand=n*rand, ### randomize the original data
         nrandroll=rollmean(nrand, 7, na.pad=T)) ### calculate 7-day rolling average of randomized data
u2_rand=u2_foldat %>% ### repeat for user2
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
  mutate(nroll=rollmean(n, 7, na.pad=T), 
         rand=rnorm(n,mean=1,sd=0.1),
         nrand=sqrt(n)*rand,
         nrandroll=rollmean(nrand, 7, na.pad=T)) ### here, square root was applied to the tally to force separation observable separation between the two

randdat=rbind(u1_rand,u2_rand) %>% ### slap the two datafgrames together
  filter(!is.na(nrandroll)) %>% ### remove NA values
  filter(MD!="01_01", ### remove certain problematic dates
         MD!="01_02",
         MD!="01_03") %>% 
  mutate(newn=predict(smooth.spline(nroll,nknots=3100))$y, ### smooth the rolling average data
         newnrand=predict(smooth.spline(nrandroll,nknots=3100))$y) ### smooth the rolling average data

randdat %>% ### call data
  ggplot()+ ### make a plot
  geom_ribbon(aes(x=MD,ymax=log2(nroll),fill=interaction(Year,user),group=interaction(Year,user)),ymin=0,alpha=0.9)+ ### make a ribbon, ie filled area 
  geom_line(aes(x=MD,y=log2(nroll),color=interaction(Year,user),group=interaction(Year,user)),alpha=0.3,stroke=1)+ ### draw line sto accentuate the ends of the ribbon 
  scale_fill_manual(values=c("#006d2c","#2ca25f","#66c2a4","#99d8c9","#ccece6", ### pick colors of ribbon fill
                             "#d0d1e6","#a6bddb","#74a9cf","#2b8cbe","#045a8d"))+
  scale_color_manual(values=c("#ccece6","#99d8c9","#66c2a4","#2ca25f","#006d2c", ### pick colors of lines, same pallete but opposite directions to increase contrast
                              "#045a8d","#2b8cbe","#74a9cf","#a6bddb","#d0d1e6"))+
  theme_classic()+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+ ### adjust overall aesthetics
  ylab("Number of followers")+xlab("follower creation date") ### label axes

Scene=randdat %>% ### repeat, but make more mysterious
    ggplot()+
    geom_ribbon(aes(x=MD,ymax=log2(newn),fill=interaction(Year,user),group=interaction(Year,user)),ymin=0,alpha=0.9)+
    geom_line(aes(x=MD,y=log2(newn),color=interaction(Year,user),group=interaction(Year,user)),alpha=0.5,lwd=0.3)+
    scale_fill_manual(values=c("#006d2c","#2ca25f","#66c2a4","#99d8c9","#ccece6",
                             "#d0d1e6","#a6bddb","#74a9cf","#2b8cbe","#045a8d"))+
    scale_color_manual(values=c("#ccece6","#99d8c9","#66c2a4","#2ca25f","#006d2c",
                                "#045a8d","#2b8cbe","#74a9cf","#a6bddb","#d0d1e6"))+
    theme_void()+theme(legend.position="none")

#ggsave("Birdsong.png",Scene,units="cm",width=40,height=20,bg="transparent") ### save the more artistic version
#ggsave("Birdsong.svg",Scene,units="cm",width=40,height=20,bg="transparent")
#ggsave("Birdsong.pdf",Scene,units="cm",width=40,height=20,bg="transparent")
