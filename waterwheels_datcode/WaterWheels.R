library(tidyverse)
library(readr)
library(ggsci)
library(gganimate)
theme_set(theme_void())

darb=read_tsv("WaterWheels.tsv") %>% 
  rename(Gage_ft=`110364_00065`,
         Gage_qual=`110364_00065_cd`,
         Discharge_ft3ps=`110365_00060`,
         Discharge_qual=`110365_00060_cd`) %>% 
  mutate(Gage_m=0.3048*Gage_ft,
         Discharge_m3ps=0.0283168*Discharge_ft3ps) %>% 
  filter(grepl("A",Gage_qual) & !grepl("e",Gage_qual)  | grepl("A",Discharge_qual) & !grepl("e",Discharge_qual) ) %>% 
  separate(datetime,c("date","time"),sep=" ",remove=F) %>% 
  separate(date,c("year","month","day"),sep="-",remove=F) %>% 
  separate(time,c("hr","min","sec"),sep=":",remove=F) 
  

darb_daymean_disc=darb %>% 
  filter(!is.na(Discharge_m3ps)) %>% 
  group_by(date,year,month,day) %>% 
  summarise(discharge_m3ps_day_mean=mean(Discharge_m3ps),
            discharge_m3ps_day_std=sd(Discharge_m3ps),
            discharge_m3ps_day_n=sum(!is.na(Discharge_m3ps))) %>% 
  ungroup() %>% 
  mutate(discharge_m3ps_day_ste=discharge_m3ps_day_std/sqrt(discharge_m3ps_day_n),
         discharge_m3ps_day_95CI=qt(1 - (0.05 / 2), discharge_m3ps_day_n - 1) * discharge_m3ps_day_ste,
         discharge_m3ps_day_chng=discharge_m3ps_day_mean-lag(discharge_m3ps_day_mean))

darb_monmean_disc=darb %>% 
  filter(!is.na(Discharge_m3ps)) %>% 
  group_by(year,month) %>% 
  summarise(discharge_m3ps_mon_mean=mean(Discharge_m3ps),
            discharge_m3ps_mon_std=sd(Discharge_m3ps),
            discharge_m3ps_mon_n=sum(!is.na(Discharge_m3ps))) %>% 
  ungroup() %>% 
  mutate(discharge_m3ps_mon_ste=discharge_m3ps_mon_std/sqrt(discharge_m3ps_mon_n),
         discharge_m3ps_mon_95CI=qt(1 - (0.05 / 2), discharge_m3ps_mon_n - 1) * discharge_m3ps_mon_ste,
         discharge_m3ps_mon_chng=discharge_m3ps_mon_mean-lag(discharge_m3ps_mon_mean))

darb_yermean_disc=darb %>% 
  filter(!is.na(Discharge_m3ps)) %>% 
  group_by(year) %>% 
  summarise(discharge_m3ps_yer_mean=mean(Discharge_m3ps),
            discharge_m3ps_yer_std=sd(Discharge_m3ps),
            discharge_m3ps_yer_n=sum(!is.na(Discharge_m3ps))) %>% 
  ungroup() %>% 
  mutate(discharge_m3ps_yer_ste=discharge_m3ps_yer_std/sqrt(discharge_m3ps_yer_n),
         discharge_m3ps_yer_95CI=qt(1 - (0.05 / 2), discharge_m3ps_yer_n - 1) * discharge_m3ps_yer_ste,
         discharge_m3ps_yer_chng=discharge_m3ps_yer_mean-lag(discharge_m3ps_yer_mean))


darb_mrgmean1=full_join(darb_yermean_disc,darb_monmean_disc) %>% 
  mutate(rand=rnorm(discharge_m3ps_mon_n,mean=1,sd=0.3))
darb_mrgmean=full_join(darb_mrgmean1 %>% 
                         select(-rand),darb_daymean_disc) %>% 
  mutate(rand=rnorm(discharge_m3ps_day_n,mean=1,sd=0.3))

darb_mrgmean1 %>% 
  filter(year>1999) %>% 
  ggplot() +
  geom_point(aes(x=month,y=year,color=discharge_m3ps_mon_mean*rand))+
  scale_color_viridis_c()+theme_classic()

darb_mrgmean %>% 
  filter(year>1999) %>% 
  ggplot() +
  geom_spoke(aes(x=day,y=month,color=log10(discharge_m3ps_day_mean),angle=(discharge_m3ps_mon_chng/max(discharge_m3ps_mon_chng))*180*rand,group=year),radius=0.5,size=1)+
  scale_color_viridis_c()+theme_classic()
darb_mrgmean %>% 
  filter(year>1999) %>% 
  ggplot() +
  geom_spoke(aes(x=day,y=month,angle=discharge_m3ps_day_mean*180,color=discharge_m3ps_mon_chng,group=year),radius=0.5,size=1)+
  scale_color_viridis_c()+theme_classic()
darb_mrgmean %>% 
  filter(year>1999) %>% 
  ggplot() +
  geom_spoke(aes(x=day,y=month,angle=discharge_m3ps_mon_chng*180*rand,color=sqrt(discharge_m3ps_mon_ste),group=year),radius=0.5,size=1)+
  scale_color_viridis_c()+theme_classic()

darb_mrgmean %>% 
  filter(year>=1990,year<=2020) %>% 
  group_by(year,month) %>% 
  mutate(scaled=scale(discharge_m3ps_day_ste)) %>% 
  filter(!is.na(scaled)) %>% 
  summarise(min=min(scaled),max=max(scaled),men=mean(scaled),med=median(scaled)) %>% 
  arrange(desc(max))
P=darb_mrgmean %>% 
  filter(year>=1990,year<=2020) %>% 
  group_by(year,month) %>% 
  mutate(scaled=scale(discharge_m3ps_day_ste)) %>% 
  #arrange(desc(year)) %>% 
  ggplot() +
  #geom_spoke(aes(x=day,y=month,color=log10(discharge_m3ps_day_mean),angle=scale(discharge_m3ps_day_chng)*180,group=year),radius=0.5,size=1)+
  #geom_spoke(aes(x=day,y=month,color=log10(discharge_m3ps_day_mean),angle=scale(discharge_m3ps_day_chng)*-180,group=year),radius=0.5,size=1)+
  #geom_spoke(aes(x=day,y=month,color=as.numeric(as.character(year)),angle=scale(discharge_m3ps_day_chng)*180,group=year),radius=0.5,size=1)+
  #geom_spoke(aes(x=day,y=month,color=as.numeric(as.character(year)),angle=scale(discharge_m3ps_day_chng)*-180,group=year),radius=0.5,size=1)+
  geom_spoke(aes(x=day,y=month,color=as.numeric(as.character(year)),angle=scaled*30,group=year),radius=0.5,size=1,alpha=0.3,lineend="round")+
  #geom_spoke(aes(x=day,y=month,color=as.numeric(as.character(year)),angle=cumsum(scale(discharge_m3ps_day_chng))*-180,group=year),radius=0.5,size=1)+
  #scale_color_viridis_c(option="mako",direction=1)+theme(legend.position="none")
  #scale_color_distiller(palette="Oranges",direction=-1)+theme(legend.position="none")
  scale_color_distiller(palette="YlOrBr",direction=-1)+theme(legend.position="none")
  #scale_color_material("yellow")+theme(legend.position="none")
ggsave("test.png",P,units="cm",width=15,height=10,bg="white")


P +
  transition_states(year)+
  shadow_mark()
P +
  transition_time(as.numeric(as.character(year)))
P +
  transition_reveal(color)

animate(P + transition_states(year,transition_length = 1, state_length = 2)+shadow_mark(),
        nframes=160, fps=20, width=600, height=400, end_pause=0, renderer = gifski_renderer())
anim_save("test.gif")
