#pkgs
library(tidyverse)
library(ggsci)
library(ggspatial)
require(maps)
require(rworldmap)
set.seed(42)

#import data
#e_obs = read_csv("Obsv_230804_gbif-observations-dwca/FamilyEricaceae_awkifcolumn37equalsEricaceae_wcolids_observations.csv") %>% 
#  separate(dateIdentified, c("date", "time"), sep = " ", remove = F) %>% 
#  separate(date, c("YEAR", "MONTH", "DAY"), sep = "-", remove = F) %>% 
#  mutate(LAT = round(decimalLatitude, digits = 2), 
#         LON = round(decimalLongitude, digits = 2))
##exported for ease of github
write_tsv(e_obs %>% 
            filter(countryCode  ==  "US") %>% 
            select_at(vars(-contains(c("recorded", "identified", "Remark")))),
          "EricaSmear.tsv")
e_obs = read_tsv("EricaSmear.tsv")

#e_tax = read_csv("Taxa_230804_inaturalist-taxonomy.dwca/FamilyEricaceae_awkifcolumn9equalsEricaceae_wcolids_taxa.csv")
#obs and tax data originating from: https://www.inaturalist.org/pages/developers
#write_tsv(e_tax, "EricaSmear_tax.tsv")
e_tax = read_tsv("EricaSmear_tax.tsv")

#set up manual data for massaging imported datasets
mycohets = c("Allotropa", 
           "Cheilotheca", 
           "Hemitomes", 
           #"Hypopitys", #all Hypopitys have been reclassified to Monotropa
           "Monotropa", 
           "Monotropastrum", 
           "Monotropsis", 
           "Pityopus", 
           "Pleuricospora", 
           "Pterospora", 
           "Sarcodes", 
           "Orthilia", #mixotrophic, ~50% from mycohet lifestyle
           "Pyrola")  #facultative/variable within genus?
nonmycohets = c("Chimaphila", 
              "Moneses")
#source Wikipedia
#started with: https://en.wikipedia.org/wiki/Monotropoideae
#then checked each genera's page  
# Note, there are other groups, e.g. orchids and Hyobanche and eevn more, that have at least a partially mycoheterotrophic lifestyle.

colors_for_genera = c("#fff7bc", 
                    "#fee391", 
                    "#fec44f", 
                    "#fe9929", 
                    "#ec7014", 
                    "#cc4c02", 
                    "#993404", 
                    "#662506", 
                    "#cb181d", 
                    "#a50f15", 
                    "#addd8e", 
                    "#78c679", 
                    "#67a9cf", 
                    "#3690c0", 
                    "gray50")

#massage and make summaries of imported data
mycohet_e_tax = e_tax %>% 
  mutate(GENUS = ifelse(genus %in% mycohets | genus %in% nonmycohets, genus, "other"), 
         GENUS = factor(GENUS, levels = c(mycohets, nonmycohets, "other")), 
         SUBFAM = ifelse(genus %in% mycohets | genus %in% nonmycohets, "Monotropoideae/Pyroloideae", "other subfam"), 
         SUBFAM = factor(SUBFAM, levels = c("Monotropoideae/Pyroloideae", "other subfam")), 
         METAB = ifelse(genus == "Orthilia" | genus  == "Pyrola", "mixed", ifelse(genus %in% mycohets, "mycohet", "photosynth")), 
         METAB = factor(METAB, levels = c("mycohet", "mixed", "photosynth")))

mycohet_e_obs = full_join(e_obs %>% 
                          rename(obsid = id) %>% 
                          select(-modified, -references) %>% 
                          mutate(taxonID = as.character(taxonID)), 
                        mycohet_e_tax %>% 
                          rename(taxid = id) %>% 
                          select(-modified, -references) %>% 
                          mutate(taxonID = sub(".*\\/", "", taxonID))) %>% 
  arrange(GENUS) %>% 
  mutate(scientificName = factor(scientificName, levels = unique(scientificName)))

mycohet_e_obs_smry = mycohet_e_obs %>% 
  filter(countryCode == "US") %>% 
  filter(!grepl("other", SUBFAM), LAT <= 50, LON >= -140) %>% 
  group_by(GENUS, LON, LAT) %>% 
  summarise(count = n())

#Visualize
mycohet_e_obs_smry %>% 
  arrange(count) %>% 
  ggplot()+
  geom_point(aes(x = LON, y = LAT, alpha = count, size = count + 1, color = GENUS), position = position_jitter(width = 10, height = 10), stroke = 0)+
  coord_quickmap()+
  #scale_fill_manual(values = colors_for_genera, drop = F)+
  scale_color_manual(values = colors_for_genera, drop = F)+
  scale_alpha_continuous(trans = "log2") + 
  scale_size_area(max_size = 10, trans = "log10") + 
  theme_void() + theme(legend.position = "none")

#ggsave("EricaSmear.png", units = "cm", width = 10, height = 8)
