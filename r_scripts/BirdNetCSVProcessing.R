##Data processing
# Load Packages -----------------------------------------------------------

library(vegan)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)

#load functions for this project
source('r_scripts/utils.R')

# Read in files -----------------------------------------------------------
#load site info
site_info = read.csv('data/costa_rica_site_info.csv')
names(site_info)
site_info$loc_site = gsub('_Site_', '_Audio-', site_info$Site)

#reformat cluster names to match the audio data files
site_info <- site_info %>%
  mutate(loc_site = str_replace_all(loc_site,
                                    c("Rancho" = "Rancho_bajo",
                                      "LaBalsa" = "La_Balsa",
                                      "Sendero"= "Golfito",
                                      "Mangrove" = "Mangroves",
                                      "SQ260" = "SQ260_1",
                                      "Tarde" = "La Tarde",
                                      "Nuevo" = "Rio_Nuevo",
                                      "Gamba" = 'La Gamba',
                                      "Palma" = "La Palma",
                                      "LosPlanes" = "Los planes",
                                      "Miramar" = "Mirenmar",
                                      "Elsi" = "Elsi_Croc",
                                      'LaReserva' = 'Indigenous_reserve',
                                      'Corcovado 1' = 'Corcovado')))

#Read in files from the simulation and ground truth, format and save as large database
###Read in groundtruth
gt = read.csv('data/groundtruth_detections_80thresh_PB.csv')
gt$site_name = gsub(" ", "-", gt$site_name)
gt$loc_site = paste0(gt$loc_name, '_', gt$site_name)
gt = filter(gt, loc_name != 'Luna')

#format some site names to match site info csv
gt = gt %>%
  mutate(site_name = case_match(site_name,
                                "Audio-1-92" ~ 'Audio-1',
                                "Audio-4-375" ~ 'Audio-4',
                                "Audio-5-58" ~ 'Audio-5',
                                "Audio-6-673" ~ 'Audio-6',
                                "Danta-14-grassland" ~ 'Audio-14',
                                .default = site_name))
gt$loc_site = paste0(gt$loc_name, '_', gt$site_name)
names(gt)[5] = 'sampled_specs'
#remove sites that are not in site info
gtsites = unique(gt$loc_site) #279
sites = unique(site_info$loc_site) #341
gtonly = gtsites[!(gtsites %in% sites)] #which sites are in gt but not site_info
#no gps coordinates for elsi_croc_audioCB - remove this
gt = gt[!(gt$loc_site %in% gtonly), ] 

#the species detections are in a single column, in a list format. Some rows are just '[]', where no species were detected. Here, I reformat so that there is one row per species detection.
#format species detections
gt$sampled_specs = gsub("\\[|]", '', gt$sampled_specs)
#The max number of detections in one sampling period is 20 - so here we make 20 columns with one detection per column.
gt = separate(gt, sampled_specs, into = paste('detection_', 1:20, sep = ''),
              sep = ',', extra = 'merge', fill = 'right', remove = FALSE)
#convert to longform - select sites with detections 
gt = pivot_longer(gt, cols = c(6:25))
#don't want to delete all the sites.
#remove NAs - a blank cell means there was a recording, but no detection. NAs are extra cells created through the separate function.
gt = gt[!is.na(gt$value),]
#no detection, but sampled, = NA
gt$value[gt$value==""] = NA
#remove 'sampled_specs' and 'name' cols
gt = gt[,c(1,2,3,4,6,8)]
#tidy up species names
#remove extra quotation marks
gt$value = gsub(" [']", "", gt$value)
gt$value = gsub("[']", "", gt$value)

#summarise the total number of species detected at each cluster
speciesgt = gt %>%
  group_by(loc_name) %>%
  summarise(species_number = n_distinct(value, na.rm = TRUE))

#total number of bird detections:
nrow(gt[!is.na(gt$value),]) #1819
#number of sites where birds are detected
unique(gt$loc_site[!is.na(gt$value)]) #126
unique(gt$loc_name[!is.na(gt$value)]) 


#add land use
gt_lu = left_join(gt, site_info[,c('loc_site', 'Habitat' )], by = 'loc_site')
write.csv(gt_lu, 'data/groundtruth_formatted.csv', row.names = FALSE)
valid_specs = unique(gt$value[!is.na(gt$value)])

# Read in experiments and modify CSVs to be more manageable ---------
#1. Random
Rndpath = 'data/Iterations_proportional_set/Random/'
RndSummary = create_summary_df(file_path = Rndpath, valid_specs)
#this produces a csv of 50 iterations, each goes through 1,2,4,8 samplers, and each row is a site with total species detections, 

#for each site, i want average species richness, plus max and min
RndRoute_avg = RndSummary %>%
  group_by(group, loc_name, experiment) %>%
  summarise(avg_SpeciesNumber = mean(species_number), max_SpeciesNumber = max(species_number),
            min_SpeciesNumber = min(species_number), 
            sd_SpeciesNumner = sd(species_number),
            len = n())

names(RndRoute_avg)
RndRoute_avg$CI95_low = RndRoute_avg$avg_SpeciesNumber - (qt(0.975,df=49)*RndRoute_avg$sd_SpeciesNumner/sqrt(50))
RndRoute_avg$CI95_up = RndRoute_avg$avg_SpeciesNumber + (qt(0.975,df=49)*RndRoute_avg$sd_SpeciesNumner/sqrt(50))

#2. Nearest Neighbour routed
NNpath = 'data/Iterations_proportional_set/NNClusRoute/'
time1= Sys.time()
NNSummary = create_summary_df(file_path = NNpath, valid_specs)
time = Sys.time() - time1

#summarise data per site
#for each site, i want average species richness, plus max and min
NN_avg = NNSummary %>%
  group_by(group, loc_name, experiment) %>%
  summarise(avg_SpeciesNumber = mean(species_number),
            max_SpeciesNumber = max(species_number),
            min_SpeciesNumber = min(species_number),
            sd_SpeciesNumner = sd(species_number),
            len = n())

NN_avg$CI95_low = NN_avg$avg_SpeciesNumber - (qt(0.975,df=49)*NN_avg$sd_SpeciesNumner/sqrt(50))
NN_avg$CI95_up = NN_avg$avg_SpeciesNumber + (qt(0.975,df=49)*NN_avg$sd_SpeciesNumner/sqrt(50))

#3. Adaptive sampling - explorative
Adplorpath = 'data/Iterations_proportional_set/Adaptivedets_factor0_3/'
AdplorSummary = create_summary_df(file_path = Adplorpath, valid_specs)
#this produces a csv of 50 iterations, each goes through 1,2,4,8 samplers, and each row is a site with total species detections, 

#summarise data per site
#for each site, i want average species richness, plus max and min
Adplor_avg = AdplorSummary %>%
  group_by(group, loc_name, experiment) %>%
  summarise(avg_SpeciesNumber = mean(species_number),
            max_SpeciesNumber = max(species_number),
            min_SpeciesNumber = min(species_number),
            sd_SpeciesNumner = sd(species_number),
            len = n())

Adplor_avg$CI95_low = Adplor_avg$avg_SpeciesNumber - (qt(0.975,df=49)*Adplor_avg$sd_SpeciesNumner/sqrt(50))
Adplor_avg$CI95_up = Adplor_avg$avg_SpeciesNumber + (qt(0.975,df=49)*Adplor_avg$sd_SpeciesNumner/sqrt(50))

# Adaptive - eploitative --------------------------------------------------
Adploipath = 'data/Iterations_proportional_set/Adaptivedets_factor100/'
AdploiSummary = create_summary_df(file_path = Adploipath, valid_specs)

#summarise data per site
#for each site, i want average species richness, plus max and min
Adploi_avg = AdploiSummary %>%
  group_by(group, loc_name, experiment) %>%
  summarise(avg_SpeciesNumber = mean(species_number),
            max_SpeciesNumber = max(species_number),
            min_SpeciesNumber = min(species_number),
            sd_SpeciesNumner = sd(species_number),
            len = n())

Adploi_avg$CI95_low = Adploi_avg$avg_SpeciesNumber - (qt(0.975,df=49)*Adploi_avg$sd_SpeciesNumner/sqrt(50))
Adploi_avg$CI95_up = Adploi_avg$avg_SpeciesNumber + (qt(0.975,df=49)*Adploi_avg$sd_SpeciesNumner/sqrt(50))


#save output
outputdir = 'data/simulationresults_prop_set/'
if (!dir.exists(outputdir)) dir.create(outputdir)
write.csv(RndSummary, paste0(outputdir, 'RandomSummary.csv'), row.names = FALSE)
write.csv(RndRoute_avg,  paste0(outputdir, 'Random_summary_averaged.csv'), row.names = FALSE)
write.csv(NNSummary,  paste0(outputdir, 'NNClusSummary.csv'), row.names = FALSE)
write.csv(NN_avg,  paste0(outputdir, 'NNClusSummary_averaged.csv'), row.names = FALSE)
write.csv(AdSummary,  paste0(outputdir, 'AdaptSummary.csv'), row.names = FALSE)
write.csv(Ad_avg,  paste0(outputdir, 'AdaptSummary_averaged.csv'), row.names = FALSE)
write.csv(AdplorSummary, paste0(outputdir, 'AdaptexploreSummary.csv'), row.names = FALSE)
write.csv(Adplor_avg, paste0(outputdir, 'AdaptexploreSummary_averaged.csv'), row.names = FALSE)
write.csv(AdploiSummary, paste0(outputdir, 'AdaptexploitSummary.csv'), row.names = FALSE)
write.csv(Adploi_avg, paste0(outputdir, 'AdaptexploitSummary_averaged.csv'), row.names = FALSE)

save(gt_lu, RndSummary, RndRoute_avg, NNSummary, NN_avg, AdploiSummary, Adploi_avg, AdplorSummary, Adplor_avg, file =  paste0(outputdir, 'allsummaries.RData'))

#Quick plot of summaries

RndRoute_avg$facet_label = paste0("N Samplers: ",RndRoute_avg$experiment)
ggplot() +
  geom_col(data = speciesgt, aes(x = reorder(loc_name, -species_number), y = species_number, fill = 'Missed')) +
  geom_col(data = RndRoute_avg, aes(x = loc_name, y = avg_SpeciesNumber, fill = "Sampled")) +
  geom_errorbar(data = RndRoute_avg, aes(x = loc_name, ymin = CI95_low, ymax = CI95_up, width = 0.4)) +
  facet_wrap(~facet_label, scales = 'free_x') +
  theme_bw()+
  theme(axis.text.x = element_text(angle=60, hjust = 1),
        legend.title = element_blank(),
        legend.text = element_text(family="Arial", color = "black", size = 16,face="bold"),
        strip.background = element_blank()) +
  scale_fill_manual(values = c('Missed' = '#D55E00', "Sampled" = '#009E73')) +
  ylab('Species Richness') +
  xlab('Location')

# routed sampling -------------------------------

NN_avg$facet_label = paste0("N Samplers: ",NN_avg$experiment)
ggplot() +
  geom_col(data = speciesgt, aes(x = reorder(loc_name, -species_number), y = species_number, fill = 'Missed')) +
  geom_col(data = NN_avg, aes(x = loc_name, y = avg_SpeciesNumber, fill = "Sampled")) +
  geom_errorbar(data = NN_avg, aes(x = loc_name, ymin = CI95_low, ymax = CI95_up, width = 0.4)) +
  facet_wrap(~facet_label, scales = 'free_x') +
  theme_bw()+
  theme(axis.text.x = element_text(angle=60, hjust = 1),
        legend.title = element_blank(),
        legend.text = element_text(family="Arial", color = "black", size = 16,face="bold"),
        strip.background = element_blank()) +
  scale_fill_manual(values = c('Missed' = '#D55E00', "Sampled" = '#009E73')) +
  ylab('Species Richness') +
  xlab('Location')


# Adaptive sampling - avian -----------------------------------------------
Ad_avg$facet_label = paste0("N Samplers: ",Ad_avg$experiment)
ggplot() +
  geom_col(data = speciesgt, aes(x = reorder(loc_name, -species_number), y = species_number, fill = 'Missed')) +
  geom_col(data = Ad_avg, aes(x = loc_name, y = avg_SpeciesNumber, fill = "Sampled")) +
  geom_errorbar(data = Ad_avg, aes(x = loc_name, ymin = CI95_low, ymax = CI95_up, width = 0.4)) +
  facet_wrap(~facet_label, scales = 'free_x') +
  theme_bw()+
  theme(axis.text.x = element_text(angle=60, hjust = 1),
        legend.title = element_blank(),
        legend.text = element_text(family="Arial", color = "black", size = 16,face="bold"),
        strip.background = element_blank()) +
  scale_fill_manual(values = c('Missed' = '#D55E00', "Sampled" = '#009E73')) +
  ylab('Species Richness') +
  xlab('Location')