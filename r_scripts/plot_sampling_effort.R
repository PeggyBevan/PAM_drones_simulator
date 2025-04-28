#this script uses functions from the BirdNet CSV processing script - mod_sim_csv and format_dets

# packages ----------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(patchwork)

source('r_scripts/utils.R')

#functions needed:
#mod_sim_csv
#format_dets

# Load Site Info & Format -------------------------------------------------
site_info = read.csv('data/costa_rica_site_info.csv')
names(site_info)
site_info$loc_site = gsub('_Site_', '_Audio-', site_info$Site)
site_info$cluster = sub("\\_.*", "", site_info$loc_site)
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
                                      'LaReserva' = 'Indigenous_reserve'
                                      )),
         cluster = case_match(cluster,
                              "Rancho" ~ "Rancho_bajo",
                              "LaBalsa" ~ "La_Balsa",
                              "Sendero"~ "Golfito",
                              "Mangrove" ~ "Mangroves",
                              "SQ260" ~ "SQ260_1",
                              "Tarde" ~ "La Tarde",
                              "Nuevo" ~ "Rio_Nuevo",
                              "Gamba" ~ 'La Gamba',
                              "Palma" ~ "La Palma",
                              "LosPlanes" ~ "Los planes",
                              "Miramar" ~ "Mirenmar",
                              "Elsi" ~ "Elsi_Croc",
                              'LaReserva' ~ 'Indigenous_reserve',
                              .default = cluster
                              ))

unique(site_info$cluster)
#remove clusters not used in avian dataset
tokeep = c("SQ259","Rancho_bajo","Buho","SQ282","Golfito","Tigre",
           "Corcovado","La_Balsa","Ronnie" ,"Mangroves","Piro",
           "Danta","Sabalo","Lomas","Cira","Teak","SQ260_1","Naim",
           "Bonito","La Tarde","Rio_Nuevo","Nicuesa","Drake",
           "SQ283","Alto","La Gamba","SQ258","La Palma","Los planes",
           "SQ235", "Mirenmar","Marvin","Derek","Elsi_Croc"
           ,"Indigenous_reserve")

site_info = site_info[site_info$cluster %in% tokeep,]
site_info$site_name = sub('.*_', '', site_info$loc_site)
site_info$site_name = sub('-', " ", site_info$site_name)

# Site statistics ---------------------------------------------------------

#number of clusters
length(unique(site_info$cluster)) #35
#average number of sites per cluster
mean(table(site_info$cluster)) #8.71
#min and max
range(table(site_info$cluster)) #4, 28
#number of sites
unique(site_info$loc_site) %>% length()


# Plot distribution of land-use types in each cluster ---------------------

#order sites by proportion of old growth
luorder = site_info %>%
  group_by(cluster, Habitat) %>%
  summarise('n' = n()) %>%
  pivot_wider(names_from = 'Habitat', values_from = 'n') %>%
  arrange(desc(`Old Growth`), desc(Secondary))

site_info$cluster = factor(site_info$cluster, levels = luorder$cluster)
ggplot(site_info, aes(x = cluster, fill = Habitat,)) +
  geom_bar(position='fill') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.05,hjust = 1))

#get long/lat for each cluster
site_info = site_info %>%
  group_by(cluster) %>%
  mutate(Long_cluster = mean(Long),
         Lat_cluster = mean(Lat))

ggplot(site_info, aes(x = Long, y = Lat, colour = cluster)) +
  geom_label(aes(x = Long_cluster, y = Lat_cluster, label = site_info$cluster), check_overlap = T) +
  geom_point() +
  theme_classic()


# Measuring sampling evenness over time -----------------------------------

#plot temporal coverage of random vs routed
#for every time period, want to know how many times a site was sampled.

#y axis will be sites, x axis will be time of day
#geom_matrix - colour will be frequency
#dataframe - for each site, group by date

gt_lu = read.csv('data/groundtruth_formatted.csv')
gt_lu$datetime2 = as.POSIXct(gt_lu$datetime, format = '%Y-%m-%d %H:%M:%OS', tz = "America/Costa_Rica")
gt_lu$time = format(gt_lu$datetime2, format = '%H:%M:%OS')
gt_lu$date = format(gt_lu$datetime2, format = '%Y-%m-%d')

#remove duplicates 
gt_lu = gt_lu[!(duplicated(gt_lu[c(1,2,3,8)])),]

#for each cluster, work out number of days sampled
gt_lu = gt_lu %>%
  group_by(loc_name) %>%
  mutate(effort_days = n_distinct(date))
#average temporal freq: how many times was a site surveyed at a particular hour?
samp_freq = gt_lu %>%
  group_by(loc_name, site_name, time, effort_days) %>%
  summarise(freq = n())

samp_freq$loc_site = paste0(samp_freq$loc_name,'_', samp_freq$site_name)

#remove 9 & 16:00 as these are not part of original regime
samp_freq = filter(samp_freq, !time %in% c('09:00:00', '16:00:00'))

gt_sampfreq = ggplot(samp_freq, aes(x = time, y = loc_name)) +
  geom_tile(aes(fill = freq)) +
  theme(axis.text.x = element_text(angle = 45))
gt_sampfreq
#for most clusters, there is even sampling across time. 
# sampling evenness loop --------------------------------------------------
options(dplyr.summarise.inform = FALSE)

#For each iteration, measure sampling evenness per prop_sampler. 
Rndpath = 'data/Iterations_proportional_set/Random/'
Rnd_even = avg_evenness(Rndpath, gt_lu)
Rnd_even$method = 'Random'

NNpath = 'data/Iterations_proportional_set/NNClusRoute/'
NNeven = avg_evenness(NNpath, gt_lu)
NNeven$method = 'Routed'

Adplorpath = 'data/Iterations_proportional_set/Adaptivedets_factor0_3/'
Adploreven = avg_evenness(Adplorpath, gt_lu)
Adploreven$method = 'Adaptive_Explorative'

Adploipath = 'data/Iterations_proportional_set/Adaptivedets_factor100/'
Adploieven = avg_evenness(Adploipath, gt_lu)
Adploieven$method = 'Adaptive_Exploitative'

#now we have sampling evenness for each cluster, for each proportion of samplers, for each method. 
#Start by looking at just random
#get average evenness for each cluster per prop. 

#group three methods together
all_even = rbind(Rnd_even, NNeven, Adploieven, Adploreven)

even_avg = all_even %>%
  group_by(method, prop_samplers) %>%
  summarise(mean_evenness = mean(evenness, na.rm = TRUE),
            sd_evenness = sd(evenness, na.rm = TRUE),
            mean_visits = mean(total_visits),
            sd_visits = sd(total_visits),
            len = n())

#even_avg$method[even_avg$method == 'NearestNeighbour'] = 'Routed'
even_avg$method = factor(even_avg$method, levels = c('Random', 'Routed', 'NearestNeighbour', 'Adaptive_Explorative', 'Adaptive_Exploitative'))
#now make a plot
p1 = ggplot(even_avg, aes(x = prop_samplers, y = mean_evenness, ymin = mean_evenness - sd_evenness, ymax = mean_evenness + sd_evenness, colour = method)) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_errorbar(position = position_dodge(width = 0.1), width = 0.1) +
  #geom_line() +
  labs(x = 'Sampling Intensity', y = 'Sampling Evenness', colour = 'Sampling Method') +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  #facet_wrap(~method) +
  theme_bw() +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#CC79A7", '#0010B7')) +
  theme(text = element_text(size = 14)) +
    ggtitle('a')
p1

# Travel distance ---------------------------------------------------------

#loop for collection travel and charger time

get_travel_time = function(file_path) {
  files = list.files(file_path)
  sum_dist_full = NULL
  for (file in 1:length(files)) {
    simdf = read.csv(paste0(file_path, files[file]))
    print(paste0('File Loaded ', files[file]))
    #find total distance traveled by each sampler
    sum_dist = simdf %>%
      group_by(prop_samplers, loc_name,sampler_num) %>%
      summarise(dist = sum(total_dist),
                charger_visits = sum(charger_visits))
    sum_dist$iteration = file
    sum_dist_full = rbind(sum_dist_full, sum_dist)
  }
  sum_dist_full$method = basename(file_path)
  return(sum_dist_full)
}


#get number of sites per cluster from groundtruth 
distbysite = site_info %>%
  group_by(cluster) %>%
  summarise(n_sites = n_distinct(site_name))

Rndpath = 'data/Iterations_proportional_set/Random/'
rnd_sum_dist = get_travel_time(Rndpath)
NNpath = "data/Iterations_proportional_set/NNClusRoute/"
nn_sum_dist = get_travel_time(NNpath)
Adplorpath = "data/Iterations_proportional_set/Adaptivedets_factor0_3/"
Adplor_sum_dist = get_travel_time(Adplorpath)
Adploipath = "data/Iterations_proportional_set/Adaptivedets_factor100/"
Adploi_sum_dist = get_travel_time(Adploipath)

all_sumdist = rbind(rnd_sum_dist, nn_sum_dist, Adplor_sum_dist, Adploi_sum_dist)

all_sumdist$method = c(rep('Random', nrow(rnd_sum_dist)), rep('Routed', nrow(nn_sum_dist)), rep('Adaptive_Explorative', nrow(Adplor_sum_dist)), rep('Adaptive_Exploitative', nrow(Adploi_sum_dist)))

all_sumdist2 = all_sumdist %>%
  #get total distance for all samplers at a cluster
  group_by(method, iteration, prop_samplers, loc_name) %>%
  summarise(total_dist = sum(dist),
            total_charger_visits = sum(charger_visits)) %>%
  mutate(loc_name = sub('.pickle', '', loc_name)) %>%
  left_join(distbysite, by = c('loc_name'='cluster')) %>%
  mutate(dist_persite = total_dist/n_sites,
         charges_per_site = total_charger_visits/n_sites)

#for the entire dataset, I want the average distance travelled at each cluster.
#to account for variation in site number and possible distances, normalise total distances by site number. 
distavg = all_sumdist2 %>%
  group_by(method, prop_samplers) %>%
  summarise(mean_distbysite = mean(dist_persite),
            sd_distbysite = sd(dist_persite),
            mean_chargebysite = mean(charges_per_site),
            sd_chargebysite = sd(charges_per_site),
            n = n()
            ) %>%
  mutate(dist_upper_ci =  mean_distbysite + (1.96 * sd_distbysite/sqrt(n)),
         dist_lower_ci =  mean_distbysite - (1.96 * sd_distbysite/sqrt(n)),
         charge_upper_ci =  mean_chargebysite + (1.96 * sd_chargebysite/sqrt(n)),
         charge_lower_ci =  mean_chargebysite - (1.96 * sd_chargebysite/sqrt(n)))

distavg$method = factor(distavg$method, levels = c('Random', 'Routed', 'Adaptive_Explorative', 'Adaptive_Exploitative'))
#distavg$method[distavg$method == 'NNClusRoute'] = 'Routed'

p2 = ggplot(distavg, aes(x = prop_samplers, y = mean_distbysite, colour = method, ymin = dist_lower_ci, ymax = dist_upper_ci)) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_errorbar(position = position_dodge(width = 0.1), width = 0.1) +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  labs(y = 'Average travel time (m) (normalised)', x = 'Sampling Intensity', colour = 'Sampling Method') +
  theme_bw() +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#CC79A7", '#0010B7')) +
  theme(text = element_text(size = 14)) +
  ggtitle('b')
p2

p3 = ggplot(distavg, aes(x = prop_samplers, y = mean_chargebysite, colour = method, ymin = charge_lower_ci, ymax = charge_upper_ci)) +
  geom_point(position = position_dodge(width = 0.1)) +
  geom_errorbar(position = position_dodge(width = 0.1), width = 0.1) +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  labs(y = 'Average Charger visits (normalised)', x = 'Sampling Intensity', colour = 'Sampling Method') +
  theme_bw() +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#CC79A7", '#0010B7')) +
  theme(text = element_text(size = 14)) +
  ggtitle('c')
p3
p1+p2 + p3+ plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom')


