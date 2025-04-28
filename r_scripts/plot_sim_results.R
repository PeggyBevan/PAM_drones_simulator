#plotting
library(lme4)
library(sjPlot)
library(ggeffects)
library(ggplot2)
library(effects)
library(patchwork)

#Plots to make
#Species richness per cluster bar chart - colour by experiment with error bars

#plot 

#number of species detected, coloured by sampling technique
load('data/simulationresults_prop_set/allsummaries.RData')

#plot random, routed together to show impact of sampling intensity

# RndRoute_avg$simulation_type = 'Random'
# NN_avg$simulation_type = 'NearestNeighbour'
# Ad_avg$simulation_type = 'Adaptive'
# 
# 
# all_avg = rbind(RndRoute_avg, NN_avg, Ad_avg)
# all_avg_species = filter(all_avg, group=='Cluster')
# all_avg_lu = filter(all_avg, group == 'LandUse')

RndSummary$simulation_type = 'Random'
NNSummary$simulation_type = 'Routed'
AdplorSummary$simulation_type = 'Adaptive Explorative'
AdploiSummary$simulation_type = 'Adaptive Exploitative'
names(AdploiSummary) = names(RndSummary)
names(AdplorSummary) = names(RndSummary)
all_summary = rbind(RndSummary, NNSummary, AdplorSummary, AdploiSummary)


# Results Section 1: Impact of sampling completeness on diversity  --------

# Start by showing how sampling completeness impacts data collection - random & routed --------

Summary = rbind(RndSummary, NNSummary)
SummaryS = Summary[Summary$group == 'Cluster',]
#species richness should be converted to species detection accuracy

#match species richness in gt_lu with species richness in RndSummaryS
gt_SR = gt_lu %>%
  group_by(loc_name) %>%
  filter(!is.na(value)) %>%
  summarise(GTSpeciesRichness = n_distinct(value))
#this deletes levels of loc_name that only contain NAs

gt_lu_SR = gt_lu %>%
  group_by(Habitat) %>%
  filter(!is.na(value)) %>%
  summarise(GTSpeciesRichness = n_distinct(value))

SummaryS = left_join(SummaryS, gt_SR, by = 'loc_name')
#give values missing from gt_SR a value of 0
SummaryS$GTSpeciesRichness[is.na(SummaryS$GTSpeciesRichness)] = 0

#total species detected per experiment
richness = SummaryS %>%
  group_by(simulation_type, iteration, experiment) %>%
  summarise(across(4:22, sum))

#create binary species presence/absence
richness[,c(4:22)] = ifelse(richness[,4:22] > 0, 1, 0)
richness$species_richness = rowSums(richness[,c(4:22)])

richness$experiment = as.numeric(richness$experiment)
richness$iteration = as.factor(richness$iteration)
richness$simulation_type = factor(richness$simulation_type, levels = c('Random', 'Routed'))

#model this
mod_rnd = glmer(species_richness ~ experiment*simulation_type + (1|iteration), data = richness, family = 'poisson')
summary(mod_rnd)

agg_SR = richness %>%
  group_by(simulation_type, experiment) %>%
  summarise(avg_species_richness = mean(species_richness),
            species_sd = sd(species_richness)) %>%
  mutate(
    lower_ci = avg_species_richness - 1.96 * species_sd/sqrt(50),
    upper_ci = avg_species_richness + 1.96 * species_sd/sqrt(50)
  )

fig2a = ggplot(agg_SR, aes(x = experiment, y = avg_species_richness, colour = simulation_type)) +
  #geom_col(colour = 'black', fill = 'grey', alpha = 0.8) +
  geom_point(size = 2, position = position_dodge(width = 0.15)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), position = position_dodge(width = 0.15), width = 0.1) +
  geom_point(data = richness, aes(x = experiment, y = species_richness), alpha = 0.6, shape = 1, size = 1, position = position_jitterdodge(dodge.width = 0.15, jitter.height = 0.1)) + 
  theme_bw() +
  labs(y = 'BirdNET Species Detected', x = 'Sampling Intensity', colour = 'Sampling Strategy') +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  scale_y_continuous(expand = c(0,0.01), limits = c(8,20), breaks = seq(0,20,2)) +
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = 'top') +
  scale_color_manual(values=c("#E69F00", "#56B4E9")) +
  geom_hline(yintercept = 19, linetype = 'dashed', colour = 'black') +
  annotate("text", x = 0.2, y = 19.2, label = "Maximum", hjust = 0.4, vjust = 0, size = 4) +
  ggtitle('b')
fig2a

#Species richness per land-use type
Summarylu = rbind(RndSummary, NNSummary)
Summarylu = Summarylu[Summarylu$group == 'LandUse',]

Summarylu = left_join(Summarylu, gt_lu_SR, by = c('loc_name'='Habitat'))
#give values missing from gt_SR a value of 0
Summarylu$GTSpeciesRichness[is.na(Summarylu$GTSpeciesRichness)] = 0

#get mean and CIs for each land use type
Summarylu_avg = Summarylu %>%
  group_by(simulation_type, experiment, loc_name) %>%
  summarise(mean_species_detected = mean(species_number),
            sd_species_detected = sd(species_number),
            GTSpeciesRichness = mean(GTSpeciesRichness),
            len = n())

Summarylu_avg$CI95_low_species = Summarylu_avg$mean_species_detected - (qt(0.975,df=49)*Summarylu_avg$sd_species_detected/sqrt(50))
Summarylu_avg$CI95_up_species = Summarylu_avg$mean_species_detected + (qt(0.975,df=49)*Summarylu_avg$sd_species_detected/sqrt(50))

Summarylu_avg$simulation_type = factor(Summarylu_avg$simulation_type, levels = c('Random', 'Routed'))

Summarylu_avg$loc_name = factor(Summarylu_avg$loc_name, levels = c('Grassland', 'Secondary', 'Old Growth', 'Palm', 'Mangrove', 'Teak'))

sampint = c("#FFAABB","#77AADD",  "#99DDFF","#DD77AA","#009E73")  
#bar plot
fig2b = ggplot(data = Summarylu_avg) +
  #geom_col(data = gt_lu_SR, aes(x = reorder(Habitat, -GTSpeciesRichness), y = GTSpeciesRichness), fill = 'grey') +
  geom_col(data = Summarylu_avg, aes(x = loc_name, y = mean_species_detected, fill = experiment), position = 'dodge') +
  scale_fill_manual(values = sampint) +
  geom_errorbar(data = Summarylu_avg, aes(x = loc_name, ymin = CI95_low_species, ymax = CI95_up_species, group = experiment), position = 'dodge') +
  facet_wrap(~simulation_type) +
  theme_bw() +
  labs(x = 'Land-use Type', y = 'BirdNET Species Detected', fill = 'Sampling Intensity') +
  theme(axis.text.x = element_text(angle=60, hjust = 1),
        legend.title = element_text(family = "Arial", color = "black", size = 12),
        legend.text = element_text(family="Arial", color = "black", size = 12), 
        legend.position = 'top',
        strip.background = element_blank(),
        text = element_text(size = 12),
        axis.text = element_text(size = 12)) +
  scale_y_continuous(expand = c(0,0.01), limits = c(0,20)) + ggtitle('a')

fig2b

fig2b+ fig2a  +plot_layout(ncol = 1) &
  theme(legend.position = 'top')

# Results Section 2: Adaptive sampling ------------------------------------

#all data = allsummary
#total species detected per experiment
richness1 = all_summary %>%
  filter(group == 'Cluster') %>%
  group_by(simulation_type, iteration, experiment) %>%
  summarise(across(4:22, sum))

#create binary species presence/absence
richness1[,c(4:22)] = ifelse(richness1[,4:22] > 0, 1, 0)
richness1$species_richness = rowSums(richness1[,c(4:22)])

richness1$experiment = as.numeric(richness1$experiment)
richness1$iteration = as.factor(richness1$iteration)
richness1$simulation_type = factor(richness1$simulation_type, levels = c('Random', 'Routed', 'Adaptive', 'Adaptive Explorative', 'Adaptive Exploitative'))
#save this df for the model
richness_mod = richness1
#plot random vs adaptive - not including routed here because we already showed theres not much difference.
richness1 = filter(richness1, simulation_type %in% c('Random', 'Adaptive Explorative', 'Adaptive Exploitative'))

agg1_SR = richness1 %>%
  group_by(simulation_type, experiment) %>%
  summarise(avg_species_richness = mean(species_richness),
            species_sd = sd(species_richness)) %>%
  mutate(
    lower_ci = avg_species_richness - 1.96 * species_sd/sqrt(50),
    upper_ci = avg_species_richness + 1.96 * species_sd/sqrt(50)
  )

ggplot(agg1_SR, aes(x = experiment, y = avg_species_richness, colour = simulation_type)) +
  #geom_col(colour = 'black', fill = 'grey', alpha = 0.8) +
  geom_hline(yintercept = 19, linetype = 'dashed', colour = 'black')+
  geom_point(data = richness1, aes(x = experiment, y = species_richness), alpha = 0.3, shape = 1, size = 1, position = position_jitterdodge(dodge.width = 0.15, jitter.height = 0.1, jitter.width = 0.05)) + 
  geom_point(size = 2, position = position_dodge(width = 0.15)) +
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), position = position_dodge(width = 0.15), width = 0.1) +
  theme_bw() +
  labs(y = 'BirdNet Species Detected', x = 'Sampling Intensity', colour = 'Sampling Strategy') +
  scale_x_continuous(breaks = seq(0,1,0.2)) +
  scale_y_continuous(expand = c(0,0.01), limits = c(8,20), breaks = seq(0,20,2)) +
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  scale_color_manual(values=c("#E69F00", "#CC79A7", '#0010B7')) +
  annotate("text", x = 0.2, y = 19.2, label = "Maximum", hjust = 0.2, vjust = 0, size = 4)


#statistical test: impact of adaptive sampling on species richness
mod_ad = glmer(species_richness ~ experiment*simulation_type + (1|iteration), data = richness_mod, family = 'poisson')
summary(mod_ad)

