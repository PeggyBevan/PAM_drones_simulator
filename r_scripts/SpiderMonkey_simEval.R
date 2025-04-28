#load spider monkey data and perform analyses, comparing groundtruth data with simulation data.
#first, do all analyses of groundtruth.
#next, do random vs routed analyses, compared to groundtruth.

library(dplyr)
library(lme4)
library(sp)
library(spdep)
library(car)
library(DHARMa)
library(lmtest)
library(data.table)
library(ggeffects)
library(sjPlot)
library(ROCR)
library(tidyr)
library(jtools) #summ

options(scipen = 999) #no e numbers

#spidermonkeys



# load data ---------------------------------------------------------------
#groundtruth
gtSM = read.csv('data/SpiderMonkeys_Processeddata.csv')
#remove unneccessary columns
#remove row number column, month, and all environmental variables except habitat and forest cover 100
#secondary 200km, house 1000
#gtSM = gtSM[,c(2:12,15,17,19,20,26,27,71)]

#site info
site_info = read.csv('data/costa_rica_site_info.csv')

#simulated data
#random
file_path = 'data/Iterations_proportional_set/RandomSM/'
file_list = list.files(file_path)
RndSM = NULL
for (i in 1:length(file_list)) {
  file = read.csv(paste0(file_path, file_list[i]))
  RndSM = rbind(RndSM, file)  
}
RndSM$simulation_type = 'Random'
#random
file_path = 'data/Iterations_proportional_set/NNClusRouteSM/'
file_list = list.files(file_path)
NNSM = NULL
for (i in 1:length(file_list)) {
  file = read.csv(paste0(file_path, file_list[i]))
  NNSM = rbind(NNSM, file)  
}
NNSM$simulation_type = 'NearestNeighbour'

simSM = rbind(RndSM, NNSM)

#remove duplicates - where simulation type, iteration, , prop_samplers, site_name, trajectory dates are the same
simSM = simSM[!duplicated(simSM[,c(17,6,11,15,16)]),]

#there is no issue with replication or removing fake recording periods because of the difference in data structure
#compared to the avian dataset.

#the cluster names also match those in site info.
length(unique(gtSM$Cluster))
unique(gt$Cluster)
unique(site_info$Cluster)

# Groundtruth analysis ----------------------------------------------------

#stats to get. 

## Min forest cover for spider monkeys (raw) -------------------------------

min_gt = min(gtSM$ForestCover100[gtSM$Presence==1])
min_gt #0.79
#which sites?
unique(gtSM$Site[gtSM$Presence==1 & gtSM$ForestCover100==min_gt])
#just Neuevo Site 3.
#Nuevo is an area of 100% forest cover where no spider monkeys were detected. but they were detected in one area of 79% forest cover within this (strange!) - this shows why we should model the threshold rather than use the absolute value.


## Threshold for spider monkey occupancy (GLM) -----------------------------------
#Thresholded level of forest cover
#coefficient from logistic regression model.

#Lawson et al. 2023 - they tried an occupancy model but found that the occupancy estimates did not differ
# from naive occupancy models, and the framework could not cope with complicated variables.
#So, they used a logistic regression glm. 

#in the paper, they tested different spatial scales of forest cover etc.
#they found 200m forest cover radius had the strongest response, secondary roads 200m
#primary roads, 1000, buildings 1000m

## Adding autocovariate from Lawson et al. ----------------------------------------------------
# define cell coordinates 
coords <- as.matrix(cbind(gtSM$Lat, gtSM$Long))
# construct autcovariate - increase neigbourhood dist (nbs) by increments of 0.1 till no cells have zero neighbours
ac <- autocov_dist(as.numeric(gtSM$Presence), coords, nbs = 2.1, longlat = TRUE)
# combine with cell coordinates
AC <- data.frame(ac = ac, x = gtSM$Lat, y = gtSM$Long)
gtSM<-cbind(gtSM, ac)

# GLM logistic regression model -------------------------------------------
#the best fitting model was forest cover / secondary road/ building. (with month as random effect?)
#month was does not change model fit, so was not included in model. 

#Secondary 200 = Secondary roads within 200m
#House 1000 = proportion of buildings within 1000m
gtglm1 <- glm(Presence ~ ForestCover100 + Secondary200km + House1000km +ac, family = binomial,
              data = gtSM)

gtglm2 <- glm(Presence ~ ForestCover100 + ac , family = binomial,
              data = gtSM)

summary(gtglm1)
summary(gtglm2)
AIC(gtglm1, gtglm2)

#We will keep all variables to match the analysis of Lawson et al.,

ggpredict(gtglm2, c('ForestCover100')) %>% 
  plot() + 
  labs(x = 'Forest Cover',
       y = 'Occupancy',
       title = '') +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

prediction1 = predict(gtglm2, type = 'response')
pred = prediction(prediction1, gtSM$Presence)
perf = performance(pred, 'sens', 'spec')

plot(perf, avg = 'threshold', colorize = TRUE, lwd = 3, print.cutoffs.at=seq(0,1,by=0.1))

#find the optimal cutoff point - the maximum sum of sensitivity and specificity
cutoffgt = perf@alpha.values[[1]][which.max(perf@x.values[[1]]+perf@y.values[[1]])]

ndata = data.frame(ForestCover100 = seq(0,1,0.001), ac = mean(gtSM$ac))
testpred = predict(gtglm2, ndata, type = 'response')
#find the value of forest cover where probability = cutoffgt
threshold = testpred[testpred > cutoffgt]
for_cut = ndata$ForestCover100[as.numeric(names(threshold)[1])]

#so you need at least 90% forest cover for spider monkeys.
#visually check that these thresholds make sense:
F3A = ggpredict(gtglm2, c('ForestCover100')) %>% 
  plot() + 
  labs(x = 'Forest Cover',
       y = 'Occupancy Probability',
       title = '') +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  #geom_hline(yintercept = cutoffgt) +
  geom_vline(xintercept = for_cut, linetype = 'dashed', colour = 'red', alpha = 0.8) +
  geom_point(data = gtSM, aes(x = ForestCover100, y = Presence)) +
  annotate('text', label = for_cut, x = 0.9, y = 0.75, vjust = -0.5, angle = 90) + 
  ggtitle('a')


#plot coefficient
forcoefgt = coef(gtglm1)[2]
forCIlowgt = confint(gtglm1)[2]
forCIupgt = confint(gtglm1)[7]
#test VIF
forcoefs = data.frame(method = 'Groundtruth', Predictor = 'ForestCover100',  coefficient = forcoefgt[[1]], lowCI = forCIlowgt, upCI = forCIupgt)
ggplot(forcoefs, aes(x = method, y = coefficient, ymin = lowCI, ymax = upCI)) +
  geom_point() +
  geom_errorbar() +
  ylim(0,25) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ylab('Forest Cover Coefficient')

#takes in glm model and dataset model is from
#assumes model only has forest cover and ac as variables.
getcutoff = function(model, dataset){
  prediction1 = predict(model, type = 'response')
  pred = prediction(prediction1, dataset$Presence)
  perf = performance(pred, 'sens', 'spec')
  
  cutoffgt = perf@alpha.values[[1]][which.max(perf@x.values[[1]]+perf@y.values[[1]])]
  
  ndata = data.frame(ForestCover100 = seq(0,1,0.001), ac = mean(dataset$ac))
  
  testpred = predict(model, ndata, type = 'response', re.form = NA)
  #find the value of forest cover where probability = cutoffgt
  threshold = testpred[testpred > cutoffgt]
  for_cut = ndata$ForestCover100[as.numeric(names(threshold)[1])]
  return(for_cut)
}

gtcutoff = getcutoff(gtglm2, gtSM)

# Random vs Routed analysis --------------------------------------------------------

#What format do we need the data in?
#We need to group the data by site, and then sum the number of detections for each site.
#We then need to add the covariates to the data.
sim_grpd = simSM %>%
  group_by(simulation_type, iteration,prop_samplers,loc_name, site_name, trajectory_dts) %>%
  summarise(total_dets = sum(detected, na.rm = T)) %>%
  mutate(Presence = ifelse(total_dets > 0, 1, 0))

#fill in missing days
days = unique(sim_grpd$trajectory_dts)
sites = unique(sim_grpd$site_name)

complete = crossing(simulation_type = unique(sim_grpd$simulation_type), iteration = unique(sim_grpd$iteration), prop_samplers = unique(sim_grpd$prop_samplers), site_name = sites, trajectory_dts = days)
#join loc_name with matching site_name values
complete = left_join(complete, site_info[,c(2,5)], by = c('site_name' = 'Site'))
names(complete)[6] = 'loc_name'

sim = sim_grpd %>%
  full_join(complete, by = c('simulation_type', 'iteration', 'prop_samplers', 'loc_name', 'site_name', 'trajectory_dts'))

#remove total dets column
sim = sim[,-7]
sim2 = pivot_wider(sim, names_from = 'trajectory_dts', values_from = 'Presence')

#create binary presence/absence variable per site
sim2$Presence = ifelse(rowSums(sim2[,c(6:12)], na.rm = T) > 0, 1, 0)

#add covariates - 'ForestCover100', 'secondary200km', 'house1000km' 'ac'
sim2 = left_join(sim2, gtSM[,c(2,71,26,27,94)], by = c('site_name' = 'Site'))



## Min forest cover for spider monkeys (raw) -------------------------------

#for each simulation type and prop_samplers and iteration, work out the min forest cover where a monkey is detected
sim2$simulation_type = factor(sim2$simulation_type, levels = c('Random', 'NearestNeighbour'))

min_sim = sim2 %>%
  group_by(simulation_type, iteration, prop_samplers) %>%
  summarise(min_fc = min(ForestCover100[Presence==1]))

ggplot(data = min_sim, aes(x = prop_samplers, y = min_fc, color = prop_samplers)) +
  geom_boxplot(aes(group = prop_samplers)) +
  labs(x = 'Sampling Intensity',
       y = 'Proportion of forest cover',
       color = 'Simulation Type') +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 16))+
  facet_wrap(~simulation_type) +
  geom_hline(yintercept = 0.79, linetype = 'dashed')+
  scale_x_continuous(breaks = seq(0,1,0.2))

## Threshold for spider monkey occupancy (GLM) -----------------------------------

#groundtruth coefs are in forcoefs
forcoefs$prop_samplers = 'gt'
forcoefs = forcoefs[,c(1,6,2,3,4,5)]
glmcoefs = data.frame(method = c(rep('Random',5), rep('NearestNeighbour',5)),prop_samplers = rep(c('0.2','0.4','0.6','0.8', '1'), 2), Predictor = 'ForestCover100', coefficient = NA, lowCI = NA, upCI = NA)
forcoefs= rbind(forcoefs, glmcoefs)

#for each method, this is the model formula:
simglm1 = glmer(Presence ~ ForestCover100 + Secondary200km 
              + House1000 + ac + (1|iteration),
              family = binomial, 
              data = sim2[sim2$simulation_type==simtype,])
simglm2 <- glmer(Presence ~ ForestCover100 + ac 
               + (1|iteration),
               family = binomial,
             data = sim2[sim2$simulation_type==simtype,])

#run one model, get the cutoff value
#another way of doing this is running this for each iteration, and getting mean cutoff with confidence intervals.
#get cutoffs for each method at each value of prop_samplers
#try this as one large model and then as iterative models to see if there is a difference

#one large model
simcoefs = NULL
for (simtype in unique(sim2$simulation_type)) {
  for (prop in unique(sim2$prop_samplers)) {
    df = filter(sim2, simulation_type == simtype & prop_samplers == prop)
    simglm1 = glmer(Presence ~ ForestCover100 + Secondary200km 
                    + House1000 + ac + (1|iteration),
                    family = binomial, 
                    data = df)
    simglm2 <- glmer(Presence ~ ForestCover100 + ac 
                     + (1|iteration),
                     family = binomial,
                     data = df)
    
    routecutoff = getcutoff(simglm2, df)
    #save coefficient and 95% confidence intervals
    s = summ(simglm2, confint = TRUE)
    sd = as.data.frame(s$coeftable)
    sd$prop_samplers = prop
    sd$simulation_type = simtype
    sd$cutoff = routecutoff
    simcoefs = rbind(simcoefs, sd)
  }
}

ggplot(simcoefs, aes(x = prop_samplers, y = cutoff, color = simulation_type)) +
  geom_boxplot(aes(group = prop_samplers)) +
  facet_wrap(~simulation_type) +
  geom_hline(yintercept = for_cut, linetype = 'dashed')


#PROGRESS 3/12: RE-RUN THIS BIG LOOP AND GET WEIGHTED MEAN OF ALL THE COEFFICIENTS :)
#add groundtruth data to sim_grpd
gtSM$simulation_type = 'Groundtruth'
gtSM$iteration = 1
gtSM$prop_samplers = 'GT'
#order to match sim_grpd
gtSM2 = gtSM[,c(18:20,12,1,2:9,11,13:17)]
names(gtSM2) = names(sim2)

sim2 = rbind(gtSM2, sim2)

#now try one model per iteration
simcoefs_i = NULL
for (simtype in unique(sim2$simulation_type)) {
  df1 = filter(sim2, simulation_type == simtype)
  for (iter in unique(df1$iteration)) {
    df2 = filter(df1, iteration == iter)
    for (prop in unique(df2$prop_samplers)) {
      df = filter(df2, prop_samplers == prop)
      print(paste0(simtype, ' ', iter, ' ', prop))
      simglm1 = glm(Presence ~ ForestCover100 + Secondary200km 
                      + House1000 + ac,
                      family = binomial, 
                      data = df)
      simglm2 <- glm(Presence ~ ForestCover100 + ac,
                       family = binomial,
                       data = df)
      
      routecutoff = getcutoff(simglm2, df)
      #save coefficient and 95% confidence intervals
      s = summ(simglm2, confint = TRUE)
      s2 = summ(simglm2)
      sd = as.data.frame(s$coeftable)
      sd$Predictor = rownames(sd)
      sd2 = as.data.frame(s2$coeftable)
      sd = left_join(sd, sd2[,c(1,2)], by = 'Est.')
      sd$prop_samplers = prop
      sd$simulation_type = simtype
      sd$iteration = iter
      sd$cutoff = routecutoff
      simcoefs_i = rbind(simcoefs_i, sd)
    }
  }
}

#make rownames into a column

simcoefs_i$simulation_type = factor(simcoefs_i$simulation_type, levels = c('Groundtruth', 'Random', 'NearestNeighbour', 'Routed'))
simcoefs_i$simulation_type[simcoefs_i$simulation_type=='NearestNeighbour'] = 'Routed'

simcoefs_i$prop_samplers = factor(simcoefs_i$prop_samplers)
F3B = ggplot(simcoefs_i[simcoefs_i$Predictor=='ForestCover100' & simcoefs_i$simulation_type != 'Groundtruth',], aes(x = prop_samplers, y = cutoff, fill = simulation_type)) +
  geom_point(aes(colour = simulation_type), position = position_jitterdodge(), alpha = 0.5)+
  geom_boxplot(alpha = 0.3) +
  geom_hline(yintercept = for_cut, linetype = 'dashed', colour = 'red', alpha =0.8) +
  labs(x = 'Sampling Intensity',
       y = 'Minimum forest cover treshold',
       fill = 'Sampling Strategy',
       colour = '') +
  scale_fill_manual(values=c("#E69F00", "#56B4E9")) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9"), guide = 'none') +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12)) +
  annotate('text', label = for_cut, x = 6, y =for_cut, vjust = 0.5) +
  coord_cartesian(xlim = c(1,5), clip = 'off')+
  ggtitle('b')

#fig3 
F3A + F3B + plot_layout(ncol = 1)


#plot responses to forest cover
#create weighted average of all coefficients. 
#each iteration gets a weight, which is the inverse of the standard error squared.
#the weighted mean is the estimate * weight, divided by the sum of all weights
#the 95% confidence interval is calculated from a weighted standard error
wm = simcoefs_i %>%
  filter(Predictor == 'ForestCover100') %>%
  group_by(simulation_type, prop_samplers) %>%
  summarize(
    weighted_mean = sum((1 / S.E.^2) * Est.) / sum(1 / S.E.^2),
    weighted_se = sqrt(1 / sum(1 / S.E.^2)),
    .groups = 'drop'
  ) %>%
  mutate(
    lower_ci = weighted_mean - 1.96 * weighted_se,
    upper_ci = weighted_mean + 1.96 * weighted_se
  )


ggplot(wm, aes(x = prop_samplers , y = weighted_mean, ymin = lower_ci, ymax = upper_ci, color = simulation_type)) +
  geom_point(position = position_dodge(width= 0.5)) +
  geom_errorbar(position = position_dodge(width = 0.5), width = 0.2) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  labs(x = 'Sampling Intensity',
       y = 'Forest Cover Coefficient (weighted mean)',
       color = 'Simulation Type') +
  theme_bw() +
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 16)) +
  scale_color_manual(values = c(c("black","#E69F00", "#56B4E9")))

#plotting the weighted mean is very similar to the  glm that has iteration as a random effect!


#next, plot the curves for each model as an example. 

#select iteration one of each simulation type, and plot curve with cutoff point
samplermods = NULL
for (i in 1:n_distinct(sim2$prop_samplers)) {
  prop = unique(sim2$prop_samplers)[i]
  df = subset(sim2, prop_samplers==i)
  glmac <- glm(Presence ~ ForestCover100^2+ac, family = binomial,
                  data = df)
  summary(glmac)
  samplermods[[i]] = glmac
  
  
  #plot coefficient
  forcoef = coef(glmac)[2]
  forCImin = confint(glmac)[2]
  forCImax = confint(glmac)[5]
  
  glmcoefs[glmcoefs$n_samplers==i,2] = forcoef
  glmcoefs[glmcoefs$n_samplers==i,3] = forCImin
  glmcoefs[glmcoefs$n_samplers==i,4] = forCImax
}

coefs = ggplot(glmcoefs, aes(x = n_samplers, y = for_coef, ymin = CI_low, ymax = CI_up)) +
  geom_point() +
  geom_errorbar() +
  ylim(0,25) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  ylab('Forest Cover Coefficient')

s1_cut = getcutoff(samplermods[[1]], sm_grpd[sm_grpd$n_samplers==1,])

s1 = ggpredict(samplermods[[1]], c('ForestCover100')) %>%
  plot() +
  geom_point(data = sm_grpd[sm_grpd$n_samplers==1,], aes(x = ForestCover100, y = Presence)) +
  labs(x = 'Forest Cover',
       y = 'Occupancy',
       title = paste0('Number of samplers: ', 1)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  geom_vline(xintercept = s1_cut, linetype = 'dashed') +
  annotate('text', label = s1_cut, x = s1_cut, y = 0.85, vjust = -0.5, angle = 90)
s1
s2_cut = getcutoff(samplermods[[2]], sm_grpd[sm_grpd$n_samplers==2,])

s2 = ggpredict(samplermods[[2]], c('ForestCover100')) %>%
  plot() +
  geom_point(data = sm_grpd[sm_grpd$n_samplers==2,], aes(x = ForestCover100, y = Presence)) +
  labs(x = 'Forest Cover',
       y = 'Occupancy',
       title = paste0('Number of samplers: ', 2)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) + 
  geom_vline(xintercept = s2_cut, linetype = 'dashed') +
  annotate('text', label = s2_cut, x= s2_cut, vjust = -0.5, angle = 90, y = 0.85)
s2

s4_cut = getcutoff(samplermods[[4]], sm_grpd[sm_grpd$n_samplers==4,])

s4 = ggpredict(samplermods[[4]], c('ForestCover100')) %>%
  plot() +
  geom_point(data = sm_grpd[sm_grpd$n_samplers==4,], aes(x = ForestCover100, y = Presence)) +
  labs(x = 'Forest Cover',
       y = 'Occupancy',
       title = paste0('Number of samplers: ', 4)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) + 
  geom_vline(xintercept = s4_cut, linetype = 'dashed')  +
  annotate('text', label = s4_cut, x = s4_cut, y = 0.85, vjust = -0.5, angle = 90)
s4

s8_cut = getcutoff(samplermods[[8]], sm_grpd[sm_grpd$n_samplers==8,])

s8 = ggpredict(samplermods[[8]], c('ForestCover100')) %>%
  plot() +
  geom_point(data = sm_grpd[sm_grpd$n_samplers==8,], aes(x = ForestCover100, y = Presence)) +
  labs(x = 'Forest Cover',
       y = 'Occupancy',
       title = paste0('Number of samplers: ', 8)) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  geom_vline(xintercept = s8_cut, linetype = 'dashed')  +
  annotate('text', label = s8_cut, x = s8_cut, y = 0.85, vjust = -0.5, angle = 90)
s8

gt_cut = getcutoff(fcglmac, gt)
gt_plt = ggpredict(fcglmac, c('ForestCover100')) %>%
  plot() +
  geom_point(data = gt, aes(x = ForestCover100, y = Presence)) +
  labs(x = 'Forest Cover',
       y = 'Occupancy',
       title = paste0('Groundtruth')) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) +
  geom_vline(xintercept = gt_cut, linetype = 'dashed')  +
  annotate('text', label = gt_cut, x = gt_cut, y = 0.85, vjust = -0.5, angle = 90)
  
  
gt_plt

library(patchwork)

gt_plt + s1 + s2 + s4 + s8 + coefs

s1 + s2 + s4 + s8 + plot_layout(nrow = 1)


min(sm_grpd$ForestCover100[sm_grpd$n_samplers==1 & sm_grpd$Presence==1])
min(sm_grpd$ForestCover100[sm_grpd$n_samplers==2 & sm_grpd$Presence==1])
min(sm_grpd$ForestCover100[sm_grpd$n_samplers==4 & sm_grpd$Presence==1])
min(sm_grpd$ForestCover100[sm_grpd$n_samplers==8 & sm_grpd$Presence==1])
