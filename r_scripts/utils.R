#functions for PAM_drones_simulator

# Functions ---------------------------------------------------------------

format_dets = function(df) {
  #convert sampled species list into usable format
  #outputs one row per detection
  df$sampled_specs = gsub("\\[|]", '', df$sampled_specs)
  df = separate(df, sampled_specs, into = paste('detection_', 1:20, sep = ''),
                sep = ',', extra = 'merge', fill = 'right', remove = FALSE)
  #convert to longform - select sites with detections 
  df = pivot_longer(df, cols = c(7:26))
  #don't want to delete all the sites.
  #remove NAs
  df = df[!is.na(df$value),]
  #no detection, but sampled, = NA
  df$value[df$value==""] = NA
  #remove 'sampled_specs' and 'name' cols
  df = df[,c(1,2,3,4,5,7,8,10)]
  #tidy up species names
  #remove extra quotation marks
  df$value = gsub(" [']", "", df$value)
  df$value = gsub("[']", "", df$value)
  return(df)
}

mod_sim_csv = function(simdf, gt) {
  #assumes that simdf has '.pickle' in loc_name
  #and column names "n_samplers", 'prop_samplers', 'total_dist', "loc_name","site_name","datetime","sampled_specs", "loc_site"
  simdf$loc_name = gsub('.pickle', '', simdf$loc_name)
  simdf$site_name = gsub(" ", "-", simdf$site_name)
  simdf$loc_site = paste0(simdf$loc_name, '_', simdf$site_name)
  simdf = simdf %>%
    mutate(site_name = case_match(site_name,
                                  "Audio-1-92" ~ 'Audio-1',
                                  "Audio-4-375" ~ 'Audio-4',
                                  "Audio-5-58" ~ 'Audio-5',
                                  "Audio-6-673" ~ 'Audio-6',
                                  "Danta-14-grassland" ~ 'Audio-14',
                                  .default = site_name))
  
  #update loc_site values
  simdf$loc_site = paste0(simdf$loc_name, '_', simdf$site_name)
  simdf = simdf[,c('n_samplers', 'prop_samplers', 'loc_name', 'site_name', 'trajectory_dts','sampled_specs', "loc_site", 'total_dist')]
  #change detections column name to match simdf
  names(simdf)[5] = 'datetime'
  
  #remove duplicates - where n_samplers, loc_site and datetime are the same.
  #table(simdf$n_samplers, simdf$loc_name)
  simdf = simdf %>%
    group_by(n_samplers) %>%
    distinct(.keep_all = TRUE)
  #remove rows where not recording in groundtruth
  return(simdf)
}

get_speciesRichness = function(det_mats) {
  #now with species number per site
  speciesn = NULL
  NSamplers = length(det_mats)
  for (n in 1:NSamplers) {
    species = as.data.frame(specnumber(det_mats[[n]]))
    simpsons = as.data.frame(diversity(det_mats[[n]], index = 'simpson'))
    species$loc_name = rownames(species)
    species = cbind(species, simpsons)
    dets = as.data.frame(det_mats[[n]])
    dets$loc_name = rownames(dets)
    #dets = pivot_wider(dets, names_from = Var2, values_from = Freq)
    species1 = left_join(species, dets, by = 'loc_name')
    species1$experiment = names(det_mats)[n]
    
    speciesn = rbind(speciesn, species1)
  }
  names(speciesn)[1] = 'species_number'
  names(speciesn)[3] = 'simpsons_div'
  return(speciesn)
}

create_summary_df = function(file_path, valid_specs = valid_specs) {
  ###
  #'read in results from each iteration of drone simulation, and convert to a summary table
  #'  one row per site per experiment per iteration.
  #'  plus one row per land-use type per experiment per iteration
  #'  extract species richness, simpsons index and individual species detections per site.
  ###
  summary_df = NULL
  files = list.files(file_path)
  for (i in 1:length(files)) {
    simdf = read.csv(paste0(file_path, files[i]))
    print(paste0('File Loaded ', files[i]))
    #format detections csv to remove repeats
    simdf = mod_sim_csv(simdf, gt)
    #format dections to make one row per species
    simdf = format_dets(simdf)
    #remove recording periods that do not match with groundtruth
    #semi_join returns rows from simdf that have a match in gt
    simdf <- simdf %>%
      semi_join(gt, by = c("datetime", "loc_site"))
    #add land use data - left_join adds habitat column, matched by loc_site
    simdf = left_join(simdf, site_info[,c('loc_site', 'Habitat' )], by = 'loc_site')
    #create matrix of species detections per cluster
    det_mats_loc = NULL
    det_mats_lu = NULL
    for (samp in 1:length(unique(simdf$prop_samplers))) {
      sampler_prop = unique(simdf$prop_samplers)[samp]
      subset = simdf[simdf$prop_samplers==sampler_prop,]
      #save matrix by cluster/loc_name
      #convert matrix values to numeric, whilst keeping rownames and colnames.
      mat = apply(as.matrix(table(subset$loc_name, subset$value)), c(1,2), as.numeric)
      #save matrix by land-use type
      mat_lu = apply(as.matrix(table(subset$Habitat, subset$value)), c(1,2), as.numeric)
      #make sure it has 19 columns and add any missing species columns
      for (spec in valid_specs) {
        if (!spec %in% colnames(mat)) {
          mat <- cbind(mat, spec = 0)
          colnames(mat)[ncol(mat)] <- spec }
        if (!spec %in% colnames(mat_lu)) {
          mat_lu <- cbind(mat_lu, spec = 0)
          colnames(mat_lu)[ncol(mat_lu)] <- spec }}
      
      det_mats_loc[[samp]] = mat
      det_mats_lu[[samp]] = mat_lu
      rm(subset)
    }
    names(det_mats_loc) = unique(simdf$prop_samplers)
    names(det_mats_lu) = unique(simdf$prop_samplers)
    #get metrics for each site, save this. 
    #should be approx 200 rows per iteration instead of 60,000
    speciesnum = get_speciesRichness(det_mats_loc)
    speciesnum$group = 'Cluster'
    speciesn_LU = get_speciesRichness(det_mats_lu)
    speciesn_LU$group = 'LandUse'
    speciesnum = rbind(speciesnum, speciesn_LU)
    speciesnum$iteration = i
    summary_df = rbind(summary_df, speciesnum)
  }
  return(summary_df)
}

avg_evenness = function(file_path, gt) {
  ### Takes in list of files and groundtruth observations, 
  ## formats to one row per detection and calculates sampling frequency per hour
  full_even = NULL
  files = list.files(file_path)
  for (file in 1:length(files)) {
    simdf = read.csv(paste0(file_path, files[file]))
    print(paste0('File Loaded ', files[file]))
    #create loc_site combo
    simdf = mod_sim_csv(simdf, gt)
    simdf = format_dets(simdf)
    #remove duplicates
    simdf = simdf[!duplicated(simdf[c(1,2,3,4,5,6)]),]
    
    #remove recording periods that do not match with groundtruth
    #semi_join returns rows from simdf that have a match in gt
    simdf <- simdf %>%
      semi_join(gt, by = c("datetime", "loc_site"))
    
    #simplify to sampling frequency
    simdf$datetime2 = as.POSIXct(simdf$datetime, 
                                 format = '%Y-%m-%d %H:%M:%OS',
                                 tz = "America/Costa_Rica")
    simdf$time = format(simdf$datetime2, format = '%H:%M:%OS')
    #simdf$date = format(simdf$datetime2, format = '%Y-%m-%d')
    
    #calculate sampling freq per hour per site
    samp_freq = simdf %>%
      group_by(loc_name, loc_site, n_samplers, prop_samplers, time) %>%
      summarise(freq = n())
    
    #remove 9 & 16:00 as these are not part of original regime
    samp_freq = filter(samp_freq, !time %in% c('09:00:00', '16:00:00'))
    
    evenness = samp_freq %>%
      group_by(prop_samplers, loc_name) %>%
      summarise(total_visits = sum(freq),
                # Proportion of visits for each site-time combination
                proportions = list(freq / sum(freq)),
                # Shannon Index H
                H = -sum((freq / sum(freq)) * log(freq / sum(freq)),
                         na.rm = TRUE),
                # Maximum Shannon Index H_max
                H_max = log(n()),
                # Evenness
                evenness = H / H_max
      ) %>%
      ungroup()
    evenness$iteration = file
    full_even = rbind(full_even, evenness)
  } #end file loop
  return(full_even)
}

