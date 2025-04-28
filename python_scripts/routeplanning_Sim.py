#testing out route planning
import os
import pickle
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
#import utils
from utils import get_sites_with_gps, read_site_file, calculate_distance_matrix, find_nearest_neighbor, cluster_sites, find_charger_site
def is_det_valid(det, filter_specs):
    return det['confidence'] > 0.8 and det['scientific_name'] in filter_specs
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
import configs

SAMPLER_SPEED = configs.SAMPLER_SPEED          # in meters per second
SAMPLER_BATTERY_S = configs.SAMPLER_BATTERY_S  # minutes in seconds
SKIP_LOCS = configs.SKIP_LOCS #skipping luna due to ambiguous site names
VALID_SPECS = configs.VALID_SPECS
SAVE_TRAJECTORY_FIGS = configs.SAVE_TRAJECTORY_FIGS
SAVE_SPECIESSAMPLED_FIGS = configs.SAVE_SPECIESSAMPLED_FIGS
comb_dets_dir = configs.comb_dets_dir
TAKEOFF_PEN = configs.TAKEOFF_PEN

#move working directory up one folder
os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


all_comb_d_fs = os.listdir(comb_dets_dir)
sub_comb_d_fs = all_comb_d_fs[:1]

all_site_infos = get_sites_with_gps(verbose=True)
all_site_info_names = np.asarray(['{} {}'.format(s['loc_name'], s['site_name']) for s in all_site_infos])

#set up figure
if SAVE_SPECIESSAMPLED_FIGS:
    fig = plt.figure(figsize=(17,12))
    axs = fig.subplots(2, 2, sharex='all', sharey='all')
    axs = np.ravel(axs)

np.random.seed(42)
for iter in range(1,51):
    all_exps = []
    exps_list = []
    print('Iteration {}'.format(iter))
    for comb_d_f in all_comb_d_fs:
    #for comb_d_f in sub_comb_d_fs:
        if comb_d_f.split('.')[0] in SKIP_LOCS: continue
        
        locs_sampled = [] #create list of sites sampled
        loc_dets = [] #create list of ?

        all_missed_specs = []
        all_common_specs = []
    
        # Load detection file for this cluster
        coordinates, loc_dets_mat, loc_sites, loc_dts = read_site_file(comb_d_f, all_site_infos, all_site_info_names)
        #these have already been filtered to remove sites with no GPS
        #create distance matrix
        
        #filter out hours where no recording is happening
        loc_dets_df = pd.DataFrame(loc_dets_mat)
        loc_dets_df = loc_dets_df.dropna(how = 'all', axis=1)
        #remove rows where all values are None - the index values will stay the same even if rows are removed.
        loc_dets_df = loc_dets_df.dropna(how = 'all', axis=0)
        
        #save date values so they can be converted back to real dates
        nn_dts = loc_dets_df.columns.tolist()
        nn_sites = loc_dets_df.index.tolist()
        loc_sites = [loc_sites[i] for i in nn_sites]
        #convert back to numpy array - this still includes some 'None' values
        nn_dets_mat = loc_dets_df.values
        coordinates = [coordinates[i] for i in nn_sites]
        dist_mat = calculate_distance_matrix(coordinates)
        charger_site_ix = find_charger_site(coordinates)
        max_samplers = len(loc_sites)
        
        for n_samp_ix, prop_samplers in enumerate([0.2, 0.4, 0.6, 0.8, 1.0]):
            
            N_SAMPLERS = int(max_samplers * prop_samplers)
            if N_SAMPLERS == 0: 
                N_SAMPLERS = 1    
            #split the sites into clusters
            '''
            #need to assign each sampler to a cluster
            #cluster number is the same as sampler number (but clusters start from 0)
            #but, if sampler number is higher than cluster number, it should just stay at charger site.
            #labels is indexed from 0
            #N_Samplers starts from 1
            #charging site is indexed from 0
            #n_clusters is equal to the actual number of sites, but the highest cluster label will be n_clusters - 1
            #sampler ID is indexed from 0 - so this should line up with cluster labels.
            '''
            n_clusters = min(N_SAMPLERS, len(loc_sites))
            labels = cluster_sites(coordinates, N_SAMPLERS = n_clusters)
            # Combine labels and loc_sites into a dictionary
            # create dictionary to save sampler trajectories
            all_samplers = []
            for i in range(N_SAMPLERS):
                sampler = dict()
                sampler['loc_name'] = comb_d_f
                sampler['sampler_num'] = i
                sampler['curr_site_ix'] = charger_site_ix
                sampler['sampled_specs'] = []
                sampler['trajectory_sites'] = []
                sampler['trajectory_dts'] = []
                sampler['trajectory_dts_ix'] = []
                sampler['battery_s'] = SAMPLER_BATTERY_S
                sampler['charger_site_ix'] = charger_site_ix
                sampler['n_samplers'] = N_SAMPLERS
                sampler['prop_samplers'] = prop_samplers
                sampler['charger_visits'] = [] #this will be 1 if the drone returns to the charger because of low battery
                sampler['battery_level_s'] = []
                sampler['total_dist'] = 0
                all_samplers.append(sampler)

            #print('Running simulation at {}'.format(comb_d_f))

            # Each step of the simulation proceeds by one hour
            timestep = 0
            visited = np.zeros(dist_mat.shape[0], dtype=np.int32)   
            while timestep < len(nn_dts)-1:
                curr_hr_ix = nn_dts[timestep] #this corresponds to the original date which can be
                #fetched from the loc_dts list
                
                #in the first time step of the simulation, the drone will
                #start at a site which was recording in the real dataset    
                chosen_next_sites = []
                for smplr in all_samplers:
                    cluster_label = smplr['sampler_num']
                    
                    #which sites can this sampler go to? 
                    sites_in_cluster = np.where(labels == cluster_label)[0]  
                    #in the first round, send each sampler to it's cluster.
                    if timestep == 0:
                        starter_opts = np.where((nn_dets_mat[sites_in_cluster, 0] != None))[0]
                        if len(starter_opts) == 0:
                            smplr['curr_site_ix'] = sites_in_cluster[0]
                            #if no sites in cluster available, go to any site in cluster to wait.
                        else:
                            #choose the first site available
                            smplr['curr_site_ix'] = starter_opts[0]
                    
                    # Keep track of where each sampler has been
                    smplr['trajectory_sites'].append(smplr['curr_site_ix'])  # adds the current site (index) to trajectory
                    smplr['trajectory_dts'].append(loc_dts[curr_hr_ix])  # adds the actual date to trajectory
                    smplr['trajectory_dts_ix'].append(curr_hr_ix)
                    smplr['battery_level_s'].append(smplr['battery_s'])
                    # Make sure battery hasn't died on sampler
                    assert(smplr['battery_s'] >= 0)
                    
                    # Add any valid detections from this site/hour to the sampler
                    curr_dets = nn_dets_mat[smplr['curr_site_ix'], timestep]  # sample detections at current site index and current hour
                    # if there are detections, extract common name and add them to the list of sampled species
                    if curr_dets != None:
                        det_list = [d['common_name'] for d in curr_dets if curr_dets and is_det_valid(d, VALID_SPECS)]
                    else:
                        det_list = []
                    smplr['sampled_specs'].append(det_list)
                    #include current site in visited
                    visited[smplr['curr_site_ix']] += 1
                    
                    # Figure out where the sampler can go next within this cluster
                    # and remove places where other samplers are already going
                    next_site_opts = sites_in_cluster[np.where((nn_dets_mat[sites_in_cluster, timestep+1] != None))[0]]
                    # return sites that are not already chosen by another sampler
                    next_site_opts = list(filter(lambda a: a not in chosen_next_sites, next_site_opts))
                    #exclude current site from options so it doesn't get stuck at one site
                    next_site_opts = [i for i in next_site_opts if i != smplr['curr_site_ix']]

                    # If at site with charger then recharge battery (assume recording still occurs while recharging)
                    if smplr['curr_site_ix'] == smplr['charger_site_ix']:
                        smplr['battery_s'] = SAMPLER_BATTERY_S
                    
                    #if there are options in next_site_opts which are false in 'visited', 
                    #we want to favour these.
                    next_site_opts_lessvisited = [i for i in next_site_opts if visited[i] < visited[smplr['curr_site_ix']]]

                    #print('Choosing next site')
                    # Finding next site: nearest neighbor:
                    #if there are available sites to visit, that haven't already been visited
                    if len(next_site_opts_lessvisited) > 0:
                        nearest, trave_to_neighbor = find_nearest_neighbor(next_site_opts = next_site_opts_lessvisited, 
                                                                        dist_mat = dist_mat,
                                                                        curr_site_ix = smplr['curr_site_ix'],
                                                                        curr_battery = smplr['battery_s'],
                                                                        charger_site_ix = smplr['charger_site_ix'])                      
                    elif len(next_site_opts) > 0:
                    #if there are sites available that have already been visited
                        nearest, trave_to_neighbor = find_nearest_neighbor(next_site_opts = next_site_opts, 
                                                                        dist_mat = dist_mat,
                                                                        curr_site_ix = smplr['curr_site_ix'], 
                                                                        curr_battery = smplr['battery_s'], 
                                                                        charger_site_ix = smplr['charger_site_ix'])                        
                    else: 
                        nearest = smplr['curr_site_ix']
                        trave_to_neighbor = 0
                    #if nearest = None, it means that the sampler can't make it to the next site & needs to return to charger.
                    #but, if len(next_site_opts) is 0, it just means there are no sites available to visit, and the sampler can stay at current site.
                    #so nearest = current site.
                    
                    #if no potential next sites, return to charger
                    if nearest is None:
                        smplr['curr_site_ix'] = smplr['charger_site_ix']
                        smplr['charger_visits'].append(1)
                        smplr['total_dist'] += dist_mat[smplr['curr_site_ix'], smplr['charger_site_ix']]
                        #include battery consumption to return to charger
                        # smplr['battery_s'] -= (dist_mat[smplr['curr_site_ix'], smplr['charger_site_ix']] / SAMPLER_SPEED) + TAKEOFF_PEN
                        #print('no next site available, staying put')
                    else:
                        #if nearest is not none, change current site to 'nearest'
                        smplr['curr_site_ix'] = nearest
                        smplr['battery_s'] -= trave_to_neighbor
                        chosen_next_sites.append(nearest)
                        smplr['charger_visits'].append(0)
                        smplr['total_dist'] += trave_to_neighbor
                        #visited[nearest] += 1   
                    #if there are no next sites, then stay at the same site.
                    #current site stays the same, battery stays the same
                    
                    
                timestep +=1
            #end while loop

            #Save sampling data for this sampler group & change to pd format
            samp_list = [] #turn each sampler trajectory for this cluster into a dataframe
            for samp in all_samplers:
                dets_df = pd.DataFrame(samp)
                site_list = [loc_sites[i] for i in samp['trajectory_sites']] #convert indexes to site names
                dets_df['site_name'] = site_list
                samp_list.append(dets_df)
            samp_df = pd.concat(samp_list)   
            
            locs_sampled.append(all_samplers) #save entire dictionary
            loc_dets.append(samp_df) #save dataframe of samplers to list
            # Find all species captured by the samplers across all sites
            all_sampled_specs = []
            for sampler_ix, sampler in enumerate(all_samplers):
                all_sampled_specs.extend(spec for sublist in sampler['sampled_specs'] for spec in sublist)
            unq_sampler_specs = np.unique(all_sampled_specs)

            # Find all species that would've been captured by 100% sampling
            # need to work out how to do this without nn_ixs
            # all_nn_dets = loc_dets_mat[nn_ixs]
            nn_ixs = np.where(nn_dets_mat != None)
            dets_mat2 = nn_dets_mat[nn_ixs]
            all_valid_dets = []
            for nn_dets in dets_mat2:
                all_valid_dets.extend([d for d in nn_dets if is_det_valid(d, VALID_SPECS)])
            all_dets_specs = [d['common_name'] for d in all_valid_dets]
            unq_all_specs = np.unique(all_dets_specs)

            # Figure out which species were captured or missed by the samplers
            all_missed_specs.extend(list(set(unq_all_specs) - set(unq_sampler_specs)))
            all_common_specs.extend(list(set(unq_all_specs).intersection(unq_sampler_specs)))

            #print('Sampled {}/{} species'.format(len(np.unique(all_sampled_specs)), len(np.unique(all_dets_specs))))

            if SAVE_TRAJECTORY_FIGS:
                fig_save_dir = os.path.join('figs', 'NNClusRoute_sim-n_samplers-{}'.format(N_SAMPLERS))
                if not os.path.exists(fig_save_dir): os.makedirs(fig_save_dir)
                #create matrix of off/on
                loc_on_off_mat = np.array([[1 if cell is not None else 0 for cell in row] for row in loc_dets_mat])
                #crop this to same shape as nn_dets_mat
                samp_start_ix = nn_dts[0]
                samp_end_ix = nn_dts[len(nn_dts)-1]
                loc_on_off_mat = loc_on_off_mat[:, samp_start_ix:samp_end_ix+1]

                tot_hrs = int((loc_dts[nn_dts[len(nn_dts)-1]] - loc_dts[nn_dts[0]]).total_seconds() // 3600)             
                xlab_hrs = np.linspace(0, tot_hrs, num=10)
                xlab_tds = [timedelta(seconds=hr*3600) for hr in xlab_hrs]
                xlab_dts = [loc_dts[nn_dts[0]] + td for td in xlab_tds]
                xlabs = [dt.strftime('%Y-%m-%d') for dt in xlab_dts]

                fig2 = plt.figure(figsize=(15,5))    
                axs1 = fig2.subplots(1, 2, width_ratios=[3, 1])
                axs1 = np.ravel(axs1)
                plt.sca(axs1[0])
                plt.matshow(loc_on_off_mat, aspect='auto', fignum=0, cmap = 'Greys', alpha = 0.5)
                for sampler_ix, sampler in enumerate(all_samplers):
                    traj_dts = [samp - samp_start_ix for samp in sampler['trajectory_dts_ix']]
                    plt.plot(traj_dts, sampler['trajectory_sites'], label='Sampler {}'.format(sampler_ix))
                    plt.scatter(traj_dts, sampler['trajectory_sites'])
                    for dt, site, visit in zip(traj_dts, sampler['trajectory_sites'], sampler['charger_visits']):
                         if visit == 1:
                            plt.scatter(dt+1, site+1, color='red')
            
                plt.title('Sampler trajectories')
                plt.xlabel('Time')
                plt.ylabel('Sites')
                plt.gca().set_xticks(xlab_hrs)
                plt.gca().set_xticks(xlab_hrs)
                plt.gca().set_xticklabels(xlabs, fontsize=10)
                plt.legend(bbox_to_anchor=(-0.1, 1), loc='upper right')
                plt.yticks(range(len(loc_sites)), loc_sites)

                plt.sca(axs1[1])
                #add second fig with map
                lats = [coord[0] for coord in coordinates]
                lons = [coord[1] for coord in coordinates]
                plt.scatter(lons, lats, c=labels, cmap='viridis')
                plt.scatter(coordinates[charger_site_ix][1], coordinates[charger_site_ix][0], color='red')
                for i, (lon, lat) in enumerate(zip(lons, lats)):
                    plt.text(lon, lat, str(i), fontsize=8)

                plt.savefig(os.path.join(fig_save_dir, 'sampler_trajectories_{}.png'.format(comb_d_f.split('.')[0])))
                plt.clf()
            #end trajectory figs
        #end comb_d_f loop
        
        #save output from this iteration - this site
        all_exps.append(locs_sampled)
        locs_pd = pd.concat(loc_dets)
        exps_list.append(locs_pd)
        
        if SAVE_SPECIESSAMPLED_FIGS: #THIS WON'T WORK FOR PROPORTIONAL SAMPLER METHOD
            missed_specs, missed_spec_counts = np.unique(all_missed_specs, return_counts=True)
            sampled_specs, sampled_spec_counts = np.unique(all_common_specs, return_counts=True)

            specs = np.unique(np.hstack((all_missed_specs,all_common_specs)))
            specs_num_missed = []
            specs_num_sampled = []
            for s in specs:
                num_sampled = sampled_spec_counts[np.where((sampled_specs == s))[0]]
                num_missed = missed_spec_counts[np.where((missed_specs == s))[0]]

                num_sampled = num_sampled[0] if len(num_sampled) > 0 else 0
                num_missed = num_missed[0] if len(num_missed) > 0 else 0
                print('{}: {} sampled, {} missed'.format(s, num_sampled, num_missed))

                specs_num_sampled.append(num_sampled)
                specs_num_missed.append(num_missed)

            sort_ix = np.argsort([a + b for a, b in zip(specs_num_missed, specs_num_sampled)])[::-1]
            specs = specs[sort_ix]
            specs_num_missed = np.asarray(specs_num_missed)[sort_ix]
            specs_num_sampled = np.asarray(specs_num_sampled)[sort_ix]

            plt.figure(fig.number)
            plt.sca(axs[n_samp_ix])
            plt.bar(specs, specs_num_sampled, label='Sampled')
            plt.bar(specs, specs_num_missed, bottom=specs_num_sampled, label='Missed')
            plt.legend()
            plt.xticks(rotation=90)
            plt.xlabel('Species')
            plt.ylabel('Locations species detected')
            plt.title('{} samplers'.format(N_SAMPLERS))
        #end save_species_sampled_figs

    #end Cluster loop
    if SAVE_SPECIESSAMPLED_FIGS: #this will no longer work for proportional sampling
        plt.suptitle('Species occurrence detectability')
        plt.tight_layout()
        plt.savefig(os.path.join('figs', 'NNClusRoute_specs_samplability.png'))
    
    exps_pd = pd.concat(exps_list)
    exps_pd['iteration'] = iter
    csv_save_dir = os.path.join('data', 'Iterations_proportional_set', 'NNClusRoute')
    if not os.path.exists(csv_save_dir): os.makedirs(csv_save_dir)
    exps_pd.to_csv(os.path.join(csv_save_dir, 'NNClusRoute_Iter{:02}.csv'.format(iter)), index=False)
#end iter loop

'''    
#looking at cluster size for each iteration
#only look at where n_samplers = 2
#what do we consider problematic? 
    #if the cluster size ratio is consistently skewed (>50% of iterations), this will impact
    #data collection by the samplers.
    #what counts as skewed? if there are 5 sites, a split of 3/2 is ok. but a split of 4/1 is not. 
    #if there are 10 sites, a split of 6/4 is ok. but a split of 7/3 is not. 
    #if there are 20 sites, a split of 11/9 is ok, a split of 12/8 is acceptable, even 13/7. 12/6 is too much.
    #so the maximum split we want is 60/40. This is 2/5.
    
labels_clusters_2s = [dic for dic in labels_clusters if dic['n_samplers'] == 2]    
clust_size = []
for dic in labels_clusters_2s:
    n_sites = len(dic['labels'])
    first_clust_size = len(np.where(dic['labels'] == 0)[0])
    second_clust_size = len(np.where(dic['labels'] == 1)[0])
    #is the min cluster size less than the ideal cluster size?
    if min(first_clust_size, second_clust_size) < (n_sites*0.3):
        crossed = 1
    else: crossed = 0
    clust = {'loc_name': dic['loc_name'],
             'iteration': dic['iteration'],
             'crossed': crossed
    }
    clust_size.append(clust)


#get number of instances where crossed = 1
sum_crossed = sum([dic['crossed'] for dic in clust_size])
prop_crossed = (sum_crossed/len(clust_size))*100

#30% of the time, the cluster size is less than the ideal cluster size.

#isolate the instances where the cluster size is less than the ideal cluster size
crossed_exps = [dic for dic in clust_size if dic['crossed'] == 1]
#list the unique loc_names in crossed_exps
unique_crossed_locs = np.unique([dic['loc_name'] for dic in crossed_exps])
#create a frequency table of how many tiems each loc_name is in crossed_exps
freq_crossed_locs = {loc: len([dic for dic in crossed_exps if dic['loc_name'] == loc]) for loc in unique_crossed_locs}
#i guess we care if there is an uneven number of samples in more than 50% of iterations. 
#single out loc_names where the frequency is greater than 25
freq_crossed_locs = {k: v for k, v in freq_crossed_locs.items() if v > 25}

#Where is this a problem?
# Buho.pickle'): 29  - 4/2 split. just over half the time, not such a big issue
#  np.str_('Cira.pickle'): 34, - There are only 4 sites, so sometimes it will be a 3/1 split
#  np.str_('Danta.pickle'): 49,  - This is quite problematic, 4/10 split nearly 100% of iterations
#  np.str_('Piro.pickle'): 31, 5/2 split, just over half the time.
#  np.str_('SQ259.pickle'): 28, only just over half the time, not a major difference
#  np.str_('SQ283.pickle'): 49, - 5/2 split, nearly 100% of iterations
#  np.str_('Sabalo.pickle'): 30, - most of the time its a 50/50 split
#  np.str_('Teak.pickle'): 49} - there are two clusters split 2/6 geographically.

#Danta, SQ283 and Teak are the most problematic.
'''