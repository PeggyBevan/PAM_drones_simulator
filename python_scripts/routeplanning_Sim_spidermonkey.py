#simulation with clustered route planning, using spider monkey data
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

# Load groundtruth data
spiders = pd.read_csv('data/SpiderMonkeys_Processeddata.csv')

clusters = spiders['Cluster'].unique()
#print(clusters)

np.random.seed(42)
for iter in range(1,51):
    all_exps = []
    exps_list = []
    print('Iteration {}'.format(iter))
    for cluster in clusters:
        if cluster in SKIP_LOCS:
            continue
        locs_sampled = []
        loc_dets = []
        
        detections = spiders[spiders['Cluster'] == cluster]
        cols_to_select = detections.columns[1:9].tolist() + ['Lat', 'Long']
        dets_sub = detections[cols_to_select].copy()
        dets_mat = dets_sub.iloc[:,1:8].to_numpy()

        days = cols_to_select[1:8]
        sites = dets_sub['Site'].tolist()
        loc_site_coords = dets_sub[['Lat', 'Long']].to_numpy()

        #find the charger site
        dist_mat = calculate_distance_matrix(loc_site_coords)
        #set charger site as the one closest to the middle
        charger_site_ix = find_charger_site(loc_site_coords)
        max_samplers = len(sites)
        for n_samp_ix, prop_samplers in enumerate([0.2, 0.4, 0.6, 0.8, 1.0]):
            
            N_SAMPLERS = int(max_samplers * prop_samplers)
            if N_SAMPLERS == 0: 
                N_SAMPLERS = 1
            #split the sites into clusters
            n_clusters = min(N_SAMPLERS, len(sites))
            labels = cluster_sites(loc_site_coords, N_SAMPLERS = n_clusters)

            #set up samplers at their first location
            all_samplers = []
            for i in range(N_SAMPLERS):
                sampler = dict()
                sampler['loc_name'] = cluster
                sampler['sampler_num'] = i
                sampler['curr_site_ix'] = charger_site_ix
                sampler['detected'] = []
                sampler['trajectory_sites'] = []
                sampler['trajectory_dts'] = []
                sampler['trajectory_dts_ix'] = []
                sampler['battery_s'] = SAMPLER_BATTERY_S
                sampler['charger_site_ix'] = charger_site_ix
                sampler['n_samplers'] = N_SAMPLERS
                sampler['prop_samplers'] = prop_samplers
                sampler['charger_visits'] = []
                sampler['battery_level_s'] = []
                sampler['total_dist'] = 0
                all_samplers.append(sampler)

            timestep = 0
            visited = np.zeros(dist_mat.shape[0], dtype=np.int32)  # keep track of how many times each site has been visited
            while timestep < len(days):
                curr_day = days[timestep] #this corresponds to the original date which can be
                #fetched from the loc_dts list

                # Keep track of which sites will be populated next so two samplers don't go to same place
                chosen_next_sites = []
                for smplr in all_samplers:
                    cluster_label = smplr['sampler_num']
                    
                    #which sites can this sampler go to? 
                    sites_in_cluster = np.where(labels == cluster_label)[0]  
                    #in the first round, send each sampler to it's cluster.
                    if timestep == 0:
                        starter_opts = sites_in_cluster
                        #choose the first site available
                        smplr['curr_site_ix'] = starter_opts[0]
                    
                    #is the charger site recording in the first hour?
                    #in this dataset, there are very few 'Nones' and there is always recording on day 1
                    #keep track of where the sampler has been:
                    smplr['trajectory_sites'].append(smplr['curr_site_ix']) #adds the current site to trajectory
                    smplr['trajectory_dts'].append(curr_day) #adds the actual date to trajectory
                    smplr['trajectory_dts_ix'].append(timestep) #adds the actual date to trajectory
                    smplr['battery_level_s'].append(smplr['battery_s'])
                    
                    # Make sure battery hasn't died on sampler
                    assert(smplr['battery_s'] >= 0)
                    # If at site with charger then recharge battery (assume recording still occurs while recharging)
                    if smplr['curr_site_ix'] == smplr['charger_site_ix']:
                        smplr['battery_s'] = SAMPLER_BATTERY_S
                    # Add if monkey's were detected on this day at this site
                    curr_dets = dets_mat[smplr['curr_site_ix'], timestep] #sample detections at current site index and current hour
                    #if there was a spider monkey detection, add to the 'detected' list
                    smplr['detected'].append(curr_dets)
                    #include current site in visited
                    visited[smplr['curr_site_ix']] += 1
                    
                    # Figure out where the sampler can go next within this cluster
                        # and remove places where other samplers are already going
                    if timestep < len(days)-1:
                        next_site_opts = sites_in_cluster[np.where(dets_mat[sites_in_cluster, timestep+1] != "NA")[0]]
                        # return sites that are not already chosen by another sampler
                        next_site_opts = list(filter(lambda a: a not in chosen_next_sites, next_site_opts))
                        #exclude current site from options so it doesn't get stuck at one site
                        next_site_opts = [i for i in next_site_opts if i != smplr['curr_site_ix']]
                    else:
                        next_site_opts = []
                    
                    next_site_opts_lessvisited = [i for i in next_site_opts if visited[i] < visited[smplr['curr_site_ix']]]
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
                        smplr['total_dist'] += trave_to_neighbor
                        smplr['charger_visits'].append(0)
                        chosen_next_sites.append(nearest)
                        #visited[nearest] += 1   
                    #if there are no next sites, then stay at the same site.
                    #current site stays the same, battery stays the same
                timestep += 1
            #end while loop
            #Save sampling data for this location   
            #change to pd format
            samp_list = [] #turn each sampler trajectory for this cluster into a dataframe
            for samp in all_samplers:
                dets_df = pd.DataFrame(samp)
                site_list = [sites[i] for i in samp['trajectory_sites']] #convert indexes to site names
                dets_df['site_name'] = site_list
                samp_list.append(dets_df)
            samp_df = pd.concat(samp_list)   
                
            locs_sampled.append(all_samplers) #save entire dictionary
            loc_dets.append(samp_df) #save dataframe of samplers to list
            if SAVE_TRAJECTORY_FIGS:
                fig_save_dir = os.path.join('figs', 'NNClusRoute_spidermonkey-n_samplers-{}'.format(N_SAMPLERS))
                if not os.path.exists(fig_save_dir): os.makedirs(fig_save_dir)
                #create matrix of off/on
                loc_on_off_mat = np.array([[1 if cell != 'NA' else 0 for cell in row] for row in dets_mat])
                
                xlab_days = np.arange(0,7)
                xlab = [day.replace('.', ' ') for day in days]
                #xlabs = [dt.strftime('%Y-%m-%d') for dt in xlab_dts]

                fig2 = plt.figure(figsize=(15,5))    
                axs1 = fig2.subplots(1, 2, width_ratios=[3, 1])
                axs1 = np.ravel(axs1)
                plt.sca(axs1[0])
                plt.matshow(loc_on_off_mat, aspect='auto', fignum=0, cmap = 'Greys', alpha = 0.5)
                for sampler_ix, sampler in enumerate(all_samplers):
                    traj_dts = [samp for samp in sampler['trajectory_dts_ix']]
                    plt.plot(traj_dts, sampler['trajectory_sites'], label='Sampler {}'.format(sampler_ix))
                    plt.scatter(traj_dts, sampler['trajectory_sites'])
            
                plt.title('Sampler trajectories')
                plt.xlabel('Time')
                plt.ylabel('Sites')
                plt.gca().set_xticks(xlab_days)
                plt.gca().set_xticks(xlab_days)
                plt.gca().set_xticklabels(xlab, fontsize=10)
                plt.legend(bbox_to_anchor=(-0.1, 1), loc='upper right')
                plt.yticks(range(len(sites)), sites)

                plt.sca(axs1[1])
                #add second fig with map
                lats = [coord[0] for coord in loc_site_coords]
                lons = [coord[1] for coord in loc_site_coords]
                plt.scatter(lons, lats, c=labels, cmap='viridis')
                plt.scatter(loc_site_coords[charger_site_ix][1], loc_site_coords[charger_site_ix][0], color='red')
                for i, (lon, lat) in enumerate(zip(lons, lats)):
                    plt.text(lon, lat, str(i), fontsize=8)

                plt.savefig(os.path.join(fig_save_dir, 'sampler_trajectories_{}.png'.format(cluster)))
                plt.clf()
                plt.close(fig2)
            #end trajectory figs
        #end cluster loop
        all_exps.append(locs_sampled)
        locs_pd = pd.concat(loc_dets)
        exps_list.append(locs_pd)
    #end N_SAMPLERS loop
    exps_pd = pd.concat(exps_list)
    exps_pd['iteration'] = iter
    csv_save_dir = os.path.join('data', 'Iterations_proportional_set', 'NNClusRouteSM')
    if not os.path.exists(csv_save_dir): os.makedirs(csv_save_dir)
    exps_pd.to_csv(os.path.join(csv_save_dir, 'NNClusRouetSM_Iter{:02}.csv'.format(iter)), index=False)
#end iter loop
                        

                        
                    
