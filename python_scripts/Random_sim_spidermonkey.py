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
for iter in range(1, 51):
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

        days = cols_to_select[1:8] #7 days
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
            while timestep < len(days):
                curr_day = days[timestep] #this corresponds to the original date which can be
                #fetched from the loc_dts list

                # Keep track of which sites will be populated next so two samplers don't go to same place
                chosen_next_sites = []
                for smplr in all_samplers:
                    #is the charger site recording in the first hour?
                    #in this dataset, there are very few 'Nones' and there is always recording on day 1
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
                    
                    # Figure out where we can go next and remove places where other samplers are already going
                    if timestep < len(days)-1: 
                        next_site_opts = np.where((dets_mat[:, timestep+1] != "NA"))[0]
                        #return sites what are not already chosen by another sampler
                        next_site_opts = list(filter(lambda a: a not in chosen_next_sites, next_site_opts))
                    else:
                        next_site_opts = []
                    
                    # If there are valid options, then choose a random next site
                    if len(next_site_opts) > 0:
                        next_site_ix = np.random.choice(next_site_opts, 1)[0]
                        next_site = sites[next_site_ix]
                        curr_site_ix = smplr['curr_site_ix']

                        # Determine if we can make it to the next site and still back to the charger
                        #current site latlong  = 
                        trave_to_next_site_s = (dist_mat[curr_site_ix, next_site_ix] / SAMPLER_SPEED) + TAKEOFF_PEN
                        trave_to_charger_from_next_site_s = (dist_mat[next_site_ix, charger_site_ix] / SAMPLER_SPEED) + TAKEOFF_PEN
                        
                        # Move the sampler to the next site and decrease battery accordingly
                        if smplr['battery_s'] >= trave_to_next_site_s + trave_to_charger_from_next_site_s:
                            smplr['battery_s'] -= trave_to_next_site_s
                            smplr['curr_site_ix'] = next_site_ix #move to next site
                            smplr['charger_visits'].append(0)
                            smplr['total_dist'] += dist_mat[curr_site_ix, next_site_ix]
                        else: #if not enough battery
                            trave_to_charger_s = (dist_mat[curr_site_ix, charger_site_ix] / SAMPLER_SPEED) + TAKEOFF_PEN
                            #at the moment, if current site is charger site, takeoff pen still given.
                            smplr['battery_s'] -= trave_to_charger_s
                            smplr['curr_site_ix'] = smplr['charger_site_ix'] #send drone back to charger.
                            smplr['charger_visits'].append(1)
                            smplr['total_dist'] += dist_mat[curr_site_ix, charger_site_ix]

                        chosen_next_sites.append(smplr['curr_site_ix'])
                    else:
                        smplr['charger_visits'].append(0)
                    #if there are no available sites, the drone will remain in position.
                    # this means that there will be some replicates of detection data in the final output
                    # which needs to be tidied before analysis.   

                timestep += 1

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
                fig_save_dir = os.path.join('figs', 'Random_sim_spidermonkey-n_samplers-{}'.format(N_SAMPLERS))
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
                    for dt, site, visit in zip(traj_dts, sampler['trajectory_sites'], sampler['charger_visits']):
                         if visit == 1:
                            plt.scatter(dt+1, site+1, color='red')
            
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
                plt.scatter(lons, lats)
                plt.scatter(loc_site_coords[charger_site_ix][1], loc_site_coords[charger_site_ix][0], color='red')
                for i, (lon, lat) in enumerate(zip(lons, lats)):
                    plt.text(lon, lat, str(i), fontsize=8)

                plt.savefig(os.path.join(fig_save_dir, 'sampler_trajectories_{}.png'.format(cluster)))
                plt.clf()
                plt.close(fig2)
            #end trajectory figs
        #end N_sampler loop
        all_exps.append(locs_sampled)
        locs_pd = pd.concat(loc_dets)
        exps_list.append(locs_pd)
    #end cluster loop
    exps_pd = pd.concat(exps_list)
    exps_pd['iteration'] = iter
    csv_save_dir = os.path.join('data', 'Iterations_proportional_set', 'RandomSM')
    if not os.path.exists(csv_save_dir): os.makedirs(csv_save_dir)
    exps_pd.to_csv(os.path.join(csv_save_dir, 'RandomSM_Iter{:02}.csv'.format(iter)), index=False)
#end iter loop