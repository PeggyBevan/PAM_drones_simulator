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
for iter in range(1, 51):
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
        #but is cropped to start and end of survey at this site.
        nn_dets_mat = loc_dets_df.values
        coordinates = [coordinates[i] for i in nn_sites]
        dist_mat = calculate_distance_matrix(coordinates)
        
        #set charger site as the one closest to the middle
        charger_site_ix = find_charger_site(coordinates)
        max_samplers = len(loc_sites)
        prev_N_SAMPLERS = 0
        for n_samp_ix, prop_samplers in enumerate([0.2, 0.4, 0.6, 0.8, 1.0]):
            
            N_SAMPLERS = int(max_samplers * prop_samplers)
            if N_SAMPLERS == 0: 
                N_SAMPLERS = 1

            #there will be some cases where N_SAmplers is the same value for dfferent prop_samplers
            #i am still re-running the simulation for these cases, for ease of data-analysis, but need to be conscious of this in data analysis.
            

            #create dictionary to save sampler trajectories
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
            #visited = np.zeros(dist_mat.shape[0], dtype=np.int32)   
            while timestep < len(nn_dts)-1:
                curr_hr_ix = nn_dts[timestep] #this corresponds to the original date which can be
                #fetched from the loc_dts list
                
                #in the first time step of the simulation, the drone will
                #start at a site which was recording in the real dataset    
                chosen_next_sites = []
                for smplr in all_samplers:
                    #is the charger site recording in the first hour?
                    #in the first time step of the simulation, the drone will
                    #start at a site which was recording in the real dataset    
                    if timestep == 0:
                        if nn_dets_mat[charger_site_ix, 0] == None: #if charging site is not recording
                            smplr['curr_site_ix'] = np.where(nn_dets_mat[:, 0] != None)[0][0] #choose the first site available
                            #otherwise, stay at charger site.
                    # Keep track of where each sampler has been
                    smplr['trajectory_sites'].append(smplr['curr_site_ix'])  # adds the current site (index) to trajectory
                    smplr['trajectory_dts'].append(loc_dts[curr_hr_ix])  # adds the actual date to trajectory
                    smplr['trajectory_dts_ix'].append(curr_hr_ix) #adds date index to trajectory
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
                    
                    # Figure out where we can go next and remove places where other samplers are already going
                    next_site_opts = np.where((nn_dets_mat[:, timestep+1] != None))[0]
                    #return sites what are not already chosen by another sampler
                    next_site_opts = list(filter(lambda a: a not in chosen_next_sites, next_site_opts))
                    if smplr['curr_site_ix'] == smplr['charger_site_ix']:
                            smplr['battery_s'] = SAMPLER_BATTERY_S
                        
                    # If there are valid options, then choose a random next site
                    if len(next_site_opts) > 0:
                        next_site_ix = np.random.choice(next_site_opts, 1)[0]
                        curr_site_ix = smplr['curr_site_ix']
                        
                        # Determine if we can make it to the next site and still back to the charger
                        trave_to_next_site_s = (dist_mat[curr_site_ix, next_site_ix] / SAMPLER_SPEED) + TAKEOFF_PEN
                        trave_to_charger_from_next_site_s = (dist_mat[next_site_ix, charger_site_ix] / SAMPLER_SPEED) + TAKEOFF_PEN
                        
                        # Move the sampler to the next site and decrease battery accordingly
                        if smplr['battery_s'] > trave_to_next_site_s + trave_to_charger_from_next_site_s:
                            smplr['battery_s'] -= trave_to_next_site_s
                            smplr['curr_site_ix'] = next_site_ix #move to next site
                            smplr['charger_visits'].append(0)
                            smplr['total_dist'] += dist_mat[curr_site_ix, next_site_ix]
                        
                        else:
                            #trave_to_charger_s = (dist_mat[curr_site_ix, charger_site_ix] / SAMPLER_SPEED) + TAKEOFF_PEN
                            #smplr['battery_s'] -= trave_to_charger_s
                            smplr['curr_site_ix'] = smplr['charger_site_ix'] #send drone back to charger.
                            smplr['charger_visits'].append(1) #mark this visit as a charging visit
                            smplr['total_dist'] += dist_mat[curr_site_ix, smplr['charger_site_ix']]

                        chosen_next_sites.append(smplr['curr_site_ix'])
                    else:
                        smplr['charger_visits'].append(0)
                #if there are no available sites, the drone will remain in position.
                # this means that there will be some replicates of detection data in the final output
                # which needs to be tidied before analysis.   
                        
                timestep +=1
            #end while loop

            #Save sampling data for this sampler group & change to pd format
            samp_list = [] #turn each sampler trajectory for this cluster/sampler combo into a dataframe
            for samp in all_samplers:
                dets_df = pd.DataFrame(samp)
                site_list = [loc_sites[i] for i in samp['trajectory_sites']] #convert indexes to site names
                dets_df['site_name'] = site_list
                #dets_df['n_samplers'] = N_SAMPLERS
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
            #need to work out how to do this without nn_ixs
            #all_nn_dets = loc_dets_mat[nn_ixs]
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
                fig_save_dir = os.path.join('figs', 'Random_sim-n_samplers-{}'.format(N_SAMPLERS))
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
                    # Determine the color based on 'charger_visit'
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
                plt.scatter(lons, lats)
                plt.scatter(coordinates[charger_site_ix][1], coordinates[charger_site_ix][0], color='red')
                for i, (lon, lat) in enumerate(zip(lons, lats)):
                    plt.text(lon, lat, str(i), fontsize=8)

                plt.savefig(os.path.join(fig_save_dir, 'sampler_trajectories_{}.png'.format(comb_d_f.split('.')[0])))
                plt.clf()
            #end trajectory figs
        #end n_sampler loop
        
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

    #end CLUSTER loop
    if SAVE_SPECIESSAMPLED_FIGS: #this will no longer work for proportional sampling
        plt.suptitle('Species occurrence detectability')
        plt.tight_layout()
        plt.savefig(os.path.join('figs', 'Random_specs_samplability.png'))
    
    exps_pd = pd.concat(exps_list)
    exps_pd['iteration'] = iter
    csv_save_dir = os.path.join('data', 'Iterations_proportional_set', 'Random')
    if not os.path.exists(csv_save_dir): os.makedirs(csv_save_dir)
    exps_pd.to_csv(os.path.join(csv_save_dir, 'Random_Iter{:02}.csv'.format(iter)), index=False)
#end iter loop