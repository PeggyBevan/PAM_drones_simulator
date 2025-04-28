import numpy as np
import os
import csv
import pickle
import geopy.distance
from sklearn.cluster import KMeans
from statistics import mean
from math import sqrt


comb_dets_dir = 'drone-data-hpc'
site_info_csv_f = 'data/costa_rica_site_info.csv'
save_dir = 'combined_data'

def get_sites_with_gps(verbose=False):
    all_sites = []
    comb_dets_dir = 'drone-data-hpc'
    for comb_d_f in os.listdir(comb_dets_dir):

        with open(os.path.join(comb_dets_dir, comb_d_f), 'rb') as f_handle:
            _, loc_sites, _ = pickle.load(f_handle)

        for loc_s in loc_sites:
            site = dict()
            site['loc_name'] = comb_d_f
            site['site_name'] = loc_s

            if len(loc_s.split(' ')) > 2: 
                sn = loc_s.split(' ')[-2]
            else:
                sn = loc_s.split(' ')[-1]
            sn = sn.split('-')[0]
            site['site_num'] = sn

            site['loc_site_str'] = '{}{}'.format(comb_d_f.split('.')[0].replace('_', ' '), site['site_num']).lower().replace(' ', '')

            all_sites.append(site)

    all_site_names = np.asarray([s['loc_site_str'] for s in all_sites])

    with open(site_info_csv_f, newline='') as csvfile:
        rdr = csv.reader(csvfile, delimiter=',')
        for row_ix, row in enumerate(rdr):
            if row_ix != 0:
                _loc_and_site = row[1]
                _loc = _loc_and_site.split('_')[0].lower()
                if _loc == 'sq260': _loc = 'sq2601'
                if _loc == 'mangrove': _loc = 'mangroves'
                if _loc == 'palma': _loc = 'lapalma'
                if _loc == 'elsi': _loc = 'elsicroc'
                if _loc == 'rancho': _loc = 'ranchobajo'
                if _loc == 'miramar': _loc = 'mirenmar'
                if _loc == 'nuevo': _loc = 'rionuevo'
                if _loc == 'gamba': _loc = 'lagamba'
                if _loc == 'tarde': _loc = 'latarde'
                if _loc == 'lareserva': _loc = 'indigenousreserve'
                if _loc == 'sendero': _loc = 'golfito'

                _site_num = _loc_and_site.split('_')[-1]
                _month = row[6]
                _elev = row[7]
                _lat = row[8]
                _long = row[9]

                match_site_ix = np.where((all_site_names == '{}{}'.format(_loc, _site_num).lower().replace(' ', '')))[0]
                if len(match_site_ix) > 0:
                    all_sites[match_site_ix[0]]['lat'] = float(_lat)
                    all_sites[match_site_ix[0]]['long'] = float(_long)
                    all_sites[match_site_ix[0]]['elev'] = float(_elev)
                    all_sites[match_site_ix[0]]['month'] = _month           
                    
    sites = []
    for s in all_sites:
        if 'lat' in s.keys(): 
            sites.append(s)
        else:
            if verbose:
                print('No GPS coords for {}'.format(s['loc_site_str']))

    if verbose:
        print('{} folders with audio, found GPS for {} sites'.format(len(all_site_names), len(sites)))
    
    return sites 

def read_site_file(comb_d_f, all_site_infos, all_site_info_names):
    """
    Read coordinates from a pickle file.
    Assumes the data is in columns labeled 'X', 'Y', and 'Demand'.
    """
    comb_d_path = os.path.join(comb_dets_dir, comb_d_f)

    with open(comb_d_path, 'rb') as f_handle:
            loc_dets_mat, loc_sites, loc_dts = pickle.load(f_handle)
    
    # Get GPS coordinates for each of the sites
    loc_site_coords = []
    keep_site_ixs = []
    
    for loc_s_ix, loc_s in enumerate(loc_sites):
        match_ixs = np.where((all_site_info_names == '{} {}'.format(comb_d_f, loc_s)))[0]
        
        if len(match_ixs) == 0:
            print('No GPS for {}: {} - excluding from simulations'.format(comb_d_f, loc_s))
            continue
        
        site_info = all_site_infos[match_ixs[0]] #extract site info for this site
        site_coords = (site_info['lat'], site_info['long'])
        loc_site_coords.append(site_coords) #create list of sites and coords
        keep_site_ixs.append(loc_s_ix) 
    
    # Only keep analysing sites we have GPS coordinates for
    loc_dets_mat = loc_dets_mat[keep_site_ixs, :]
    loc_sites = np.asarray(loc_sites)[keep_site_ixs]

    
    coordinates = loc_site_coords
    #demands = 1
    return coordinates, loc_dets_mat, loc_sites, loc_dts


def calculate_distance_matrix(coordinates):
    """
    Calculate the distance matrix between coordinates.
    """
    num_points = len(coordinates)
    dist_matrix = np.zeros((num_points, num_points))

    for i in range(num_points):
        for j in range(num_points):
            dist_matrix[i, j] = calculate_distance(coordinates, i, j)

    return dist_matrix

#change this to geopy distance?
def calculate_distance(coordinates, i, j, SAMPLER_SPEED = 15):
    """
    Calculate the Euclidean distance between two points & return in number of seconds
    """
    
    distance = (geopy.distance.distance(coordinates[i], coordinates[j]).m) / SAMPLER_SPEED
    return distance

def calculate_total_distance(route, dist_matrix):
    """
    Calculate the total distance of a given route using the distance matrix.
    """
    total_distance = 0
    num_points = len(route)

    for i in range(num_points - 1):
        current_node = route[i]
        next_node = route[i + 1]
        total_distance += dist_matrix[current_node, next_node]

    return total_distance

def find_nearest_neighbor(next_site_opts, dist_mat, curr_site_ix, curr_battery, charger_site_ix = 0, SAMPLER_SPEED = 15, TAKEOFF_PEN = 16):
    """
    Given a list of sites and distance matrix, find the nearest neighbor.
    """
    min_dist = float('inf')
    nearest = None
    for neighbor in next_site_opts:
        trave_to_neighbor = (dist_mat[curr_site_ix, neighbor] / SAMPLER_SPEED) + TAKEOFF_PEN
        trave_home_from_neighbor = (dist_mat[neighbor, charger_site_ix] / SAMPLER_SPEED) + TAKEOFF_PEN
        # if travel to neighbor + travel to charger from neighbor is less than capacity
        if trave_to_neighbor + trave_home_from_neighbor <= curr_battery and dist_mat[curr_site_ix, neighbor] < min_dist:
            # if the next stop is within capacity, and the distance is less than the current min
            # save as nearest
            nearest = neighbor
            min_dist = dist_mat[curr_site_ix, neighbor]
    return nearest, trave_to_neighbor
                        
def cluster_sites(coordinates, N_SAMPLERS) :
    """
    takes a list of coordinates and separates into n clusters.
    """
    clusters = KMeans(n_clusters=N_SAMPLERS).fit(coordinates)
    labels = clusters.labels_
    return labels

def find_charger_site(coordinates):
    """"
    Select the site closest to the middle of all sites to be charger
    """
    mean_x = np.mean([coord[0] for coord in coordinates])
    mean_y = np.mean([coord[1] for coord in coordinates])
    distances = [sqrt((x - mean_x)**2 + (y - mean_y)**2) for x, y in coordinates]
    closest_index = distances.index(min(distances))
    return closest_index

def is_det_valid(det, filter_specs):
    return det['confidence'] > 0.8 and det['scientific_name'] in filter_specs

def find_curr_dets(curr_site, timestep, dets_mat):
    """
    Find the current detections at a site & timestep
    Requires current site, timestep, and detection matrix
    """
    curr_dets = dets_mat[curr_site, timestep]  # sample detections at current site index and current hour
    # if there are detections, extract common name and add them to the list of sampled species
    if curr_dets != None:
        det_list = []
        if len(curr_dets) > 0:
            for d in curr_dets:
                if is_det_valid(d, VALID_SPECS):
                    det_list.append(d['common_name'])
                else:
                    det_list = []
        else: det_list = []
                        

if __name__ == '__main__':
    get_sites_with_gps(verbose=True)


