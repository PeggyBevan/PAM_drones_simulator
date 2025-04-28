#here, load all pickle files and combine into a csv

import os
import pickle
#from tqdm import tqdm 
import numpy as np 
import pandas as pd

MIN_DET_CONF = 0.8
VALID_SPECS = ['Ramphastos ambiguus', 'Volatinia jacarina', 'Brotogeris jugularis',
    'Pitangus sulphuratus', 'Piranga rubra', 'Tyrannus melancholicus',
    'Coereba flaveola', 'Megarynchus pitangua', 'Todirostrum cinereum',
    'Saltator maximus', 'Myiozetetes similis', 'Thraupis episcopus',
    'Melanerpes rubricapillus', 'Ramphocelus passerinii', 'Stilpnia larvata',
    'Leptotila verreauxi', 'Piaya cayana', 'Progne chalybea', 'Amazilia tzacatl']

dets_dir = 'drone-data-hpc'
det_fs = os.listdir(dets_dir)

save_dir = 'data'
if not os.path.exists(save_dir): os.makedirs(save_dir)

all_gt = []
for f_ix, f in enumerate(det_fs):
    det_path = os.path.join(dets_dir, f)
    with open(det_path, 'rb') as f_handle:
        loc_dets_mat, loc_sites, loc_dts = pickle.load(f_handle)
    #filter out times where recording wasn't happening 'not-nones'.
    nn_ixs = np.where(loc_dets_mat != None) #this is a matrix of sites and dates
    nn_site_ixs = nn_ixs[0]
    nn_hr_ixs = nn_ixs[1]
    loc_name = f.split('.')[0]

    site_list = [loc_sites[i] for i in nn_site_ixs]
    date_list = [loc_dts[i] for i in nn_hr_ixs]

    loc_dets_dicts = []
    for site, date in zip(nn_site_ixs, nn_hr_ixs):
        loc_dets_dicts.append(loc_dets_mat[site, date])
    
    loc_dets_values = []
    for d in loc_dets_dicts:
        if len(d) == 0:
            loc_dets_values.append([])
        else:
            dets = [det['common_name'] for det in d if det['confidence'] > 0.8 and det['scientific_name'] in VALID_SPECS]
            loc_dets_values.append(dets)

    loc_df = pd.DataFrame({'n_samplers': 'groundtruth', 'loc_name': loc_name, 'site_name': site_list, 'datetime': date_list, 'detections': loc_dets_values})
    all_gt.append(loc_df)

gt_df = pd.concat(all_gt)
gt_df.to_csv(os.path.join(save_dir, 'groundtruth_detections_80thresh_PB.csv'), index=False)


