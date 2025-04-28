
#configs for drone settings
SAMPLER_SPEED = 15          # cruise speed in meters per second
SAMPLER_BATTERY_S = 60*20   # minutes in seconds - limit to 20 minutes
SKIP_LOCS = ['Luna'] #skipping luna due to ambiguous site names
VALID_SPECS = ['Ramphastos ambiguus', 'Volatinia jacarina', 'Brotogeris jugularis',
    'Pitangus sulphuratus', 'Piranga rubra', 'Tyrannus melancholicus',
    'Coereba flaveola', 'Megarynchus pitangua', 'Todirostrum cinereum',
    'Saltator maximus', 'Myiozetetes similis', 'Thraupis episcopus',
    'Melanerpes rubricapillus', 'Ramphocelus passerinii', 'Stilpnia larvata',
    'Leptotila verreauxi', 'Piaya cayana', 'Progne chalybea', 'Amazilia tzacatl']
SAVE_TRAJECTORY_FIGS = False
SAVE_SPECIESSAMPLED_FIGS = False
comb_dets_dir = 'drone-data-hpc'
TAKEOFF_PEN = 26.68 #seconds


#takeoff speed is half of the cruise speed - 7.5m/s
#imaging we cruise at 50m height, it will take 50/7.5 = 6.67 seconds to reach 50m
#and it takes twice as much energy because fighting gravity.
#propellors are going twice as hard for twice as long. 
#so takeoff penalty = 6.67 * 2 = 13.34s
#so takeoff penalty = 13.34*2 = 26.68s

