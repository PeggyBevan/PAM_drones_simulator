# PAM_drones_simulator
Simulating passive acoustic monitoring using autonomous UAVS
This project uses an existing, labelled acousic dataset from costa rica to simulate surveys of autonomous UAVs.
The UAV simulations are run in python, and statistical analysis and visualisation is run in R.

Data
-
The scripts for bird analyses expect a series of pickle files in a folder called 'drone-data-hpc'. These contain arrays of birdnet detections for each site for each cluster by the hour.
The scripts for spider monkey analyses expect a csv in the 'data' folder called 'SpiderMonkeys_Processeddata.csv'


Python requirements
This code was run in Python V3.11.9
Package requirements:
OS
Pandas
Numpy
datetime
csv
pickle
geopy
sklearn
statistics
math

Supporting files
-
utils.py: functions to support the main scripts

configs.py: simulation configurations that are the same across all simulations.

**Simulation files (Python)**
-

**random_sim.py, routeplanning_sim.py, adaptive_sim.py**
These are all scripts to run the simulation of UAV-based bird surveys. They take approximately 6 minutes to run on a macbook pro with Intel i7 chip and 16GB RAM. 
Each script uses a different sampling strategy 
Each script will run simulations for 5 values of sampling intensity (0.2,0.4,0.6,0.8 & 1) and run 50 iterations of these.
The adaptive_sim script takes 3x longer as it runs the simulations with three values of weighting factor.

**Random_sim_spidermonkey.p, routeplanning_sim_spidermonkey.py**
These are all scripts to run the simulation of UAV-based spider monkey surveys. They take approximately 2 minutes to run on a macbook pro with Intel i7 chip and 16GB RAM.
As with the bird surveys, each script uses a different sampling strategy, and runs simulations with 5 values of sampling intensity over 50 iterations.

**Analysis files (R)**
- 
The following packages are required
vegan
ggplot2
tidyr
dplyr
stringr

**BirdNetCSVProcessing.R**
This script takes all detections from all UAV-based bird survey simulations and complies into one script. 
Takes approximately 10 minutes to run on a macbook pro with Intel i7 chip and 16GB RAM.

**plot_sim_results.R**
Plot figures 2 and 4 - showing output of UAV-based bird surveys.














