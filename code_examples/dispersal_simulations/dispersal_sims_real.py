# Simulates for each sigma on every map file provided - as long as this completes in 72 hours should be fine!

# We want to iterate over every file in our file list
# Parrellise by running having separate job for each sigma (6 different sigma)

import os
import sys
import math
import shutil
import logging
import csv
import sqlite3
from pycoalescence import DispersalSimulation
import time
# map_dir = "/Volumes/Seagate 3TB/Data/Maps/FragmentMaps/CentralAmericaExtractions_select"
# tmpdir = "/Volumes/Seagate 3TB/Results/FragmentedLandscapes/Dispersal/disp{}/".format(job_num)
# if not os.path.exists(tmpdir):
# 	os.mkdir(tmpdir)
# out_csv = "/Volumes/Seagate 3TB/Results/FragmentedLandscapes/Dispersal/dispersal_{}.csv".format(job_num)
# Simulates for each sigma on every map file provided - as long as this completes in 72 hours should be fine!

# We want to iterate over every file in our file list
# Parrellise by running having separate job for each sigma (6 different sigma)


map_dir = "/work/set114/Panama/Data/FragmentMaps/CentralAmericaExtractions/"
file_list = os.listdir(map_dir)
job_num = int(sys.argv[1])
sigma_list = [1, 2, 4, 8, 16, 32]
sigma_ref = job_num % 6
sigma = sigma_list[sigma_ref]
number_steps = [10**x for x in range(-10, -4, 1)]
# Okay, let's run a simulation batch on each node, so each node only runs once for each file, but runs both dispersal
# and coalescence simulations
job_ref = job_num*len(file_list)
seed = 0
tmpdir = os.environ['TMPDIR']
csv_name = os.path.join(tmpdir, "dispersal_{}.csv".format(job_num))
out_csv = "/work/set114/Results/FragmentedLandscapes/Dispersal/dispersal_real_{}.csv".format(job_num)
start = time.time()
with open(csv_name, "wb") as csvfile:
	csvwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
	csvwriter.writerow(["file", "sigma", "seed", "number_steps", "distance"])
	for file in file_list:
		print(file)
		if ".tif" in file and ".aux" not in file:
			map_file = os.path.join(map_dir, file)
			# # Now run a dispersal simulation
			out_db = os.path.join(tmpdir, "dispersal_{}_{}.db".format(seed, sigma))
			m = DispersalSimulation(dispersal_db=out_db)
			if not os.path.exists(out_db):
				try:
					m.set_map(map_file)
					m.set_simulation_parameters(number_repeats=1000, seed=1, landscape_type="tiled_fine", sigma=sigma,
												number_steps=number_steps)
					m.run_mean_distance_travelled()
					# m.test_mean_distance_travelled(number_repeats=1000, number_steps=100000, output_database=out_db,
					# 							   map_file=map_file, seed=seed, dispersal_method="normal",
					# 							   sigma=sigma, landscape_type="tiled")
				except Exception as e:
					logging.warning(str(e))
			for ref in m.get_database_references():
				csvwriter.writerow([file, sigma, seed, m.get_database_parameters()[ref]["number_steps"],
									m.get_mean_distance_travelled(parameter_reference=ref)])
			seed += 1
		current=time.time()
		if current - start > 60*60*23:
			break
shutil.move(csv_name, out_csv)

