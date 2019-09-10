# Simulates for each sigma on every map file provided - as long as this completes in 72 hours should be fine!

# We want to iterate over every file in our file list
# Parrellise by running having separate job for each sigma (6 different sigma)

import os
import sys
import math
import logging
import csv
import sqlite3

import shutil

from PyCoalescence import Map

# map_dir = "/Volumes/Seagate 3TB/Data/Maps/FragmentMaps/Random_select"
job_num = int(sys.argv[1])
# tmpdir = "/Volumes/Seagate 3TB/Results/FragmentedLandscapes/Dispersal/disp{}/".format(job_num+10)
# if not os.path.exists(tmpdir):
# 	os.mkdir(tmpdir)
# out_csv = "/Volumes/Seagate 3TB/Results/FragmentedLandscapes/Dispersal/dispersal_{}.csv".format(job_num + 10)

map_dir = "/work/set114/Panama/Data/FragmentMaps/Random/"
tmpdir = os.environ['TMPDIR']
file_list = os.listdir(map_dir)
sigma_list = [1, 2, 4, 8, 16, 32]
sigma_ref = job_num % 6
sigma = sigma_list[sigma_ref]
# Okay, let's run a simulation batch on each node, so each node only runs once for each file, but runs both dispersal
# and coalescence simulations
job_ref = job_num*len(file_list)
seed = 0
# tmpdir = os.environ['TMPDIR']
csv_name = os.path.join(tmpdir, "dispersal_{}.csv".format(job_num))
out_csv = "/work/set114/Panama/Results/FragmentedLandscapes/dispersal_random_{}.csv".format(job_num)
with open(csv_name, "w") as csvfile:
	csvwriter = csv.writer(csvfile, delimiter=' ', quotechar='|', quoting=csv.QUOTE_MINIMAL)
	csvwriter.writerow(["file", "sigma", "seed", "mean_dispersal", "mean_distance_travelled",
						"mean_distance_travelled_100", "mean_distance_travelled_1000",
						"mean_distance_travelled_10000", #"mean_distance_travelled_100000",
						"stdev_dispersal", "stdev_distance_travelled",
						"stdev_distance_travelled_100", "stdev_distance_travelled_1000"])
	for file in file_list:
		print(file)
		if ".tif" in file and ".aux" not in file:
			map_file = os.path.join(map_dir, file)
			# # Now run a dispersal simulation
			out_db = os.path.join(tmpdir, "dispersal_{}_{}.db".format(seed, sigma))
			m = Map(dispersal_db=out_db)
			if not os.path.exists(out_db):
				try:
					m.test_mean_dispersal(number_repeats=100000, output_database=out_db, map_file=map_file, seed=seed,
										  dispersal_method="normal", sigma=sigma, landscape_type="tiled")
					m.test_mean_distance_travelled(number_repeats=1000, number_steps=10, output_database=out_db,
												   map_file=map_file, seed=seed, dispersal_method="normal",
												   sigma=sigma, landscape_type="tiled")
					m.test_mean_distance_travelled(number_repeats=1000, number_steps=100, output_database=out_db,
												   map_file=map_file, seed=seed, dispersal_method="normal",
												   sigma=sigma, landscape_type="tiled")
					m.test_mean_distance_travelled(number_repeats=1000, number_steps=1000, output_database=out_db,
												   map_file=map_file, seed=seed, dispersal_method="normal",
												   sigma=sigma, landscape_type="tiled"),
					m.test_mean_distance_travelled(number_repeats=1000, number_steps=10000, output_database=out_db,
												   map_file=map_file, seed=seed, dispersal_method="normal",
												   sigma=sigma, landscape_type="tiled")
					# m.test_mean_distance_travelled(number_repeats=1000, number_steps=100000, output_database=out_db,
					# 							   map_file=map_file, seed=seed, dispersal_method="normal",
					# 							   sigma=sigma, landscape_type="tiled")
				except Exception as e:
					logging.warning(str(e))
			csvwriter.writerow([file, sigma, seed,
								m.get_mean_dispersal(parameter_reference=1),
								m.get_mean_distance_travelled(parameter_reference=2),
								m.get_mean_distance_travelled(parameter_reference=3),
								m.get_mean_distance_travelled(parameter_reference=4),
								m.get_mean_distance_travelled(parameter_reference=5),
								# m.get_mean_distance_travelled(parameter_reference=6),
								m.get_stdev_dispersal(parameter_reference=1),
								m.get_stdev_distance_travelled(parameter_reference=2),
								m.get_stdev_distance_travelled(parameter_reference=3),
								m.get_stdev_distance_travelled(parameter_reference=4)])
			seed += 1
shutil.move(csv_name, out_csv)
