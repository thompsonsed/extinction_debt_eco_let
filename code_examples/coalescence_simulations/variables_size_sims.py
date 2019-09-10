import sqlite3

import os
import pandas as pd
import sys
from pycoalescence import CoalescenceTree, Simulation

def generate_simulation_list(source_directory):
	"""
	Generates a list of all simulation objects
	:param source_directory: directory to read file types from
	:return:
	"""
	sim_list = []
	for sim_type in ["Contrived", "Real", "Random"]:
		search_path = os.path.join(source_directory, sim_type)
		file_list = [x for x in os.listdir(search_path) if ".tif" in x]
		total = 0
		tmp_list = []
		for sigma in [1, 2, 4, 8, 16, 32]:
			for seed in range(10):
				for file in file_list:
					if "10000" in file:
						continue
					if total < 10:
						tmp_list.append({"file" : os.path.join(search_path, file),
									 "sigma" : sigma,
									 "seed" : seed,
									 "type" : sim_type})
						total += 1
					else:
						total = 0
						sim_list.append(tmp_list)
						tmp_list = []
	return sim_list

if __name__ == "__main__":
	number = int(sys.argv[1])
	tmpdir = os.environ['TMPDIR']
	fragment_maps_dir = "/work/set114/Panama/Data/FragmentMaps/Select/"
	output_dir = "/work/set114/Panama/Results/FragmentedLandscapes/"
	this_simulation = generate_simulation_list(fragment_maps_dir)[number]
	output = []
	for s in this_simulation:
		file = s["file"]
		seed = s["seed"]
		sigma = s["sigma"]
		sim_type = s["type"]
		c = Simulation()
		c.set_simulation_parameters(seed=seed, job_type=number,
								output_directory=tmpdir,
								min_speciation_rate=0.0001, sigma=sigma, deme=1, sample_size=1.0,
								max_time=70 * 60 * 60, dispersal_method="normal", landscape_type="tiled_fine")
		c.set_speciation_rates([0.0001, 0.001, 0.01, 0.1])
		c.set_map_files(sample_file="null", fine_file=file)
		try:
			c.finalise_setup()
			c.run_coalescence()
		except (IOError, sqlite3.OperationalError, sqlite3.DatabaseError):
			os.remove(c.output_database)
			c.config_open = False
			c.finalise_setup()
			c.run_coalescence()
		t = CoalescenceTree(c)
		for ref in t.get_community_references():
			richness = t.get_species_richness(ref)
			output.append({"file" : file, "seed" : seed, "job_type" : number, "sigma" : sigma,
						   "speciation_rate" : t.get_community_parameters(ref)["speciation_rate"],
						   "richness" : richness, "type" : sim_type
						   })
	output_df = pd.DataFrame(output)
	output_df.to_csv(os.path.join(output_dir, "variable_size_{}.csv".format(number)))
