"""
Contains the functions for generating maps with the desired coverage from real landscapes
"""
import math
import os
from random import randint

import numpy as np
from osgeo import gdal
from pycoalescence import Map


def extract_landscape(desired_cover, size, number, output_dir, main_map=None):
	"""
	Searches across forest cover of the amazon to find regions within 1% cover of the desired cover
	proportion, then makes the forest cover exactly the desired percentage, extracts the map and saves
	it to the output directory.

	:param float desired_cover: the desired proportion forest cover
	:param int size: the size of the landscapes to output
	:param int number: the desired number of fragments to find
	:param str output_dir: the desired output location for files (stored as map_size_desiredcover_number.tif)
	:param Map main_map: if supplied, uses the Map object for data handling. This can save on re-loading the same file
						 from disk multiple times, providing major performance gains for large files or slow read
						 speeds.
	"""
	if main_map is None:
		main_map = Map(file="../../Data/Maps/FragmentMaps/CSA_extract_int.tif")
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	x_size, y_size = main_map.get_x_y()
	map_cover_list = []
	x = 0
	y = 0
	while y < y_size - size:
		while x < x_size - size:
			tmp_arr = main_map.get_cached_subset(x, y, size, size)
			tmp_sum = np.sum(np.rint(tmp_arr))
			if tmp_sum > size * size * (desired_cover - 0.01) and \
							tmp_sum < size * size * (desired_cover + 0.01):
				map_cover_list.append([x, y])
				x += size
			x += 1
			if len(map_cover_list) >= number:
				break
		x = 0
		y += size
	if len(map_cover_list) < number:
		max_num = len(map_cover_list)
		print("Not found desired number of maps: wanted {}, found {}".format(number,
																			 len(map_cover_list)))
	else:
		max_num = number
	# print(map_cover_list)
	for i in range(max_num):
		output_file = os.path.join(output_dir, "map_{}_{}_{}.tif".format(size, desired_cover, i))
		x, y = map_cover_list[i]
		this_arr = np.rint(main_map.get_cached_subset(x, y, size, size))
		# print("cached shape: {}".format(main_map._cached_array.shape))
		# print("x, y: {}, {}".format(x, y))
		# print("i: {}".format(i))
		# print(this_arr)
		# print("shape: {}".format(this_arr.shape))
		total = np.sum(np.rint(this_arr))
		desired_number = int(math.floor(size * size * desired_cover))
		diff = desired_number - total
		# print(desired_number)
		# print(diff)
		while diff > 0:
			randx = randint(1, size - 2)
			randy = randint(1, size - 2)
			# Find a fragment edge and add a pixel at a time
			if this_arr[randy, randx] == 0.0:
				if this_arr[randy + 1, randx] == 1.0 or this_arr[randy, randx + 1] == 1.0 or \
								this_arr[randy - 1, randx] == 1.0 or this_arr[randy, randx - 1] == 1.0:
					this_arr[randy, randx] = 1.0
					diff -= 1
		while diff < 0:
			# Random numbers are offsetted so we don't pick from the edges of the map!
			# Yay for an easy hack!
			randx = randint(1, size - 2)
			randy = randint(1, size - 2)
			# Find a fragment edge and remove a pixel at a time
			if this_arr[randy, randx] == 1.0:
				if this_arr[randy + 1, randx] == 0.0 or this_arr[randy, randx + 1] == 0.0 or \
								this_arr[randy - 1, randx] == 0.0 or this_arr[randy, randx - 1] == 0.0:
					this_arr[randy, randx] = 0.0
					diff += 1
		# check sum adds up
		if np.sum(np.rint(this_arr)) != desired_number:
			raise ValueError("Sum doesn't add up! {} != {}".format(np.sum(np.rint(this_arr)), desired_number))
		# Now save the file
		geotransform = (0, 1, 0, 0, 0, -1)
		# Stupid attempt to fix the no data issue for writing out numpy arrays
		for i in range(min(this_arr.shape)):
			if this_arr[i, i] == 0:
				this_arr[i, i] == 0
		output_raster = gdal.GetDriverByName('GTiff').Create(output_file,
															 size, size, 1,
															 gdal.GDT_Float32)
		output_raster.SetGeoTransform(geotransform)
		out_band = output_raster.GetRasterBand(1)
		out_band.WriteArray(this_arr)
		out_band.FlushCache()
		out_band.SetNoDataValue(-99)
		del this_arr
		del output_raster
		# break


def random_landscape(desired_cover, size, number, output_dir):
	"""
	Creates random landscapes of the prescribed size, with the correct amount of habitat, and saves them
	all to the output directory as random_size_desired_cover_number.tif
	:param desired_cover: the desired amount of habitat cover
	:param size: the desired size of the landscape
	:param number: the desired number of random maps to create
	:param output_dir: the desired output directory
	:return:
	"""
	for i in range(number):
		output_file = os.path.join(output_dir, "random_{}_{}_{}.tif".format(size, desired_cover, i))
		generate_random(size, output_file, desired_cover)
		# print("x, y: {}, {}".format(x, y))
		# arr = np.zeros([size, size])
		# total = 0
		# desired_number = int(desired_cover * size * size)
		# while total < desired_number:
		# 	randx = randint(0, size - 1)
		# 	randy = randint(0, size - 1)
		# 	if arr[randx, randy] == 0:
		# 		arr[randx, randy] = 1
		# 		total += 1
		#
		# # check sum adds up
		# if np.sum(arr) != desired_number:
		# 	raise ValueError("Sum doesn't add up!")
		# # Now save the file
		# geotransform = (0, 1, 0, 0, 0, -1)
		# # Stupid attempt to fix the no data issue for writing out numpy arrays
		# for i in range(min(arr.shape)):
		# 	if arr[i, i] == 0:
		# 		arr[i, i] == 0
		# output_file = os.path.join(output_dir, "random_{}_{}_{}.tif".format(size, desired_cover, i))
		# output_raster = gdal.GetDriverByName('GTiff').Create(output_file,
		# 													 size, size, 1,
		# 													 gdal.GDT_Float32)
		# output_raster.SetGeoTransform(geotransform)
		# out_band = output_raster.GetRasterBand(1)
		# out_band.WriteArray(arr)
		# out_band.FlushCache()
		# out_band.SetNoDataValue(-99)
		# del arr
		# del output_raster

def tile_map(input_map, output_map, number):
	"""
	Tiles the input_map the provided number of times in both the x and y directions to produce the
	output_map.
	:param input_map: the input map to read from
	:param output_map: the output map to generate
	:param number: the number of times to tiles the landscape
	"""
	map_in = Map(input_map)
	map_in.open()
	output_map = Map(output_map)
	# x_size = map_in.get_x_y()[0] * number
	# y_size = map_in.get_x_y()[1] * number
	output_map.data = np.tile(map_in.data, (number, number))
	output_map.create()


def generate_random(size, output_file, percentage_cover):
	"""
	Generates a random map of the specified size at the specified output location
	:param size: the x and y size of the map
	:param output_file: path to the output file (should be a tif file).
	:param percentage_cover: the percentage cover of the output landscape
	"""
	arr = np.zeros([size, size])
	total = 0
	while total < int(size * size * percentage_cover):
		randx = randint(0, size - 1)
		randy = randint(0, size - 1)
		if arr[randx, randy] == 0:
			arr[randx, randy] = 1
			total += 1

	# check sum adds up
	if np.sum(arr) != int(size * size * percentage_cover):
		raise ValueError("Sum doesn't add up! {} != {}".format(np.sum(arr),
															   size * size * percentage_cover))
	# Now save the file
	geotransform = (0, 1, 0, 0, 0, -1)
	# Stupid attempt to fix the no data issue for writing out numpy arrays
	for i in range(min(arr.shape)):
		if arr[i, i] == 0:
			arr[i, i] == 0
	output_raster = gdal.GetDriverByName('GTiff').Create(output_file,
														 size, size, 1,
														 gdal.GDT_Float32)
	output_raster.SetGeoTransform(geotransform)
	out_band = output_raster.GetRasterBand(1)
	out_band.WriteArray(arr)
	out_band.FlushCache()
	out_band.SetNoDataValue(-99)
	del arr
	del output_raster

def tile_from_directory(number, file_list, output_file):
	"""
	Tiles tif files from the given directory sequentially the provided number of times in the x and y dimension to
	produce the output file.

	:param number: the number of times to tile in the x and y dimension
	:param file_list: list of file paths to combine
	:param output_file: the path to the output file to be created
	"""
	if os.path.exists(output_file):
		raise IOError("File already exists at {}".format(output_file))
	output_map = Map()
	if len(file_list) == 0:
		raise IOError("Input list does not contain any files.")
	first_map = Map(file_list[0])
	size_x = first_map.get_x_y()[0]
	size_y = first_map.get_x_y()[1]
	output_map.data = np.zeros(shape=(size_y * number, size_x * number))
	start_x = 0
	start_y = 0
	count = 0
	for file in file_list:
		if count > number**2:
			break
		m = Map(file)
		m.open()
		output_map.data[start_y: start_y+size_y, start_x: start_x + size_x] = m.data
		start_x += size_x
		if start_x + size_x > size_x * number:
			start_y += size_y
			start_x = 0
		if start_y + size_y > size_y * number:
			break
		count += 1
	output_map.create(output_file)
