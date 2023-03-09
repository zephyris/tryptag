# user setting: verbose output from accessing data
print_status = True

# user setting: path in which to cache image data (up to ~8Tb without zip file, ~16Tb with)
data_cache_path = "_tryptag_cache"

# user setting: remove zip files after download (~halves data cache size)
remove_zip_files = False

# variables which will contain tryptag data
gene_list = None
zenodo_index = None

# image properties
um_per_px = 6.5 / 63 # physical pixel size / corrected magnification

# function to fetch zenodo record information, respecting request rate limits
def _fetch_zenodo_record_json(zenodo_id):
	# fetch Zenodo record JSON
	from urllib.request import urlopen
	from urllib.error import HTTPError
	import json
	import time
	url = "https://zenodo.org/api/records/"+str(zenodo_id)
	# Zenodo queries are rate limited, so request with that in mind
	response = None
	while response is None:
		try:
			response = urlopen(url)
		except HTTPError as e:
			if print_status: print("	Zenodo rate limit reached, waiting to retry")
			time.sleep(60) # testing shows the rate limiter resets after 50s, so 60s for a bit of space
	zenodo_json = json.loads(response.read().decode(response.info().get_param('charset') or 'utf-8-sig'))
	zenodo_record_id = str(zenodo_json["id"])
	return zenodo_json

# fetch gene list/metadata, recording information in gene_list and zenodo_index
def fetch_gene_list():
	global gene_list
	global zenodo_index
	# if the global variable is None then not loaded yet
	if gene_list is None:
		import urllib.request
		# fetch Zenodo record JSON, to get latest version doi
		zenodo_json = _fetch_zenodo_record_json(6862289)
		zenodo_record_id = str(zenodo_json["id"])
		if print_status: print("	fetching gene list from latest zenodo version id: "+zenodo_record_id)
		# load zenodo record index
		if zenodo_index is None:
			if print_status: print("	fetching zenodo record id index")
			# load DOI index
			zenodo_index = {}
			url = "https://zenodo.org/record/"+zenodo_record_id+"/files/plate_doi_index.tsv?download=1"
			response = urllib.request.urlopen(url)
			for line in response:
				line = line.decode(response.info().get_param('charset') or 'utf-8-sig').splitlines()[0].split("\t")
				doi_data = {}
				doi_data["master_record_id"] = line[0].split(".")[-1]
				zenodo_index[line[1]] = doi_data
		# load localisations table
		gene_list = {}
		if print_status: print("	fetching gene data table")
		url = "https://zenodo.org/record/"+zenodo_record_id+"/files/localisations.tsv?download=1"
		response = urllib.request.urlopen(url)
		for line in response:
			line = line.decode(response.info().get_param('charset') or 'utf-8-sig').splitlines()[0].split("\t")
			if line[0] != "Gene ID":
				# if not the header line, grab gene data
				gene_list[line[0]] = {}
				termini = ["c", "n"]
				offsets = [5, 16]
				for t in range(len(termini)):
					if line[offsets[t] + 0] == "cell line generated":
						terminus_data = {
							"plate": line[offsets[t] + 1].split(" ")[0],
							"well": line[offsets[t] + 1].split(" ")[1],
							"loc": line[offsets[t] + 8],
							"signl_low": line[offsets[t] + 9],
							"signal_background": line[offsets[t] + 10]
						}
						if terminus_data["plate"] == "V1115_210180807":
							terminus_data["plate"] = "V1115_20180807"
						terminus_data["zenodo_id"] = zenodo_index[terminus_data["plate"]]["master_record_id"]
						gene_list[line[0]][termini[t]] = terminus_data

# verbose progress bar for file download
_progress_bar = None
def _show_progress_bar(block_num, block_size, total_size):
	import progressbar
	global _progress_bar
	if _progress_bar is None:
		_progress_bar = progressbar.ProgressBar(maxval=total_size)
		_progress_bar.start()
	downloaded = block_num * block_size
	if downloaded < total_size:
		_progress_bar.update(downloaded)
	else:
		_progress_bar.finish()
		_progress_bar = None

# md5 hash function for checking zip integrity
def _file_md5_hash(path, blocksize = 2**20):
	import hashlib
	m = hashlib.md5()
	with open(path, "rb") as file:
		while True:
			buffer = file.read(blocksize)
			if not buffer:
				break
			m.update(buffer)
	return m.hexdigest()

# get microscopy data for a given gene_id and terminus
# updates gene_list and zenodo_index
# places data in data_cache_path
def fetch_data(gene_id, terminus):
	global data_cache_path
	terminus = terminus.lower()
	# load tryptag data, if not already
	fetch_gene_list()
	if terminus in gene_list[gene_id]:
		import os
		# check if the data cache directory exists, and make if not
		if not os.path.isdir(data_cache_path):
			if print_status: print("	making data cache directory")
			os.mkdir(data_cache_path)
		# target paths for zip file and data subdirectory
		plate = gene_list[gene_id][terminus]["plate"]
		zip_path = os.path.join(data_cache_path, plate+".zip")
		dir_path = os.path.join(data_cache_path, plate)
		if not os.path.isfile(zip_path) and not os.path.isfile(os.path.join(dir_path, "_"+plate+".zip.md5")):
			# fetch the processed microscopy data from Zenodo
			if print_status: print("	making plate data directory for: "+plate)
			import urllib.request
			import shutil
			if "record_id" not in zenodo_index[plate]:
				# if not already translated, fetch the latest from master Zenodo ID
				if print_status: print("	remapping zenodo id to latest version for id: "+gene_list[gene_id][terminus]["zenodo_id"])
				import json
				# fetch Zenodo record JSON, to get latest version doi
				zenodo_json = _fetch_zenodo_record_json(gene_list[gene_id][terminus]["zenodo_id"])
				zenodo_index[plate]["record_doi"] = zenodo_json["doi"]
				zenodo_index[plate]["record_id"] = str(zenodo_json["id"])
				for file in zenodo_json["files"]:
					if file["key"].endswith("_processed.zip"):
						zenodo_index[plate]["record_url"] = file["links"]["self"]
						zenodo_index[plate]["record_md5"] = file["checksum"].split(":")[-1]
			# download data
			zip_md5 = 0
			zip_path_temp = os.path.join(data_cache_path, "temp.zip")
			while zip_md5 != zenodo_index[plate]["record_md5"]:
				if print_status:
					print("	downloading data: "+plate+".zip")
					urllib.request.urlretrieve(zenodo_index[plate]["record_url"], zip_path_temp, _show_progress_bar)
				else:
					urllib.request.urlretrieve(zenodo_index[plate]["record_url"], zip_path_temp)
				if print_status: print("	checking md5 of zip file")
				zip_md5 = _file_md5_hash(zip_path_temp)
				if zip_md5 != zenodo_index[plate]["record_md5"]:
					if print_status: print("	md5 of downloaded zip is incorrect ("+zip_md5+"), retrying download")
			shutil.move(zip_path_temp, zip_path)
			with open(zip_path+".md5", "w") as file: file.write(zip_md5)
		if not os.path.isfile(os.path.join(dir_path, "_"+plate+".zip.md5")):
			# unzip data
			if print_status: print("	decompressing data: "+plate+".zip")
			import zipfile
			import shutil
			try:
				with zipfile.ZipFile(zip_path) as archive:
					count = 0
					for file in archive.namelist():
						# loop through all files, finding files ending with the cell roi suffix
						suffix = "_roisCells.txt"
						if file.endswith(suffix):
							count += 1
							source_path = os.path.split(file)
							# infer main tif and thresholded tif image filenames
							file_roi = file
							file_tif = file[:-len(suffix)]+".tif"
							file_thr = file[:-len(suffix)]+"_thr.tif"
							if file_roi in archive.namelist() and file_tif in archive.namelist() and file_thr in archive.namelist():
								# if all exist, decompress and move to plate directory
								archive.extract(file_roi, data_cache_path)
								archive.extract(file_tif, data_cache_path)
								archive.extract(file_thr, data_cache_path)
								target_path = os.path.join(data_cache_path, plate)
								shutil.move(os.path.join(data_cache_path, file_roi), os.path.join(target_path, os.path.split(file_roi)[-1]))
								shutil.move(os.path.join(data_cache_path, file_tif), os.path.join(target_path, os.path.split(file_tif)[-1]))
								shutil.move(os.path.join(data_cache_path, file_thr), os.path.join(target_path, os.path.split(file_thr)[-1]))
							else:
								print("")
								print("== missing file for "+os.path.split(file_tif)[-1]+" ==")
							if print_status: print(".", end = "", flush = True)
				os.rmdir(os.path.join(data_cache_path, plate, "microscopeImagesAutoprocessed"));
				# copy source zip md5 to the data directory, also marks decompression as complete
				shutil.copyfile(zip_path+".md5", os.path.join(dir_path, "_"+plate+".zip.md5"))
				if remove_zip_files:
					os.remove(zip_path)
					os.remove(zip_path+".md5")
				if print_status:
					print("")
					print("	decompressed "+str(count)+" fields of view")
			except BadZipfile:
				print ("== invalid zip file: "+plate+".zip ==")
				if remove_zip_files: os.remove(zip_path)
		# count fields of view and number of cells
		base_path = os.path.join(data_cache_path, plate)
		if "fields_count" not in gene_list[gene_id][terminus] and os.path.isdir(base_path):
			if print_status: print("	counting image data files for: "+gene_id+" "+terminus)
			import glob
			cells = []
			for i in range(len(glob.glob(os.path.join(base_path, gene_id+"_4_"+terminus.upper()+"_*_roisCells.txt")))):
				with open(os.path.join(base_path, gene_id+"_4_"+terminus.upper()+"_"+str(i + 1)+"_roisCells.txt")) as cells_file:
					lines = cells_file.readlines()
					cells.append([])
					for line in lines:
						line = line.split("\t")
						if line[0] != "cell":
							cell_data = {}
							cell_data["wand"] = (int(line[1]), int(line[2]))
							cell_data["centre"] = (float(line[3]), float(line[4]))
							cell_data["angle"] = float(line[7])
							cells[-1].append(cell_data)
			gene_list[gene_id][terminus]["fields_count"] = len(cells)
			gene_list[gene_id][terminus]["cells_per_field"] = [len(x) for x in cells]
			gene_list[gene_id][terminus]["cells"] = cells.copy()

# checks data in data_cache_path to make sure it is the latest version
def check_data():
	import os
	# load tryptag data, if not already
	fetch_gene_list()
	for plate in zenodo_index:
		zip_path = os.path.join(data_cache_path, plate+".zip")
		dir_path = os.path.join(data_cache_path, plate)
		if os.path.isdir(dir_path) or os.path.isfile(zip_path):
			zenodo_json = _fetch_zenodo_record_json(zenodo_index[plate]["master_record_id"])
			for file in zenodo_json["files"]:
				if file["key"].endswith("_processed.zip"):
					zenodo_index[plate]["record_md5"] = file["checksum"].split(":")[-1]
		if os.path.isdir(dir_path):
			md5_file = os.path.join(dir_path, "_"+plate+".zip.md5")
			if not os.path.isfile(md5_file):
				if print_status: print("	md5 file mising from image directory")
			else:
				with open(md5_file, "r") as file:
					md5 = file.read()
					if md5 != zenodo_index[plate]["record_md5"]:
						if print_status: print("	md5 file in directory does not match zenodo record id for zip")
		if os.path.isfile(zip_path):
			if remove_zip_files:
				if print_status: print("	zip found which should have been removed")
			md5_file = os.path.join(zip_path+".md5")
			with open(md5_file, "r") as file:
				md5 = file.read()
				if md5 != zenodo_index[plate]["record_md5"]:
					if print_status: print("	md5 file for zip file does not match zenodo record id for zip")

# force load of all data by iterating through every gene and terminus
def fetch_all_data():
	fetch_gene_list()
	for gene in gene_list:
		for terminus in ["c", "n"]:
			fetch_data(gene, terminus)

# global variables for caching last field of view loaded
_field_base_path_pil = None
_thresholds_pil = None
_channels_pil = None

# open field of view, returning array with one PIL image per channel and threshold image
# [phase_(gray) mng_(green) dna_(blue) phase_threshold dna_threshold] OR
# [pha, mng, dna, pth, dth]
# channel images are mode F, 32-bit float, threshold images are mode L, 8 bit (0 or 255)
def open_field_pillow(gene, terminus, field):
	from PIL import Image, ImageSequence
	import os
	terminus = terminus.lower()
	fetch_data(gene, terminus)
	global _field_base_path_pil
	global _thresholds_pil
	global _channels_pil
	field_base_path = os.path.join(data_cache_path, gene_list[gene][terminus]["plate"], gene+"_4_"+terminus.upper()+"_"+str(field + 1))
	if field_base_path != _field_base_path_pil:
		_field_base_path_pil = field_base_path
		field_image = Image.open(_field_base_path_pil+".tif")
		field_threshold = Image.open(_field_base_path_pil+"_thr.tif")
		_channels_pil = []
		for channel in ImageSequence.Iterator(field_image):
			_channels_pil.append(channel.copy().convert("F"))
		_thresholds_pil = []
		for threshold in ImageSequence.Iterator(field_threshold):
			_thresholds_pil.append(threshold.copy().convert("L"))
	return [_channels_pil[0].copy(), _channels_pil[1].copy(), _channels_pil[2].copy(), _thresholds_pil[0].copy(), _thresholds_pil[1].copy()]

def open_cell_pillow(gene, terminus, field, cell, rotate = False, width = 323):
	from PIL import Image, ImageSequence, ImageOps, ImageDraw, ImageChops, ImageMath
	if rotate:
		width_inter = width * 1.5 # greater than width * 2**0.5
		height = round(width / 2)
	else:
		width_inter = width
	cell_data = gene_list[gene][terminus]["cells"][field][cell]
	crop_centre = cell_data["centre"]
	fill_centre = cell_data["wand"]
	# open field
	channels = open_field_pillow(gene, terminus, field)
	# process phase threshold image to only have cell of interest
	# halve image values, flood fill cell of interst, subtract 127 then double image values
	# really need a image_change_values(min, max, v) function!
	channels[3] = ImageChops.multiply(channels[3], Image.new("L", channels[3].size, color = 127))
	ImageDraw.floodfill(channels[3], (fill_centre[0], fill_centre[1]), 255, thresh = 0)
	channels[3] = ImageChops.subtract(channels[3], Image.new("L", channels[3].size, color = 127))
	channels[3] = ImageMath.eval("a * 2", a = channels[3])
	# Crop (and rotate)
	cell_channels = []
	for channel in channels:
		# if crop outside of image bounds, then first increase canvas size
		offs = 0
		if crop_centre[0] - width_inter / 2 < 0 or crop_centre[1] - width_inter / 2 < 0 or crop_centre[0] + width_inter / 2 > channel.width or crop_centre[1] + width_inter / 2 > channel.height:
			channel = ImageOps.expand(channel, border = round(width_inter / 2), fill = 0)
			offs = round(width_inter / 2)
		# square crop
		channel = channel.crop((crop_centre[0] + offs - width_inter / 2, crop_centre[1] + offs - width_inter / 2, crop_centre[0] + offs + width_inter / 2, crop_centre[1] + offs + width_inter / 2))
		if rotate:
			# if rotating, rotate then crop to final dimensions
			channel = channel.rotate(-gene_list[gene][terminus]["cells"][field][cell]["angle"], resample = Image.BICUBIC)
			channel = channel.crop((width_inter / 2 - width / 2, width_inter / 2 - height / 2, width_inter / 2 + width / 2, width_inter / 2 + height / 2))
		cell_channels.append(channel)
	return cell_channels

# global variables for caching last field of view loaded
_field_base_path_sk = None
_thresholds_sk = None
_channels_sk = None

# open field of view, returning array with one scikit image/np array image per channel and threshold image
# [phase_(gray) mng_(green) dna_(blue) phase_threshold dna_threshold] AKA
# [pha, mng, dna, pth, dth]
# channel images are mode F, 32-bit float, threshold images are mode L, 8 bit (0 or 255)
def open_field(gene, terminus, field):
	import skimage
	import numpy
	import os
	terminus = terminus.lower()
	fetch_data(gene, terminus)
	global _field_base_path_sk
	global _thresholds_sk
	global _channels_sk
	field_base_path = os.path.join(data_cache_path, gene_list[gene][terminus]["plate"], gene+"_4_"+terminus.upper()+"_"+str(field + 1))
	if field_base_path != _field_base_path_sk:
		_field_base_path_sk = field_base_path
		field_image = skimage.io.imread(_field_base_path_sk+".tif")
		field_image = numpy.moveaxis(field_image, [0, 1, 2], [1, 2, 0]) # why loading in a weird order?
		field_threshold = skimage.io.imread(_field_base_path_sk+"_thr.tif")
		_channels_sk = []
		for channel in field_image:
			_channels_sk.append(channel.copy())
		_thresholds_sk = []
		for threshold in field_threshold:
			_thresholds_sk.append(threshold.copy())
	return [_channels_sk[0].astype("uint16").copy(), _channels_sk[1].astype("uint32").copy(), _channels_sk[2].astype("uint16").copy(), _thresholds_sk[0].astype("uint8").copy(), _thresholds_sk[1].astype("uint8").copy()]

def _skimage_crop(image, x, y, w, h):
	return image[int(y):int(y + h), int(x):int(x + w)]

def _skimage_change_values(image, min, max, v):
	image[image >= min & image <= max] = v
	return

def open_cell(gene, terminus, field, cell, rotate = False, width = 323):
	import skimage
	import numpy
	if rotate:
		width_inter = width * 1.5 # greater than width * 2**0.5
		height = round(width / 2)
	else:
		width_inter = width
	cell_data = gene_list[gene][terminus]["cells"][field][cell]
	crop_centre = cell_data["centre"]
	fill_centre = cell_data["wand"]
	# open field
	channels = open_field(gene, terminus, field)
	# process phase threshold image to only have cell of interest, nb. xy swapped in skimage arrays
	channels[3][channels[3] == 255] = 127
	channels[3]=skimage.morphology.flood_fill(channels[3], (fill_centre[1], fill_centre[0]), 255)
	channels[3][channels[3] == 127] = 0
	# Crop (and rotate)
	cell_channels = []
	for channel in channels:
		# if crop outside of image bounds, then first increase canvas size
		offs = 0
		half_width_inter = round(width_inter / 2)
		if crop_centre[0] - half_width_inter < 0 or crop_centre[1] - half_width_inter < 0 or crop_centre[0] + half_width_inter > channel.shape[1] or crop_centre[1] + half_width_inter > channel.shape[0]:
			channel = numpy.pad(channel, ((half_width_inter, half_width_inter), (half_width_inter, half_width_inter)), mode="median")
			offs = half_width_inter
		# square crop
		channel = _skimage_crop(channel, crop_centre[0] + offs - half_width_inter, crop_centre[1] + offs - half_width_inter, width_inter, width_inter)
		if rotate:
			# if rotating, rotate then crop to final dimensions
			channel_dtype = channel.dtype # have have to force data type and use preserve_range=True to prevent rotate from mangling the data
			channel = skimage.transform.rotate(channel, -gene_list[gene][terminus]["cells"][field][cell]["angle"], preserve_range=True).astype(channel_dtype)
			channel = _skimage_crop(channel, half_width_inter - width / 2, half_width_inter - height / 2, width, height)
		cell_channels.append(channel)
	return cell_channels