class TrypTag:
  def __init__(self):
    # user setting: verbose output from accessing data
    self.print_status = True

    # user setting: path in which to cache image data (up to ~8Tb without zip file, ~16Tb with)
    self.data_cache_path = "_tryptag_cache"

    # user setting: remove zip files after download (~doubles data cache size if False)
    self.remove_zip_files = True

    # MAGIC NUMBERS:
    # master zenodo record id
    self.master_zenodo_id = 6862289 # tryptag
    #self.master_zenodo_id = 7258722 # targeted bsf
    self._data_cache_size = 222 * float(17 << 30) + float(25 << 30) # 222 plates @ 17 GiB/plate, + 25 GiB for 1 temporary zip

    # variables which will contain tryptag data
    self.gene_list = None
    self.zenodo_index = None
    self.localisation_ontology = None

    # image properties
    self.um_per_px = 6.5 / 63 # physical pixel size / corrected magnification

    # verbose progress bar for file download
    self._progress_bar = None

    # global variables for caching last field of view loaded
    self._field_base_path_sk = None
    self._thresholds_sk = None
    self._channels_sk = None

  # function to fetch text from a zenodo url, respecting request rate limits
  def _fetch_zenodo_text(self, url):
    from urllib.request import urlopen
    from urllib.error import HTTPError
    import time
    # Zenodo queries are rate limited, so request with that in mind
    response = None
    while response is None:
      try:
        response = urlopen(url)
      except HTTPError as e:
        if self.print_status: print("  Zenodo query rate limit reached, waiting to retry")
        time.sleep(60) # testing shows the rate limiter resets after 50s, so 60s for a bit of space
    text = response.read().decode(response.info().get_param("charset") or "utf-8-sig")
    return text

  # function to fetch zenodo record information
  def _fetch_zenodo_record_json(self, zenodo_id):
    import json
    # fetch Zenodo record JSON
    url = "https://zenodo.org/api/records/"+str(zenodo_id)
    if self.print_status: print("  Fetching Zenodo record for "+str(zenodo_id)+" from: "+url)
    text = self._fetch_zenodo_text(url)
    return json.loads(text)

  # function to fetch text from a file by name in a zenodo record, using zenodo_json 
  def _fetch_zenodo_record_file(self, zenodo_json, file_name):
    from urllib.request import urlopen
    from urllib.error import HTTPError
    for file in zenodo_json["files"]:
      if file["key"] == file_name:
        url = file["links"]["self"]
        if self.print_status: print("  Fetching file "+file_name+" from: "+url)
        with urlopen(url) as response:
           text = response.read().decode(response.info().get_param("charset") or "utf-8-sig")
           return text

  # function for parsing localisation annotation strings
  def _parse_localisation_annotation(self, string):
    import re
    # Regex to extract the term and modifiers from strings such as "term[modifier1,modifier2]"
    annotation_re = re.compile(r'^(?P<term>[^[\]]+)(?:|\[(?P<modifiers>[^\]]+)\])$')
    annotation_list = []
    # Split a string into comma-separated parts, skipping commas that are inside square brackets
    for annotation_string in re.split(r',\s*(?![^\[]*\])', string):
        m = annotation_re.match(annotation_string)
        if m is None:
            # Raise error if an annotation string is malformatted
            raise ValueError()
        term, modifiers = m.groups()
        entry = {
            "term": term,
        }
        if modifiers:
            entry["modifiers"] = modifiers.split(',')
        annotation_list.append(entry)
    return annotation_list

  # fetch gene list/metadata
  # records information in self.gene_list and self.zenodo_index
  def fetch_gene_list(self):
    # if the global variable is None then not loaded yet
    if self.gene_list is None:
      import urllib.request
      # fetch Zenodo record JSON, to get latest version doi
      print("Fetching gene list from Zenodo, record ID: "+str(self.master_zenodo_id))
      zenodo_json = self._fetch_zenodo_record_json(self.master_zenodo_id)
      zenodo_record_id = str(zenodo_json["id"])
      if self.print_status: print("  Using latest Zenodo version, record ID: "+zenodo_record_id)
      # load zenodo record index
      if self.zenodo_index is None:
        # load DOI index
        # TODO: Reformat to use URL from _fetch_zenodo_record_json
        self.zenodo_index = {}
        url = "https://zenodo.org/record/"+zenodo_record_id+"/files/plate_doi_index.tsv?download=1"
        if self.print_status: print("  Fetching plate to Zenodo ID mapping from: "+url)
        response = urllib.request.urlopen(url)
        for line in response:
          line = line.decode(response.info().get_param("charset") or "utf-8-sig").splitlines()[0].split("\t")
          doi_data = {}
          doi_data["master_record_id"] = line[0].split(".")[-1]
          self.zenodo_index[line[1]] = doi_data
      # load localisations table
      # TODO: Reformat to use URL from _fetch_zenodo_record_json
      self.gene_list = {}
      url = "https://zenodo.org/record/"+zenodo_record_id+"/files/localisations.tsv?download=1"
      if self.print_status: print("  Fetching gene data table from: "+url)
      response = urllib.request.urlopen(url)
      for line in response:
        line = line.decode(response.info().get_param("charset") or "utf-8-sig").splitlines()[0].split("\t")
        if line[0] == "Gene ID":
          indices = {}
          for l in range(len(line)):
            indices[line[l]] = l
        else:
          # if not the header line, grab gene data
          self.gene_list[line[0]] = {}
          termini = ["c", "n"]
          for t in termini:
            p = t.upper() + " "
            if line[indices[p + "status"]] == "cell line generated":
              # data from columns which must be present
              terminus_data = {
                "plate": line[indices[p + "plate and well"]].split(" ")[0],
                "well": line[indices[p + "plate and well"]].split(" ")[1],
                "loc": self._parse_localisation_annotation(line[indices[p + "localisation"]]),
                "primer_f": line[indices[p + "primer F"]],
                "primer_r": line[indices[p + "primer R"]]
              }
              # additional data, optional columns
              if p+" classified as faint" in indices:
                terminus_data["signl_low"]: line[indices[p + "fainter than parental"]]
                terminus_data["signal_background"]: line[indices[p + "classified as faint"]]
              # look up zenodo id
              terminus_data["zenodo_id"] = self.zenodo_index[terminus_data["plate"]]["master_record_id"]
              self.gene_list[line[0]][t] = terminus_data

  # recursively called function to build a list of localisation terms
  # Adds a list of parents, a hierachy down to the root node
  # Adds a list of children, a simple list of all children of the current node
  def _parse_localisation_ontology(self, ontology_json, parent=["root"]):
    terms = {}
    # get list of localisations or sublocalisations to iterate through
    if "localisation" in ontology_json:
      loclist = ontology_json["localisation"]
    else:
      loclist = ontology_json["sublocalisation"]
    # iterate through list
    for l in range(len(loclist)):
      # add current entry to terms dict
      terms[loclist[l]["name"]] = loclist[l].copy()
      terms[loclist[l]["name"]]["parent"] = parent
      if "sublocalisation" in loclist[l]:
        # record children to current entry 
        terms[loclist[l]["name"]]["children"] = []
        del terms[loclist[l]["name"]]["sublocalisation"]
        for s in range(len(loclist[l]["sublocalisation"])):
          terms[loclist[l]["name"]]["children"].append(loclist[l]["sublocalisation"][s]["name"])
        # recurse through subterms
        terms.update(self._parse_localisation_ontology(loclist[l], parent=parent + [loclist[l]["name"]]))
    return terms

  # load ontologies for intelligent localisation based searching
  def fetch_ontologies(self):
    if self.localisation_ontology is None:
      import json
      zenodo_json = self._fetch_zenodo_record_json(self.master_zenodo_id)
      # load and parse to flat dict of terms with parent and children names
      self.localisation_ontology = self._parse_localisation_ontology(json.loads(self._fetch_zenodo_record_file(zenodo_json, "localisation_ontology.json")))

  # master localisation search function
  # searches for a match of each term in a localisation list (ie. annotation of a gene id/terminus localisation)
  # against a query localisation term, unless that localisation has a modifier in the exclude_modifiers list
  # if match_subterms is True, then also matches aganst parent structures of the localisation term
  def _localisation_match(self, localisations, query_term, match_subterms=True, exclude_modifiers=["weak", "<10%"]):
    self.fetch_ontologies()
    # iterate through each annotated localisation
    for l in range(len(localisations)):
      # check if the current localisation term should be exlcuded from matching based on modifiers
      modifiers_excluded = False
      print(localisations[l])
      if "modifiers" in localisations[l] and exclude_modifiers is not None:
        for modifier in exclude_modifier:
          if modifier in localisations[l]["modifiers"]:
            modifiers_excluded = True
            break
      # if not excluded based on modifiers, search for an annotation query_term match
      if modifiers_excluded == False:
        # exact match
        if localisations[l]["term"] == query_term:
          return True
        # parent match, if matching subterms
        if match_subterms and query_term in self.localisation_ontology[localisations[l]["term"]]["parent"]:
          return True
    # no matches, so return false
    return False

  # simple localisation match, forgiving of querying nonexistent data
  # checks if a query term match any of a gene id/terminus' annotations or their parents
  # does not match if the gene id/terminus' annotation for that structure has a weak or <10% modifier
  def localisation_match(self, gene_id, terminus, query_term):
    # get gene localisation
    self.fetch_gene_list()
    self.fetch_ontologies()
    if gene_id in self.gene_list:
      if terminus in self.gene_list[gene_id]:
        localisation = self.gene_list[gene_id][terminus]["loc"]
        # test for localisation match
        return self._localisation_match(localisation, query_term)
    return False

  # general progress bar function
  def _show_progress_bar(self, block_num, block_size, total_size):
    import progressbar
    if self._progress_bar is None:
      self._progress_bar = progressbar.ProgressBar(maxval=total_size)
      self._progress_bar.start()
    downloaded = block_num * block_size
    if downloaded < total_size:
      self._progress_bar.update(downloaded)
    else:
      self._progress_bar.finish()
      self._progress_bar = None

  # md5 hash function (for checking zip integrity)
  def _file_md5_hash(self, path, blocksize = 2**20):
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
  # updates self.gene_list and self.zenodo_index
  # places data in self.data_cache_path
  def fetch_data(self, gene_id, terminus):
    import os
    import shutil
    import glob.glob
    import urllib.request
    from zipfile import ZipFile
    from filelock import FileLock
    terminus = terminus.lower()
    # load tryptag data, if not already
    self.fetch_gene_list()
    if terminus in self.gene_list[gene_id]:
      # check if the data cache directory exists, and make if not
      if not os.path.isdir(self.data_cache_path):
        if self.print_status: print("Making data cache directory: "+self.data_cache_path)
        os.mkdir(self.data_cache_path)
        # check disk usage
        space_required = self._data_cache_size
        if self.remove_zip_files == True:
          space_required = space_reqired * 2 # ~ double if retaining zips
        total, used, free = shutil.disk_usage(self.data_cache_path)
        if free < space_required:
          if self.print_status: print("! Insufficient free disk space for full data cache: "+str(round(free / float(2 << 40), 2))+" / "+str(round(space_required / float(2 << 40), 2))+" TiB available !")
      # target paths for zip file and data subdirectory
      plate = self.gene_list[gene_id][terminus]["plate"]
      zip_path = os.path.join(self.data_cache_path, plate+".zip")
      dir_path = os.path.join(self.data_cache_path, plate)
      ziplock_path = os.path.join(self.data_cache_path, plate+".zip.lock")
      dirlock_path = os.path.join(self.data_cache_path, plate+".dir.lock")
      # setup plate-specific lock for threadsafe zip download
      lock = FileLock(ziplock_path)
      lock.acquire()
      # download zip
      try:
        if not os.path.isfile(zip_path) and not os.path.isfile(os.path.join(dir_path, "_"+plate+".zip.md5")):
          print("Fetching data for gene ID "+gene_id+", tagged at "+terminus+" terminus")
          # fetch the processed microscopy data from Zenodo
          if self.print_status: print("  Making plate data directory for: "+plate)
          if "record_id" not in self.zenodo_index[plate]:
            # if not already translated, fetch the latest Zenodo ID from master Zenodo ID
            # fetch Zenodo record JSON, to get latest version doi
            zenodo_json = self._fetch_zenodo_record_json(self.gene_list[gene_id][terminus]["zenodo_id"])
            self.zenodo_index[plate]["record_doi"] = zenodo_json["doi"]
            self.zenodo_index[plate]["record_id"] = str(zenodo_json["id"])
            for file in zenodo_json["files"]:
              if file["key"].endswith("_processed.zip"):
                self.zenodo_index[plate]["record_url"] = file["links"]["self"]
                self.zenodo_index[plate]["record_md5"] = file["checksum"].split(":")[-1]
          # download data
          zip_md5 = 0
          zip_path_temp = os.path.join(self.data_cache_path, plate+".zip.tmp")
          while zip_md5 != self.zenodo_index[plate]["record_md5"]:
            if self.print_status:
              print("  Downloading data from: "+self.zenodo_index[plate]["record_url"])
              urllib.request.urlretrieve(self.zenodo_index[plate]["record_url"], zip_path_temp, self._show_progress_bar)
            else:
              urllib.request.urlretrieve(self.zenodo_index[plate]["record_url"], zip_path_temp)
            if self.print_status: print("  Checking MD5 hash of: "+plate+".zip.tmp")
            zip_md5 = self._file_md5_hash(zip_path_temp)
            if zip_md5 != self.zenodo_index[plate]["record_md5"]:
              if self.print_status: print("! MD5 of downloaded zip is incorrect ("+zip_md5+"), retrying download !")
          if self.print_status: print("  Download complete")
          shutil.move(zip_path_temp, zip_path)
          with open(zip_path+".md5", "w") as file: file.write(zip_md5)
      finally:
        lock.release()
      # setup plate-specific lock for threadsafe decompression
      lock = FileLock(dirlock_path)
      lock.acquire()
      # decompress zip
      try:
        if not os.path.isfile(os.path.join(dir_path, "_"+plate+".zip.md5")):
          # unzip data
          if self.print_status: print("  Decompressing image data from: "+plate+".zip")
          try:
            with ZipFile(zip_path) as archive:
              suffix = "_roisCells.txt"
              count_checked = 0
              total_names = len(archive.namelist())
              count_decompressed = 0
              missing = []
              # do decompression
              for file in archive.namelist():
                count_checked += 1
                if self.print_status: self._show_progress_bar(count_checked, 1, total_names)
                # loop through all files, finding files ending with the cell roi suffix and not starting with control (or common misspellings)
                if file.endswith(suffix) and not (os.path.split(file)[-1].startswith("control") or os.path.split(file)[-1].startswith("ontrol") or os.path.split(file)[-1].startswith("Control")):
                  source_path = os.path.split(file)
                  # infer main tif and thresholded tif image filenames
                  file_roi = file
                  file_tif = file[:-len(suffix)]+".tif"
                  file_thr = file[:-len(suffix)]+"_thr.tif"
                  if file_roi in archive.namelist() and file_tif in archive.namelist() and file_thr in archive.namelist():
                    count_decompressed += 1
                    # if all exist, decompress and move to plate directory
                    archive.extract(file_roi, self.data_cache_path)
                    archive.extract(file_tif, self.data_cache_path)
                    archive.extract(file_thr, self.data_cache_path)
                    target_path = os.path.join(self.data_cache_path, plate)
                    shutil.move(os.path.join(self.data_cache_path, file_roi), os.path.join(target_path, os.path.split(file_roi)[-1]))
                    shutil.move(os.path.join(self.data_cache_path, file_tif), os.path.join(target_path, os.path.split(file_tif)[-1]))
                    shutil.move(os.path.join(self.data_cache_path, file_thr), os.path.join(target_path, os.path.split(file_thr)[-1]))
                  else:
                    missing.append(os.path.split(file_tif)[-1])
            if self.print_status and len(missing) > 0:
              print("! Expected files were missing in the zip for the following gene IDs !")
              print(" ".join(missing))
            # remove any subdirectories remaining from decompression
            obj = os.scandir(os.path.join(self.data_cache_path, plate))
            for entry in obj :
              if entry.is_dir():
                shutil.rmtree(os.path.join(self.data_cache_path, plate, entry.name))
            obj.close()
            # copy source zip md5 to the data directory, also marks decompression as complete
            shutil.copyfile(zip_path+".md5", os.path.join(dir_path, "_"+plate+".zip.md5"))
            if self.remove_zip_files:
              os.remove(zip_path)
              os.remove(zip_path+".md5")
            if self.print_status:
              print("  Decompressed "+str(count_decompressed)+" fields of view")
          except BadZipfile:
            print ("! Zip file invalid: "+plate+".zip !")
            if self.remove_zip_files: os.remove(zip_path)
      finally:
        lock.release()
      # count fields of view and number of cells
      base_path = os.path.join(self.data_cache_path, plate)
      if "fields_count" not in self.gene_list[gene_id][terminus] and os.path.isdir(base_path):
        if self.print_status: print("  Counting image data files for: "+gene_id+" "+terminus)
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
        self.gene_list[gene_id][terminus]["fields_count"] = len(cells)
        self.gene_list[gene_id][terminus]["cells_per_field"] = [len(x) for x in cells]
        self.gene_list[gene_id][terminus]["cells"] = cells.copy()

  # checks if cached image data is available for gene and terminus
  def check_if_cached(self, gene_id, terminus):
    # False if gene_list or zenodo_index have not been loaded yet
    if self.gene_list is None:
      return False
    if self.zenodo_index is None:
      return False
    # path for data subdirectory
    plate = self.gene_list[gene_id][terminus]["plate"]
    dir_path = os.path.join(self.data_cache_path, plate)
    # False if MD5 not yet copied to data directory
    if not os.path.isfile(os.path.join(dir_path, "_"+plate+".zip.md5")):
      return False
    # otherwise cached
    return True

  # checks data in self.data_cache_path to make sure it is the latest version
  def check_data_cache(self):
    import os
    # load tryptag data, if not already
    self.fetch_gene_list()
    print("Checking data cache for errors")
    for plate in self.zenodo_index:
      zip_path = os.path.join(self.data_cache_path, plate+".zip")
      dir_path = os.path.join(self.data_cache_path, plate)
      if os.path.isdir(dir_path) or os.path.isfile(zip_path):
        zenodo_json = self._fetch_zenodo_record_json(self.zenodo_index[plate]["master_record_id"])
        for file in zenodo_json["files"]:
          if file["key"].endswith("_processed.zip"):
            self.zenodo_index[plate]["record_md5"] = file["checksum"].split(":")[-1]
      if os.path.isdir(dir_path):
        md5_file = os.path.join(dir_path, "_"+plate+".zip.md5")
        if not os.path.isfile(md5_file):
          if self.print_status: print("  MD5 file mising from image directory:", plate)
        else:
          with open(md5_file, "r") as file:
            md5 = file.read()
            if md5 != self.zenodo_index[plate]["record_md5"]:
              if self.print_status: print("  MD5 file in directory does not match zenodo record md5:", plate)
      if os.path.isfile(zip_path):
        if self.remove_zip_files:
          if self.print_status: print("  Zip found which should have been removed:", plate)
        md5_file = os.path.join(zip_path+".md5")
        with open(md5_file, "r") as file:
          md5 = file.read()
          if md5 != self.zenodo_index[plate]["record_md5"]:
            if self.print_status: print("  MD5 file for zip file does not match zenodo record MD5:", plate)

  # force load of all data by iterating through every gene and terminus
  def fetch_all_data(self):
    self.fetch_gene_list()
    for gene in self.gene_list:
      for terminus in ["c", "n"]:
        self.fetch_data(gene, terminus)

  # open field of view, returning array with one scikit image/np array image per channel and threshold image
  # [phase_(gray) mng_(green) dna_(blue) phase_threshold dna_threshold] AKA
  # [pha, mng, dna, pth, dth]
  # channel images are mode F, 32-bit float, threshold images are mode L, 8 bit (0 or 255)
  def open_field(self, gene, terminus, field):
    import skimage.io
    import numpy
    import os
    terminus = terminus.lower()
    self.fetch_data(gene, terminus)
    field_base_path = os.path.join(self.data_cache_path, self.gene_list[gene][terminus]["plate"], gene+"_4_"+terminus.upper()+"_"+str(field + 1))
    if field_base_path != self._field_base_path_sk:
      self._field_base_path_sk = field_base_path
      field_image = skimage.io.imread(self._field_base_path_sk+".tif")
      field_image = numpy.moveaxis(field_image, [0, 1, 2], [1, 2, 0]) # why loading in a weird order?
      field_threshold = skimage.io.imread(self._field_base_path_sk+"_thr.tif")
      self._channels_sk = []
      for channel in field_image:
        self._channels_sk.append(channel.copy())
      self._thresholds_sk = []
      for threshold in field_threshold:
        self._thresholds_sk.append(threshold.copy())
    return [self._channels_sk[0].astype("uint16").copy(), self._channels_sk[1].astype("uint32").copy(), self._channels_sk[2].astype("uint16").copy(), self._thresholds_sk[0].astype("uint8").copy(), self._thresholds_sk[1].astype("uint8").copy()]

  def _skimage_crop(self, image, x, y, w, h):
    return image[int(y):int(y + h), int(x):int(x + w)]

  def _skimage_change_values(self, image, min, max, v):
    image[image >= min & image <= max] = v
    return

  # master function for opening a cell
  def _open_cell(self, gene, terminus, field, crop_centre, fill_centre, phathr = None, dnathr = None, angle = 0, rotate = False, width = 323):
    import skimage.morphology
    import skimage.transform
    import numpy
    if rotate:
      width_inter = width * 1.5 # greater than width * 2**0.5
      height = round(width / 2)
    else:
      width_inter = width
    # open field
    channels = self.open_field(gene, terminus, field)
    # replace with custom phathr and dnathr, if defined
    if not phathr is None:
      channels[3] = phathr.copy()
    if not dnathr is None:
      channels[4] = dnathr.copy()
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
      channel = self._skimage_crop(channel, crop_centre[0] + offs - half_width_inter, crop_centre[1] + offs - half_width_inter, width_inter, width_inter)
      if rotate:
        # if rotating, rotate then crop to final dimensions
        channel_dtype = channel.dtype # have have to force data type and use preserve_range=True to prevent rotate from mangling the data
        channel = skimage.transform.rotate(channel, angle, preserve_range=True).astype(channel_dtype)
        channel = self._skimage_crop(channel, half_width_inter - width / 2, half_width_inter - height / 2, width, height)
      cell_channels.append(channel)
    return cell_channels

  # open a cell, cropped from a field of view
  # uses the phase and dna threshold images from tryptag
  # cell x, y coordinate in phase threshold from tryptag
  def open_cell(self, gene, terminus, field, cell, rotate = False, width = 323):
    self.fetch_data(gene, terminus)
    cell_data = self.gene_list[gene][terminus]["cells"][field][cell]
    crop_centre = cell_data["centre"]
    fill_centre = cell_data["wand"]
    angle = cell_data["angle"]
    return self._open_cell(gene, terminus, field, crop_centre, fill_centre, phathr = None, dnathr = None, angle = angle, rotate = rotate, width = width)

  # opens a custom cell from a field of view
  # instead of using tryptag thresholded phase and dna images, user-provided phathr and dnathr (uint8, 255 = object)
  # crop_centre is the (x, y) tuple around which to crop
  # fill_centre is a (x, y) tuple of a pixel which is in the target cell object (255) in phathr
  def open_cell_custom(self, gene, terminus, field, phathr, dnathr, crop_centre, fill_centre, angle = False, rotate = False, width = 323):
    return self._open_cell(gene, terminus, field, crop_centre, fill_centre, phathr = phathr, dnathr = dnathr, angle = angle, rotate = rotate, width = width)
