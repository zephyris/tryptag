from functools import cached_property
import urllib.request
import os
import shutil
import glob
import urllib.request
from zipfile import ZipFile, BadZipFile

import numpy
import progressbar
from filelock import FileLock
import skimage.io
import skimage.morphology
import skimage.transform

class TrypTag:
  def __init__(
      self,
      verbose=True,
      data_cache_path="./_tryptag_cache",
      remove_zip_files=True,
      master_zenodo_id=6862289,
      data_cache_size=222 * 17 * float(2 << 30) + 25 * float(2 << 30), # 222 plates @ 17 GiB/plate, + 25 GiB for 1 temporary zip
      um_per_px=6.5 / 63,
  ):
    """Initialise TrypTag data access object
    
    Keyword arguments:
    verbose -- print verbose output from accessing data (default `True`)
    data_cache_path -- the directory that will hold the downloaded TrypTag data (default `./_tryptag_cache`, requires up to 8Tb)
    remove_zip_files -- whether to remove zip files after download and extraction (default `True`, doubles data cache size if `False`)
    master_zenodo_id -- Zenodo ID of the master TrypTag deposition (default `6862289`, do not change unless you know what you're doing)
    data_cache_size -- Size (in bytes) of the data cache (default chosen for TrypTag's 222 plates, roughly 3.8 Tb)
    um_per_px -- physical pixel size / corrected magnification
    """
    # user setting: verbose output from accessing data
    self.print_status = verbose

    # user setting: path in which to cache image data (up to ~8Tb without zip file, ~16Tb with)
    self.data_cache_path = data_cache_path

    # user setting: remove zip files after download (~doubles data cache size if False)
    self.remove_zip_files = remove_zip_files

    # MAGIC NUMBERS:
    # master zenodo record id
    self.master_zenodo_id = master_zenodo_id # tryptag
    #self.master_zenodo_id = 7258722 # targeted bsf
    self._data_cache_size = data_cache_size

    # image properties
    self.um_per_px = um_per_px

    # verbose progress bar for file download
    self._progress_bar = None

    # global variables for caching last field of view loaded
    self._field_base_path_sk = None
    self._thresholds_sk = None
    self._channels_sk = None

  # function to fetch text from a zenodo url, respecting request rate limits
  # TODO: Cache per session(?)
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

  # function to fetch zenodo record information, using _fetch_zenodo_text
  def _fetch_zenodo_record_json(self, zenodo_id):
    import json
    # fetch Zenodo record JSON
    url = "https://zenodo.org/api/records/"+str(zenodo_id)
    if self.print_status: print("  Fetching Zenodo record for "+str(zenodo_id)+" from: "+url)
    text = self._fetch_zenodo_text(url)
    return json.loads(text)

  # function to fetch text from a file by name in a zenodo record, using zenodo_json and _fetch_zenodo_text
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

  @cached_property
  def zenodo_record_id(self):
    zenodo_json = self._fetch_zenodo_record_json(self.master_zenodo_id)
    return str(zenodo_json["id"])

  @cached_property
  def zenodo_index(self):
    # load plate to zenodo record index
    zenodo_index = {}
    # download plate_doi_index.tsv from master record
    zenodo_json = self._fetch_zenodo_record_json(self.zenodo_record_id)
    lines = self._fetch_zenodo_record_file(zenodo_json, "plate_doi_index.tsv").splitlines()
    # parse tsv to dict indexed by plateID_YYYYMMDD
    for line in lines:
      line = line.split("\t")
      zenodo_index[line[1]] = {"master_record_id": line[0].split(".")[-1]}
    return zenodo_index

  # fetch gene list/metadata
  # records information in self.gene_list and self.zenodo_index
  @cached_property
  def gene_list(self):
    # fetch Zenodo record JSON, to get latest version doi
    if self.print_status: print("Fetching gene list from Zenodo, record ID: "+str(self.master_zenodo_id))
    if self.print_status: print("  Using latest Zenodo version, record ID: "+self.zenodo_record_id)
    # load localisations table
    # download localisations.tsv from master record
    zenodo_json = self._fetch_zenodo_record_json(self.zenodo_record_id)
    lines = self._fetch_zenodo_record_file(zenodo_json, "localisations.tsv").splitlines()
    gene_list = {}
    # parse line by line, expects the first line to be headers and grab indices
    for line in lines:
      line = line.split("\t")
      if line[0] == "Gene ID":
        indices = {}
        for l in range(len(line)):
          indices[line[l]] = l
      else:
        # if not the header line, grab gene data
        gene_list[line[0]] = {}
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
            gene_list[line[0]][t] = terminus_data
    return gene_list

  @cached_property
  def worklist_all(self):
    # all genes ids/termini, as a worklist
    return [
      {
        "gene_id": gene_id,
        "terminus": terminus,
      } for gene_id, terminus in
      ((gene_id, terminus) for gene_id, gene_entry in self.gene_list.items() for terminus in ["n", "c"] if terminus in gene_entry)
    ]

  @cached_property
  def worklist_parental(self):
    # parental dummy gene id/termini entries, as a worklist
    return [
      {
        "gene_id": gene_id,
        "terminus": terminus,
      } for gene_id, terminus in
      ((gene_id, terminus) for gene_id, gene_entry in self.gene_list.items() for terminus in ["n", "c"] if "wild-type" in gene_id and terminus in gene_entry)
    ]

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
    # if the root node, add a root pseudolocalisation
    if parent == ["root"]:
      terms["root"] = {
        "name": "root",
        "parent": [],
        "children": [x["name"] for x in ontology_json["localisation"]]
      }
    return terms

  # load ontologies for intelligent localisation based searching
  @cached_property
  def localisation_ontology(self):
    import json
    zenodo_json = self._fetch_zenodo_record_json(self.master_zenodo_id)
    # load and parse to flat dict of terms with parent and children names
    return self._parse_localisation_ontology(json.loads(self._fetch_zenodo_record_file(zenodo_json, "localisation_ontology.json")))

  # localisation match function for searches
  # searches for a match of each term in a localisation list for a gene id and terminus against
  # a query localisation term, unless that localisation has a modifier in the exclude_modifiers list
  # if match_subterms is True, then also matches aganst parent structures of the localisation term
  # if required_modifiers is not None, then it must also match all of them
  def localisation_match(self, gene_id, terminus, query_term, match_subterms=True, exclude_modifiers=["weak", "<10%"], required_modifiers=None):
    # get query localisation
    localisations = self.gene_list[gene_id][terminus]["loc"]
    # get ontology, and lead to keyerror if query_term not in ontology
    ont = self.localisation_ontology[query_term]
    # iterate through each annotated localisation
    for l in range(len(localisations)):
      # check if any necessary modifiers are present
      modifiers_included = True
      if required_modifiers is not None:
        modifiers_count = 0
        if "modifiers" in localisations[l]:
          for modifier in required_modifiers:
            if modifier in localisations[l]["modifiers"]:
              modifiers_count += 1
        if modifiers_count < len(required_modifiers):
          modifiers_included = False
      if modifiers_included == False:
        break
      # check if the current localisation term should be exlcuded from matching based on modifiers
      modifiers_excluded = False
      if "modifiers" in localisations[l] and exclude_modifiers is not None:
        for modifier in exclude_modifiers:
          if modifier in localisations[l]["modifiers"]:
            modifiers_excluded = True
            break
      # if not excluded based on modifiers, search for an annotation query_term match
      if modifiers_excluded == True:
        break
      # exact match
      if localisations[l]["term"] == query_term:
        return True
      # parent match, if matching subterms
      if match_subterms == True and query_term in self.localisation_ontology[localisations[l]["term"]]["parent"]:
        return True
    # no matches, so return false
    return False

  # get a list of gene hits from a localisation_match
  def localisation_search(self, query_term, match_subterms=True, exclude_modifiers=["weak", "<10%"], include_modifiers=None, required_modifiers=None):
    # check all against query
    hits = []
    for gene_id in self.gene_list:
      for terminus in ["n", "c"]:
        if terminus in self.gene_list[gene_id]:
          if self.localisation_match(gene_id, terminus, query_term, match_subterms=match_subterms, exclude_modifiers=exclude_modifiers, required_modifiers=required_modifiers):
            hits.append({
              "gene_id": gene_id,
              "terminus": terminus
            })
    return hits

  # general progress bar function
  def _show_progress_bar(self, block_num, block_size, total_size):
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
    terminus = terminus.lower()
    # load tryptag data, if not already
    if terminus in self.gene_list[gene_id]:
      # check if the data cache directory exists, and make if not
      if not os.path.isdir(self.data_cache_path):
        if self.print_status: print("Making data cache directory: "+self.data_cache_path)
        os.mkdir(self.data_cache_path)
        # check disk usage
        space_required = self._data_cache_size
        if self.remove_zip_files == False:
          space_required = space_required * 2 # ~ double if retaining zips
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
          os.mkdir(dir_path)
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
          except BadZipFile:
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
    # Check for attributes in `__dict__` - `functools.cached_property` writes an attribute of the same name
    if not (("gene_list" in self.__dict__) and ("zenodo_index" in self.__dict__)):
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
    for gene in self.gene_list:
      for terminus in ["c", "n"]:
        self.fetch_data(gene, terminus)

  # open field of view, returning array with one scikit image/np array image per channel and threshold image
  # [phase_(gray) mng_(green) dna_(blue) phase_threshold dna_threshold] AKA
  # [pha, mng, dna, pth, dth]
  # channel images are mode F, 32-bit float, threshold images are mode L, 8 bit (0 or 255)
  def open_field(self, gene, terminus, field):
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

  # worker for multiprocess/thread analysis
  # iterates through work_list of {"gene_id": gene_id, "terminus": terminus} entries
  # runs analysis_function on gene_id and terminus in each entry
  # returns list of result objects {"gene_id": gene_id, "terminus": terminus, "result": result_from_analysis_function}
  def _list_analysis_worker(self, work_list, analysis_function, tryptag=None, worker_index=None):
    # if passed a TrypTag instance then use a copy (for running as a spawned thread/process), otherwise use self as current_tryptag
    if tryptag is None:
      current_tryptag = self
    else:
      from copy import deepcopy
      current_tryptag = deepcopy(tryptag)
    results = []
    for i in range(len(work_list)):
      entry = work_list[i]
      if self.print_status and worker_index is not None: print("  Worker index", worker_index + 1, "processing worklist entry", i + 1, "of", len(work_list))
      result = {
        "gene_id": entry["gene_id"],
        "terminus": entry["terminus"]
      }
      current_tryptag.fetch_data(entry["gene_id"], entry["terminus"])
      current_result = analysis_function(current_tryptag, entry["gene_id"], entry["terminus"])
      result.update({
          "result": current_result
      })
      results.append(result)
    return results

  # simple handler for multiprocess/thread analysis
  # analyses work_list of {"gene_id": gene_id, "terminus": terminus} entries
  # uses analysis_function, which must:
  #   take only tryptag, gene_id and terminus as arguments
  #   not use global variables and return a single variable as a result
  # returns a list of {"gene_id": gene_id, "terminus": terminus, "result": result_from_analysis_function}
  def analyse_list(self, work_list, analysis_function, workers=None, multiprocess_mode=None):
    import concurrent.futures
    import multiprocessing
    import numpy
    if self.print_status: print("Analysing worklist")
    # deduplicate work_list
    dedup_work_list = []
    for entry in work_list:
      if entry not in dedup_work_list:
        dedup_work_list.append(entry)
    if multiprocess_mode is None:
      # run by direct call of the _list_analysis_worker function
      if self.print_status: print("  Single process")
      results = self._list_analysis_worker(work_list, analysis_function)
    else:
      # get number of workers, default to number of cpus
      if workers is None:
        workers = multiprocessing.cpu_count()
      # split dedup_work_list list into lists for each worker
      split_work_list = numpy.array_split(dedup_work_list, workers)
      if multiprocess_mode == "process":
        # setup executor as a process pool
        if self.print_status: print("  Parallel processes with", workers, "workers")
        executor = concurrent.futures.ProcessPoolExecutor(max_workers=workers)
      elif multiprocess_mode == "thread":
        # setup executor as a thread pool
        if self.print_status: print("  Parallel threads with", workers, "workers")
        executor = concurrent.futures.ThreadPoolExecutor(max_workers=workers)
      # trigger fetch of gene_list and zenodo_index, prior to copying for each thread
      self.zenodo_index
      self.gene_list
      # loop for setting up futures
      results = []
      futures = []
      for i in range(len(split_work_list)):
        # pass each split_work_list list item to a _list_analysis_worker function
        if len(split_work_list[i]) > 0:
          future = executor.submit(self._list_analysis_worker, work_list=split_work_list[i], analysis_function=analysis_function, tryptag=self, worker_index=i)
        futures.append(future)
      for future in concurrent.futures.as_completed(futures):
        # concatenate results as they are returned
        results += future.result()
      # return concatenated results
    return results

  # Gets a list of all genes and analyses
  def analyse_all(self, function):
    list = []
    for gene_id in self.gene_list:
      for terminus in ["n", "c"]:
        if terminus in self.gene_list[gene_id]:
          list.append({
              "gene_id": gene_id,
              "terminus": terminus
          })
    return (self.analyse_cells(list, function))
