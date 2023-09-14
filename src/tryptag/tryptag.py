from functools import cached_property
import urllib.request
import os
import shutil
import glob
from zipfile import ZipFile, BadZipFile
from typing import NamedTuple, Any, Callable

import numpy
from tqdm.auto import tqdm
from filelock import FileLock
import skimage.io
import skimage.morphology
import skimage.transform

class CellLine:
  """
  CellLine object holding information on a cell line dataset - life cycle stage, tagged gene and tagged terminus.
  """
  def __init__(
    self,
    gene_id: str,
    terminus: str,
    life_stage: str = None
  ):
    self.gene_id = gene_id
    self.terminus = terminus
    self.life_stage = life_stage

  def __repr__(self):
    return self.__str__()

  def __str__(self):
    return " ".join([str(x) for x in ["gene_id =", self.gene_id, "terminus =", self.terminus, "life_stage =", self.life_stage]])
  
  def __hash__(self):
    return hash((self.gene_id, self.terminus, self.life_stage))

class CellImage():
  """
  CellImage object holding information on a (cropped) section of a TrypTag image containing a specific cell.
  """
  def __init__(
    self,
    phase: Any,
    mng: Any,
    dna: Any,
    phase_mask: Any,
    dna_mask: Any,
    phase_mask_othercells: Any,
    rotated: bool,
    field_index: int = 0,
    cell_index: int = 0,
    cell_line: CellLine | None = None
  ):
    self.phase = phase
    self.mng = mng
    self.dna = dna
    self.phase_mask = phase_mask
    self.dna_mask = dna_mask
    self.phase_mask_othercells = phase_mask_othercells
    self.rotated = rotated
    self.field_index = field_index
    self.cell_index = cell_index
    self.cell_line = cell_line

  def __repr__(self):
    string = "\n".join([str(x) for x in ["phase", self.phase, "mng", self.mng, "dna", self.dna, "phase_mask", self.phase_mask, "dna_mask", self.dna_mask, "phase_mask_othercells", self.phase_mask_othercells]])
    if self.cell_line is not None:
      string += "\n" + str(self.cell_line)
    string += "\n" + " ".join([str(x) for x in ["field_index =", self.field_index, "cell_index =", self.cell_index, "rotated =", self.rotated]])
    return string

  def __str__(self):
    string = ""
    if self.cell_line is not None:
      string += str(self.cell_line)
    string += " " + " ".join([str(x) for x in ["field_index =", self.field_index, "cell_index =", self.cell_index, "rotated =", self.rotated]])
    return string

class FieldImage():
  """
  FieldImage object holding information on a TrypTag image.
  """
  def __init__(
    self,
    phase: Any,
    mng: Any,
    dna: Any,
    phase_mask: Any,
    dna_mask: Any,
    field_index: int = None,
    cell_line: Any = None
  ):
    self.phase = phase
    self.mng = mng
    self.dna = dna
    self.phase_mask = phase_mask
    self.dna_mask = dna_mask
    self.field_index = field_index
    self.cell_line = cell_line

  def __repr__(self):
    string = "\n".join([str(x) for x in ["phase", self.phase, "mng", self.mng, "dna", self.dna, "phase_mask", self.phase_mask, "dna_mask", self.dna_mask]])
    if self.cell_line is not None:
      string += "\n" + str(self.cell_line)
    string += "\n" + " ".join([str(x) for x in ["field_index =", self.field_index]])
    return string

  def __str__(self):
    string = ""
    if self.cell_line is not None:
      string += str(self.cell_line)
    string += " " + " ".join([str(x) for x in ["field_index =", self.field_index]])
    return string

class _tqdmDownload(tqdm):
  def urllib_callback(self, transferred_blocks, block_size, total_size):
    self.total = total_size
    return self.update(transferred_blocks * block_size - self.n)

class TrypTag:
  def __init__(
      self,
      verbose=True,
      data_cache_path="./_tryptag_cache",
      remove_zip_files=True,
      master_zenodo_id=6862289,
      data_cache_plates=222,
      data_cache_platesize=17 * float(2 << 30),
      data_cache_zipsize=14 * float(2 << 30),
      um_per_px=6.5 / 63,
      life_stages=["procyclic"]
  ):
    """Initialise TrypTag data access object
    
    Keyword arguments:
    verbose -- print verbose output from accessing data (default `True`)
    data_cache_path -- the directory that will hold the downloaded TrypTag data (default `./_tryptag_cache`, requires up to 8Tb)
    remove_zip_files -- whether to remove zip files after download and extraction (default `True`, doubles data cache size if `False`)
    master_zenodo_id -- Zenodo ID of the master TrypTag deposition (default `6862289`, do not change unless you know what you're doing)
    data_cache_plates -- Number of plates in data cache (default chosen for TrypTag's 222 plates)
    data_cache_platesize -- Size (in bytes) of one plate directory (default chosen for TrypTag's plates, 17 Gb)
    data_cache_zipsize -- Size (in bytes) of one plate zip file (default chosen for TrypTag's plates, 14 Gb)
    um_per_px -- physical pixel size / corrected magnification
    life_stages -- list of life cycle stages covered by this dataset, loading will default to first entry
    """
    # user setting: verbose output from accessing data
    self.print_status = verbose

    # user setting: path in which to cache image data
    self.data_cache_path = data_cache_path

    # user setting: remove zip files after download (~doubles data cache size if False)
    self.remove_zip_files = remove_zip_files

    # MAGIC NUMBERS:
    # master zenodo record id
    self.master_zenodo_id = master_zenodo_id
    # data cache information
    self._data_cache_plates = data_cache_plates
    self._data_cache_platesize = data_cache_platesize
    self._data_cache_zipsize = data_cache_zipsize
    self._data_cache_size = data_cache_plates * data_cache_platesize + data_cache_zipsize
    self._data_cache_zipsize = data_cache_plates * data_cache_zipsize

    # image properties
    self.um_per_px = um_per_px

    # standard terminus names
    self.termini = ["n", "c"]

    # global variable for caching strings requested from urls
    self._url_str_cache = {}

    # global variables for caching last field of view loaded
    self._field_base_path_sk = None
    self._thresholds_sk = None
    self._channels_sk = None

    # variables for handling of different data sources
    self.life_stages = life_stages

  def _fetch_zenodo_text(self, url: str, use_file_cache: bool = True) -> str:
    """
    Data retrieval helper function. Fetch a string from any Zenodo URL, respecting rate limits.
    Caches data as text files in `_zenodo` subdirectory of `data_cache_path`.
    Also caches data in memory in `self._url_str_cache` to boost performance.

    :param url: Zenodo URL.
    :param use_file_cache: Read from cached file, if available.
    :return: String of the contents of URL.
    """
    # check memory cache and return string if in memory
    if url in self._url_str_cache:
      if self.print_status: print("  Using memory cache")
      return self._url_str_cache[url]
    # check file cache and, if cached, read file and return string
    zenodo_cache_path = os.path.join(self.data_cache_path, "_zenodo")
    file_cache_path = os.path.join(zenodo_cache_path, url[len("https://"):].replace("/", "-")+".txt")
    if os.path.isfile(file_cache_path):
      if self.print_status: print("  Using file cache")
      with open(file_cache_path, "r") as f:
        return f.read()
    # otherwise, fetch from url
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
        time.sleep(60) # MAGIC NUMBER: testing shows the rate limiter resets after 50s, so 60s for a bit of space
    text = response.read().decode(response.info().get_param("charset") or "utf-8-sig")
    # cache in memory
    self._url_str_cache[url] = text
    # make directories for file cache
    if not os.path.isdir(self.data_cache_path):
      if self.print_status: print("Making data cache directory: "+self.data_cache_path)
      os.mkdir(self.data_cache_path)
    if not os.path.isdir(zenodo_cache_path):
      if self.print_status: print("  making zenodo cache directory: "+zenodo_cache_path)
      os.mkdir(zenodo_cache_path)
    # cache as file, using filelock to avoid double writes
    lock = FileLock(file_cache_path+".lock")
    lock.acquire()
    try:
      # write to cache
      with open(file_cache_path, "w") as f:
        f.write(text)
    finally:
      lock.release()
    # return
    return text

  def _fetch_zenodo_record_json(self, zenodo_id: str):
    """
    Data retrieval helper function. Fetch a record JSON for a Zenodo ID.

    :param zenodo_id: Zenodo ID
    :return: JSON record information for `zenodo_id`.
    """
    import json
    # fetch Zenodo record JSON
    url = "https://zenodo.org/api/records/"+str(zenodo_id)
    if self.print_status: print("  Fetching Zenodo record for "+str(zenodo_id)+" from: "+url)
    text = self._fetch_zenodo_text(url)
    return json.loads(text)

  def _fetch_zenodo_record_file(self, zenodo_json, file_name: str) -> str:
    """
    Data retrieval helper function. Fetch the contents of a text file in a Zenodo ID.

    :param zenodo_json: JSON record information for a Zenodo ID.
    :param file_name: File name.
    :return: String of the contents of `file_name`, fetched from Zenodo.
    """
    from urllib.request import urlopen
    from urllib.error import HTTPError
    for file in zenodo_json["files"]:
      if file["key"] == file_name:
        url = file["links"]["self"]
        if self.print_status: print("  Fetching file "+file_name+" from: "+url)
        return self._fetch_zenodo_text(url)

  # function for parsing localisation annotation strings
  def _parse_localisation_annotation(self, string: str) -> list:
    """
    Parse standard localisation annotation string in the form term[modifier,modifier],term,term[modifier] etc.
    :return: List of dicts containing localisation terms and modifiers in the form `{"term": term, "modifiers": [modifier, modifier]}`
    """
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
  def zenodo_record_id(self) -> str:
    """
    Master Zenodo ID from which to load data.
    """
    zenodo_json = self._fetch_zenodo_record_json(self.master_zenodo_id)
    return str(zenodo_json["id"])

  @cached_property
  def zenodo_index(self) -> list:
    """
    List of Zenodo ID to plateID_YYYYMMDD mappings, for data retrieval.
    """
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

  @cached_property
  def gene_list(self) -> list:
    """
    Master data list, organised as dicts within `gene_list[life_stage][gene_id][terminus]` structure.
    `life_stage` is the life cycle stage, for tryptag always `"procyclic"`.
    `gene_id` is a TriTrypDB gene ID, eg. `"Tb927.5.3250"`.
    `terminus` is `"n"` or `"c"`.
    """
    # fetch Zenodo record JSON, to get latest version doi
    if self.print_status: print("Fetching gene list from Zenodo, record ID: "+str(self.master_zenodo_id))
    if self.print_status: print("  Using latest Zenodo version, record ID: "+self.zenodo_record_id)
    # load localisations table
    # download localisations.tsv from master record
    zenodo_json = self._fetch_zenodo_record_json(self.zenodo_record_id)
    lines = self._fetch_zenodo_record_file(zenodo_json, "localisations.tsv").splitlines()
    # load into default life stage
    # TODO: Handle loading of multiple different life stages
    gene_list = {self.life_stages[0]: {}}
    # parse line by line, expects the first line to be headers and grab indices
    for line in lines:
      line = line.split("\t")
      if line[0] == "Gene ID":
        indices = {}
        for l in range(len(line)):
          indices[line[l]] = l
      else:
        # if not the header line, grab gene data
        gene_list[self.life_stages[0]][line[0]] = {}
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
            gene_list[self.life_stages[0]][line[0]][t] = terminus_data
    return gene_list

  def worklist_all(self, life_stage: str = None) -> list:
    """
    All `gene_id` and `terminus` combinations with data, as a `CellLine` object containing `life_stage`, `gene_id` and `terminus`.

    :param life_stage: Which life cycle stage to load data from, default is first entry in `self.life_stages`.
    """
    if life_stage is None:
      life_stage = self.life_stages[0]
    return [
      CellLine(
        life_stage=life_stage,
        gene_id=gene_id,
        terminus=terminus,
       ) for gene_id, terminus in
      ((gene_id, terminus) for gene_id, gene_entry in self.gene_list[life_stage].items() for terminus in self.termini if terminus in gene_entry)
    ]

  def worklist_parental(self, life_stage: str = None) -> list:
    """
    All dummy `gene_id` and `terminus` combinations which correspond to parental cell line samples, as a `CellLine` object containing `life_stage`, `gene_id` and `terminus`.
    """
    if life_stage is None:
      life_stage = self.life_stages[0]
    return [
      CellLine(
        life_stage=life_stage,
        gene_id=gene_id,
        terminus=terminus,
       ) for gene_id, terminus in
      ((gene_id, terminus) for gene_id, gene_entry in self.gene_list[life_stage].items() for terminus in self.termini if "wild-type" in gene_id and terminus in gene_entry)
    ]

  @cached_property
  def localisation_ontology(self) -> list:
    """
    Localisation ontology annotation terms for intelligent localisation-based searching.
    """

    def _parse_localisation_ontology(ontology_json, parent: list = ["root"]) -> dict:
      """
      Parse localisaiton ontology from hierachical object to list of objects with parent and children lists.
      Recursively called to assemble list.

      :param ontology_json: Parsed ontology JSON.
      :param parent: List of parent terms, used in recursive calls to build ontology.
      :return: List of dicts, each containing information about that localisation ontology term.
      """
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
          terms.update(_parse_localisation_ontology(loclist[l], parent=parent + [loclist[l]["name"]]))
      # if the root node, add a root pseudolocalisation
      if parent == ["root"]:
        terms["root"] = {
          "name": "root",
          "parent": [],
          "children": [x["name"] for x in ontology_json["localisation"]]
        }
      return terms

    import json
    zenodo_json = self._fetch_zenodo_record_json(self.master_zenodo_id)
    # load and parse to flat dict of terms with parent and children names
    return _parse_localisation_ontology(json.loads(self._fetch_zenodo_record_file(zenodo_json, "localisation_ontology.json")))

  def _localisation_match(self, cell_line, query_term: str, match_subterms: bool = True, exclude_modifiers: list = ["weak", "<10%"], required_modifiers: list = None) -> bool:
    """
    Test if the localisation annotation of `life_stage`, `gene_id` and `terminus` in `gene_list` match the query.

    :param cell line: `CellLine` object containing `life_stage`, `gene_id` and `terminus`.
    :param query_term: Search query annotation term from the localisation ontology.
    :param match_subterms: Whether to also match child/subterm/substructures of `query_term`, default `True`.
    :param exclude_modifiers: List of modifier terms none of which can be matched, default `"weak"` and `"<10%"`.
    :param required_modifiers: List of modifier terms all of which must be matched, default `None`.
    :return: Whether or not this is a localisation match.
    """
    # get query localisation
    localisations = self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["loc"]
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

  def localisation_search(self, query_term: str, life_stage: str = None, match_subterms: bool = True, exclude_modifiers: list = ["weak", "<10%"], required_modifiers: list = None) -> list:
    """
    Get a worklist of `gene_id` and `terminus` hits where any of the localisation annotations match the query.

    :param life_stage: Life cycle stage, default to first entry in `self.life_stages`.
    :param query_term: Search query annotation term from the localisation ontology.
    :param match_subterms: Whether to also match child/subterm/substructures of `query_term`, default `True`.
    :param exclude_modifiers: List of modifier terms none of which can be matched, default `"weak"` and `"<10%"`.
    :param required_modifiers: List of modifier terms all of which must be matched, default `None`.
    :return: List of `CellLine` objects of the hits, containing `life_stage`, `gene_id` and `terminus`.
    """
    # determine life stage
    if life_stage is None:
      life_stage = self.life_stages[0]
    # check all against query
    hits = []
    for cell_line in self.worklist_all(life_stage):
      if self._localisation_match(cell_line,  query_term, match_subterms=match_subterms, exclude_modifiers=exclude_modifiers, required_modifiers=required_modifiers):
        hits.append(cell_line)
    return hits

  def gene_id_search(self, gene_id_list: list, life_stage: str = None) -> list:
    """
    Use a list of gene ids to build a worklist of all tagged termini for that gene.

    :param gene_id_list: List of strings which are gene IDs to search for.
    :param query_term: Search query annotation term from the localisation ontology.
    """
    # determine life stage
    if life_stage is None:
      life_stage = self.life_stages[0]
    # build list
    hits = []
    for cell_line in self.worklist_all(life_stage):
      if cell_line.gene_id in gene_id_list:
        hits.append(cell_line)
    return hits

  def check_data_cache_usage(self, exact: bool = False):
    """
    Check disk usage of the data cache vs. free space.
    Prints information about disk usage for data cache, free space and whether the full cache will fit.

    :param exact: If `True`, then calculate usage exactly. Otherwise, approximate estimate from number of directories.
    """
    if exact:
      # full exact check, sum all file sizes in all directories
      str_prefix = ""
      sum_dir = 0
      sum_zip = 0
      for file in os.listdir(self.data_cache_path):
        if os.path.isfile(os.path.join(self.data_cache_path, file)):
          sum_zip += os.path.getsize(os.path.join(self.data_cache_path, file))
        elif os.path.isdir(os.path.join(self.data_cache_path, file)):
          for subfile in os.listdir(os.path.join(self.data_cache_path, file)):
            sum_dir += os.path.getsize(os.path.join(self.data_cache_path, file, subfile))
    else:
      # quick approximation, count directories and multiply by expected size
      str_prefix = "~"
      sum_dir = 0
      sum_zip = 0
      for file in os.listdir(self.data_cache_path):
        if file.endswith(".zip"):
          sum_zip += self._data_cache_zipsize
        if os.path.isdir(os.path.join(self.data_cache_path, file)):
          sum_dir += self._data_cache_platesize
    if self.print_status:
      print("Current data cache:")
      print("  Directory usage: "+str_prefix+str(round(sum_dir / float(2 << 40), 4))+" TiB")
      print("  Zip file usage: "+str_prefix+str(round(sum_zip / float(2 << 40), 4))+" TiB")
    # check disk space
    space_required = self._data_cache_size
    if self.remove_zip_files == False:
      space_required += self._data_cache_zipsize # add if retaining zips
    # subtract used space
    space_required -= sum_dir
    if self.remove_zip_files == False:
      space_required -= sum_zip
    total, used, free = shutil.disk_usage(self.data_cache_path)
    # warn if not enough space
    if free < space_required:
      if self.print_status: print("! Insufficient free disk space for full data cache: "+str(round(free / float(2 << 40), 4))+" / "+str(round(space_required / float(2 << 40), 4))+" TiB available !")

  def _count_cells(self, cell_line):
    """
    Counts the number of fields of view and cells in them for for the requested cell_line (`terminus`, `gene_id`, `terminus`) from the data cache.
    Records the data in `self.gene_list`, as `fields_count` `int`, `cells_per_field` `list` and `cells` `list` of `list`.

    :param cell line: `CellLine` object containing `life_stage`, `gene_id` and `terminus`.
    """
    base_path = os.path.join(self.data_cache_path, self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["plate"])
    if "cells" not in self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus] and os.path.isdir(base_path):
      cells = []
      for i in range(len(glob.glob(os.path.join(base_path, cell_line.gene_id+"_4_"+cell_line.terminus.upper()+"_*_roisCells.txt")))):
        with open(os.path.join(base_path, cell_line.gene_id+"_4_"+cell_line.terminus.upper()+"_"+str(i + 1)+"_roisCells.txt")) as cells_file:
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
      self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["fields_count"] = len(cells)
      self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["cells_per_field"] = [len(x) for x in cells]
      self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["cells"] = cells.copy()
      if self.print_status: print("  Counted", self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["fields_count"], "image data files with", sum(self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["cells_per_field"]), "cells for", cell_line.life_stage, cell_line.gene_id, cell_line.terminus)

  def fetch_data(self, cell_line):
    """
    Downloads and caches microscopy data for the plate containing a `gene_id` and `terminus` combination
    Places data in `self.data_cache_path`
    Updates self.gene_list[life_stage][gene_id][terminus] with information about number of fields of view and cells

    :param cell line: `CellLine` object containing `life_stage`, `gene_id` and `terminus`.
    :return: List of dicts of all cells for this `life_stage`, `gene_id` and `terminus` in the form `{"field_index": field_index, "cell_index": cell_index}`
    """

    def _file_md5_hash(path: str, blocksize:int = 2**20) -> str:
      """
      Calculates MD5 hash of the file at `path`.
      """
      import hashlib
      m = hashlib.md5()
      with open(path, "rb") as file:
        while True:
          buffer = file.read(blocksize)
          if not buffer:
            break
          m.update(buffer)
      return m.hexdigest()

    # check if the data cache directory exists, and make if not
    if not os.path.isdir(self.data_cache_path):
      if self.print_status: print("Making data cache directory: "+self.data_cache_path)
      os.mkdir(self.data_cache_path)
    # target paths for zip file and data subdirectory
    plate = self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["plate"]
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
        self.check_data_cache_usage()
        print("Fetching data for "+cell_line.life_stage+" form, gene ID "+cell_line.gene_id+" tagged at "+cell_line.terminus+" terminus")
        # fetch the processed microscopy data from Zenodo
        if self.print_status: print("  Making plate data directory for: "+plate)
        if "record_id" not in self.zenodo_index[plate]:
          # if not already translated, fetch the latest Zenodo ID from master Zenodo ID
          # fetch Zenodo record JSON, to get latest version doi
          zenodo_json = self._fetch_zenodo_record_json(self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["zenodo_id"])
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
            with _tqdmDownload(unit='B', unit_scale=True, unit_divisor=1024, miniters=1, desc=plate) as progress:
              urllib.request.urlretrieve(self.zenodo_index[plate]["record_url"], zip_path_temp, progress.urllib_callback)
          else:
            urllib.request.urlretrieve(self.zenodo_index[plate]["record_url"], zip_path_temp)
          if self.print_status: print("  Checking MD5 hash of: "+plate+".zip.tmp")
          zip_md5 = _file_md5_hash(zip_path_temp)
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
              count_decompressed = 0
              missing = []
              # do decompression
              for file in tqdm(archive.namelist(), unit=' files', miniters=1, desc=plate, disable=not self.print_status, smoothing=0):
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
      self._count_cells(cell_line)

  def cell_list(self, cell_line) -> list:
    """
    Returns a list fields of view and cells in them for for the requested cell_line (`terminus`, `gene_id`, `terminus`).

    :param cell line: `CellLine` object containing `life_stage`, `gene_id` and `terminus`.
    :return: List of dicts of all cells for this `life_stage`, `gene_id` and `terminus` in the form `{"field_index": field_index, "cell_index": cell_index}`
    """
    # fetch data
    self.fetch_data(cell_line)
    # return field/cell list
    fieldcell_list = []
    for field_index, cell_list in enumerate(self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["cells"]):
      for cell_index, cell in enumerate(cell_list):
        fieldcell_list.append({"field_index": field_index, "cell_index": cell_index})
    return fieldcell_list

  def check_if_cached(self, cell_line) -> bool:
    """
    Checks if data is cached for a given `gene_id` and `terminus`

    :param cell line: `CellLine` object containing `life_stage`, `gene_id` and `terminus`.
    :return: If the data is already cached
    """
    # path for data subdirectory
    plate = self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["plate"]
    dir_path = os.path.join(self.data_cache_path, plate)
    # False if MD5 not yet copied to plate data directory
    if not os.path.isfile(os.path.join(dir_path, "_"+plate+".zip.md5")):
      return False
    # otherwise cached
    return True

  def check_data_cache_integrity(self):
    """
    Checks all microscopy in the data cache, reporting any missmatches of MD5 hash per data plate.
    Any MD5 hash missmatches imply a new version of data for that plate may be available.
    Checks the number of images available per cell line for obvious missing data.
    Does not do a full check and may well throw `FileNotFoundError` if data cache is very malformed.
    Temporarily overrides `self.print_status` to `False` so only its specific output is shown.
    """
    import os
    # temporarily silence verbose output
    original_print_status = self.print_status
    self.print_status = False
    print("Checking data cache for errors")
    for plate in self.zenodo_index:
      zip_path = os.path.join(self.data_cache_path, plate+".zip")
      dir_path = os.path.join(self.data_cache_path, plate)
      # if either zip or dir exist
      if os.path.isdir(dir_path) or os.path.isfile(zip_path):
        # get md5 of source file from zenodo entry latest version
        zenodo_json = self._fetch_zenodo_record_json(self.zenodo_index[plate]["master_record_id"])
        for file in zenodo_json["files"]:
          if file["key"].endswith("_processed.zip"):
            self.zenodo_index[plate]["record_md5"] = file["checksum"].split(":")[-1]
      # first, check md5s
      if os.path.isdir(dir_path):
        # check dir saved md5 of source zip vs zenodo latest
        md5_file = os.path.join(dir_path, "_"+plate+".zip.md5")
        if not os.path.isfile(md5_file):
          # missing md5, likely partial decompression
          print("  MD5 file mising from image directory:", plate)
        else:
          with open(md5_file, "r") as file:
            # md5 missmatch check
            md5 = file.read()
            if md5 != self.zenodo_index[plate]["record_md5"]:
              print("  MD5 file in directory does not match zenodo record md5:", plate)
      if os.path.isfile(zip_path):
        # check zip saved md5 vs zenodo latest
        if self.remove_zip_files:
          # zip present, but should have been removed
          print("  Zip found which should have been removed:", plate)
        md5_file = os.path.join(zip_path+".md5")
        if not os.path.isfile(md5_file):
          # missing md5, likely incomplete zip download
          print("  MD5 file mising for zip file:", plate)
        else:
          with open(md5_file, "r") as file:
            # md5 missmatch check
            md5 = file.read()
            if md5 != self.zenodo_index[plate]["record_md5"]:
              print("  MD5 file for zip file does not match zenodo record MD5:", plate)
      # then, check if there are images for every gene in the plate
      if os.path.isdir(dir_path):
        # for all life_stage/gene_id/terminus
        for life_stage in self.life_stages:
          for cell_line in self.worklist_all(life_stage):
            if self.gene_list[life_stage][cell_line.gene_id][cell_line.terminus]["plate"] == plate:
              if "cells" in self.gene_list[life_stage][cell_line.gene_id][cell_line.terminus]:
                # ensure cells are counted
                self.cell_list(cell_line)
                if len(self.gene_list[life_stage][cell_line.gene_id][cell_line.terminus]["cells"]) == 0:
                  print("  No images found, but images expected, for", life_stage, gene_id, terminus, "in", plate)
        # restore verbose output state
    self.print_status = original_print_status

  def fetch_all_data(self):
    """
    Fetches all microscopy data and stores it in the data cache.
    """
    # for all life_stage/gene_id/terminus
    for life_stage in self.life_stages:
      for cell_line in self.worklist_all(life_stage):
        self.fetch_data(cell_line)

  def _open_field(self, cell_line, field_index: int = 0, custom_field_image = None) -> list:
    """
    Returns field of view image data.

    :param cell line: `CellLine` object containing `life_stage`, `gene_id` and `terminus`.
    :param field_index: Index of the field of view. If not set, then `0`.
    :param custom_field_image: `FieldImage` object containing custom field images to use. Images can be skimage image or `None`. Entries of None will use tryptag default. If not set or `None`, then use all tryptag defaults.
    :return: List with one `skimage` image per image channel and threshold image. List is in the order `[phase, mng, dna, phase_mask, dna_mask]`.
    """
    # check for none life_stage, replace with default
    if cell_line.life_stage is None:
      cell_line.life_stage = self.life_stages[0]
    # ensure data is fetched
    self.fetch_data(cell_line)
    # determine base path for files
    field_base_path = os.path.join(self.data_cache_path, self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["plate"], cell_line.gene_id+"_4_"+cell_line.terminus.upper()+"_"+str(field_index + 1))
    if field_base_path != self._field_base_path_sk:
      # if not the last path opened, open field and threshold images and store base path in self
      self._field_base_path_sk = field_base_path
      field_image = skimage.io.imread(self._field_base_path_sk+".tif")
      field_image = numpy.moveaxis(field_image, [0, 1, 2], [1, 2, 0]) # MAGIC NUMBER: correct loading in a weird order
      field_threshold = skimage.io.imread(self._field_base_path_sk+"_thr.tif")
      # 'clean' the field_threshold images, ensure all pixels less than 255 are set to 0
      for image in field_threshold:
        image[image < 255] = 0
      # store a copy of the images in 
      self._channels_sk = []
      for channel in field_image:
        self._channels_sk.append(channel.copy())
      self._thresholds_sk = []
      for threshold in field_threshold:
        self._thresholds_sk.append(threshold.copy())
    # setup output, copying images as downstream usage may modify
    images = []
    images.append(self._channels_sk[0].astype("uint16").copy())
    images.append(self._channels_sk[1].astype("uint32").copy())
    images.append(self._channels_sk[2].astype("uint16").copy())
    images.append(self._thresholds_sk[0].astype("uint8").copy())
    images.append(self._thresholds_sk[1].astype("uint8").copy())
    if custom_field_image is not None:
      # if not None, overwrite with custom data
      if custom_field_image.phase is not None:
        images[0] = custom_field_image.phase.copy()
      if custom_field_image.mng is not None:
        images[1] = custom_field_image.mng.copy()
      if custom_field_image.dna is not None:
        images[2] = custom_field_image.dna.copy()
      if custom_field_image.phase_mask is not None:
        images[3] = custom_field_image.phase_mask.copy()
      if custom_field_image.dna_mask is not None:
        images[4] = custom_field_image.dna_mask.copy()
    return images

  def open_field(self, cell_line, field_index: int = 0) -> list:
    """
    Opens a field of view from a `gene_id`, `terminus` and `field_index`.

    :param cell line: `CellLine` object containing `life_stage`, `gene_id` and `terminus`.
    :param field_index: Index of the field of view. If not set, then `0`.
    :return: FieldImage object, containing the image channels as attributes `phase`, `mng`, `dna`, and the thresholds `phase_mask` and `dna_mask`.
    """
    [phase, mng, dna, phase_mask, dna_mask] = self._open_field(cell_line, field_index)
    return FieldImage(
      phase=phase,
      mng=mng,
      dna=dna,
      phase_mask=phase_mask,
      dna_mask=dna_mask,
      cell_line=cell_line,
      field_index=field_index,
    )

  def _open_cell(self, cell_line, field_index: int, crop_centre: tuple, fill_centre: tuple, custom_field_image = None, rotate: bool = False, angle: float = 0, width: int = 323) -> list:
    """
    Returns cell image data, returning a list with one `skimage` image per image channel and threshold image.

    :param cell line: `CellLine` object containing `life_stage`, `gene_id` and `terminus`.
    :param field_index: Index of the field of view. If not set, then `0`.
    :param crop_centre: `(x, y)` tuple around which to crop. Ignored if `width <= 0`.
    :param fill_centre: `(x, y)` at which to do a flood fill to select the target cell object in phase_mask.
    :param custom_field_image: `FieldImage` object containing custom field images to use. Images can be skimage image or `None`. Entries of None will use tryptag default. If not set or `None`, then use all tryptag defaults.
    :param rotate: Whether or not to rotate the cell, default `False`. Ignored if `width <= 0`
    :param angle: Angle in degrees to rotate cell clockwise. Default 0.
    :param width: If positive, width of cropped cell image centred on `crop_centre`. If negative, padding around the `phase_mask`. Default, `323`.
    :return: List with one `skimage` image per image channel and threshold image. List is in the order `[phase, mng, dna, phase_mask, dna_mask, phase_mask_othercells]`.
    """

    # define image crop function
    def _skimage_crop(image, x, y, w, h):
      return image[int(y):int(y + h), int(x):int(x + w)]

    # open field, passing custom images
    channels = self._open_field(cell_line, field_index, custom_field_image = custom_field_image)
    # process phase threshold image to split cell of interest from neighbouring cells, nb. xy swapped in skimage arrays
    # replace existing mask pixels with value 127
    channels[3][channels[3] == 255] = 127
    # flood fill cell of interest to 255
    channels[3]=skimage.morphology.flood_fill(channels[3], (fill_centre[1], fill_centre[0]), 255)
    # append a copy (channels index 5), pixels if equal 127 ie. other cells
    channels.append(255*(channels[3] == 127))
    # filter channels index 3 for current cell only, set pixels equal to 127 to zero
    channels[3][channels[3] == 127] = 0

    # crop (and potentially rotate) to get cell image
    cell_channels = []
    if width < 0:
      # padding mode
      label_image = skimage.measure.label(channels[3])
      (ymin, xmin, ymax, xmax) = skimage.measure.regionprops(label_image)[0]["bbox"]
      # do actual cropping
      for channel in channels:
        # negative width is padding
        padding = -width
        # if crop outside of image bounds, then first increase canvas size
        offs = 0
        if ymin - padding < 0 or xmin - padding < 0 or ymax + padding > channel.shape[0] or xmax + padding > channel.shape[1]:
          channel = numpy.pad(channel, ((padding, padding), (padding, padding)), mode="median")
          offs = padding
        # crop
        cell_channels.append(channel[ymin - padding + offs:ymax + padding + offs, xmin - padding + offs:xmax + padding + offs])
    elif width > 0:
      # fixed width mode
      if rotate:
        width_inter = width * 1.5 # greater than width * 2**0.5
        height = round(width / 2)
      else:
        width_inter = width
      # do actual cropping
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
          channel = skimage.transform.rotate(channel, -angle, preserve_range=True).astype(channel_dtype)
          channel = _skimage_crop(channel, half_width_inter - width / 2, half_width_inter - height / 2, width, height)
        cell_channels.append(channel)
    else:
      raise ValueError("`width` must be a nonzero integer, positive for fixed image width in pixels, negative for padding around phase_mask in pixels.")
    # downstream analysis (including tryptools) allowed to assume cell mask does not touch image edge, therefore set border pixels to 0 (may clip large cells)
    cell_channels[3][:, [0, -1]] = 0
    cell_channels[3][[0, -1], :] = 0
    return cell_channels

  # open a cell, cropped from a field of view
  # uses the phase and dna threshold images from tryptag
  # cell x, y coordinate in phase threshold from tryptag
  def open_cell(self, cell_line, field_index: int = 0, cell_index: int = 0, width: int = 323, rotate: bool = False) -> list:
    """
    Opens a cell from a `gene_id`, `terminus`, `field_index` and `cell_index`.
    Use `open_cell_custom` if you would like to inject custom images, coordinates and/or angle.

    :param cell line: `CellLine` object containing `life_stage`, `gene_id` and `terminus`.
    :param field_index: Index of the field of view. If not set, then `0`.
    :param cell_index: Index of the cell in the field of view. If not set, then `0`.
    :param width: If positive, width of cropped cell image (may clip very large cells). If negative, padding for crop around the `phase_mask`. Default, `323`.
    :param rotate: Whether or not to rotate the cell. Default `False`. Set to `False` if `width < 0` (padded crop mode).
    :return: CellImage object, containing the image channels as attributes `phase`, `mng`, `dna`, and the thresholds `phase_mask`, `dna_mask` and `phase_mask_othercells`.
    """
    if cell_line.life_stage is None:
      cell_line.life_stage = self.life_stages[0]
    self.fetch_data(cell_line)
    cell_data = self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["cells"][field_index][cell_index]
    crop_centre = cell_data["centre"]
    fill_centre = cell_data["wand"]
    angle = cell_data["angle"]
    if width < 0:
      rotate = False
    phase, mng, dna, phase_mask, dna_mask, phase_mask_othercells = self._open_cell(cell_line, field_index, crop_centre, fill_centre, angle = angle, rotate = rotate, width = width)
    return CellImage(
      phase=phase,
      mng=mng,
      dna=dna,
      phase_mask=phase_mask,
      dna_mask=dna_mask,
      phase_mask_othercells=phase_mask_othercells,
      cell_line=cell_line,
      field_index=field_index,
      cell_index=cell_index,
      rotated=rotate,
    )

  def open_cell_custom(self, cell_line, field_index: int = 0, cell_index: int = None, custom_field_image = None, fill_centre: tuple = None, crop_centre: tuple = None, rotate: bool = False, angle: float = None, width: int = 323) -> list:
    """
    Advanced customisable open cell, opens a cell from a `gene_id`, `terminus` and `field_index`, but with customisable images, cell coordinates and/or angle.
    This allows use of custom pth image and cell coordinates, default cell coordinates but prefiltered mng image, etc.

    :param cell line: `CellLine` object containing `life_stage`, `gene_id` and `terminus`.
    :param field_index: Index of the field of view. If not set, then `0`.
    :param cell_index: Index of the cell in the field of view, used for finding default `fill_centre` and `crop_centre` coordinates. If not set or `None`, then use `fill_centre` coordinate.
    :param custom_field_image: `FieldImage` object containing custom field images to use. Images can be skimage image or `None`. Entries of None will use tryptag default. If not set or `None`, then use all tryptag defaults.
    :param fill_centre: `(x, y)` tuple of a pixel which is in the target cell object (pixel value 255) in pth image.
    :param crop_centre: `(x, y)` tuple around which to crop, otherwise crop around `fill_centre`.
    :param rotate: Whether or not to rotate the cell. Default `False`. Set to `False` if `width < 0` (padded crop mode).
    :param angle: Angle in degrees to rotate cell clockwise. If not set or `None`, tryptag default.
    :param width: If positive, width of cropped cell image (may clip very large cells). If negative, padding for crop around the `phase_mask`. Default, `323`.
    :return: List with one `skimage` image per image channel and threshold image. List is in the order `[phase_(gray), mng_(green), dna_(blue), phase_threshold, dna_threshold]`, often referred to as `[pth, mng, dna, pth, dth]`.
    """
    # if cell_index is set, then use tryptag defaults unless overridden
    if cell_index is not None:
      cell_data = self.gene_list[cell_line.life_stage][cell_line.gene_id][cell_line.terminus]["cells"][field_index][cell_index]
      if crop_centre is None:
        crop_centre = cell_data["centre"]
      if fill_centre is None:
        fill_centre = cell_data["wand"]
      if angle is None:
        angle = cell_data["angle"]
    else:
      # crop_centre must be set, if not then throw error
      if fill_centre is None:
        raise ValueError("If `cell_index` is not set or `None` then `fill_centre` must be set and not be `None`")
      # if custom fill_centre, check for crop_centre otherwise default to fill_centre
      if crop_centre is None:
        crop_centre = fill_centre
    if width < 0:
      rotate = False
    [phase, mng, dna, phase_mask, dna_mask, phase_mask_othercells] = self._open_cell(cell_line, field_index, crop_centre, fill_centre, custom_field_image, angle = angle, rotate = rotate, width = width)
    return CellImage(
      phase=phase,
      mng=mng,
      dna=dna,
      phase_mask=phase_mask,
      dna_mask=dna_mask,
      phase_mask_othercells=phase_mask_othercells,
      cell_line=cell_line,
      field_index=field_index,
      cell_index=None,
      rotated=rotate,
    )

  def _list_analysis_worker(self, cell_line: CellLine, analysis_function: callable, threading_mode: bool) -> list:
    """
    Worker for multiprocess/thread parallel analysis of a `work_list`.

    :param cell_line: `CellLine` object defining the `life_stage`, `gene_id` and `terminus` combination to analyse.
    :param analysis_function: Function name to use for analysis. `analysis_function` should take exactly two arguments, `tryptag` (`TrypTag` instance) and `cell_line` (`CellLine` object) in this order.
    :param threading_mode: True if running in threading mode. This triggers a deep copy of `self` to avoid thread safety issues.
    :return: List of dicts in the form `{"life_stage": life_stage, "gene_id": gene_id, "terminus": terminus, "result": analysis_function_return}`.
    """
    # Deep copy tryptag object if running in threading mode (to avoid thread safety issues)
    if threading_mode:
      from copy import deepcopy
      self = deepcopy(self)

    if self.print_status: print(f"  Starting to process cell line { cell_line }")

    self.fetch_data(cell_line)
    result = {
      "cell_line": cell_line,
      "result": analysis_function(self, cell_line),
    }
    return result

  def analyse_list(self, work_list, analysis_function, workers=None, multiprocess_mode="process"):
    """
    Simple handler for parallel analysis of a `work_list`.

    :param work_list: List of CellLine objects, defining each `life_stage, `gene_id` and `terminus` combination to analyse.
    :param analysis_function: Function name to use for analysis. `analysis_function` should take exactly two arguments, `tryptag` (`TrypTag` instance) and `cell_line` (`CellLine` object) in this order.
    :param workers: Number of threads/processes to spawn, default is number of CPUs.
    :param multiprocess_mode: `"process"` for parallel processes, `"thread"` for parallel threads or `None` for no parallel processing (directly calls `analysis_function`).
    :return: List of dicts in the form `{"life_stage": life_stage, "gene_id": gene_id, "terminus": terminus, "result": analysis_function_return}`. These may be in a different order to `work_list`.
    """
    import concurrent.futures
    import multiprocessing
    import numpy
    
    if self.print_status: print("Analysing worklist")

    # deduplicate work_list
    dedup_work_list = set(work_list)
    for entry in dedup_work_list:
      if entry.life_stage is None:
        entry.life_stage = self.life_stages[0]

    # get number of workers, default to number of cpus
    if workers is None:
      workers = multiprocessing.cpu_count()

    if multiprocess_mode is None:
      # run in a single thread, still use the ThreadPoolExecutor since that's equivalent
      if self.print_status: print("  Single process")
      Executor = concurrent.futures.ThreadPoolExecutor
      workers = 1
    elif multiprocess_mode == "process":
      # setup executor as a process pool
      if self.print_status: print("  Parallel processes with", workers, "workers")
      Executor = concurrent.futures.ProcessPoolExecutor
    elif multiprocess_mode == "thread":
      # setup executor as a thread pool
      if self.print_status: print("  Parallel threads with", workers, "workers")
      Executor = concurrent.futures.ThreadPoolExecutor
    else:
      raise ValueError(f"Unknown multiprocess_mode '{multiprocess_mode}")

    # trigger fetch of gene_list and zenodo_index, prior to copying for each thread
    self.zenodo_index
    self.gene_list
    with Executor(workers) as executor:
      futures = [executor.submit(self._list_analysis_worker, cell_line=cell_line, analysis_function=analysis_function, threading_mode=multiprocess_mode=="thread") for cell_line in dedup_work_list]
      results = [future.result() for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), smoothing=0)]
    return results

class BSFTag(TrypTag):
  def __init__(self, *args, **kwargs):
    if 'master_zenodo_id' in kwargs:
      raise ValueError('`master_zenodo_id` is not a valid argument')
    super().__init__(self, *args, **kwargs, master_zenodo_id=7258722, data_cache_path="./_bsftag_cache", data_cache_plates=4, life_stages=["bloodstream"])
