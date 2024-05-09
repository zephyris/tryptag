import numpy
import scipy
import skimage
import math

from matplotlib import patches


def cell_signal_analysis(cell_image) -> dict:
  """
  Simple `mng` signal intensity within a cell using `phase_mask`.
  Statistics from median background-corrected `mng` signal.

  :param cell_image: `CellImage` object.
  :return: A dict of cell `mng` signal information.
  """
  mng = cell_image.mng - numpy.median(cell_image.mng)
  lab = skimage.measure.label(cell_image.phase_mask)
  props_table = skimage.measure.regionprops_table(lab, mng, properties=("area", "intensity_max", "intensity_mean"))
  return {
    "cell_area": props_table["area"][0],
    "mng_mean": props_table["intensity_mean"][0],
    "mng_sum": props_table["intensity_mean"][0] * props_table["area"][0],
    "mng_max": props_table["intensity_max"][0]
  }

def _mask_pruned_skeleton(thr, prefilter_radius, prune_length):
  """
  Skeletonisation with prefilter and pruning.

  :param cell_image: An `skimage` image / `numpy` `ndarray` of the mask for skeletonisation.
  :param prefilter_radius: Gaussian blur radius for pre-filtering.
  :param prune_length: Shortest branch length to retain, will remove the entire skeleton if a single line shorter than this.
  :return: An `skimage` image / `numpy` `ndarray` of the skeleton and an ndarray for the distance transform of the thresholded image, which is used for calculating the corrected cell length.
  """
  # pre-filter with gaussian blur
  thr_fil = skimage.filters.gaussian(thr, sigma=prefilter_radius) > 0.5
  # make skeleton
  skeleton = skimage.morphology.skeletonize(thr_fil)
  skeleton = skeleton.astype(numpy.uint8)
  # make a neighbour count skeleton, 1 = terminus, 2 = arm, >2 = branch point
  neighbours = scipy.ndimage.convolve(skeleton, [[1, 1, 1], [1, 0, 1], [1, 1, 1]]) * skeleton
  # filter for 1 neigbour only, ie terminus image, and use to list termini
  termini = neighbours.copy()
  termini[termini > 1] = 0
  termini_y, termini_x = skimage.morphology.local_maxima(termini, indices=True, allow_borders=False)
  # prune skeleton
  for t in range(len(termini_x)):
    length = 0
    cx, cy = termini_x[t], termini_y[t]
    v = neighbours[cy, cx]
    while length < prune_length + 2 and v > 0 and v < 3:
      v = 0
      # mark visited pixels with 2, if removeable (not a branch)
      if neighbours[cy, cx] < 3:
        skeleton[cy, cx] = 2
      # for all neighbours...
      for a in range(-1, 2):
        for b in range(-1, 2):
          # if a skeleton pixel, step in that direction
          if (a != 0 or b != 0) and skeleton[cy + b, cx + a] == 1:
            length += 1
            v = neighbours[cy, cx]
            cy += b
            cx += a
            # break inner loop on match
            break
        # break outer loop with inner
        else:
          continue
        break
    # if short enough then prune by replacing visited pixels (2) with 0
    if length < prune_length:
      skeleton[skeleton == 2] = 0
    else:
      skeleton[skeleton == 2] = 1
  # reskeletonise, to handle messy branch points left over
  skeleton = skimage.morphology.medial_axis(skeleton, return_distance=False).astype(numpy.uint8)
  # use thresholded image to get distance of medial axis to edge of the mask
  _, distances = skimage.morphology.medial_axis(thr, return_distance=True)
  return skeleton, distances

def cell_kn_analysis(cell_image, min_area=17, kn_threshold_area=250):
  """
  Classifies and measures DNA signal in a trypanosome kinetoplast and nucleus from `phase_mask`, `dna_mask` and `dna` images.
  Returns centroid, area, sum median background-corrected `dna` signal and the minor and major axis length of an ellipse fitted to the dna_mask, as well as the orientation of this ellipse.
  Gives particularly informative anterior-posterior morphometry when a clean midline is found.
  Largely based on: doi:10.1186/1741-7007-10-1

  :param cell_image: `CellImage` object.
  :return: A dict of cell K/N information.
  """
  # background correct dna using median
  dna = cell_image.dna - numpy.median(cell_image.dna)
  # label objects and measure signal intensity and location
  dna_lab = skimage.measure.label(cell_image.dna_mask)
  pth_props_table = skimage.measure.regionprops_table(dna_lab, cell_image.phase_mask, properties=("intensity_max", "area_convex"))
  dna_props_table = skimage.measure.regionprops_table(dna_lab, dna, properties=("intensity_mean", "intensity_max", "centroid_weighted", "minor_axis_length", "major_axis_length", "orientation"))
  dna_objects = []
  # filter dna objects
  for i in range(0, dna_lab.max()):
    # if labelled dth object overlaps cell object in pth
    if pth_props_table["intensity_max"][i] == 255 and pth_props_table["area_convex"][i] > min_area: # MAGIG NUMBER: Minimum kinetoplast area
      # get stats
      dna_objects.append({
        "centroid": {
          "x": dna_props_table["centroid_weighted-0"][i],
          "y": dna_props_table["centroid_weighted-1"][i]
        },
        "area": pth_props_table["area_convex"][i],
        "dna_sum": pth_props_table["area_convex"][i] * dna_props_table["intensity_mean"][i],
        "dna_max": dna_props_table["intensity_max"][i],
        "major_axis_length": dna_props_table["major_axis_length"][i],
        "minor_axis_length": dna_props_table["minor_axis_length"][i],
        "orientation": -dna_props_table["orientation"][i] * (180 / numpy.pi) + 90, # conversion from skimage to matplotlib
      })
  # classify dna objects as k/n
  # sort by area, classify smallest ceil(count / 2) as k
  #   ie. k = n for even, k = n + 1 for odd
  dna_objects.sort(key=lambda x: x["area"])
  count_k = 0
  for o in range(math.ceil(len(dna_objects) / 2)):
    # unless too large
    if dna_objects[o]["area"] < kn_threshold_area: # MAGIC NUMBER: Maximum area for kinetoplast
      dna_objects[o]["type"] = "k"
      count_k += 1
  for object in dna_objects:
    if "type" not in object:
      object["type"] = "n"
  count_n = len(dna_objects) - count_k
  count_kn = str(count_k)+"K"+str(count_n)+"N"
  # return XKXN string, list of kinetoplast objects and list of nucleus objects
  return {
    "count_kn": count_kn,
    "count_k": count_k,
    "count_n": count_n,
    "objects_k": [x for x in dna_objects if x["type"] == "k"],
    "objects_n": [x for x in dna_objects if x["type"] == "n"]
  }

def cell_midline_analysis(cell_image, prefilter_radius=2, min_length_pruning=15):
  """
  Analyse a trypanosome cell shape from `phase_mask` by skeletonisation.

  :param cell_image: `CellImage` object.
  :return: A dict of cell skeleton information.
  """
  pth_skeleton, distances = _mask_pruned_skeleton(cell_image.phase_mask, prefilter_radius, min_length_pruning) # MAGIC NUMBERS: Radius for prefiltering, length for pruning branches
  neighbours = scipy.ndimage.convolve(pth_skeleton, [[1, 1, 1], [1, 0, 1], [1, 1, 1]]) * pth_skeleton
  termini_count = numpy.count_nonzero(neighbours == 1)
  midline_count = numpy.count_nonzero(neighbours == 2)
  branches_count = numpy.count_nonzero(neighbours > 2)
  morphology = {
    "termini": termini_count,
    "midlines": midline_count,
    "branches": branches_count,
  }
  # trace, if a single line (two termini, zero branches)
  if termini_count == 2 and branches_count == 0:
    termini = neighbours.copy()
    termini[termini > 1] = 0
    termini_y, termini_x = skimage.morphology.local_maxima(termini, indices=True, allow_borders=False)
    # trace from index 0
    midline = [[termini_y[0], termini_x[0]]]
    v = pth_skeleton[midline[-1][0], midline[-1][1]]
    while v > 0:
      v = 0
      # mark visited pixels by setting to 0
      pth_skeleton[midline[-1][0], midline[-1][1]] = 0
      # for all neighbours...
      for a in range(-1, 2): # a is delta in x
        for b in range(-1, 2): # b is delta in y
          # if a skeleton pixel, step in that direction
          if pth_skeleton[midline[-1][0] + b, midline[-1][1] + a] == 1:
            midline.append([midline[-1][0] + b, midline[-1][1] + a])
            v = pth_skeleton[midline[-1][0], midline[-1][1]]
            # break inner loop on match
            break
        # break outer loop with inner
        else:
          continue
        break
    morphology.update({
      "midline": midline,
      "midline_pixels": len(midline),
      "mat_distances": distances
    })
  return morphology

def cell_morphology_analysis(cell_image) -> dict:
  """
  Analyses a trypanosome cell morphology from `phase_mask`, `dna_mask` and `dna` images.
  Combines `cell_midline_analysis` and `kn_analysis` to give a morphometric analysis.
  Gives particularly informative anterior-posterior morphometry when a clean midline is found.
  Largely based on: doi:10.1186/1741-7007-10-1

  :param cell_image: `CellImage` object.
  :return: A dict of cell morphology information.
  """
  midline_analysis = cell_midline_analysis(cell_image)
  kn_analysis = cell_kn_analysis(cell_image)
  dna_objects = kn_analysis["objects_k"] + kn_analysis["objects_n"]
  # get k/n positions along midline, if a single midline identified in midline_analysis
  if "midline" in midline_analysis:
    # determine nearest midline point to dna centroid, ie. position along midline, and sort
    for object in dna_objects:
      object["midline_index"] = scipy.spatial.distance.cdist([[object["centroid"]["x"], object["centroid"]["y"]]], midline_analysis["midline"]).argmin()
    dna_objects.sort(key=lambda x: x["midline_index"])
    # split back to kinetoplast and nucleus, update kn_analysis
    kn_analysis.update({
      "objects_k": [x for x in dna_objects if x["type"] == "k"],
      "objects_n": [x for x in dna_objects if x["type"] == "n"]
    })
    # if at least one kinetoplast, orient cell with anterior at the most terminus-proximal kinetoplast
    if len(kn_analysis["objects_k"]) > 0:
      # check orientation using position of kinetoplasts along cell midline from both ends
      is_reversed = False
      min_k_offs = len(midline_analysis["midline"])
      for object in kn_analysis["objects_k"]:
        k_offs_1 = object["midline_index"]
        k_offs_2 = len(midline_analysis["midline"]) - object["midline_index"]
        # if closest k to cell end, so far
        if k_offs_1 < min_k_offs or k_offs_2 < min_k_offs:
          # if closest to end (rather than start) of midline, then reversed
          if k_offs_2 < k_offs_1:
            is_reversed = True
            min_k_offs = k_offs_2
          else:
            min_k_offs = k_offs_1
      if is_reversed == True:
        # reverse midline
        midline_analysis["midline"].reverse()
        # recalculate midline indexes
        for object in dna_objects:
          object["midline_index"] = len(midline_analysis["midline"]) - object["midline_index"] - 1
        # reverse DNA objects
        dna_objects.reverse()
        kn_analysis["objects_k"].reverse()
        kn_analysis["objects_n"].reverse()
      # add cell anterior and posterior coordinates to morphology object
      midline_analysis.update({
        "anterior": midline_analysis["midline"][-1],
        "posterior": midline_analysis["midline"][0]
      })
      # add ordered kn representation
      kn_analysis.update({
        "kn_ordered": "".join([x["type"].upper() for x in dna_objects])
      })
      # calculate distance along midline
      distance = [0]
      root2 = 2 ** 0.5
      for i in range(1, len(midline_analysis["midline"])):
        # step length 1 for orthogonal adjacent, root2 for diagonal adjacent
        if abs(midline_analysis["midline"][i][0] - midline_analysis["midline"][i-1][0]) == 1 and abs(midline_analysis["midline"][i][1] - midline_analysis["midline"][i-1][1]) == 1:
          distance.append(distance[-1] + root2)
        else:
          distance.append(distance[-1] + 1)

      distances = midline_analysis["mat_distances"]
      start = midline_analysis["anterior"]
      end = midline_analysis["posterior"]
      distance_ant = distances[start[0], start[1]]
      distance_post = distances[end[0], end[1]]

      midline_analysis.update({
        "distance": distance,
        "length": distance[-1] + distance_ant + distance_post
      })
  midline_analysis.update(kn_analysis)
  return midline_analysis

def plot_morphology_analysis(cell_image, result):
  """
  Makes a `matplotlib` `pyplot` summarising morphometric analysis from `cell_kn_analysis`, `cell_midline_analysis` or `cell_morphology_analysis`.
  The returned `pyplot` object can be displayed with `.show()` method, saved with `.savefig()`, etc.

  :param cell_image: `CellImage` object.
  :param result: Dict result from `cell_kn_analysis`, `cell_midline_analysis` or `cell_morphology_analysis`.
  :return: A `matplotlib` `pyplot` object.
  """
  from skimage import data, io
  from matplotlib import cm
  from matplotlib import pyplot as plt
  # use the dna image masked with phase_mask and transparently masked with dna_mask as the background
  im = plt.imshow(cell_image.dna * (cell_image.phase_mask / 255) * (0.5 + (cell_image.dna_mask / 512)), cmap=cm.gray)
  # kinetoplast positions (cyan)
  if "objects_k" in result:
    for obj in result["objects_k"]:
      ellipse = patches.Ellipse((obj["centroid"]["y"], obj["centroid"]["x"]),
                                width=obj["major_axis_length"],
                                height=obj["minor_axis_length"],
                                angle=obj["orientation"],
                                facecolor=(0.0, 0.4, 0.8, 0.5))
      plt.gca().add_patch(ellipse)
  # nucleus positions (green)
  if "objects_n" in result:
    for obj in result["objects_n"]:
      ellipse = patches.Ellipse((obj["centroid"]["y"], obj["centroid"]["x"]),
                                width=obj["major_axis_length"],
                                height=obj["minor_axis_length"],
                                angle=obj["orientation"],
                                facecolor=(0.0, 0.8, 0.4, 0.5))

    plt.gca().add_patch(ellipse)
  # cell midline (gray)
  if "midline" in result:
    plt.plot([a[1] for a in result["midline"]], [a[0] for a in result["midline"]], color=(0.8, 0.8, 0.8, 0.2))
  # anterior and posterior points (magenta, yellow)
  if "anterior" in result:
    plt.scatter([result["anterior"][1]], [result["anterior"][0]], s=10, marker="o", color=(1.0, 0.0, 1.0, 0.5))
    plt.scatter([result["posterior"][1]], [result["posterior"][0]], s=10, marker="o", color=(1.0, 1.0, 0.0, 0.5))
  # nearest point on midline to centroid of kinetoplasts and nuclei (dark blue, green)
  if "midline" in result:
    if "objects_k" in result:
      plt.scatter([result["midline"][a["midline_index"]][1] for a in result["objects_k"]], [result["midline"][a["midline_index"]][0] for a in result["objects_k"]], s=10, marker="o", color=(0.0, 0.0, 0.8, 0.8))
    if "objects_n" in result:
      plt.scatter([result["midline"][a["midline_index"]][1] for a in result["objects_n"]], [result["midline"][a["midline_index"]][0] for a in result["objects_n"]], s=10, marker="o", color=(0.0, 0.8, 0.0, 0.8))
  return plt
