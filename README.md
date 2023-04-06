# tryptag
`tryptag` is a python module for accessing and handling [TrypTag](http://tryptag.org) genome-wide protein localisation project data.
Its primary intended use is for easy access to image data for automated image analysis.

## Installation
First, make sure `python`, `git` and `pip` are installed, on linux/mac:

```shell
sudo apt install python3 python3-pip git
```

Or, on windows:

```bash
winget install Python.Python.3.0; winget install Git.Git
python3 -m pip install -update pip
```

Next, install `tryptag` using `pip`. This requires `git`:

```shell
pip install git+https://github.com/zephyris/tryptag
```

`tryptag` requires several python modules: `numpy` `scikit-image`, `progressbar2` and `filelock`. These are automatically installed when using `pip`.

To reinstall and upgrade use `pip`:
```shell
pip install --upgrade --force-reinstall git+https://github.com/zephyris/tryptag
```

To uninstall also use `pip`:

```shell
pip uninstall tryptag
```

## Quickstart guide

To use the `tryptag` module, import the `TrypTag` class and set up an instance (normally called `tryptag`):

```python
from tryptag import TrypTag
tryptag = TrypTag()
```

Microscopy data is multiple fields of view per cell line.
It can be accessed using instances of `CellLine`, a simple class defining cell line life cycle stage, gene id (as used on [TriTrypDB](http://tritrypdb.org)), tagging terminus (`n` or `c`). There are multiple fields of view, accessed by field_index:

```python
from tryptag import CellLine
cell_line = CellLine("Tb927.9.8570", "n")
field_index = 2
field_image = tryptag.open_field(cell_line, field_index)
```

This returns an instance of a `FieldImage` object, containing the phase, mNG, DNA stain, phase mask and DNA mask images.

The cells in the phase threshold image are indexed and can be opened individually. To open a specific cell in the field of view:

```python
cell_index = 14
cell_image = tryptag.open_cell(cell_line, field_index, cell_index)
```

Similar to `open_field`, `open_cell` returns a `CellImage` object.

Images within a `FieldImage` or `CellImage` object can be accessed using dot notation: `field_image.phase`, `.mng`, `.dna`, `.phase_mask` and `.dna_mask`.
All images are `numpy` `ndarray` objects, as used by `scikit-image`:

```python
from skimage import io
io.imshow(cell_image.phase)
io.show()
```

Bear in mind that accessing a nonexistant gene id, tagging terminus, field or cell will give `KeyError` errors. For example:

```python
for field_index in range(7):
    try:
        field_image = tryptag.open_field(cell_line, field_index)
        # do your analysis here
    except:
        print("Field not found, field_index:", field_index)
```

## Full guide

### Localisation searches

`tryptag` understands the localisation annotation ontology and provides a tool for intelligent localisation searches. First setup `tryptag`:

```python
from tryptag import TrypTag
tryptag = TrypTag()
```

You can search by any of the localisation annotation terms:

```python
results = tryptag.localisation_search("nucleoplasm")
```

This returns a list of `CellLine` objects including gene id and terminus, which you can access using dot notation: `cell_line.gene_id`, `.terminus`.

### Opening images

The primary intended use of the `tryptag` module is for easy access of specific field of view and cell images for automated image analysis. See quickstart.

`open_field` and `open_cell` return a `FieldImage` or `CellImage` object containing these images which can be accessed using the appropriate dot notation. Microscopy data is in three image channels and two thresholded images:

Image channels:

1. `.phase` Phase contrast (transmitted light, overall cell morphology) uint16
2. `.mng` mNG fluorescence (green fluorescence, from the tagged protein) uint32
3. `.dna` DNA stain fluorescence (blue fluorescence, using Hoechst 33342) uint16

Thresholded images:

1. `.phase_mask` Thresholded phase contrast (cells) uint8, 255 = object, current cell of interest
2. `.dna_mask`Thresholded DNA stain (nuclei and kinetoplasts - mitochondrial DNA organelles) uint8, 255 = object, kinetoplasts or nuclei

`CellImage` objects additionally contain `.phase_mask_othercells` which is a mask of every _other_ cell in the view.

### Image analysis

`tryptag` includes `tryptools` which provides some useful tools for image analysis of _Trypanosoma brucei_ cells. First import `TrypTag` and `tryptools` and set up `tryptag`.

```python
from tryptag import TrypTag, tryptools
tryptag = TrypTag()
```

The `tryptools` methods take a `CellImage` object as an input and return various automated image analysis data.

```python
cell_image = tryptag.open_cell(CellLine(life_stage="procyclic", gene_id="Tb927.9.8570", terminus="n"))
morphology_result = tryptools.cell_morphometry_analysis(cell_image)
```

### High throughput

`tryptag` makes it easy to apply an analysis to many cell lines. First import and set up `tryptag`:

```python
from tryptag import TrypTag, tryptools
tryptag = TrypTag()
```

Define your analysis function you'd like to apply to each cell line. This example analyses mNG signal in each individual cell:

```python
def analysis_function(tryptag, cell_line):
    result = []
    fieldcell_list = tryptag.cell_list(cell_line)
    for entry in fieldcell_list:
        tryptools.cell_signal_analysis(tryptag.open_cell(cell_line, entry["field_index"], entry["cell_index"]))
    return result
```

Run the analysis using the built-in multiprocess analysis tool. This example applies this function to _all_ cell lines. Automated iteration through the entire ~5,000,000 cell dataset:

```python
worklist = tryptag.worklist_all()
results = tryptag.analyse_list(worklist, analysis_function)
```

You can also use the output of `tryptag.localisation_search` for `worklist`, or `tryptag.worklist_parental` for data from untagged parental cells.

### Microscopy data

You can use `tryptag` to download the microscopy data. In python, import the module and set up a `TrypTag` instance:

```python
from tryptag import TrypTag, CellLine
tryptag = TrypTag()
```

You can trigger download of microscopy data for a specific gene id and tagging terminus using `fetch_data`.
This looks up in which tagging plate correspond to the most recent replicate of this life cycle stage, gene ID and terminus tagging attempt, and the URL at which to find this data.
It then downloads and decompresses the data to the `data_cache_path` directory.

```python
tryptag.fetch_data(CellLine(life_stage="procyclic", gene_id="Tb927.7.1920", terminus="n"))
```

This will take a long time; to get image data for a single gene the data for an entire plate needs to be downloaded. This is typically ~10 to 20 Gb.

Look through the data cache directory and you will find the microscopy data, in one subdirectory per tagging plate and named by gene id and tagging terminus.

### General tips

Make sure you have enough free disk space to download, decompress and cache the image data. This is up to ~40 Gb for a single plate and ~4 Tb for the entire dataset.
The default cache location is `_tryptag_data` within the current working directory. You can change this to a relative or absolute path to any directory you wish - we recommend a scratch drive with sufficient space.
Make sure this is set at the start of every script:

Linux/Mac:
```python
from tryptag import TrypTag
tryptag = TrypTag(data_cache_path = "\mnt\z\my\scratch\directory")
```
or Windows:
```python
tryptag = TrypTag(data_cache_path = "Z:/my/scratch/directory")
```

Do not delete or move files from `data_cache_path`. `tryptag` does not check the plate subdirectories for integrity. You can, however, safely delete a plate subdirectory.
Interrupt of either data download or zip decompression _should_ behave gracefully, leaving partial data but not preventing later automatic re-download and/or re-decompression.

If multiple scripts using the same `data_cache_path` simultaneously try to download a plate it _should_ be handled gracefully.
One script should download and decompress the image data, while the others (silently) wait until it the image data is available.
However, for large-scale analyses, it is more robust to ensure all data is already cached. You can easily download all image data (this will probably take more than one week!):

```python
tryptag.fetch_all_data()
```

The TrypTag data may have minor errors which will be corrected over time. `fetch_all_data` always fetches the latest localisation listing from [Zenodo](https://zenodo.org/record/6862289).
Cached image data may be an older version. `tryptag` records the MD5 hash of the source zip files. If the data source (Zenodo depositions) are updated, the MD5 hash will change.
Cached data inconsistent MD5 hashes can be checked and reported (but currently not corrected) using:

```python
tryptag.check_data_cache_integrity()
```

`tryptag` gives quite verbose information about what it is currently doing to fetch data. To silence this output:

```python
from tryptag import TrypTag
tryptag = TrypTag(verbose=False)
```

Internally, most tryptag data is held in a dict of dicts variable called `gene_list`, which you can explore and access directly for advanced usage.
This gets populated with information about number of fields of view, cell locations, etc. when a method like `cell_list`, `open_cell` or `open_field` requests microscopy data.

```python
from tryptag import TrypTag, CellLine
tryptag = TrypTag()
life_stage, gene_id, terminus = "procyclic", "Tb927.9.8570", "n"
localisation = tryptag.gene_list[life_stage][gene_id][terminus]["loc"]
tryptag.cell_list(CellLine(life_stage=life_stage, gene_id=gene_id, terminus=terminus))
cell_information = tryptag.gene_list[life_stage][gene_id][terminus]["cells"]
```

# Citing

If you use the TrypTag data resource, please cite Billington _et al._ 2023 Nature Microbiology [doi:10.1038/s41564-022-01295-6](https://doi.org/10.1038/s41564-022-01295-6).
We recommend including this citation in the results or methods if TrypTag was used as part of a discovery process.
If directly using TrypTag images, please also indicate in the figure legend or similar which images are from TrypTag.

If you use the `tryptag` module to access or analyse TrypTag data, please also cite this [Github repository](https://github.com/zephyris/tryptag) and the master TrypTag Zenodo deposition [doi:10.5281/zenodo.6862289](https://doi.org/10.5281/zenodo.6862289).

You may also find the following papers of use:
1. Dean _et al._ 2019 Trends Parasitol. [doi:10.1016/j.pt.2016.10.009](https://doi.org/10.1016/j.pt.2016.10.009) Project announcement, with original aims and experimental strategy.
2. Halliday _et al._ 2019 Mol. Biochem. Parasitol. [doi:10.1016/j.molbiopara.2018.12.003](https://doi.org/10.1016/j.molbiopara.2018.12.003) Describes the localisation ontology, with landmark protein examples.
