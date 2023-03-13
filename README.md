# tryptag
`tryptag` is a python module for accessing and handling [TrypTag](http://tryptag.org) genome-wide protein localisation project data.
Its primary intended use is for easy access to image data for automated image analysis.

## Installation
Install using `pip`:

```
pip install git+https://github.com/zephyris/tryptag
```

To uninstall also use `pip`:

```
pip uninstall tryptag
```

`tryptag` requires several python modules: `numpy` `scikit-image`, `progressbar2` and `filelock`. These should be installed using `pip`.

## Quickstart guide

To use the `tryptag` module, import the `TrypTag` class and set up an instance (normally called `tryptag`):

```
from tryptag import TrypTag
tryptag = TrypTag()
```

Microscopy data is multiple fields of view per cell line.
To open a specific field of view of a tagged cell line, access it by gene id (as used on [TriTrypDB](http://tritrypdb.org)), tagging terminus (`n` or `c`) and field index (integer):

```
[phase, mng, dna, phasethreshold, dnathreshold] = tryptag.open_field(gene_id, terminus, field)
```

The cells in the phase threshold image are indexed and can be opened individually.
To open a specific cell in a field:

```
[phase, mng, dna, phasethreshold, dnathreshold] = tryptag.open_cell(gene_id, terminus, field, cell)
```

Opened images are `numpy` `ndarray` objects, as used by `scikit-image`.

Bear in mind that accessing a nonexistant gene id, tagging terminus, field or cell will give `KeyError` errors. For example:

```
gene_id = "Tb927.7.1920"
terminus = "n"
tryptag.fetch_data(gene_id, terminus)
for field in range(5):
    try:
        [pha, mng, dna, pth, dth] = tryptag.open_field(gene_id, terminus, field)
        # do your analysis here
    except:
        print("Field not found")
```

You can access the `tryptag.gene_list` dict for more intelligent iteration. To iterate through all cells for all gene ids and termini that exist, ie. automated analysis of the entire ~5,000,000 cell dataset:

```
tryptag.fetch_gene_list()
termini = ["n", "c"]
for gene_id in tryptag.gene_list:
    for terminus in termini:
        if terminus in tryptag.gene_list[gene_id]:
            tryptag.fetch_data(gene_id, terminus)
            for field in range(len(tryptag.gene_list[gene_id][terminus]["cells"])):
                for cell in range(len(tryptag.gene_list[gene_id][terminus]["cells"][field])):
                    [pha, mng, dna, pth, dth] = tryptag.open_cell(gene_id, terminus, field, cell)
                    # do your analysis here
```

## Full guide

### Tagging metadata
You can use `tryptag` to access data about tagged cell lines. In python, import the module and set up a `TrypTag` instance:

```
from tryptag import TrypTag
tryptag = TrypTag()
```

To view data about a tagged cell line by gene id, for example `Tb927.7.1920`:

```
tryptag.fetch_gene_list()
print(tryptag.gene_list["Tb927.7.1920"])
```

This will print an object containing information about the tagging of this gene.
If tagging data for a gene id does not exist then this will give a `KeyError` error. You can handle this by using:

```
gene_id = "Tb927.not.arealgeneid"
if gene_id in tryptag.gene_list:
    print(tryptag.gene_list[gene_id])
else:
    print("gene id not found")
```

or

```
gene_id = "Tb927.not.arealgeneid"
try:
    print(tryptag.gene_list[gene_id])
except:
    print("gene id not found")
```

Tagging terminus-specific information can be accessed similarly, using `n` or `c` for the terminus. Again, trying to access data for a terminus which lacks data will give a `KeyError` error:

```
gene_id = "Tb927.7.1920"
terminus = "n"
if gene_id in tryptag.gene_list:
    if terminus in tryptag.gene_list[gene_id]:
        print(tryptag.gene_list[gene_id][terminus])
    else:
        print("data for terminus not found")
else:
    print("gene id not found")
```

### Microscopy data
You can use `tryptag` to download the microscopy data. In python, import the module and set up a `TrypTag` instance:

```
from tryptag import TrypTag
tryptag = TrypTag()
```

You can fetch the data for a specific gene id and tagging terminus using `fetch_data`.
This looks up in which tagging plate correspond to the most recent replicate of this tagging, and the url at which to find this data.
It then downloads and decompresses the data to the `data_cache_path` directory.

This will take a long time, to get image data for a single gene the data for an entire plate needs to be downloaded. This is typically ~10 to 20 Gb.

```
tryptag.fetch_data("Tb927.7.1920", "n")
```

This will give a `KeyError` error if there is no data for that terminus. To fetch, for example, image data for a list of gene ids of interest, you could use:

```
gene_ids = ["Tb927.7.1920", "Tb927.1.2670", "Tb927.11.1150"]
termini = ["n", "c"]
for gene_id in gene_ids: 
    if gene_id in tryptag.gene_list:
        for terminus in termini:
            if terminus in tryptag.gene_list[gene_id]:
                tryptag.fetch_data(gene_id, terminus)
```

Look through the data cache directory and you will find the microscopy data, in one subdirectory per tagging plate and named by gene id and tagging terminus.

### Image analysis

The primary intended use of the `tryptag` module is for easy access of specific field of view and cell images for automated image analysis. See quickstart.

Microscopy data is in three image channels and two thresholded images:

Image channels:

1. Phase contrast (transmitted light, overall cell morphology) uint16
2. mNG fluorescence (green fluorescence, from the tagged protein) uint32
3. DNA stain fluorescence (blue fluorescence, using Hoechst 33342) uint16

Thresholded images:

1. Thresholded phase contrast (cells) uint8, 255 = object
2. Thresholded DNA stain (nuclei and kinetoplasts - mitochondrial DNA organelles) uint8, 255 = object

`open_field` and `open_cell` both return an array of images in this order.

### General tips

Make sure you have enough free disk space to download, decompress and cache the image data. This is ~40 Gb for a single plate and ~6 Tb for the entire dataset.
The default cache location is `_tryptag_data` within the current working directory. You can change this to a relative or absolute path to any directory you wish - we recommend a scratch drive with sufficient space.
Make sure this is set at the start of every script:

Linux/Mac:
```
from tryptag import TrypTag
tryptag = TrypTag()
tryptag.data_cache_path = "\mnt\z\my\scratch\directory"
```
or Windows:
```
tryptag.data_cache_path = "Z:/my/scratch/directory"
```

Do not delete or move files from `data_cache_path`. `tryptag` does not check the plate subdirectories for integrity. You can, however, safely delete a plate subdirectory.
Interrupt of either data download or zip decompression _should_ behave gracefully, leaving partial data but not preventing later automatic re-download and/or re-decompression.

If multiple scripts using the same `data_cache_path` simultaneously try to download a plate it _should_ be handled gracefully.
One script should download and decompress the image data, while the others (silently) wait until it the image data is available.
However, for large-scale analyses, it is more robust to ensure all data is already cached. You can use `fetch_all_data` to download all image data (this will probably take more than one week!):

```
tryptag.fetch_all_data()
```

The TrypTag data may have minor errors which will be corrected over time. `fetch_all_data` always fetches the latest localisation listing from [Zenodo](https://zenodo.org/record/6862289).
Cached image data may be an older version. `tryptag` records the MD5 hash of the source zip files. If the data source (Zenodo depositions) are updated, the MD5 hash will change.
Cached data inconsistent MD5 hashes can be checked using:

```
tryptag.check_data_cache()
```

`tryptag` gives quite verbose information about what it is currently doing to fetch data. To silence this output set `print_status` at the start of the script before running any `tryptag` functions:

```
from tryptag import TrypTag
tryptag = TrypTag()
tryptag.print_status = False
# do your analysis here
```

# Citing

If you use the TrypTag data resource, please cite Billington _et al._ 2023 Nature Microbiology [doi:10.1038/s41564-022-01295-6](https://doi.org/10.1038/s41564-022-01295-6).
We recommend including this citation in the results or methods if TrypTag was used as part of a discovery process.
If directly using TrypTag images, please also indicate in the figure legend or similar which images are from TrypTag.

If you use the `tryptag` module to access or analyse TrypTag data, please also cite this [Github repository](https://github.com/zephyris/tryptag) and the master TrypTag Zenodo deposition [doi:10.5281/zenodo.6862289](https://doi.org/10.5281/zenodo.6862289).

You may also find the following papers of use:
1. Dean _et al._ 2019 Trends Parasitol. [doi:10.1016/j.pt.2016.10.009](https://doi.org/10.1016/j.pt.2016.10.009) Project announcement, with original aims and experimental strategy.
2. Halliday _et al._ 2019 Mol. Biochem. Parasitol. [doi:10.1016/j.molbiopara.2018.12.003](https://doi.org/10.1016/j.molbiopara.2018.12.003) Describes the localisation ontology, with landmark protein examples.
