# tryptag
`tryptag` is a python module for accessing and handling [TrypTag](http://tryptag.org) genome-wide protein localisation project data.
Its primary intended use is for easy access to data for automated image analysis.

## Installation
Install (or install an updated version) using `pip`:

```
pip install git+https://github.com/zephyris/tryptag
```

To uninstall also use `pip`:

```
pip uninstall tryptag
```

`tryptag` requires several python modules: `numpy` `scikit-image` and `progressbar`. These should be installed using `pip`:

```
pip install numpy scikit-image progressbar
```

To use the `tryptag` module, import the `TrypTag` class and set up an instance (normally called `tryptag`):

```
from tryptag import TrypTag
tryptag = TrypTag()
```

## Quickstart guide
Microscopy data is multiple fields of view per cell line.
To open a specific field of view of a tagged cell line, access it by gene id (as used on [TriTrypDB](http://tritrypdb.org)), tagging terminus (`n` or `c`) and field index (integer):

```
phase, mng, dna, phasethreshold, dnathreshold = tryptag.open_field(gene_id, terminus, field)
```

The cells in the phase threshold image are indexed and can be opened individually.
To open a specific cell:

```
phase, mng, dna, phasethreshold, dnathreshold = tryptag.open_cell(gene_id, terminus, field, cell)
```

Opened images are `numpy` `ndarray` objects, as used by `scikit-image`.

Bear in mind that accessing a nonexistant gene id, tagging terminus, field and cell will give `KeyError` errors. For example:

```
gene_id = "Tb927.not.arealgeneid"
terminus = "n"
tryptag.fetch_data(gene_id, terminus)
try:
    for field in range(5):
        [pha, mng, dna, pth, dth] = tryptag.open_field(gene_id, terminus, field)
        # do your analysis here
except:
    print("gene id, terminus or field not found")
```

You can access the `tryptag.gene_list` dict for more intelligent iteration. To iterate through all cells for all gene ids and termini that exist:
```
tryptag.fetch_gene_list()
termini = ["n", "c"]
for gene_id in tryptag.gene_list:
    for terminus in termini:
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
gene_id = "Tb927.7.1920"
if gene_id in tryptag.gene_list:
    print(tryptag.gene_list[gene_id])
else:
    print("gene id not found")
```

or

```
gene_id = "Tb927.7.1920"
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

This will take a long time, to get image data for a single gene the data for an entire plate needs to be downloaded. This is typically ~20 Tb.

```
fetch_data("Tb927.7.1920", "n")
```

This will give a `KeyError` error if there is no data for that terminus. To fetch, for example, image data for a list of gene ids of interest, you could use:

```
gene_ids = ["Tb927.7.1920", "Tb927.1.2670", "Tb927.11.1150"]
termini = ["n", "c"]
for gene_id in gene_id: 
    if gene_id in tryptag.gene_list:
        for terminus in termini:
            if terminus in tryptag.gene_list[gene_id]:
                fetch_data(gene_id, terminus)
```

Look through the data cache directory and you will find your desired microscopy data, in one subdirectory per tagging plate and named by gene id and tagging terminus.

### Image analysis

The primary intended use of the `tryptag` module is for easy access of specific field of view and cell images for automated image analysis. See quickstart.

### General tips

Make sure you have enough free disk space to download, decompress and cache the image data. This is ~40 Gb for a single plate and up to ~8 Tb for the entire dataset.
The default cache location is `_tryptag_data` within the current working directory. You can change this to any directory you wish, we recommend a scratch drive with sufficient space:

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

Do not remove files from `data_cache_path`, `tryptag` does not check if files have been removed. You can, however, delete a whole plate directory.

Do not run multiple scripts using the same `data_cache_path` as they may try to write to the same file at the same time. The exception is if all data is already cached. You can run a script to download everything:

```
from tryptag import TrypTag
tryptag = TrypTag()
tryptag.fetch_all_data()
```

The TrypTag data may have minor errors which will be corrected over time. `fetch_all_data` always fetches the latest localisation listing.
Cached image data may be an older version. `tryptag` records the MD5 hash of the source zip files. If the data source (Zenodo depositions) are updated, the MD5 hash will change.
Cached data inconsistent MD5 hashes can be checked using:

```
from tryptag import TrypTag
tryptag = TrypTag()
tryptag.check_data_cache()
```

`tryptag` gives quite verbose information about what it is currenly doing to fetch data. To silence this output set `print_status` at the start of the script before running any `tryptag` functions:

```
from tryptag import TrypTag
tryptag = TrypTag()
tryptag.print_status = False
# do your analysis here
```

