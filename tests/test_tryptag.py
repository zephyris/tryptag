from tryptag import TrypTag, tryptools
from tryptag.bia import BioimageArchive
from tryptag.cache import Cache, IncompatibleCacheError
from tryptag.images import CellImage
from tryptag.datasource import CellLine
from tryptag.zenodo import Zenodo

import pytest


def test_data_origin_check(tmp_path):
    cache = Cache(tmp_path)
    BioimageArchive(cache)

    cache = Cache(tmp_path)
    with pytest.raises(IncompatibleCacheError):
        Zenodo(cache)

    cache.delete()

    cache = Cache(tmp_path)
    Zenodo(cache)

    cache = Cache(tmp_path)
    with pytest.raises(IncompatibleCacheError):
        BioimageArchive(cache)

    cache.delete()

@pytest.fixture(
    scope="session",
    params=["zenodo", "bia"],
)
def tt_instance(tmp_path_factory, request):
    datasource_name = request.param
    if datasource_name == "zenodo":
        datasource = Zenodo
    elif datasource_name == "bia":
        datasource = BioimageArchive
    else:
        raise ValueError("Unknown data source name")
    cache = Cache(tmp_path_factory.mktemp(datasource_name))
    tt = TrypTag(
        datasource=datasource(cache=cache),
    )
    yield tt

    # Teardown
    cache.delete()


def test_localisation_search_simple(tt_instance: TrypTag):
    query = "nucleoplasm"
    hits = tt_instance.localisation_search(query)
    assert len(hits) == 521


def test_localisation_search_advanced(tt_instance: TrypTag):
    query = "spindle"
    match_subterms = True
    exclude_modifiers = {"weak", "<10%", "25%"}
    required_modifiers = {"cell cycle dependent"}

    hits = tt_instance.localisation_search(
        query,
        match_subterms=match_subterms,
        required_modifiers=required_modifiers,
        exclude_modifiers=exclude_modifiers,
    )

    assert len(hits) == 1


def test_localisation_ontology(tt_instance: TrypTag):
    ontology = tt_instance.localisation_ontology
    assert len(ontology) == 50
    assert ontology["nuclear lumen"]["parent"][-1] == "nucleus"


def test_cellimage(tt_instance: TrypTag):
    cell_line = tt_instance.gene_list["wild-type.2.2h"]["N"]
    cell_image = tt_instance.open_cell(cell_line, 2, 8)
    assert isinstance(cell_image, CellImage)


@pytest.fixture(scope="session")
def parental_wl_results():
    return sorted([
        {
            'cell_line': CellLine("wild-type.1.2h", "N"),
            'result': {
                '2K1N': 32,
                '1K1N': 108,
                '0K1N': 2,
                '1K2N': 11,
                '2K2N': 45,
                '2K3N': 6,
                '3K3N': 3,
                '4K4N': 1,
                '3K2N': 6,
                '4K3N': 2
            }
        },
        {
            'cell_line': CellLine("wild-type.2.0h", "N"),
            'result': {
                '1K1N': 117,
                '3K2N': 14,
                '2K1N': 35,
                '2K2N': 28,
                '1K2N': 15,
                '4K3N': 7,
                '0K1N': 1,
                '4K4N': 1,
                '3K3N': 6,
                '1K0N': 2,
                '2K3N': 1}
        },
        {
            'cell_line': CellLine("wild-type.2.2h", "N"),
            'result': {
                '0K0N': 1,
                '1K1N': 101,
                '2K2N': 39,
                '2K1N': 46,
                '3K2N': 30,
                '4K3N': 6,
                '1K2N': 6,
                '3K3N': 17,
                '5K5N': 1,
                '2K3N': 2,
                '0K1N': 1,
                '4K4N': 2,
                '3K4N': 1,
                '1K0N': 1
            }
        },
        {
            'cell_line': CellLine("wild-type.1.0h", "N"),
            'result': {
                '1K1N': 133,
                '2K1N': 28,
                '2K2N': 44,
                '3K2N': 14,
                '3K3N': 15,
                '4K3N': 10,
                '1K2N': 13,
                '2K3N': 8,
                '4K4N': 5,
                '3K4N': 3,
                '1K3N': 1,
                '0K1N': 2,
                '0K2N': 1,
                '5K5N': 1,
                '5K4N': 1
            }
        }
    ], key=lambda e: (e["cell_line"].gene_id, e["cell_line"].terminus))


def analyse(tryptag, cell_line):
    result = {}
    fieldcell_list = tryptag.cell_list(cell_line)
    for fieldcell in fieldcell_list:
        cell_image = tryptag.open_cell(
            cell_line, fieldcell.field.index, fieldcell.index)
        kn_result = tryptools.cell_kn_analysis(cell_image)
        if kn_result["count_kn"] not in result:
            result[kn_result["count_kn"]] = 0
        result[kn_result["count_kn"]] += 1
    return result


@pytest.mark.parametrize(
    "mode,workers",
    [
        (None, None),
        ("process", None),
        ("process", 2),
        ("thread", None),
        ("thread", 2),
    ]
)
def test_analyse_worklist(
    tt_instance: TrypTag,
    parental_wl_results,
    mode,
    workers,
):
    worklist = tt_instance.worklist_parental()
    # Sort results, otherwise the comparison fails.
    results = sorted(
        tt_instance.analyse_list(
            worklist,
            analyse,
            multiprocess_mode=mode,
            workers=workers,
        ),
        key=lambda e: (e["cell_line"].gene_id, e["cell_line"].terminus)
    )

    assert results == parental_wl_results
