from tryptag import TrypTag
from tryptag.images import CellImage


def test_high_level():
    tt = TrypTag()
    cell_line = tt.gene_list["wild-type.2.2h"]["N"]
    cell_image = tt.open_cell(cell_line, 2, 8)
    assert isinstance(cell_image, CellImage)
