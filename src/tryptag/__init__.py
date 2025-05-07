import logging
import sys
from .tryptag import TrypTag, FieldImage, CellImage, CellLine
__all__ = ["TrypTag", "FieldImage", "CellImage", "CellLine"]

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
