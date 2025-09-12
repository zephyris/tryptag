from __future__ import annotations
import io
import itertools
import logging
import weakref

import numpy
import skimage
from PIL import Image
import tifffile

from .cache import FileTypes
from .datasource import Field, Cell, CellLine, DataSource

logger = logging.getLogger("tryptag.images")


def autocontrast(image: numpy.ndarray) -> numpy.ndarray[numpy.uint8]:
    """
    Quick and dirty autocontrast.

    This function takes the given `image`, scales it such that the minumum
    value becomes zero and the maximum value becomes 255 and returns a uint8
    image.

    :param image: numpy.ndarray, input image
    :returns: numpy.ndarray type uint8
    """
    cmin = image.min()
    cmax = image.max()
    return (
        (255 * ((image - cmin) / (cmax - cmin))).round(0).astype(numpy.uint8))


def _channel_to_png_bytes(channel: numpy.ndarray):
    image = Image.fromarray(channel)
    bytes_obj = io.BytesIO()
    with bytes_obj:
        image.save(bytes_obj, format="PNG")
        return bytes_obj.getvalue()


class Channel(numpy.ndarray):
    def _repr_png_(self):
        return _channel_to_png_bytes(
            autocontrast(self))


class CellImage():
    """
    CellImage object holding information on a (cropped) section of a TrypTag
    image containing a specific cell.
    """
    IMAGE_MEMBERS = [
        "phase",
        "mng",
        "dna",
        "phase_mask",
        "phase_mask_othercells"
        "dna_mask",
    ]

    phase: Channel | None = None
    mng: Channel | None = None
    dna: Channel | None = None
    phase_mask: Channel | None = None
    phase_mask_othercells: Channel | None = None
    dna_mask: Channel | None = None

    def __init__(
        self,
        cell: Cell,
        rotated: bool,
        width: int = 323,
        fill_centre: tuple[int, int] | None = None,
        crop_centre: tuple[int, int] | None = None,
        angle: float | None = None,
        custom_field_image: FieldImage | None = None
    ):
        """
        Initialise a CellImage object.

        :param cell: the Cell metadata object
        :param rotated: bool, whether the `CellImage` should be rotated such
            that the cell is displayed in the canonical way.
        :param width: int, width the image should be cropped to
        :param fill_centre: xy coordinate of the flood fill centre for the
            mask. This needs to be within the cell.
        :param crop_centre: xy coordinate of the centre of the resulting image
        :param angle: float, angle the cell subtends with the x axis
        :param custom_field_image: `FieldImage` or None, if given, the raw
            image data is taken from this `FieldImage` rather than loaded from
            the `Field`. Allows pre-processing of the images.
        """
        logger.debug(f"Creating cell image for {cell}.")

        self.rotated = rotated

        field_image = FieldImage.from_field(
            cell.field,
            custom_field_image=custom_field_image
        )
        # Hold strong reference to FieldImage as long as CellImage exists.
        # This way FieldImage._CACHE will store a weak reference to it
        # and subsequent CellImages from the same FieldImage can use
        # a cached version.
        self.field_image = field_image

        self.cell = cell
        if fill_centre is None:
            fill_centre = cell.wand
        if crop_centre is None:
            crop_centre = (int(cell.centre[0]), int(cell.centre[1]))
        if angle is None:
            angle = cell.angle

        # define image crop function
        def _skimage_crop(image, x, y, w, h):
            return image[int(y):int(y + h), int(x):int(x + w)]

        phase_mask = field_image.phase_mask.copy()
        # process phase threshold image to split cell of interest from
        # neighbouring cells, nb. xy swapped in skimage arrays
        # replace existing mask pixels with value 127
        phase_mask[phase_mask == 255] = 127
        # flood fill cell of interest to 255
        phase_mask = skimage.morphology.flood_fill(
            phase_mask,
            (fill_centre[1], fill_centre[0]),
            255,
        )
        # append a copy (channels index 5), pixels if equal 127 ie. other cells
        phase_mask_othercells = 255*(phase_mask == 127)
        # Remove other cells from phase_mask
        phase_mask[phase_mask < 255] = 0

        channels_iterator = itertools.chain(
            field_image.iter_images(copy=True),
            [
                ("phase_mask", phase_mask),  # Overwrite the original phase mask
                ("phase_mask_othercells", phase_mask_othercells)
            ]
        )
        # crop (and potentially rotate) to get cell image
        cell_channels: dict[str, numpy.ndarray] = {}
        if width < 0:
            # padding mode
            label_image = skimage.measure.label(phase_mask)
            (ymin, xmin, ymax, xmax) = skimage.measure.regionprops(
                label_image)[0]["bbox"]
            # do actual cropping
            for ch_name, channel in channels_iterator:
                # negative width is padding
                padding = -width
                # if crop outside of image bounds, then first increase canvas
                # size
                offs = 0
                if (
                    ymin - padding < 0 or
                    xmin - padding < 0 or
                    ymax + padding > channel.shape[0] or
                    xmax + padding > channel.shape[1]
                ):
                    channel = numpy.pad(
                        channel,
                        ((padding, padding), (padding, padding)),
                        mode="median"
                    )
                    offs = padding
                # crop
                cell_channels[ch_name] = (
                    channel[ymin - padding + offs:ymax + padding + offs,
                            xmin - padding + offs:xmax + padding + offs])
        elif width > 0:
            # fixed width mode
            if rotated:
                width_inter = width * 1.5  # greater than width * 2**0.5
                height = round(width / 2)
            else:
                width_inter = width
            # do actual cropping
            for ch_name, channel in channels_iterator:
                # if crop outside of image bounds, then first increase canvas
                # size
                offs = 0
                half_width_inter = round(width_inter / 2)
                if (
                    crop_centre[0] - half_width_inter < 0 or
                    crop_centre[1] - half_width_inter < 0 or
                    crop_centre[0] + half_width_inter > channel.shape[1] or
                    crop_centre[1] + half_width_inter > channel.shape[0]
                ):
                    channel = numpy.pad(
                        channel,
                        (
                            (half_width_inter, half_width_inter),
                            (half_width_inter, half_width_inter)
                        ),
                        mode="median"
                    )
                    offs = half_width_inter
                # square crop
                channel = _skimage_crop(
                    channel,
                    crop_centre[0] + offs - half_width_inter,
                    crop_centre[1] + offs - half_width_inter,
                    width_inter,
                    width_inter
                )
                if rotated:
                    # if rotating, rotate then crop to final dimensions
                    # have have to force data type and use preserve_range=True
                    # to prevent rotate from mangling the data
                    channel_dtype = channel.dtype
                    channel = skimage.transform.rotate(
                        channel,
                        -cell.angle,
                        preserve_range=True
                    ).astype(channel_dtype)
                    channel = _skimage_crop(
                        channel,
                        half_width_inter - width / 2,
                        half_width_inter - height / 2,
                        width,
                        height
                    )
                cell_channels[ch_name] = channel
        else:
            raise ValueError(
                "`width` must be a nonzero integer, positive for fixed image "
                "width in pixels, negative for padding around phase_mask in "
                "pixels.")
        # downstream analysis (including tryptools) allowed to assume cell
        # mask does not touch image edge, therefore set border pixels to 0
        # (may clip large cells)
        cell_channels["phase_mask"][:, [0, -1]] = 0
        cell_channels["phase_mask"][[0, -1], :] = 0

        for channel, image in cell_channels.items():
            setattr(self, channel, image.view(Channel))

    def __repr__(self):
        string = "\n".join([str(x) for x in [
            "phase", self.phase,
            "mng", self.mng,
            "dna", self.dna,
            "phase_mask", self.phase_mask,
            "dna_mask", self.dna_mask,
            "phase_mask_othercells",
            self.phase_mask_othercells]])
        if self.cell.field.cell_line is not None:
            string += "\n" + str(self.cell.field.cell_line)
        string += "\n" + " ".join([str(x) for x in [
            "field_index =", self.cell.field.index,
            "cell_index =", self.cell.index,
            "rotated =", self.rotated]])
        return string

    def __str__(self):
        string = ""
        if self.cell.field.cell_line is not None:
            string += str(self.cell.field.cell_line)
        string += " " + " ".join([str(x) for x in [
            "field_index =", self.cell.field.index,
            "cell_index =", self.cell.index,
            "rotated =", self.rotated]])
        return string

    def _repr_png_(self):
        phase = autocontrast(self.phase)
        mng = autocontrast(self.mng)
        dna = autocontrast(self.dna)

        output = numpy.empty(
            (phase.shape[0], phase.shape[1], 3), dtype=numpy.uint16)
        output[:, :, :] = phase[:, :, numpy.newaxis]
        output[:, :, 1] += mng
        output[:, :, [0, 2]] += dna[:, :, numpy.newaxis]
        output = numpy.clip(output, 0, 255).astype(numpy.uint8)

        return _channel_to_png_bytes(output)


class FieldImage():
    """
    FieldImage object holding a TrypTag image belonging to a CellLine.
    """
    IMAGE_MEMBERS = [
        "phase",
        "mng",
        "dna",
        "phase_mask",
        "dna_mask",
    ]

    phase: Channel | None = None
    mng: Channel | None = None
    dna: Channel | None = None
    phase_mask: Channel | None = None
    dna_mask: Channel | None = None

    cell_line: CellLine | None = None
    field: Field | None = None
    _field_index: int | None = None
    custom_field_image: FieldImage | None = None

    _CACHE: weakref.WeakValueDictionary[
        tuple[Field, FieldImage | None], FieldImage
    ] = weakref.WeakValueDictionary()

    def __init__(
        self,
        phase: numpy.ndarray | None = None,
        mng: numpy.ndarray | None = None,
        dna: numpy.ndarray | None = None,
        phase_mask: numpy.ndarray | None = None,
        dna_mask: numpy.ndarray | None = None,
        cell_line: CellLine | None = None,
        field_index: int | None = None,
        phase_contrast: (int, int) | None = None,
        mng_contrast: (int, int) | None = None,
        dna_contrast: (int, int) | None = None
    ):
        """
        Initialise a new FieldImage object.

        Note that this does not load data from the data source. All parameters
        are optional. This can be used to create a custom `FieldImage`.

        :param phase: phase microscopy channel
        :param mng: green (tag) fluorescence channel
        :param dna: magenta (DNA stain) fluorescence channel
        :param phase_mask: thresholded mask image for cell bodies
        :param dna_mask: thresholded mask image for DNA (nucleus and
            kinetoplast)
        :param cell_line: CellLine object this `FieldImage` belongs to
        :param field_index: int, index of the `Field` in the `CellLine`
        """
        self.phase = phase.view(Channel) if phase is not None else None
        self.mng = mng.view(Channel) if mng is not None else None
        self.dna = dna.view(Channel) if dna is not None else None
        self.phase_mask = (
            phase_mask.view(Channel) if phase_mask is not None else None)
        self.dna_mask = (
            dna_mask.view(Channel) if dna_mask is not None else None)
        if (cell_line is None) != (field_index is None):
            raise ValueError("need to specify both cell_line and field_index")
        elif cell_line is not None:
            self.cell_line = cell_line
            self._field_index = field_index

            if cell_line.initialised:
                self.field = cell_line.fields[field_index]

    @property
    def field_index(self):
        """Index of the `Field` in the `CellLine`."""
        if self.field is not None:
            return self.field.index
        return self._field_index

    def update(self, other: FieldImage):
        """Update this FieldImage with non-None image planes in `other`."""
        for member in self.IMAGE_MEMBERS:
            other_img: numpy.ndarray | None = getattr(other, member)
            if other_img is not None:
                setattr(self, member, other_img.copy().view(Channel))

    def iter_images(self, copy: bool = False):
        """Generator over image members. Yields a tuple (name, image) per
        iteration. Makes a copy of the images first if `copy=True`.
        """
        for member in self.IMAGE_MEMBERS:
            image: numpy.ndarray | None = getattr(self, member)
            if copy and image is not None:
                image = image.copy()
            yield (member, image)

    def __repr__(self):
        string = "\n".join([str(x) for x in [
            "phase", self.phase,
            "mng", self.mng,
            "dna", self.dna,
            "phase_mask", self.phase_mask,
            "dna_mask", self.dna_mask]])
        if self.field is not None:
            string += "\n" + str(self.field.cell_line)
            string += "\n" + " ".join([str(x) for x in [
                "field_index =", self.field.index]])
        return string

    def __str__(self):
        string = ""
        if self.field is not None:
            string += str(self.field.cell_line)
            string += " " + " ".join([str(x) for x in [
                "field_index =", self.field.index]])
        return string

    def _open_image(
        self,
        datasource: DataSource,
    ):
        if self.field is None:
            raise ValueError("self.field is None")
        img_path = datasource.load_plate_file(
            self.field.cell_line.plate,
            self.field.filename(FileTypes.IMAGE),
            return_file_object=False,
        )
        tif = tifffile.TiffFile(img_path)
        try:
            channel_contrast = tif.imagej_metadata["Ranges"]
            contrast = {
                "phase": (channel_contrast[0], channel_contrast[1]),
                "mng": (channel_contrast[2], channel_contrast[3]),
                "dna": (channel_contrast[4], channel_contrast[5]),
            }
        except:
            contrast = None
        field_image = [
            tif.asarray(0),
            tif.asarray(1),
            tif.asarray(2),
        ]
        return field_image, contrast

    def _open_thresholded(
        self,
        datasource: DataSource,
    ):
        if self.field is None:
            raise ValueError("self.field is None")
        thr_path = datasource.load_plate_file(
            self.field.cell_line.plate,
            self.field.filename(FileTypes.THRESHOLDED),
            return_file_object=False,
        )
        field_threshold = skimage.io.imread(thr_path)
        # 'clean' the field_threshold images, ensure all pixels less than
        # 255 are set to 0
        for image in field_threshold:
            image[image < 255] = 0
        return field_threshold

    def _process(
            self,
            datasource: DataSource,
            custom_field_image: FieldImage | None = None
    ):
        image, contrast = self._open_image(datasource)
        thresholded = self._open_thresholded(datasource)

        # setup output, copying images as downstream usage may modify
        self.phase = image[0].astype("uint16", copy=True).view(Channel)
        self.mng = image[1].astype("uint32", copy=True).view(Channel)
        self.dna = image[2].astype("uint16", copy=True).view(Channel)
        self.phase_mask = thresholded[0].astype("uint8", copy=True).view(Channel)
        self.dna_mask = thresholded[1].astype("uint8", copy=True).view(Channel)
        self.phase_contrast = contrast.get("phase") if contrast else None
        self.mng_contrast = contrast.get("mng") if contrast else None
        self.dna_contrast = contrast.get("dna") if contrast else None

        if custom_field_image is not None:
            self.update(custom_field_image)

        # Hold a strong reference to this object in the class.
        # Makes sure that FieldImage._CACHE has the weak reference
        # for subsequent accesses to this Field Image.
        FieldImage._last_field_image = self

    @staticmethod
    def from_field(
        field: Field,
        custom_field_image: FieldImage | None = None,
    ):
        """
        Create a `FieldImage` from a given `Field` object.

        Uses channels from `custom_field_image` if given. This allows custom
        image manipulation before processing.

        :param field: `Field` metadata object
        :param custom_field_image: `FieldImage` with any custom channels
        :return: new `FieldImage` object
        """
        logger.debug(f"Loading image for field {field}...")
        try:
            field_image = FieldImage._CACHE[(field, custom_field_image)]
            logger.debug("...from cache.")
        except KeyError:
            logger.debug("...from file.")
            field_image = FieldImage()
            field_image.field = field
            field_image._process(
                field.datasource,
                custom_field_image=custom_field_image
            )

            FieldImage._CACHE[(field, custom_field_image)] = field_image
        return field_image

    def _repr_png_(self):
        phase = autocontrast(self.phase)
        mng = autocontrast(self.mng)
        dna = autocontrast(self.dna)

        output = numpy.empty(
            (phase.shape[0], phase.shape[1], 3), dtype=numpy.uint16)
        output[:, :, :] = phase[:, :, numpy.newaxis]
        output[:, :, 1] += mng
        output[:, :, [0, 2]] += dna[:, :, numpy.newaxis]
        output = numpy.clip(output, 0, 255).astype(numpy.uint8)

        return _channel_to_png_bytes(output)
