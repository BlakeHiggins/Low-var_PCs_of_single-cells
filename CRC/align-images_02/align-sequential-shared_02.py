### Script for Analysis of CyCIF CRC dataset for "Are low-variance Principal Components of single-cell omics data biologically meaningful, or just noise?"
### Paper source:  https://www.cell.com/cell/fulltext/S0092-8674(22)01571-9
### Inputs: WD-76845-044.ome.tif WD-76845-045.ome.tif WD-76845-046.ome.tif WD-76845-047.ome.tif
### outputs: one single-channel TIFF per section per channel, all on a common grid
### paste(), crop_like() and dtype_convert() are copied from Ashlar (https://github.com/labsyspharm/ashlar, ashlar/utils.py),
### Copyright (c) 2022 President and Fellows of Harvard College, MIT License.
###The GravityRegistration setup is adapted from pirt's examples/reg_pairwise.py
###(https://github.com/almarklein/pirt), Copyright (c) 2014 Almar Klein, BSD-3-Clause.
# pylab
import copy
import pathlib
import warnings

import numpy as np
import pirt
import scipy.ndimage
import skimage.util
import tifffile
import visvis as vv
from matplotlib import pyplot as plt
from skimage.exposure import rescale_intensity
from skimage.feature import match_template
from skimage.morphology import convex_hull_image
from skimage.transform import downscale_local_mean, rotate


def crop_like(img, target):
    if img.shape[0] > target.shape[0]:
        img = img[: target.shape[0], :]
    if img.shape[1] > target.shape[1]:
        img = img[:, : target.shape[1]]
    return img


def dtype_convert(img, dtype):
    """Convert an image to the requested data-type.

    This is just a wrapper around skimage.util.dtype.convert that silences its
    FutureWarning, as Ashlar pins skimage to a version before that planned
    deprecation.

    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r".*scikit-image 1\.0", FutureWarning)
        return skimage.util.dtype.convert(img, dtype)


def paste(target, img, pos, func=None):
    """Composite img into target."""
    pos = np.array(pos)
    # Bail out if destination region is out of bounds.
    if np.any(pos >= target.shape[:2]) or np.any(pos + img.shape[:2] < 0):
        return
    pos_f, pos_i = np.modf(pos)
    yi, xi = pos_i.astype("i8")
    # Clip img to the edges of the mosaic.
    if yi < 0:
        img = img[-yi:]
        yi = 0
    if xi < 0:
        img = img[:, -xi:]
        xi = 0
    target_slice = target[yi : yi + img.shape[0], xi : xi + img.shape[1]]
    img = crop_like(img, target_slice)
    # Skip expensive sub-pixel shift if fractional position is zero.
    if pos_f.any():
        if img.ndim == 2:
            img = scipy.ndimage.shift(img, pos_f)
        else:
            for c in range(img.shape[2]):
                img[..., c] = scipy.ndimage.shift(img[..., c], pos_f)
        # For any axis where there is a non-zero subpixel shift, crop out the
        # last row or column of pixels on the "losing" side. These pixels will
        # be darker than normal and will introduce artifacts in most blending
        # modes.
        y1 = None if pos_f[0] <= 0 else 1
        y2 = None if pos_f[0] >= 0 else -1
        x1 = None if pos_f[1] <= 0 else 1
        x2 = None if pos_f[1] >= 0 else -1
        img = img[y1:y2, x1:x2]
        target_slice = target_slice[y1:y2, x1:x2]
        # Exit if image area is zero after subpixel shift.
        if not np.all(img.shape):
            return
    if np.issubdtype(img.dtype, np.floating):
        np.clip(img, 0, 1, img)
    img = dtype_convert(img, target.dtype)
    if func is None:
        target_slice[:] = img
    elif isinstance(func, np.ufunc):
        func(target_slice, img, out=target_slice)
    else:
        target_slice[:] = func(target_slice, img)


def bounding_box(img):
    rows = np.any(img, axis=1)
    cols = np.any(img, axis=0)
    rmin, rmax = np.where(rows)[0][[0, -1]]
    cmin, cmax = np.where(cols)[0][[0, -1]]

    return rmin, rmax, cmin, cmax


def translation_offset(img1, img2, debug=True):
    convex_hull2 = convex_hull_image(img2 > np.percentile(img2, 99))
    row_s, row_e, col_s, col_e = bounding_box(convex_hull2)
    result = match_template(img1, img2[row_s:row_e, col_s:col_e])
    row, col = np.unravel_index(np.argmax(result), result.shape)

    if debug:
        fig = plt.figure()
        ax1 = plt.subplot(1, 3, 1)
        ax2 = plt.subplot(1, 3, 2)
        ax3 = plt.subplot(1, 3, 3, sharex=ax2, sharey=ax2)

        ax1.imshow(img2[row_s:row_e, col_s:col_e], cmap=plt.cm.gray)
        ax1.set_axis_off()
        ax1.set_title("template")

        ax2.imshow(img1, cmap=plt.cm.gray)
        ax2.set_axis_off()
        ax2.set_title("image")
        # highlight matched region
        box_h, box_w = row_e - row_s, col_e - col_s
        rect = plt.Rectangle((col, row), box_w, box_h, edgecolor="r", facecolor="none")
        ax2.add_patch(rect)

        ax3.imshow(result)
        ax3.set_axis_off()
        ax3.set_title("match_template\nresult")
        # highlight matched region
        ax3.autoscale(False)
        ax3.plot(
            col, row, "o", markeredgecolor="r", markerfacecolor="none", markersize=10
        )

    return row - row_s, col - col_s


def match_pad(imgs, offsets=None):
    shape = np.array([i.shape for i in imgs]).max(axis=0)
    if not offsets:
        offsets = [(0, 0)] * len(imgs)

    padded_imgs = [
        np.full(shape, np.percentile(rescale_intensity(i), 25), dtype=i.dtype)
        for i in imgs
    ]
    for idx, (i, o) in enumerate(zip(imgs, offsets)):
        paste(padded_imgs[idx], rescale_intensity(i), o)
    return padded_imgs


# ---------------------------------------------------------------------------- #
#                               I/O configuration                              #
# ---------------------------------------------------------------------------- #
# Directory containing the four OME-TIFFs downloaded in step 01
img_dir = pathlib.Path(r"E:\htan-tnp-sardana")

file_settings = [
    {
        "location": img_dir / "WD-76845-044.ome.tif",
        "registration channel": 8,
        "pixel size": 0.65,
        "rotation": 0,
        "number of channels": 40,
    },
    {
        "location": img_dir / "WD-76845-045.ome.tif",
        "registration channel": 12,
        "pixel size": 0.65,
        "rotation": 6.7,
        "number of channels": 57,
    },
    {
        "location": img_dir / "WD-76845-046.ome.tif",
        "registration channel": 8,
        "pixel size": 0.65,
        "rotation": 352.9,
        "number of channels": 57,
    },
    {
        "location": img_dir / "WD-76845-047.ome.tif",
        "registration channel": 12,
        "pixel size": 0.65,
        "rotation": 5,
        "number of channels": 57,
    },
]

out_dir = img_dir / "aligned"
out_dir.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------- #
#                            initial rigid alignment                           #
# ---------------------------------------------------------------------------- #
_imgs = []
for i in file_settings:
    _img = tifffile.imread(i["location"], key=i["registration channel"])
    _img = downscale_local_mean(_img, (4, 4))
    _imgs.append(_img)


imgs = [copy.copy(_imgs[0])]
for idx, i in enumerate(_imgs[1:], 1):
    img = rotate(i, angle=file_settings[idx]["rotation"])
    imgs.append(img)


img_offsets = [(0, 0)]
# commented out to speed up
# for idx, i in enumerate(imgs[1:], 1):
#     offset = translation_offset(imgs[0], i, debug=True)
#     img_offsets.append(offset)

# the computed offset using the translation_offset function
# img_offsets = [(0, 0), (-48, -68), (556, 273), (336, 68)]

# The correct offset of image 1 is found using ashlar's thumbnail approach,
#  where the image is further 20-times downscaled
img_offsets = [(0, 0), (36, -248), (556, 273), (336, 68)]

padded_imgs = match_pad(imgs, img_offsets)


# ---------------------------------------------------------------------------- #
#                        gravity registration using pirt                       #
# ---------------------------------------------------------------------------- #

# Get figure
vv.closeAll()
fig = vv.figure(1)
vv.clf()
fig.position = 200, 100, 900, 500

cropped_thumb_padded_imgs = [i[:, :] for i in padded_imgs]

# Init registration
reg = pirt.GravityRegistration(*cropped_thumb_padded_imgs)

if isinstance(reg, pirt.DiffeomorphicDemonsRegistration):
    reg.params.speed_factor = 2
    reg.params.noise_factor = 0.1
    reg.params.mapping = "backward"
    reg.params.scale_sampling = 16
    reg.params.final_grid_sampling = 35
if isinstance(reg, pirt.GravityRegistration):
    reg.params.scale_levels = 5
    reg.params.scale_sampling = 10
    reg.params.mapping = "backward"
    reg.params.deform_wise = "groupwise"
    reg.params.deform_limit = 1
    reg.params.final_scale = 1
    reg.params.final_grid_sampling = 20
    reg.params.grid_sampling_factor = 0.5  # !! important especially for Laplace !!
    reg.params.frozenedge = True
    reg.params.mass_transforms = 2
    reg.params.speed_factor = 1.0

# Register, pass the figure so that the algorithm can show progress
reg.register(1, fig)


# ---------------------------------------------------------------------------- #
#                              warp image channels                             #
# ---------------------------------------------------------------------------- #
"""
load image
rotate image
shift and pad image
apply transform
crop image
save image
"""
ref_image = tifffile.imread(
    file_settings[0]["location"], key=file_settings[0]["registration channel"]
)


for idx, setting in enumerate(file_settings[1:], start=1):
    print("idx:", idx)

    grav_deform = reg.get_final_deform(idx, 0, "backward")
    deform_shape = np.array(grav_deform.field_shape, dtype=int)

    for page in range(setting["number of channels"]):
        print("page:", page)
        img_to_transform = tifffile.imread(setting["location"], key=page)
        img_ori_max = img_to_transform.max()
        img_to_transform = rotate(img_to_transform, angle=setting["rotation"])

        # The field was computed using a downscale_local_mean(img, (4, 4)) image
        img_shifted = np.full(deform_shape * 4, 0, dtype=np.float64)
        paste(img_shifted, img_to_transform, np.array(img_offsets[idx], dtype=int) * 4)

        del img_to_transform

        img_shifted = grav_deform.resize_field(
            tuple(deform_shape * 4)
        ).apply_deformation(img_shifted)

        img_shifted = crop_like(img_shifted, ref_image)
        img_shifted = rescale_intensity(
            img_shifted,
            out_range=(0, img_ori_max),
        ).astype(np.uint16)

        ori_filename = pathlib.Path(setting["location"]).stem.replace(".ome", "")
        tifffile.imwrite(
            out_dir / f"{ori_filename}-ch_{page:02d}.tif",
            img_shifted,
            compression="zlib",
            tile=(1024, 1024),
        )

# ----------- output single channel ref image for pyramid building ----------- #
setting = file_settings[0]
for channel in range(setting["number of channels"]):
    ref_image_i = tifffile.imread(setting["location"], key=channel)
    ref_filename = pathlib.Path(setting["location"]).stem.replace(".ome", "")

    tifffile.imwrite(
        out_dir / f"{ref_filename}-ch_{channel:02d}.tif",
        ref_image_i,
        compression="zlib",
        tile=(1024, 1024),
    )
