# Non-rigid alignment of serial sections

`02_align-sequential-shared.py` registers the four sections onto `WD-76845-044`
and writes the warped channels.

## Environment

```bash
pixi install
pixi run 02_align-sequential-shared.py
```

`pixi.lock` is solved for `win-64` only. To run elsewhere, add your platform to
`platforms` in `pixi.toml` and re-solve — the dependencies themselves are not
Windows-specific.

`pirt.GravityRegistration.register()` draws its progress into a `visvis` Qt
window, so **a display is required**; the script will not run headless.

`numpy` is pinned below 1.24 in `pixi.toml`; loosen it only after checking
that `pirt` still works.

## What the script does

1. Read the registration channel of each section and downscale it 4× with
   `downscale_local_mean`.
2. Rotate each thumbnail by its `rotation` angle, then paste each into a
   common-shape canvas at its `img_offsets` position.
3. Run `pirt.GravityRegistration` groupwise over the four thumbnails.
4. For each non-reference section and each channel: read the full-resolution
   channel, apply the same rotation and offset, apply the deformation field
   resized 4× back to full resolution, crop to the reference's shape, and write
   `<section>-ch_<NN>.tif`.
5. Write the reference section's channels unchanged, so all four sections end
   up on the same grid.

Output goes to `out_dir` as one single-channel, zlib-compressed, tiled TIFF per
section per channel (211 files for this dataset).

## Configuration

The only edit needed is `img_dir` at the top of the script — the directory
holding the four OME-TIFFs from step 01. The committed value is a Windows
`E:\...` path. Output is written to `img_dir / "aligned"`. The remaining
per-section values in `file_settings` are the ones this alignment was produced
with and should be left as they are:

- **`registration channel`** — zero-based page index of the DNA channel used
  for registration: 8 for sections 044 and 046, 12 for 045 and 047.
- **`number of channels`** — 40 for section 044, 57 for the rest.
- **`rotation`** — degrees counter-clockwise about the image centre, as
  `skimage.transform.rotate` defines it.
- **`img_offsets`** — `(row, col)` translation relative to the reference, **in
  4×-downscaled pixels**; the script multiplies by 4 when warping at full
  resolution.

`translation_offset()` estimates an offset by template matching the convex hull
of one thumbnail against another. The call is commented out (line ~201) so that
reruns skip it, and its results are recorded in the script:
`[(0, 0), (-48, -68), (556, 273), (336, 68)]`. The offsets actually used
replace section 045's entry with `(36, -248)`, obtained from
[Ashlar](https://github.com/labsyspharm/ashlar)'s thumbnail approach on a
further 20×-downscaled image, because the template-matching result for that
section was wrong.
