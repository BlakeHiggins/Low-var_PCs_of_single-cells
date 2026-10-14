### Script for analysis of the 3D IMC breast carcinoma dataset for "Are low-variance Principal Components of single-cell omics data biologically meaningful, or just noise?"
### Adds the 3D centroid of each cell mask to the original authors' single-cell intensity table.
### Inputs:  full_mask_final_segmentation_hwatershed_500.00_90%.tif, model201710_mean_intensities.csv
### Output:  model201710_cells_with_positions.csv
### Data source: https://www.nature.com/articles/s43018-021-00301-w (data at https://zenodo.org/records/5782846)
### adapts centroid and count merging logic from authors' analysis code BodenmillerGroup/3D_IMC_publication (https://github.com/BodenmillerGroup/3D_IMC_publication/tree/82f43d46713a461b5a1cec3e70b5bd5eb81830ca/Python/scripts/single_cell_analysis)

### import libraries
import tifffile
import numpy as np
import pandas as pd
from skimage.measure import regionprops

### Load the segmentation mask
mask = tifffile.imread("Data/full_mask_final_segmentation_hwatershed_500.00_90%.tif")

props = regionprops(mask)
### extract centroids and lables from mask 
spatial_df = pd.DataFrame({
    "label":        [p.label for p in props],
    "centroid_z":   [p.centroid[0] for p in props],
    "centroid_y":   [p.centroid[1] for p in props],
    "centroid_x":   [p.centroid[2] for p in props]
})
### load in counts 
counts = pd.read_csv("Data/model201710_mean_intensities.csv")


counts_with_labels = counts
#### merge centroids and counts 
merged = spatial_df.merge(counts_with_labels, on="label")
### save out 
merged.to_csv("Data/model201710_cells_with_positions.csv", index=False)
print("Done.")