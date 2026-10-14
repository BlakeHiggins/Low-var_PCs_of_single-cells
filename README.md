## Are low-variance Principal Components of single-cell omics data biologically meaningful, or just noise?



### Workflow  
This repository contains the analysis code for "Are low-variance Principal Components of single-cell omics data biologically meaningful, or just noise?" with instructions for retrieving the public data sets it uses.


```
└── Low-variance_Principal_Components_of_single-cell_omics_data/
    ├── README.md
    ├── CRC/
    │   ├── crc_analysis_03.Rmd
    │   ├── align-images_02/
    │   │   ├── align-sequential-shared_02.py
    │   │   ├── align-sequential-shared_02_README.md
    │   │   ├── pixi.lock
    │   │   └── pixi.toml
    │   └── download-images_01/
    │       ├── .Rhistory
    │       ├── download-htan-files_01.md
    │       └── gen3_manifest.json
    ├── IMC/
    │   ├── IMC_make_csv_01.py
    │   └── IMC_pipeline_02.Rmd
    ├── seqFISH/
    │   └── Seqfish_pipeline_01.Rmd
    └── Tonsil/
        ├── 01_affine_registration.m
        ├── Flat-field_correction_03.Rmd
        ├── Reorder_02.Rmd
        ├── tonsil_pipeline_04.Rmd
        └── Data/
            └── affine_sweep_results.mat

```


#### **Download Public data** 

#####  **CyCIF tonsil dataset** 
Lin et al., 2018  
Paper source : https://elifesciences.org/articles/31657#content  

Single Cell Data: [elife-31657-fig6-data1-v3.zip](https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMzE2NTcvZWxpZmUtMzE2NTctZmlnNi1kYXRhMS12My56aXA-/elife-31657-fig6-data1-v3.zip?_hash=KqaC6E4lyB8cBihZ1quzv%2FgyoFJXnlDlhLEKvo1UijU%3D) 

Tonsil Image data: [TonsilA](https://omero.hms.harvard.edu/pathviewer/viewer/#?slide=514684&q_dataset=3174&tabs=SD&cgs=3509&x=4448&y=4256&t=0&z=0&zoom=3.482&rotate=0)  [TonsilB](https://omero.hms.harvard.edu/pathviewer/viewer/#?slide=514687&q_dataset=3175&tabs=SD&cgs=3619&x=4448&y=4256&t=0&z=0&zoom=3.482&rotate=0)  
Unzip and place dataTONSILA.csv and dataTONSILB.csv in Tonsil/Data/. Export each image from the viewer as TIFF and save it in the same folder as TONSILA.ome.tif.tif / TONSILB.ome.tif  

##### **seqFISH mouse hippocampus**   
Shah et al., 2016  
Paper source : https://www.cell.com/neuron/fulltext/S0896-6273(16)30702-4  

Single Cell Data: Table S6 [mmc6.xlsx](https://www.cell.com/cms/10.1016/j.neuron.2016.10.001/attachment/759be4dc-04a6-4a58-b6f6-9b52be2802db/mmc6.xlsx)  

place mmc6.xlsx in seqFISH/Data/

##### **IMC**   
Kuett et al., 2022  
Paper source : https://www.nature.com/articles/s43018-021-00301-w  

Single Cell Data : [MainHer2BreastCancerModel.zip](https://zenodo.org/records/4752030)  

unzip and place the mask and model201710_mean_intensities.csv in IMC/Data/  

##### **CRC sections**  
Lin et al., 2023  
Paper source : https://www.cell.com/cell/fulltext/S0092-8674(22)01571-9  


Download images : see 01-download-images/01_download-htan-files.md

Antibody panels: Table S3 [mmc3.xlsx](https://www.cell.com/cms/10.1016/j.cell.2022.12.028/attachment/2e9e88ba-5a67-4ee1-bb45-dcbb923bc7c0/mmc3.xlsx)

Place mmc3.xlsx in CRC/Data/ 

Cell data : BLAKE ASK SHU 


#### **Scripts for Analysis**

##### **CyCIF tonsil dataset**
Run in order 01 → 04. All inputs should be placed in Tonsil/Data/ and outputs will be stored in Tonsil/Data/

|Script|Description|Input|Output|
|--------|-------------|-------|--------|
|affine_registration_01.m|Loads the saved affine transform (Tonsil A image → Tonsil B image; the optimizer sweep that produced it is kept, commented out) and applies it to Tonsil A cell coordinates.|TONSILA.ome.tif.tif TONSILB.ome.tif dataTONSILA.csv dataTONSILB.csv affine_sweep_results.mat (included)| dataTONSILA_registered.csv|
|Reorder_02.Rmd|Reorders the Tonsil B channel columns so that cycle N corresponds to the same antibody in both sections.|dataTONSILB.csv|dataTONSILB_reordered.csv|
|Flat-field_correction_03.Rmd|Removes Tonsil A cells outside the Tonsil B field of view and cells near tile edges; log-transforms marker intensities and applies a flat-field correction.| dataTONSILA_registered.csv dataTONSILB_reordered.csv| TONSILA_aligned_corrected.csv TONSILB_reordered_corrected.csv|
|tonsil_pipeline_04.Rmd|PCA on 30 channels; spatial correlation by neighbour rank within and between sections (1000 nearest neighbours); exponential fits; predicted noise fraction per PC from repeat-stained markers; Figure 1 panels.| TONSILA_aligned_corrected.csv TONSILB_reordered_corrected.csv |Figure 1 panels|



##### **seqFISH mouse hippocampus** 
All inputs should be placed in seqFISH/Data/ and outputs will be stored in seqFISH/Data/

|Script|Description|Input|Output|
|--------|-------------|-------|--------|
|Seqfish_pipeline_01.Rmd| PCA on the 214 barcoded genes; spatial correlation by nearest-neighbour rank within each field of view; exponential decay fits; predicted noise fraction per PC; Figure 2 panels| mmc6.xlsx| Figure 2 panels| 


##### **IMC** 
Run in order 01 → 02. All inputs should be placed in IMC/Data/ and outputs will be stored in IMC/Data/

| Script | Description | Input | Output |
|--------|-------------|-------|--------|
|IMC_make_csv_01.py| Creates the CSV used in IMC_pipeline_02.Rmd from downloaded data|full_mask_final_segmentation_hwatershed_500.00_90%.tif model201710_mean_intensities.csv| model201710_cells_with_positions.csv|
|IMC_pipeline_02.Rmd| PCA on 24 markers; spatial correlation by neighbour rank within slices and between slices at every gap (300 nearest neighbours); amplitude A and length scale L from the within-slice curves; mean absolute error between within-slice and between-slice curves; figure panels.|model201710_cells_with_positions.csv|Figure 3 panels|



##### **CRC sections**
?Run in order 01 → 04?. Scripts 01 and 02 work on the downloaded images in their own folder. For script 04, place mmc3.xlsx in CRC/Data/ and the two cell tables in CRC/Data/Raw_CRC_CSV/; figure panels are written to CRC/Data/Plots/

|Script | Description | Input | Output |
|--------|-------------|-------|--------|
|01-download-images/01_download-htan-files.md|Instructions for downloading the four section images with the Gen3 client (needs an NCI Data Commons API key).| gen3_manifest.json (included) |WD-76845-044.ome.tif WD-76845-045.ome.tif WD-76845-046.ome.tif WD-76845-047.ome.tif|
|02-align-images/02_align-sequential-shared.py|Registers the four sections onto WD-76845-044 and writes the warped channels. Python, run with pixi; see the README in that folder.| the four OME-TIFFs |one single-channel TIFF per section per channel, all on a common grid|
|03_segmentation.Rmd | segmentation of cells from tifs|aligned tiffs|Sardana_CRC_46.csv Sardana_CRC_47.csv  |
|04_crc_analysis.Rmd|Quality control; PCA on 24 tumor markers (section 47); mean marker intensity of the 40 nearest non-epithelial cells (section 46); CCA; permutation test; Multivariate regression; figure panels.|mmc3.xlsx, Sardana_CRC_46.csv, Sardana_CRC_47.csv| Figure 4 panels|

#### **Requirements**

R 4.5.3 (run on Windows 11)

| Package | Version | Used by |
|---------|---------|---------|
| tidyverse | 2.0.0 | tonsil, IMC, CRC |
| dplyr | 1.2.0 | all |
| ggplot2 | 4.0.2 | all |
| stringr | 1.6.0 | seqFISH, CRC |
| purrr | 1.2.1 | seqFISH |
| readxl | 1.5.0 | seqFISH, CRC |
| cowplot | 1.2.0 | all |
| ggrepel | 0.9.8 | seqFISH |
| patchwork | 1.3.2 | seqFISH |
| ggnewscale | 0.5.2 | seqFISH, IMC |
| RANN | 2.6.2 | all |
| scales | 1.4.0 | tonsil, seqFISH, IMC |
| svglite | 2.2.2 | all (writing .svg files) |
| zoo | 1.8-15 | CRC |
| pheatmap | 1.0.13 | CRC |
| reshape2 | 1.4.5 | CRC |

MATLAB R2026a (version 26.1) with Image Processing Toolbox 26.1, for affine_registration_01.m (run on Windows 11)  

Python 3.12.10, for IMC_make_csv_01.py: tifffile 2026.9.20, pandas 3.0.6, scikit-image 0.26.0, numpy 2.5.3 (run on Windows 11)    

Python, for the CRC alignment: the environment is defined in 02-align-images/pixi.toml and pixi.lock. Install pixi, then run `pixi install` in that folder. The lock file is solved for 64-bit Windows; see the README in that folder.

Gen3 client and an NCI Data Commons API key, for downloading the CRC images.

