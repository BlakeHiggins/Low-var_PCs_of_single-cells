### Script for Analysis of CyCIF CRC dataset for "Are low-variance Principal Components of single-cell omics data biologically meaningful, or just noise?"
### Paper source:  https://www.cell.com/cell/fulltext/S0092-8674(22)01571-9
### Inputs: gen3_manifest.json
### outputs: WD-76845-044.ome.tif WD-76845-045.ome.tif WD-76845-046.ome.tif WD-76845-047.ome.tif


# Download Selected Files

- Available through NCI CRDC General Commons (CRDC-GC)

Follow the instructions below on how to access data from each of these sources. Further details are avaliable in the [HTAN Manual](https://docs.humantumoratlas.org).

* * *

#### **Files Available through NCI CRDC General Commons (CRDC-GC)**

The files listed below are available without additional access requirements.

```txt
 WD-76845-047.ome.tif
 WD-76845-046.ome.tif
 WD-76845-045.ome.tif
 WD-76845-044.ome.tif
```

Use the `gen3_manifest.json` included in this folder.

1. Install the [Gen3 Client](https://gen3.org/resources/user/gen3-client/#1-installation-instructions)
2. Get your [NCI Data Commons Framework Services API Key](https://nci-crdc.datacommons.io/identity)
3. Configure a Gen3 profile by running:

   ```bash
   gen3-client configure \
       --profile=htan \
       --cred=~/.gen3/credentials.json \
       --apiendpoint=https://nci-crdc.datacommons.io
   ```


Run the following `gen3` command.

```bash
gen3-client download-multiple \
    --profile=htan \
    --manifest=gen3_manifest.json \
    --download-path=my_htan_dir
```
