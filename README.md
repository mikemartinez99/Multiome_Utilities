# Multiome Utilities
Custom functions used for the pre-processing / analysis of scATAC, scRNA, and/or 10X multiome data.

# Table of Contents
- [scATAC](#scatac)
- [scRNA](#scrna)
- [Multiome](#multiome)

# scATAC
Functions for various preprocessing tasks. These functions can be easily automated and applied across multiple samples in parallel using the `mapply` function. 

**run_atac_qc_filtered**
The `run_atac_qc_filtered` creates a `ChromatinAssay` from the `Signac` package and adds it to a `Seurat Object` from `Seurat`

Arguments:

**1.** samples: A sample ID

**2.** files_dir: the path to files directory which will then help construct the path to the filtered h5 matrix

**.3.** annotations: A `GRanges` annotation object




