# ENACT: End-to-End Analysis and Cell Type Annotation for Visium High Definition (HD) Slides

>[!NOTE]
>This is the official repo for [ENACT](https://www.biorxiv.org/content/10.1101/2024.10.17.618905v1). The manuscript can be accessed through [BiorXiv](https://www.biorxiv.org/content/10.1101/2024.10.17.618905v1).

Spatial transcriptomics (ST) enables the study of gene expression within its spatial context in histopathology samples. To date, a limiting factor has been the resolution of sequencing based ST products. The introduction of the Visium High Definition (HD) technology opens the door to cell resolution ST studies. However, challenges remain in the ability to accurately map transcripts to cells and in cell type assignment based on spot data.

ENACT is the first tissue-agnostic pipeline that integrates advanced cell segmentation with Visium HD transcriptomics data to infer cell types across whole tissue sections. Our pipeline incorporates novel bin-to-cell assignment methods, enhancing the accuracy of single-cell transcript estimates. Validated on diverse synthetic and real datasets, our approach demonstrates high effectiveness at predicting cell types and scalability, offering a robust solution for spatially resolved transcriptomics analysis.

This repository has the code for inferring cell types from the sub-cellular transcript counts provided by VisiumHD.

This can be achieved through the following steps:

1. **Cell segmentation**: segment high resolution image using NN-based image segmentation networks such as Stardist.
2. **Bin-to-cell assignment**: Obtain cell-wise transcript counts by aggregating the VisiumHD bins that are associated with each cell
3. **Cell type inference**: Use the cell-wise transcript counts to infer the cell labels/ phenotypes using methods used for single-cell RNA seq analysis ([CellAsign](https://www.nature.com/articles/s41592-019-0529-1#:~:text=CellAssign%20uses%20a%20probabilistic%20model%20to%20assign%20single) or [CellTypist](https://pubmed.ncbi.nlm.nih.gov/35549406/#:~:text=To%20systematically%20resolve%20immune%20cell%20heterogeneity%20across%20tissues,) or [Sargent](https://www.sciencedirect.com/science/article/pii/S2215016123001966#:~:text=We%20present%20Sargent,%20a%20transformation-free,%20cluster-free,%20single-cell%20annotation) if installed) or novel approaches, and use comprehensive cell marker databases ([Panglao](https://panglaodb.se/index.html) or [CellMarker](http://xteam.xbio.top/CellMarker/) can be used as reference).

>[!NOTE]
>At this time, Sargent is currently not available in GitHub. For information on how to access [Sargent](https://doi.org/10.1016/j.mex.2023.102196) (doi: https://doi.org/10.1016/j.mex.2023.102196), please contact the paper's corresponding authors (nima.nouri@sanofi.com). We provide the results obtained by Sargent in [ENACT's Zenodo page](https://zenodo.org/records/13887921) under the following folders:
>- ENACT_supporting_files/public_data/human_colorectal/paper_results/chunks/naive/sargent_results/
>- ENACT_supporting_files/public_data/human_colorectal/paper_results/chunks/weighted_by_area/sargent_results/
>- ENACT_supporting_files/public_data/human_colorectal/paper_results/chunks/weighted_by_transcript/sargent_results/
>- ENACT_supporting_files/public_data/human_colorectal/paper_results/chunks/weighted_by_cluster/sargent_results/

<!-- 
<div style="text-align: center;">
  <img src="figs/pipelineflow.png" alt="ENACT"/>
</div> -->
![plot](figs/pipelineflow.png)

## Index of Instructions:
- [System Requirements](#system-requirements)
- [Install ENACT from Source](#install-enact-from-source)
- [Install ENACT with Pip](#install-enact-with-pip)
- [Input Files for ENACT](#input-files-for-enact)
- [Defining ENACT Configurations](#defining-enact-configurations)
- [Output Files for ENACT](#output-files-for-enact)
- [Running ENACT from Notebook](#running-enact-from-notebook)
- [Running ENACT from Terminal](#running-enact-from-terminal)
- [Running Instructions](#running-instructions)
- [Reproducing Paper Results](#reproducing-paper-results)
- [Creating Synthetic VisiumHD Datasets](#creating-synthetic-visiumhd-datasets)
- [Citing ENACT](#citing-enact)

## System Requirements
ENACT was tested with the following specifications:
* Hardware Requirements: 32 CPU, 64GB RAM, 100 GB (harddisk and memory requirements may vary depending on whole slide image size; if the weight of the wsi is small the memory requirements can be significantly decreased)

* Software: Python 3.9, (Optional) GPU (CUDA 11)

## Install ENACT from Source 
### Step 1: Clone ENACT repository
```
git clone https://github.com/Sanofi-OneAI/oneai-dda-spatialtr-enact.git
cd oneai-dda-spatialtr-enact
```
### Step 2: Setup Python environment
Start by defining the location and the name of the Conda environment in the `Makefile`:
```
ENV_DIR := /home/oneai/envs/   <---- Conda environment location
PY_ENV_NAME := enact_py_env    <---- Conda environment name
```
Next, run the following Make command to create a Conda environment with all of ENACT's dependencies
```
make setup_py_env
```

## Install ENACT with Pip
ENACT can be installed from [Pypi](https://pypi.org/project/enact-SO/) using:
```
pip install enact-SO
```

## Input Files for ENACT
ENACT requires only three files, which can be obtained from SpaceRanger’s outputs for each experiment:

1. **Whole resolution tissue image**. This will be segmented to obtain the cell boundaries that will be used to aggregate the transcript counts.
2. **tissue_positions.parquet**. This is the file that specifies the *2um* Visium HD bin locations relative to the full resolution image.
3. **filtered_feature_bc_matrix.h5**. This is the .h5 file with the *2um* Visium HD bin counts.

## Defining ENACT Configurations
ENACT users can choose to specify the configurations via one of two ways:

1. Passing them within the class constructor:
```
  from enact.pipeline import ENACT

  so_hd = ENACT(
      cache_dir="/home/oneai/test_cache",
      wsi_path="Visium_HD_Human_Colon_Cancer_tissue_image.btf",
      visiumhd_h5_path="binned_outputs/square_002um/filtered_feature_bc_matrix.h5",
      tissue_positions_path="binned_outputs/square_002um/spatial/tissue_positions.parquet",
      )
```
<details>
  <summary><strong>Full list of ENACT parameters (click to expand)</strong></summary>

  ## Parameters

  - **cache_dir (str)**:  
    Directory to cache ENACT results. This must be specified by the user.

  - **wsi_path (str)**:  
    Path to the Whole Slide Image (WSI) file. This must be provided by the user.

  - **visiumhd_h5_path (str)**:  
    Path to the Visium HD h5 file containing spatial transcriptomics data. This 
    must be provided by the user.

  - **tissue_positions_path (str)**:  
    Path to the tissue positions file that contains spatial locations of barcodes. 
    This must be provided by the user.

  - **analysis_name (str)**:  
    Name of the analysis, used for output directories and results.  
    *Default*: `"enact_demo"`.

  - **seg_method (str)**:  
    Cell segmentation method.  
    *Default*: `"stardist"`.  
    *Options*: `["stardist"]`.

  - **patch_size (int)**:  
    Size of patches (in pixels) to process the image. Use a smaller patch size to 
    reduce memory requirements.  
    *Default*: `4000`.

  - **use_hvg (bool)**:  
    Whether to use highly variable genes (HVG) during the analysis.  
    *Default*: `True`.  
    *Options*: `[True]`.

  - **n_hvg (int)**:  
    Number of highly variable genes to use if `use_hvg` is `True`.  
    *Default*: `1000`.

  - **n_clusters (int)**:  
    Number of clusters. Used only if `bin_to_cell_method` is `"weighted_by_cluster"`.  
    *Default*: `4`.

  - **bin_representation (str)**:  
    Representation type for VisiumHD bins.  
    *Default*: `"polygon"`.  
    *Options*: `["polygon"]`.

  - **bin_to_cell_method (str)**:  
    Method to assign bins to cells.  
    *Default*: `"weighted_by_cluster"`.  
    *Options*: `["naive", "weighted_by_area", "weighted_by_gene", "weighted_by_cluster"]`.

  - **cell_annotation_method (str)**:  
    Method for annotating cell types.  
    *Default*: `"celltypist"`.  
    *Options*: `["celltypist", "sargent" (if installed), "cellassign"]`.

  - **cell_typist_model (str)**:  
    Path to the pre-trained CellTypist model for cell type annotation. Only used if 
    `cell_annotation_method` is `"celltypist"`.  
    Refer to [CellTypist Models](https://www.celltypist.org/models) for a list of 
    available models.  
    *Default*: `""` (empty string).

  - **run_synthetic (bool)**:  
    Whether to run synthetic data generation for testing purposes.  
    *Default*: `False`.

  - **segmentation (bool)**:  
    Flag to run the image segmentation step.  
    *Default*: `True`.

  - **bin_to_geodataframes (bool)**:  
    Flag to convert the bins to GeoDataFrames.  
    *Default*: `True`.

  - **bin_to_cell_assignment (bool)**:  
    Flag to run bin-to-cell assignment.  
    *Default*: `True`.

  - **cell_type_annotation (bool)**:  
    Flag to run cell type annotation.  
    *Default*: `True`.

  - **cell_markers (dict)**:  
    A dictionary of cell markers used for annotation. Only used if `cell_annotation_method` 
    is one of `["sargent", "cellassign"]`.

  - **chunks_to_run (list)**:  
    Specific chunks of data to run the analysis on, typically for debugging.  
    *Default*: `[]` (runs all chunks).

  - **configs_dict (dict)**:  
    Dictionary containing ENACT configuration parameters. If provided, the values 
    in `configs_dict` will override any corresponding parameters passed directly 
    to the class constructor. This is useful for running ENACT with a predefined 
    configuration for convenience and consistency.  
    *Default*: `{}` (uses the parameters specified in the class constructor).

</details>

2. Specifying configurations in a `yaml` file: (sample file located under `config/configs.yaml`):
```yaml
    analysis_name: <analysis-name>                              <---- custom name for analysis. Will create a folder with that name to store the results
    run_synthetic: False                                        <---- True if you want to run bin to cell assignment on synthetic dataset, False otherwise
    cache_dir: <path-to-store-enact-outputs>                    <---- path to store pipeline outputs
    paths:
        wsi_path: <path-to-whole-slide-image>                   <---- path to whole slide image
        visiumhd_h5_path: <path-to-counts-file>                 <---- location of the 2um x 2um gene by bin file (filtered_feature_bc_matrix.h5) from 10X Genomics 
        tissue_positions_path: <path-to-tissue-positions>       <---- location of the tissue of the tissue_positions.parquet file from 10X genomicsgenomics
    steps:
        segmentation: True                                      <---- True if you want to run segmentation
        bin_to_geodataframes: True                              <---- True to convert bin to geodataframes
        bin_to_cell_assignment: True                            <---- True to bin-to-cell assignment
        cell_type_annotation: True                              <---- True to run cell type annotation
    params:
      bin_to_cell_method: "weighted_by_cluster"                 <---- bin-to-cell assignment method. Pick one of ["naive", "weighted_by_area", "weighted_by_gene", "weighted_by_cluster"]
      cell_annotation_method: "celltypist"                      <---- cell annotation method. Pick one of ["cellassign", "celltypist"]
      cell_typist_model: "Human_Colorectal_Cancer.pkl"          <---- CellTypist model weights to use. Update based on organ of interest if using cell_annotation_method is set to "celltypist"
      seg_method: "stardist"                                    <---- cell segmentation method. Stardist is the only option for now
      patch_size: 4000                                          <---- defines the patch size. The whole resolution image will be broken into patches of this size. Reduce if you run into memory issues
      use_hvg: True                                             <---- True only run analysis on top n highly variable genes. Setting it to False runs ENACT on all genes in the counts file
      n_hvg: 1000                                               <---- number of highly variable genes to use 
      n_clusters: 4                                             <---- number of cell clusters to use for the "weighted_by_cluster" method. Default is 4.
  cell_markers:                                                 <---- cell-gene markers to use for cell annotation. Only applicable if params/cell_annotation_method is "cellassign" or "sargent"
      Epithelial: ["CDH1","EPCAM","CLDN1","CD2"]
      Enterocytes: ["CD55", "ELF3", "PLIN2", "GSTM3", "KLF5", "CBR1", "APOA1", "CA1", "PDHA1", "EHF"]
      Goblet cells: ["MANF", "KRT7", "AQP3", "AGR2", "BACE2", "TFF3", "PHGR1", "MUC4", "MUC13", "GUCA2A"]
```

## Output Files for ENACT
ENACT outputs all its results under the `cache` directory which gets automatically created at run time:
```
.
└── cache/
    └── <anaylsis_name> /
        ├── chunks/
        │   ├── bins_gdf/
        │   │   └── patch_<patch_id>.csv
        │   ├── cells_gdf/
        │   │   └── patch_<patch_id>.csv
        │   └── <bin_to_cell_method>/
        │       ├── bin_to_cell_assign/
        │       │   └── patch_<patch_id>.csv
        │       ├── cell_ix_lookup/
        │       │   └── patch_<patch_id>.csv
        │       └── <cell_annotation_method>_results/
        │           ├── cells_adata.csv
        │           └── merged_results.csv
        └── cells_df.csv
```
ENACT breaks down the whole resolution image into "chunks" (or patches) of size `patch_size`. Results are provided per-chunk under the `chunks` directory.
* `bins_gdf`:Folder containing GeoPandas dataframes representing the 2um Visium HD bins within a given patch
* `cells_gdf`: Folder containing GeoPandas dataframes representing cells segmented in the tissue
* `<bin_to_cell_method>/bin_to_cell_assign`: Folder contains dataframes with the transcripts assigned to each cells
* `<bin_to_cell_method>/cell_ix_lookup`: Folder contains dataframes defining the indices and coordinates of the cells
* `<bin_to_cell_method>/<cell_annotation_method>_results/cells_adata.csv`: Anndata object containing the results from ENACT (cell coordinates, cell types, transcript counts)
* <`bin_to_cell_method>/<cell_annotation_method>_results/merged_results.csv`: Dataframe (.csv) containing the results from ENACT (cell coordinates, cell types)

## Running ENACT from Notebook
The [demo notebook](ENACT_demo.ipynb) provides a step-by-step guide on how to install and run ENACT on VisiumHD public data using notebook.

## Running ENACT from Terminal
This section provides a guide for running ENACT on the [Human Colorectal Cancer sample](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc) provided on 10X Genomics' website.
### Step 1: Install ENACT from Source 
Refer to [Install ENACT from Source](#install-enact-from-source)

### Step 2: Download the necessary files from the 10X Genomics website:

1.  Whole slide image: full resolution tissue image
```
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_tissue_image.btf
```

2. Visium HD output file. The transcript counts are provided in a .tar.gz file that needs to be extracted:
```
curl -O https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_binned_outputs.tar.gz
tar -xvzf Visium_HD_Human_Colon_Cancer_binned_outputs.tar.gz
```
Locate the following two files from the extracted outputs file.
```
.
└── binned_outputs/
    └── square_002um/
        ├── filtered_feature_bc_matrix.h5   <---- Transcript counts file (2um resolution)
        └── spatial/
            └── tissue_positions.parquet    <---- Bin locations relative to the full resolution image
```

### Step 3: Update input file locations and parameters under `config/configs.yaml`

Refer to [Running Instructions](#running-instructions) for a full list of ENACT parameters to change.

Below is a sample configuration file to use to run ENACT on the Human Colorectal cancer sample:

```yaml
analysis_name: "colon-demo"
run_synthetic: False # True if you want to run bin to cell assignment on synthetic dataset, False otherwise.
cache_dir: "cache/ENACT_outputs"                                                                          # Change according to your desired output location
paths:  
  wsi_path: "<path_to_data>/Visium_HD_Human_Colon_Cancer_tissue_image.btf"                                # whole slide image path
  visiumhd_h5_path: "<path_to_data>/binned_outputs/square_002um/filtered_feature_bc_matrix.h5"            # location of the 2um x 2um gene by bin file (filtered_feature_bc_matrix.h5) from 10X Genomics.   
  tissue_positions_path: "<path_to_data>/binned_outputs/square_002um/spatial/tissue_positions.parquet"    # location of the tissue of the tissue_positions.parquet file from 10X genomics
steps:
  segmentation: True # True if you want to run segmentation
  bin_to_geodataframes: True # True to convert bin to geodataframes
  bin_to_cell_assignment: True # True to assign cells to bins
  cell_type_annotation: True # True to run cell type annotation
params:
  seg_method: "stardist" # Stardist is the only option for now
  patch_size: 4000 # Defines the patch size. The whole resolution image will be broken into patches of this size
  bin_representation: "polygon"  # or point TODO: Remove support for anything else
  bin_to_cell_method: "weighted_by_cluster" # or naive
  cell_annotation_method: "celltypist"
  cell_typist_model: "Human_Colorectal_Cancer.pkl"
  use_hvg: True # Only run analysis on highly variable genes + cell markers specified
  n_hvg: 1000 # Number of highly variable genes to use
  n_clusters: 4 
  chunks_to_run: []
cell_markers:
  # Human Colon
  Epithelial: ["CDH1","EPCAM","CLDN1","CD2"]
  Enterocytes: ["CD55", "ELF3", "PLIN2", "GSTM3", "KLF5", "CBR1", "APOA1", "CA1", "PDHA1", "EHF"]
  Goblet cells: ["MANF", "KRT7", "AQP3", "AGR2", "BACE2", "TFF3", "PHGR1", "MUC4", "MUC13", "GUCA2A"]
  Enteroendocrine cells: ["NUCB2", "FABP5", "CPE", "ALCAM", "GCG", "SST", "CHGB", "IAPP", "CHGA", "ENPP2"]
  Crypt cells: ["HOPX", "SLC12A2", "MSI1", "SMOC2", "OLFM4", "ASCL2", "PROM1", "BMI1", "EPHB2", "LRIG1"]
  Endothelial: ["PECAM1","CD34","KDR","CDH5","PROM1","PDPN","TEK","FLT1","VCAM1","PTPRC","VWF","ENG","MCAM","ICAM1","FLT4"]     
  Fibroblast: ["COL1A1","COL3A1","COL5A2","PDGFRA","ACTA2","TCF21","FN"]
  Smooth muscle cell: ["BGN","MYL9","MYLK","FHL2","ITGA1","ACTA2","EHD2","OGN","SNCG","FABP4"]
  B cells: ["CD74", "HMGA1", "CD52", "PTPRC", "HLA-DRA", "CD24", "CXCR4", "SPCS3", "LTB", "IGKC"]
  T cells: ["JUNB", "S100A4", "CD52", "PFN1P1", "CD81", "EEF1B2P3", "CXCR4", "CREM", "IL32", "TGIF1"]
  NK cells: ["S100A4", "IL32", "CXCR4", "FHL2", "IL2RG", "CD69", "CD7", "NKG7", "CD2", "HOPX"]

```

## Running Instructions
This section provides a guide on running ENACT on your own data
### Step 1: Install ENACT from Source 
Refer to [Install ENACT from Source](#install-enact-from-source)

### Step 2: Define the Location of ENACT's Required Files
Define the locations of ENACT's required files in the `config/configs.yaml` file. Refer to [Input Files for ENACT](#input-files-for-enact)
```yaml
    analysis_name: <analysis-name>                              <---- custom name for analysis. Will create a folder with that name to store the results
    cache_dir: <path-to-store-enact-outputs>                    <---- path to store pipeline outputs
    paths:
        wsi_path: <path-to-whole-slide-image>                   <---- path to whole slide image
        visiumhd_h5_path: <path-to-counts-file>                 <---- location of the 2um x 2um gene by bin file (filtered_feature_bc_matrix.h5) from 10X Genomics. 
        tissue_positions_path: <path-to-tissue-positions>       <---- location of the tissue of the tissue_positions.parquet file from 10X genomics
```

### Step 3: Define ENACT configurations
Define the following core parameters in the `config/configs.yaml` file:
```yaml
    params:
      bin_to_cell_method: "weighted_by_cluster"                 <---- bin-to-cell assignment method. Pick one of ["naive", "weighted_by_area", "weighted_by_gene", "weighted_by_cluster"]
      cell_annotation_method: "celltypist"                      <---- cell annotation method. Pick one of ["cellassign", "celltypist", "sargent" (if installed)]
      cell_typist_model: "Human_Colorectal_Cancer.pkl"          <---- CellTypist model weights to use. Update based on organ of interest if using cell_annotation_method is set to
```
Refer to [Defining ENACT Configurations](#defining-enact-configurations) for a full list of parameters to configure. If using CellTypist, set `cell_typist_model` to one of the following models based on the organ and species under study: [CellTypist models](https://www.celltypist.org/models#:~:text=CellTypist%20was%20first%20developed%20as%20a%20platform%20for). 

### Step 4: Define Cell Gene Markers (Only applies for cell_annotation_method is "cellassign" or "sargent")
Define the cell gene markers in `config/configs.yaml` file. Those can be expert annotated or obtained from open-source databases such as [Panglao](https://panglaodb.se/index.html) or [CellMarker](http://xteam.xbio.top/CellMarker/). Example cell markers for human colorectal cancer samples:
```yaml
  cell_markers:
    Epithelial: ["CDH1","EPCAM","CLDN1","CD2"]
    Enterocytes: ["CD55", "ELF3", "PLIN2", "GSTM3", "KLF5", "CBR1", "APOA1", "CA1", "PDHA1", "EHF"]
    Goblet cells: ["MANF", "KRT7", "AQP3", "AGR2", "BACE2", "TFF3", "PHGR1", "MUC4", "MUC13", "GUCA2A"]
    Enteroendocrine cells: ["NUCB2", "FABP5", "CPE", "ALCAM", "GCG", "SST", "CHGB", "IAPP", "CHGA", "ENPP2"]
    Crypt cells: ["HOPX", "SLC12A2", "MSI1", "SMOC2", "OLFM4", "ASCL2", "PROM1", "BMI1", "EPHB2", "LRIG1"]
    Endothelial: ["PECAM1","CD34","KDR","CDH5","PROM1","PDPN","TEK","FLT1","VCAM1","PTPRC","VWF","ENG","MCAM","ICAM1","FLT4"]     
    Fibroblast: ["COL1A1","COL3A1","COL5A2","PDGFRA","ACTA2","TCF21","FN"]
    Smooth muscle cell: ["BGN","MYL9","MYLK","FHL2","ITGA1","ACTA2","EHD2","OGN","SNCG","FABP4"]
    B cells: ["CD74", "HMGA1", "CD52", "PTPRC", "HLA-DRA", "CD24", "CXCR4", "SPCS3", "LTB", "IGKC"]
    T cells: ["JUNB", "S100A4", "CD52", "PFN1P1", "CD81", "EEF1B2P3", "CXCR4", "CREM", "IL32", "TGIF1"]
    NK cells: ["S100A4", "IL32", "CXCR4", "FHL2", "IL2RG", "CD69", "CD7", "NKG7", "CD2", "HOPX"]
```
### Step 5: Run ENACT
```
make run_enact
```

## Reproducing Paper Results
This section provides a guide on how to reproduce the ENACT paper results on the [10X Genomics Human Colorectal Cancer VisumHD sample](https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc). 
Here, ENACT is run on various combinations of bin-to-cell assignment methods and cell annotation algorithms.

### Step 1: Install ENACT from Source 
Refer to [Install ENACT from Source](#install-enact-from-source)

### Step 2: Run ENACT on combinations of bin-to-cell assignment methods and cell annotation algorithms
3. Run the following command which will download all the supplementary file from [ENACT's Zenodo page](https://zenodo.org/records/13887921) and programmatically run ENACT with various combinations of bin-to-cell assignment methods and cell annotation algorithms:
```
make reproduce_results
```

## Creating Synthetic VisiumHD Datasets

1. To create synthetic VisiumHD dataset from Xenium or seqFISH+ data, run and follow the instructions of the notebooks in [src/synthetic_data](src/synthetic_data).

2. To run the ENACT pipeline with the synthetic data, set the following parameters in the `config/configs.yaml` file: 

```yaml
run_synthetic: True                                        <---- True if you want to run bin to cell assignment on synthetic dataset, False otherwise.
```
    
3. Run ENACT:
```
make run_enact
```

## Citing ENACT
If you use this repository or its tools in your research, please cite the following:
```
@article {Kamel2024.10.17.618905,
	author = {Kamel, Mena and Song, Yiwen and Solbas, Ana and Villordo, Sergio and Sarangi, Amrut and Senin, Pavel and Mathew, Sunaal and Ayestas, Luis Cano and Wang, Seqian and Classe, Marion and Bar-Joseph, Ziv and Planas, Albert Pla},
	title = {ENACT: End-to-End Analysis of Visium High Definition (HD) Data},
	elocation-id = {2024.10.17.618905},
	year = {2024},
	doi = {10.1101/2024.10.17.618905},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {Motivation: Spatial transcriptomics (ST) enables the study of gene expression within its spatial context in histopathology samples. To date, a limiting factor has been the resolution of sequencing based ST products. The introduction of the Visium High Definition (HD) technology opens the door to cell resolution ST studies. However, challenges remain in the ability to accurately map transcripts to cells and in assigning cell types based on the transcript data. Results: We developed ENACT, the first tissue-agnostic pipeline that integrates advanced cell segmentation with Visium HD transcriptomics data to infer cell types across whole tissue sections. Our pipeline incorporates novel bin-to-cell assignment methods, enhancing the accuracy of single-cell transcript estimates. Validated on diverse synthetic and real datasets, our approach is both scalable and effective offering a robust solution for spatially resolved transcriptomics analysis. Availability and implementation: ENACT source code is available at https://github.com/Sanofi-Public/enact-pipeline. Experimental data is available at https://zenodo.org/records/13887921. Supplementary information: Supplementary data are available at BiorXiv online.Competing Interest StatementThe authors have declared no competing interest.},
	URL = {https://www.biorxiv.org/content/early/2024/10/20/2024.10.17.618905},
	eprint = {https://www.biorxiv.org/content/early/2024/10/20/2024.10.17.618905.full.pdf},
	journal = {bioRxiv}
}
```
