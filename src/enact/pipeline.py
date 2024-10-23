"""Class for defining methods for VisiumHD pipeline
"""

import os
from csbdeep.utils import normalize
import geopandas as gpd
import numpy as np
import pandas as pd
from PIL import Image
import scanpy as sc
from scipy import sparse
import shapely
from shapely.geometry import Polygon, Point
from shapely import wkt
from stardist.models import StarDist2D
import tifffile as tifi
from tqdm import tqdm
import yaml
import anndata
import scvi
import seaborn as sns
from scvi.external import CellAssign
import logging
import ssl
import argparse

Image.MAX_IMAGE_PIXELS = None
from .utils.logging import get_logger
from .assignment_methods.naive import naive_assignment
from .assignment_methods.weight_by_area import weight_by_area_assignment
from .assignment_methods.weight_by_gene import (
    weight_by_gene_assignment,
    weight_by_cluster_assignment,
)

class ENACT:
    """Class for methods for the ENACT pipeline"""

    def __init__(
        self,
        cache_dir="",
        wsi_path="",
        visiumhd_h5_path="",
        tissue_positions_path="",
        analysis_name="enact_demo",
        seg_method="stardist",
        patch_size=4000,
        use_hvg=True,
        n_hvg=1000,
        n_clusters=4,
        bin_representation="polygon",
        bin_to_cell_method="weighted_by_cluster",
        cell_annotation_method="celltypist",
        cell_typist_model="",
        run_synthetic=False,
        segmentation=True,
        bin_to_geodataframes=True,
        bin_to_cell_assignment=True,
        cell_type_annotation=True,
        cell_markers={},
        chunks_to_run=[],
        configs_dict={}
    ):
        """
Initialize the class with the following parameters:

Args:
    cache_dir (str): Directory to cache ENACT results. This must be populated 
        by the user.
    wsi_path (str): Path to the Whole Slide Image (WSI) file. This must be 
        populated by the user.
    visiumhd_h5_path (str): Path to the Visium HD h5 file containing spatial 
        transcriptomics data. This must be populated by the user.
    tissue_positions_path (str): Path to the tissue positions file that 
        contains spatial locations of barcodes. This must be populated by the 
        user.
    analysis_name (str): Name of the analysis, used for output directories and 
        results. Default is "enact_demo".
    seg_method (str): Cell segmentation method. Default is "stardist". 
        Options: ["stardist"].
    patch_size (int): Size of patches (in pixels) to process the image. Use a 
        smaller patch size to reduce memory requirements. Default is 4000.
    use_hvg (bool): Whether to use highly variable genes (HVG) during the 
        analysis. Default is True. Options: [True].
    n_hvg (int): Number of highly variable genes to use if `use_hvg` is True. 
        Default is 1000.
    n_clusters (int): Number of clusters. Only used if `bin_to_cell_method` is 
        "weighted_by_cluster". Default is 4.
    bin_representation (str): Representation type for VisiumHD bins. Default is 
        "polygon". Options: ["polygon"].
    bin_to_cell_method (str): Method to assign bins to cells. Default is 
        "weighted_by_cluster". Options: ["naive", "weighted_by_area", 
        "weighted_by_gene", "weighted_by_cluster"].
    cell_annotation_method (str): Method for annotating cell types. Default is 
        "celltypist". Options: ["celltypist", "sargent" (if installed), 
        "cellassign"].
    cell_typist_model (str): Path to the pre-trained CellTypist model for cell 
        type annotation. Only used if `cell_annotation_method` is 
        "celltypist". Refer to https://www.celltypist.org/models for a full 
        list of models. Default is an empty string.
    run_synthetic (bool): Whether to run synthetic data generation for testing 
        purposes. Default is False.
    segmentation (bool): Flag to run the image segmentation step. Default is 
        True.
    bin_to_geodataframes (bool): Flag to convert the bins to GeoDataFrames. 
        Default is True.
    bin_to_cell_assignment (bool): Flag to run bin-to-cell assignment. Default 
        is True.
    cell_type_annotation (bool): Flag to run cell type annotation. Default is 
        True.
    cell_markers (dict): A dictionary of cell markers used for annotation. Only 
        used if `cell_annotation_method` is one of ["sargent", "cellassign"].
    chunks_to_run (list): Specific chunks of data to run the analysis on for 
        debugging purposes. Default is an empty list (runs all chunks).
    configs_dict (dict): Dictionary containing ENACT configuration parameters. 
        If provided, the values in `configs_dict` will override any 
        corresponding parameters passed directly to the class constructor. This 
        is useful for running ENACT with a predefined configuration for 
        convenience and consistency. Default is an empty dictionary (i.e., 
        using the parameters defined in the class constructor).
"""
  
        # Todo: add class documentation
        user_configs = {
            "analysis_name": analysis_name,
            "run_synthetic": run_synthetic,
            "cache_dir": cache_dir,
            "paths": {
                "wsi_path": wsi_path,
                "visiumhd_h5_path": visiumhd_h5_path,
                "tissue_positions_path": tissue_positions_path,
            },
            "steps": {
                "segmentation": segmentation,
                "bin_to_geodataframes": bin_to_geodataframes,
                "bin_to_cell_assignment": bin_to_cell_assignment,
                "cell_type_annotation": cell_type_annotation,
            },
            "params": {
                "seg_method": seg_method,
                "patch_size": patch_size,
                "bin_representation":bin_representation,
                "bin_to_cell_method": bin_to_cell_method,
                "cell_annotation_method": cell_annotation_method,
                "cell_typist_model": cell_typist_model,
                "use_hvg": use_hvg,
                "n_hvg": n_hvg,
                "n_clusters": n_clusters,
                "chunks_to_run": chunks_to_run,
            },
            "cell_markers": cell_markers,
        }
        self.configs = user_configs
        if configs_dict != {}:
            # If user specifies a configs_dict -> it overwrites all other user-specified parameters
            self.overwrite_configs(configs_dict)

        if self.configs["cache_dir"] == "":
            raise ValueError(f"Error: Please provide a value for 'cache_dir'.")
        
        if self.configs["params"]["cell_annotation_method"] == "celltypist":
            if self.configs["params"]["cell_typist_model"] == "":
                raise ValueError(
                    f"Error: Please provide a value for 'cell_typist_model'. "
                    "Refer to https://www.celltypist.org/models for a full list of models."
                )
        
        if self.configs["params"]["cell_annotation_method"] in ["sargent", "cellassign"]:
            if self.configs["cell_markers"] == {}:
                raise ValueError(
                    f"Error: Please provide a value for 'cell_markers'."
                )

        # Load input files
        core_paths = ["wsi_path", "visiumhd_h5_path", "tissue_positions_path"]
        for core_path in core_paths:
            if self.configs["paths"][core_path] == "":
                raise ValueError(
                    f"Error: Please provide a value for '{core_path}'. "
                )
        self.initiate_instance_variables()
        self.load_configs()
    
    def overwrite_configs(self, configs_dict):
        """
        Function overwrites the configurations with the content in the configs_dict
        """
        for key, value in configs_dict.items():
            if key in ["cell_markers"] or not isinstance(value, dict):
                self.configs[key] = value
            else:
                for sub_key, sub_value in value.items():
                    self.configs[key][sub_key] = sub_value

    def initiate_instance_variables(self):
        """
        Creates instance variables for all the keys defined in the config file
        """
        run_details = "ENACT running with the following configurations: \n"
        kwargs = {}
        for key, value in self.configs.items():
            if key in ["cell_markers"] or not isinstance(value, dict):
                setattr(self, key, value)
                run_details+= f" {key}: {value}\n"
                kwargs[key] = value
            else:
                for sub_key, sub_value in value.items():
                    setattr(self, sub_key, sub_value)
                    run_details+= f" {sub_key}: {sub_value}\n"
                    kwargs[sub_key] = sub_value
        self.run_details = run_details
        self.kwargs = kwargs

    def load_configs(self):
        """Loading the configuations and parameters"""
        # Generating paths
        self.cache_dir = os.path.join(self.cache_dir, self.analysis_name)
        self.cells_df_path = os.path.join(self.cache_dir, "cells_df.csv")
        self.cells_layer_path = os.path.join(self.cache_dir, "cells_layer.png")
        self.cell_chunks_dir = os.path.join(self.cache_dir, "chunks", "cells_gdf")
        self.bin_chunks_dir = os.path.join(self.cache_dir, "chunks", "bins_gdf")
        self.bin_assign_dir = os.path.join(
            self.cache_dir, "chunks", self.bin_to_cell_method, "bin_to_cell_assign"
        )
        self.cell_ix_lookup_dir = os.path.join(
            self.cache_dir, "chunks", self.bin_to_cell_method, "cell_ix_lookup"
        )
        os.makedirs(self.cache_dir, exist_ok=True)
        self.cellannotation_results_dir = os.path.join(
            self.cache_dir,
            "chunks",
            self.bin_to_cell_method,
            f"{self.cell_annotation_method}_results",
        )
        os.makedirs(self.cellannotation_results_dir, exist_ok=True)
        os.makedirs(self.cell_chunks_dir, exist_ok=True)
        self.logger = get_logger("ENACT", self.cache_dir)
        self.logger.info(f"<initiate_instance_variables> {self.run_details}")

    def load_image(self, file_path=None):
        """Load image from given file path
        Arguments:
            file_path {string} : path to the file that we are trying to load
        Returns:
            np.array -- loaded image as numpy array
        """
        if file_path == None:
            file_path = self.wsi_path
        img_arr = tifi.imread(file_path)
        crop_bounds = self.get_image_crop_bounds()
        x_min, y_min, x_max, y_max = crop_bounds
        img_arr = img_arr[y_min:y_max, x_min:x_max, :]
        self.logger.info("<load_image> Successfully loaded image!")
        return img_arr, crop_bounds

    def get_image_crop_bounds(self, file_path=None):
        """Get the crop location of the image to adjust the coordinates accordingly

        Args:
            file_path (_type_): _description_

        Returns:
            _type_: _description_
        """
        if file_path == None:
            file_path = self.wsi_path
        tissue_pos_list = pd.read_parquet(self.tissue_positions_path)

        # Cleaning up, removing negative coords,removing out of tissue bins
        tissue_pos_list_filt = tissue_pos_list[tissue_pos_list.in_tissue == 1]
        tissue_pos_list_filt["pxl_row_in_fullres"] = tissue_pos_list_filt[
            "pxl_row_in_fullres"
        ].astype(int)
        tissue_pos_list_filt["pxl_col_in_fullres"] = tissue_pos_list_filt[
            "pxl_col_in_fullres"
        ].astype(int)
        tissue_pos_list_filt = tissue_pos_list_filt.loc[
            (tissue_pos_list_filt.pxl_row_in_fullres >= 0)
            & (tissue_pos_list_filt.pxl_col_in_fullres >= 0)
        ]
        x_min = tissue_pos_list_filt["pxl_col_in_fullres"].min()
        y_min = tissue_pos_list_filt["pxl_row_in_fullres"].min()
        x_max = tissue_pos_list_filt["pxl_col_in_fullres"].max()
        y_max = tissue_pos_list_filt["pxl_row_in_fullres"].max()
        return (x_min, y_min, x_max, y_max)

    def normalize_image(self, image, min_percentile=5, max_percentile=95):
        """_summary_

        Args:
            image (_type_): _description_
            min_percentile (int, optional): _description_. Defaults to 5.
            max_percentile (int, optional): _description_. Defaults to 95.

        Returns:
            _type_: _description_
        """
        # Adjust min_percentile and max_percentile as needed
        image_norm = normalize(image, min_percentile, max_percentile)
        self.logger.info("<normalize_image> Successfully normalized image!")
        return image_norm

    def segment_cells(self, image, prob_thresh=0.005):
        """_summary_

        Args:
            image (_type_): _description_
            prob_thresh (float, optional): _description_. Defaults to 0.005.

        Returns:
            _type_: _description_
        """
        if self.seg_method == "stardist":
            # Adjust nms_thresh and prob_thresh as needed
            # ssl._create_default_https_context = ssl._create_unverified_context
            self.stardist_model = StarDist2D.from_pretrained("2D_versatile_he")
            labels, polys = self.stardist_model.predict_instances_big(
                image,
                axes="YXC",
                block_size=4096,
                prob_thresh=prob_thresh,
                nms_thresh=0.001,
                min_overlap=128,
                context=128,
                normalizer=None,
                n_tiles=(4, 4, 1),
            )
            self.logger.info("<run_segmentation> Successfully segmented cells!")
            return labels, polys
        else:
            self.logger.warning("<run_segmentation> Invalid cell segmentation model!")
            return None, None

    def convert_stardist_output_to_gdf(self, cell_polys, save_path=None):
        """Convert stardist output to geopandas dataframe

        Args:
            cell_polys (_type_): _description_
            save_path (_type_, optional): _description_. Defaults to None.
        """
        if save_path == None:
            save_path = self.cells_df_path
        # Creating a list to store Polygon geometries
        geometries = []
        centroids = []
        cell_x, cell_y = [], []

        # Iterating through each nuclei in the 'polys' DataFrame
        for nuclei in range(len(cell_polys["coord"])):
            # Extracting coordinates for the current nuclei and converting them to (y, x) format
            coords = [
                (y, x)
                for x, y in zip(
                    cell_polys["coord"][nuclei][0], cell_polys["coord"][nuclei][1]
                )
            ]
            # Creating a Polygon geometry from the coordinates
            poly = Polygon(coords)
            centroid = list(poly.centroid.coords)[0]
            centroids.append(centroid)
            geometries.append(poly)
            cell_x.append(centroid[0])
            cell_y.append(centroid[1])

        # Creating a GeoDataFrame using the Polygon geometries
        gdf = gpd.GeoDataFrame(geometry=geometries)
        gdf["id"] = [f"ID_{i+1}" for i, _ in enumerate(gdf.index)]
        gdf["cell_x"] = cell_x
        gdf["cell_y"] = cell_y
        gdf["centroid"] = centroids
        # Save results to disk
        gdf.to_csv(save_path)
        return gdf

    def split_df_to_chunks(self, df, x_col, y_col, output_dir):
        """
        Break the cells df into files, of size patch_size x patch_size

        Args:
            df (_type_): _description_
        """
        os.makedirs(output_dir, exist_ok=True)
        i = 0
        # Need to divide into chunks of patch_size pixels by patch_size pixels
        df[["patch_x", "patch_y"]] = (df[[x_col, y_col]] / self.patch_size).astype(int)
        df["patch_id"] = df["patch_x"].astype(str) + "_" + df["patch_y"].astype(str)
        self.logger.info(
            f"<split_df_to_chunks> Splitting into chunks. output_dir: {output_dir}"
        )
        unique_patches = df.patch_id.unique()
        for patch_id in tqdm(unique_patches, total=len(unique_patches)):
            patch_cells = df[df.patch_id == patch_id]
            if len(patch_cells) == 0:
                continue
            patch_cells.to_csv(os.path.join(output_dir, f"patch_{patch_id}.csv"))

    def load_visiumhd_dataset(self, crop_bounds):
        """Loads the VisiumHD dataset and adjusts the
        coordinates to the cropped image

        Args:
            crop_bounds (tuple): crop bounds

        Returns:
            AnnData: AnnData object with the VisiumHD data
            int: bin size in pixels
        """

        # Accounting for crop bounds
        if crop_bounds is not None:
            x1, y1, _, _ = crop_bounds
        else:
            x1, y1 = (0, 0)
        # Load Visium HD data
        adata = sc.read_10x_h5(self.visiumhd_h5_path)

        # Load the Spatial Coordinates
        df_tissue_positions = pd.read_parquet(self.tissue_positions_path)

        # Set the index of the dataframe to the barcodes
        df_tissue_positions = df_tissue_positions.set_index("barcode")

        # Create an index in the dataframe to check joins
        df_tissue_positions["index"] = df_tissue_positions.index

        # *Important step*: Representing coords in the cropped WSI frame
        df_tissue_positions["pxl_row_in_fullres"] = (
            df_tissue_positions["pxl_row_in_fullres"] - y1
        )
        df_tissue_positions["pxl_col_in_fullres"] = (
            df_tissue_positions["pxl_col_in_fullres"] - x1
        )
        # Adding the tissue positions to the meta data
        adata.obs = pd.merge(
            adata.obs, df_tissue_positions, left_index=True, right_index=True
        )

        first_row = df_tissue_positions[
            (df_tissue_positions["array_row"] == 0)
            & (df_tissue_positions["array_col"] == 0)
        ]["pxl_col_in_fullres"]
        second_row = df_tissue_positions[
            (df_tissue_positions["array_row"] == 0)
            & (df_tissue_positions["array_col"] == 1)
        ]["pxl_col_in_fullres"]
        bin_size = second_row[0] - first_row[0]
        if self.configs["params"]["use_hvg"]:
            # Keeping the top n highly variable genes + the user requested cell markers
            n_genes = self.configs["params"]["n_hvg"]
            # Normalizing to median total counts
            adata_norm = adata.copy()
            sc.pp.normalize_total(adata_norm)
            # Logarithmize the data
            sc.pp.log1p(adata_norm)
            sc.pp.highly_variable_genes(adata_norm, n_top_genes=n_genes)

            hvg_mask = adata_norm.var["highly_variable"]
            cell_markers = [
                item
                for sublist in self.configs["cell_markers"].values()
                for item in sublist
            ]
            missing_markers = set(cell_markers) - set(hvg_mask.index)
            self.logger.info(
                f"<load_visiumhd_dataset> Missing the following markers: {missing_markers}"
            )
            available_markers = list(set(cell_markers) & set(hvg_mask.index))
            hvg_mask.loc[available_markers] = True
            adata = adata[:, hvg_mask]
        return adata, bin_size

    def generate_bin_polys(self, bins_df, x_col, y_col, bin_size):
        """Represents the bins as Shapely polygons

        Args:
            bins_df (pd.DataFrame): bins dataframe
            x_col (str): column with the bin centre x-coordinate
            y_col (str): column with the bin centre y-coordinate
            bin_size (int): bin size in pixels

        Returns:
            list: list of Shapely polygons
        """
        geometry = []
        # Generates Shapely polygons to represent each bin
        if self.bin_representation == "point":
            # Geometry column is just the centre (x, y) for a VisiumHD bin
            geometry = [Point(xy) for xy in zip(bins_df[x_col], bins_df[y_col])]
        elif self.bin_representation == "polygon":
            self.logger.info(
                f"<generate_bin_polys> Generating bin polygons. num_bins: {len(bins_df)}"
            )
            half_bin_size = bin_size / 2
            bbox_coords = pd.DataFrame(
                {
                    "min_x": bins_df[x_col] - half_bin_size,
                    "min_y": bins_df[y_col] - half_bin_size,
                    "max_x": bins_df[x_col] + half_bin_size,
                    "max_y": bins_df[y_col] + half_bin_size,
                }
            )
            geometry = [
                shapely.geometry.box(min_x, min_y, max_x, max_y)
                for min_x, min_y, max_x, max_y in tqdm(
                    zip(
                        bbox_coords["min_x"],
                        bbox_coords["min_y"],
                        bbox_coords["max_x"],
                        bbox_coords["max_y"],
                    ),
                    total=len(bins_df),
                )
            ]
        else:
            self.logger.warning("<generate_bin_polys> Select a valid mode!")
        return geometry

    def convert_adata_to_cell_by_gene(self, adata):
        """Converts the AnnData object from bin-by-gene to
        cell-by-gene AnnData object.

        Args:
            adata (AnnData): bin-by-gene AnnData

        Returns:
            AnnData: cell-by-gene AnnData
        """
        # Group the data by unique cell IDs
        groupby_object = adata.obs.groupby(["id"], observed=True)

        # Extract the gene expression counts from the AnnData object
        counts = adata.X.tocsr()

        # Obtain the number of unique nuclei and the number of genes in the expression data
        N_groups = groupby_object.ngroups
        N_genes = counts.shape[1]

        # Initialize a sparse matrix to store the summed gene counts for each nucleus
        summed_counts = sparse.lil_matrix((N_groups, N_genes))

        # Lists to store the IDs of polygons and the current row index
        polygon_id = []
        row = 0
        # Iterate over each unique polygon to calculate the sum of gene counts.
        for polygons, idx_ in groupby_object.indices.items():
            summed_counts[row] = counts[idx_].sum(0)
            row += 1
            polygon_id.append(polygons)
        # Create an AnnData object from the summed count matrix
        summed_counts = summed_counts.tocsr()
        cell_by_gene_adata = anndata.AnnData(
            X=summed_counts,
            obs=pd.DataFrame(polygon_id, columns=["id"], index=polygon_id),
            var=adata.var,
        )
        return cell_by_gene_adata

    def generate_bins_gdf(self, adata, bin_size):
        """Convert the bins Anndata object to a geodataframe

        Args:
            adata (_type_): _description_

        Returns:
            _type_: _description_
        """
        bin_coords_df = adata.obs.copy()
        geometry = self.generate_bin_polys(
            bins_df=bin_coords_df,
            x_col="pxl_col_in_fullres",
            y_col="pxl_row_in_fullres",
            bin_size=bin_size,
        )
        bins_gdf = gpd.GeoDataFrame(bin_coords_df, geometry=geometry)
        return bins_gdf

    def assign_bins_to_cells(self, adata):
        """Assigns bins to cells based on method requested by the user

        Args:
            adata (_type_): _description_
        """
        os.makedirs(self.bin_assign_dir, exist_ok=True)
        os.makedirs(self.cell_ix_lookup_dir, exist_ok=True)
        if self.configs["params"]["chunks_to_run"]:
            chunk_list = self.configs["params"]["chunks_to_run"]
        else:
            chunk_list = os.listdir(self.cell_chunks_dir)
        self.logger.info(
            f"<assign_bins_to_cells> Assigning bins to cells using {self.bin_to_cell_method} method"
        )
        for chunk in tqdm(chunk_list, total=len(chunk_list)):
            if os.path.exists(os.path.join(self.cell_ix_lookup_dir, chunk)):
                continue
            if chunk in [".ipynb_checkpoints"]:
                continue
            if not os.path.exists(os.path.join(self.bin_chunks_dir, chunk)):
                continue

            # Loading the cells geodataframe
            cell_gdf_chunk_path = os.path.join(self.cell_chunks_dir, chunk)
            cell_gdf_chunk = gpd.GeoDataFrame(pd.read_csv(cell_gdf_chunk_path))
            cell_gdf_chunk["geometry"] = cell_gdf_chunk["geometry"].apply(wkt.loads)
            cell_gdf_chunk.set_geometry("geometry", inplace=True)
            # Loading the bins geodataframe
            bin_gdf_chunk_path = os.path.join(self.bin_chunks_dir, chunk)
            bin_gdf_chunk = gpd.GeoDataFrame(pd.read_csv(bin_gdf_chunk_path))
            bin_gdf_chunk["geometry"] = bin_gdf_chunk["geometry"].apply(wkt.loads)
            bin_gdf_chunk.set_geometry("geometry", inplace=True)

            # Perform a spatial join to check which coordinates are in a cell nucleus
            result_spatial_join = gpd.sjoin(
                bin_gdf_chunk,
                cell_gdf_chunk[["geometry", "id", "cell_x", "cell_y"]],
                how="left",
                predicate="intersects",
            )

            # Only keeping the bins that overlap with a cell
            result_spatial_join = result_spatial_join[
                ~result_spatial_join["index_right"].isna()
            ]

            # Getting unique bins and overlapping bins
            barcodes_in_overlaping_polygons = pd.unique(
                result_spatial_join[result_spatial_join.duplicated(subset=["index"])][
                    "index"
                ]
            )
            result_spatial_join["unique_bin"] = ~result_spatial_join["index"].isin(
                barcodes_in_overlaping_polygons
            )
            # Filter the adata object to contain only the barcodes in result_spatial_join
            # shape: (#bins_overlap x #genes)
            expanded_adata = adata[result_spatial_join["index"]]
            # Adding the cell ids to the anndata object (the cell that the bin is assigned to)
            # Can have duplicate bins (i.e. "expanded") if a bin is assigned to more than one cell
            expanded_adata.obs["id"] = result_spatial_join["id"].tolist()

            # Reshape the anndata object to (#cells x #genes)
            filtered_result_spatial_join = result_spatial_join[
                result_spatial_join["unique_bin"]
            ]
            filtered_adata = adata[filtered_result_spatial_join["index"]]
            filtered_adata.obs["id"] = filtered_result_spatial_join["id"].tolist()

            unfilt_result_spatial_join = result_spatial_join.copy()
            self.logger.info("<assign_bins_to_cells> done spatial join")

            if result_spatial_join.empty:
                self.logger.info(
                    "result_spatial_join is empty, skipping bin-to-cell assignment."
                )
                continue

            elif self.bin_to_cell_method == "naive":
                result_spatial_join = naive_assignment(result_spatial_join)
                expanded_adata = filtered_adata.copy()

            elif self.bin_to_cell_method == "weighted_by_area":
                result_spatial_join, expanded_adata = weight_by_area_assignment(
                    result_spatial_join, expanded_adata, cell_gdf_chunk
                )

            elif self.bin_to_cell_method == "weighted_by_gene":
                unique_cell_by_gene_adata = self.convert_adata_to_cell_by_gene(
                    filtered_adata
                )
                result_spatial_join, expanded_adata = weight_by_gene_assignment(
                    result_spatial_join, expanded_adata, unique_cell_by_gene_adata
                )

            elif self.bin_to_cell_method == "weighted_by_cluster":
                unique_cell_by_gene_adata = self.convert_adata_to_cell_by_gene(
                    filtered_adata
                )
                result_spatial_join, expanded_adata = weight_by_cluster_assignment(
                    result_spatial_join,
                    expanded_adata,
                    unique_cell_by_gene_adata,
                    n_clusters=self.configs["params"]["n_clusters"],
                )
            else:
                self.logger.info("ERROR", self.bin_to_cell_method)
            self.logger.info("<assign_bins_to_cells> convert_adata_to_cell_by_gene")
            cell_by_gene_adata = self.convert_adata_to_cell_by_gene(expanded_adata)
            del expanded_adata

            # Save the gene to cell assignment results to a .csv file
            chunk_gene_to_cell_assign_df = pd.DataFrame(
                cell_by_gene_adata.X.toarray(),
                columns=cell_by_gene_adata.var_names,
            )

            chunk_gene_to_cell_assign_df = chunk_gene_to_cell_assign_df.loc[
                :, ~chunk_gene_to_cell_assign_df.columns.duplicated()
            ].copy()
            # Saving counts to cache
            chunk_gene_to_cell_assign_df.to_csv(
                os.path.join(self.bin_assign_dir, chunk)
            )

            # Getting number of bins shared between cells
            overlaps_df = (
                unfilt_result_spatial_join.groupby(["id", "unique_bin"])
                .count()["in_tissue"]
                .reset_index()
            )
            overlaps_df = overlaps_df.pivot(
                index="id", columns="unique_bin", values="in_tissue"
            ).fillna(0)
            try:
                overlaps_df.columns = ["num_shared_bins", "num_unique_bins"]
            except:
                overlaps_df.columns = ["num_unique_bins"]
                overlaps_df["num_shared_bins"] = 0
            cell_gdf_chunk = cell_gdf_chunk.merge(
                overlaps_df, how="left", left_on="id", right_index=True
            )
            cell_gdf_chunk[["num_shared_bins", "num_unique_bins"]] = cell_gdf_chunk[
                ["num_shared_bins", "num_unique_bins"]
            ].fillna(0)
            # Save index lookup to store x and y values and cell index
            index_lookup_df = cell_by_gene_adata.obs.merge(
                cell_gdf_chunk, how="left", left_index=True, right_on="id"
            )[
                ["cell_x", "cell_y", "num_shared_bins", "num_unique_bins", "id"]
            ].reset_index(
                drop=True
            )
            index_lookup_df["num_transcripts"] = chunk_gene_to_cell_assign_df.sum(
                axis=1
            )
            index_lookup_df["chunk_name"] = chunk
            index_lookup_df.to_csv(os.path.join(self.cell_ix_lookup_dir, chunk))
            self.logger.info(
                f"{self.bin_to_cell_method} mean count per cell: {chunk_gene_to_cell_assign_df.sum(axis=1).mean()}"
            )

    def assign_bins_to_cells_synthetic(self):
        """Assigns bins to cells based on method requested by the user

        Args:
            adata (_type_): _description_
        """
        os.makedirs(self.bin_assign_dir, exist_ok=True)
        os.makedirs(self.cell_ix_lookup_dir, exist_ok=True)
        if self.configs["params"]["chunks_to_run"]:
            chunk_list = self.configs["params"]["chunks_to_run"]
        else:
            chunk_list = os.listdir(self.cell_chunks_dir)

        self.logger.info(
            f"<assign_bins_to_cells_synthetic> Assigning bins to cells using {self.bin_to_cell_method} method"
        )
        for chunk in tqdm(chunk_list, total=len(chunk_list)):
            if os.path.exists(os.path.join(self.cell_ix_lookup_dir, chunk)):
                continue
            if chunk in [".ipynb_checkpoints"]:
                continue
            # Loading the cells geodataframe
            cell_gdf_chunk_path = os.path.join(self.cell_chunks_dir, chunk)
            cell_gdf_chunk = gpd.GeoDataFrame(pd.read_csv(cell_gdf_chunk_path))
            cell_gdf_chunk["geometry"] = cell_gdf_chunk["geometry"].apply(wkt.loads)
            cell_gdf_chunk.set_geometry("geometry", inplace=True)
            cell_gdf_chunk["geometry"] = cell_gdf_chunk["geometry"].buffer(0)
            # Loading the bins geodataframe
            bin_gdf_chunk_path = os.path.join(self.bin_chunks_dir, chunk)
            bin_gdf_chunk = gpd.GeoDataFrame(pd.read_csv(bin_gdf_chunk_path))
            bin_gdf_chunk["geometry"] = bin_gdf_chunk["geometry"].apply(wkt.loads)
            bin_gdf_chunk.set_geometry("geometry", inplace=True)

            # Perform a spatial join to check which coordinates are in a cell nucleus
            result_spatial_join = gpd.sjoin(
                bin_gdf_chunk[["geometry", "assigned_bin_id", "row", "column"]],
                cell_gdf_chunk[["geometry", "cell_id"]],
                how="left",
                predicate="intersects",
            )

            # Only keeping the bins that overlap with a cell
            result_spatial_join = result_spatial_join[
                ~result_spatial_join["index_right"].isna()
            ]

            # Getting unique bins and overlapping bins
            barcodes_in_overlaping_polygons = pd.unique(
                result_spatial_join[
                    result_spatial_join.duplicated(subset=["assigned_bin_id"])
                ]["assigned_bin_id"]
            )
            result_spatial_join["unique_bin"] = ~result_spatial_join[
                "assigned_bin_id"
            ].isin(barcodes_in_overlaping_polygons)
            bin_gdf_chunk = bin_gdf_chunk.set_index("assigned_bin_id")
            adata = anndata.AnnData(
                bin_gdf_chunk.drop(columns=["geometry", "row", "column"])
            )

            # Filter the adata object to contain only the barcodes in result_spatial_join
            # shape: (#bins_overlap x #genes)
            expanded_adata = adata[result_spatial_join["assigned_bin_id"]]
            # Adding the cell ids to the anndata object (the cell that the bin is assigned to)
            # Can have duplicate bins (i.e. "expanded") if a bin is assigned to more than one cell
            expanded_adata.obs["id"] = result_spatial_join["cell_id"].tolist()
            expanded_adata.obs["index"] = result_spatial_join[
                "assigned_bin_id"
            ].tolist()

            # Reshape the anndata object to (#cells x #genes)
            filtered_result_spatial_join = result_spatial_join[
                result_spatial_join["unique_bin"]
            ]
            filtered_adata = adata[filtered_result_spatial_join["assigned_bin_id"]]
            filtered_adata.obs["id"] = filtered_result_spatial_join["cell_id"].tolist()
            if not sparse.issparse(filtered_adata.X):
                filtered_adata.X = sparse.csr_matrix(filtered_adata.X)
            if not sparse.issparse(expanded_adata.X):
                expanded_adata.X = sparse.csr_matrix(expanded_adata.X)
            self.logger.info("<assign_bins_to_cells> done spatial join")

            cell_gdf_chunk.rename(columns={"cell_id": "id"}, inplace=True)
            result_spatial_join.rename(
                columns={"assigned_bin_id": "index"}, inplace=True
            )
            result_spatial_join.rename(columns={"cell_id": "id"}, inplace=True)
            if self.bin_to_cell_method == "naive":
                result_spatial_join = naive_assignment(result_spatial_join)
                expanded_adata = filtered_adata.copy()

            elif self.bin_to_cell_method == "weighted_by_area":
                # cell_gdf_chunk = cell_gdf_chunk.set_index('index_right')
                result_spatial_join, expanded_adata = weight_by_area_assignment(
                    result_spatial_join, expanded_adata, cell_gdf_chunk
                )

            elif self.bin_to_cell_method == "weighted_by_gene":
                unique_cell_by_gene_adata = self.convert_adata_to_cell_by_gene(
                    filtered_adata
                )
                result_spatial_join, expanded_adata = weight_by_gene_assignment(
                    result_spatial_join, expanded_adata, unique_cell_by_gene_adata
                )

            elif self.bin_to_cell_method == "weighted_by_cluster":
                unique_cell_by_gene_adata = self.convert_adata_to_cell_by_gene(
                    filtered_adata
                )
                result_spatial_join, expanded_adata = weight_by_cluster_assignment(
                    result_spatial_join,
                    expanded_adata,
                    unique_cell_by_gene_adata,
                    n_clusters=self.configs["params"]["n_clusters"],
                )

            self.logger.info("<assign_bins_to_cells> convert_adata_to_cell_by_gene")

            if not sparse.issparse(expanded_adata.X):
                expanded_adata.X = sparse.csr_matrix(expanded_adata.X)
            cell_by_gene_adata = self.convert_adata_to_cell_by_gene(expanded_adata)
            del expanded_adata

            # Save the gene to cell assignment results to a .csv file
            chunk_gene_to_cell_assign_df = pd.DataFrame(
                cell_by_gene_adata.X.toarray(),
                columns=cell_by_gene_adata.var_names,
            )
            chunk_gene_to_cell_assign_df.insert(
                0, "id", cell_by_gene_adata.obs["id"].values
            )

            chunk_gene_to_cell_assign_df = chunk_gene_to_cell_assign_df.loc[
                :, ~chunk_gene_to_cell_assign_df.columns.duplicated()
            ].copy()

            # Saving counts to cache
            chunk_gene_to_cell_assign_df.to_csv(
                os.path.join(self.bin_assign_dir, chunk)
            )

            self.logger.info(f"{chunk} finished")

    def merge_files(
        self, input_folder, output_file_name="merged_results.csv", save=True
    ):
        """Merges all files in a specified input folder into a single output file.

        Args:
            input_folder (str): The path to the folder containing the input files to be merged.
            output_file_name (str): The name of the output file.
        """
        # List to store the DataFrames
        csv_list = []
        output_file = os.path.join(input_folder, output_file_name)

        if self.configs["params"]["chunks_to_run"]:
            chunk_list = self.configs["params"]["chunks_to_run"]
        else:
            chunk_list = os.listdir(input_folder)
        # Loop through all files in the directory
        for filename in chunk_list:
            if filename in ["annotated.csv", ".ipynb_checkpoints"]:
                continue
            if "merged" in filename:
                continue

            # Read each .csv file and append it to the list
            file_path = os.path.join(input_folder, filename)
            df = pd.read_csv(file_path)
            csv_list.append(df)

        # Concatenate all DataFrames in the list
        concatenated_df = pd.concat(csv_list, ignore_index=True)

        if save:
            # Save the concatenated DataFrame to the output file
            concatenated_df.to_csv(output_file, index=False)
            self.logger.info(
                f"<merge_files> files have been merged and saved to {output_file}"
            )
        return concatenated_df

    def run_cell_type_annotation(self):
        """Runs cell type annotation"""
        ann_method = self.configs["params"]["cell_annotation_method"]
        if ann_method == "sargent":
            self.logger.info(
                f"<run_cell_type_annotation> Will launch Sargent separately. "
                "Please ensure Sargent is installed."
            )
        elif ann_method == "cellassign":
            from .cellassign import CellAssignPipeline

            cellassign_obj = CellAssignPipeline(**self.kwargs)
            cellassign_obj.format_markers_to_df()
            cellassign_obj.run_cell_assign()
            self.logger.info(
                f"<run_cell_type_annotation> Successfully ran CellAssign on Data."
            )

        elif ann_method == "celltypist":
            from .celltypist import CellTypistPipeline
            celltypist_obj = CellTypistPipeline(**self.kwargs)

            celltypist_obj.run_cell_typist()
            self.logger.info(
                f"<run_cell_type_annotation> Successfully ran CellTypist on Data."
            )
        else:
            self.logger.info(
                "<run_cell_type_annotation> Please select a valid cell annotation "
                "method. options=['cellassign', 'sargent']"
            )

    def package_results(self):
        """Packages the results of the pipeline"""
        from .package_results import PackageResults

        pack_obj = PackageResults(**self.kwargs)
        ann_method = self.configs["params"]["cell_annotation_method"]
        if ann_method == "sargent":
            results_df, cell_by_gene_df = pack_obj.merge_sargent_output_files()
            adata = pack_obj.df_to_adata(results_df, cell_by_gene_df)
            pack_obj.save_adata(adata)
            self.logger.info("<package_results> Packaged Sargent results")
        elif ann_method == "cellassign":
            cell_by_gene_df = pack_obj.merge_cellassign_output_files()
            results_df = pd.read_csv(
                os.path.join(self.cellannotation_results_dir, "merged_results.csv")
            )
            adata = pack_obj.df_to_adata(results_df, cell_by_gene_df)
            pack_obj.save_adata(adata)
            self.logger.info("<package_results> Packaged CellAssign results")
        elif ann_method == "celltypist":
            cell_by_gene_df = pack_obj.merge_cellassign_output_files()
            results_df = pd.read_csv(
                os.path.join(self.cellannotation_results_dir, "merged_results.csv")
            )
            adata = pack_obj.df_to_adata(results_df, cell_by_gene_df)
            pack_obj.save_adata(adata)
            self.logger.info("<package_results> Packaged CellTypist results")
        else:
            self.logger.info(
                f"<package_results> Please select a valid cell annotation method"
            )

    def run_enact(self):
        """Runs ENACT given the user-specified configs"""
        if not self.run_synthetic:
            # Loading image and getting shape and cropping boundaries (if applicable)
            wsi, crop_bounds = self.load_image()

            # Run cell segmentation
            if self.segmentation:
                wsi_norm = self.normalize_image(
                    image=wsi, min_percentile=5, max_percentile=95
                )
                cell_labels, cell_polys = self.segment_cells(image=wsi_norm)
                cells_gdf = self.convert_stardist_output_to_gdf(
                    cell_polys=cell_polys, save_path=None
                )
                # Split the cells geodataframe to chunks
                self.split_df_to_chunks(
                    df=cells_gdf,
                    x_col="cell_x",
                    y_col="cell_y",
                    output_dir=self.cell_chunks_dir,
                )
                del cells_gdf, wsi  # Clearing memory
            else:
                # Load the segmentation results from cache
                cells_gdf = pd.read_csv(self.cells_df_path)
                cells_gdf["geometry"] = cells_gdf["geometry"].apply(wkt.loads)

            # Loading the VisiumHD reads
            if self.bin_to_geodataframes:
                bins_adata, bin_size = self.load_visiumhd_dataset(crop_bounds)
                # Convert VisiumHD reads to geodataframe objects
                bins_gdf = self.generate_bins_gdf(bins_adata, bin_size)
                # Splitting the bins geodataframe object
                self.split_df_to_chunks(
                    df=bins_gdf,
                    x_col="pxl_col_in_fullres",
                    y_col="pxl_row_in_fullres",
                    output_dir=self.bin_chunks_dir,
                )
                del bins_gdf

            # Run bin-to-cell assignment
            if self.bin_to_cell_assignment:
                bins_adata, bin_size = self.load_visiumhd_dataset(crop_bounds)
                self.assign_bins_to_cells(bins_adata)

            # Run cell type annotation
            if self.cell_type_annotation:
                self.run_cell_type_annotation()
                self.package_results()

        else:
            # Generating synthetic data
            if self.analysis_name in ["xenium", "xenium_nuclei"]:
                cells_gdf = pd.read_csv(self.cells_df_path)
                self.split_df_to_chunks(
                    df=cells_gdf,
                    x_col="cell_x",
                    y_col="cell_y",
                    output_dir=self.cell_chunks_dir,
                )
            self.assign_bins_to_cells_synthetic()


if __name__ == "__main__":
    # Creating ENACT object
    parser = argparse.ArgumentParser(description="Specify ENACT config file location.")
    parser.add_argument(
        "--configs_path", type=str, required=False, help="Config file location"
    )
    args = parser.parse_args()
    if not args.configs_path:
        configs_path = "config/configs.yaml"
    else:
        configs_path = args.configs_path
    print(f"<ENACT> Loading configurations from {configs_path}")
    with open(configs_path, "r") as stream:
        configs = yaml.safe_load(stream)
    so_hd = ENACT(configs_dict=configs)
    so_hd.run_enact()
