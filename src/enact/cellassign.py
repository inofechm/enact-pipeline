"""Class for defining methods to package pipeline outputs into AnnData objects
"""

import os
import pandas as pd
import anndata
import scanpy as sc
import scvi
import seaborn as sns
from scvi.external import CellAssign
import numpy as np
import torch

from .pipeline import ENACT

seed = 42


class CellAssignPipeline(ENACT):
    """Class for running CellAssign algorithm"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.configs = configs

    def format_markers_to_df(self):
        """Method to format marker genes to a pandas dataframe
        num gene x num cell_types
        """
        markers_dict = self.configs["cell_markers"]
        genes_set = set([item for sublist in markers_dict.values() for item in sublist])
        markers_df = pd.DataFrame(columns=markers_dict.keys(), index=sorted(genes_set))
        markers_df = markers_df.fillna(0)
        for cell_type, gene_markers in markers_dict.items():
            markers_df.loc[gene_markers, cell_type] = 1
        self.markers_df = markers_df

    def run_cell_assign(self):
        """Runs CellAssign"""
        bin_assign_df = self.merge_files(self.bin_assign_dir, save=False)
        cell_lookup_df = self.merge_files(self.cell_ix_lookup_dir, save=False)

        spatial_cols = ["cell_x", "cell_y"]
        stat_columns = ["num_shared_bins", "num_unique_bins", "num_transcripts"]
        cell_lookup_df.loc[:, "id"] = cell_lookup_df["id"].astype(str)
        cell_lookup_df = cell_lookup_df.set_index("id")
        cell_lookup_df["num_transcripts"] = cell_lookup_df["num_transcripts"].fillna(0)

        bin_assign_df.index = cell_lookup_df.index
        bin_assign_df = bin_assign_df.drop(columns=["Unnamed: 0"])
        bin_assign_df = bin_assign_df.fillna(0).astype(int)
        adata = anndata.AnnData(bin_assign_df.astype(int))

        adata.obsm["spatial"] = cell_lookup_df[spatial_cols].astype(int)
        adata.obsm["stats"] = cell_lookup_df[stat_columns].astype(int)

        lib_size = adata.X.sum(1)
        adata.obs["size_factor"] = lib_size / np.mean(lib_size)
        adata.obs["lib_size"] = lib_size

        marker_gene_mat = self.markers_df.copy()
        marker_gene_mat = marker_gene_mat.loc[
            sorted(list(set(self.markers_df.index) & set(bin_assign_df.columns)))
        ]
        bdata = adata[:, marker_gene_mat.index].copy()

        torch.manual_seed(seed)
        scvi.external.CellAssign.setup_anndata(bdata, size_factor_key="size_factor")
        model = CellAssign(bdata, marker_gene_mat, random_b_g_0=False)
        model.train()
        predictions = model.predict()

        bdata.obs["cell_type"] = predictions.idxmax(axis=1).values
        bdata.obs[adata.obsm["spatial"].columns] = adata.obsm["spatial"]
        bdata.obs[adata.obsm["stats"].columns] = adata.obsm["stats"]
        bdata.obs["chunk_name"] = cell_lookup_df["chunk_name"]
        bdata.obs.to_csv(
            os.path.join(self.cellannotation_results_dir, "merged_results.csv")
        )
        print(
            f"saved to : {os.path.join(self.cellannotation_results_dir, 'merged_results.csv')}"
        )


if __name__ == "__main__":
    # Creating CellAssignPipeline object
    cell_assign = CellAssignPipeline(configs_path="config/configs.yaml")
    cell_assign.format_markers_to_df()
