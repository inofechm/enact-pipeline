"""Class for defining methods to package pipeline outputs into AnnData objects
"""

import os
import pandas as pd
import anndata
import scanpy as sc
import seaborn as sns
import numpy as np

## Attempt to import celltypist, and prompt installation if not found
import celltypist
from celltypist import models

from .pipeline import ENACT


class CellTypistPipeline(ENACT):
    """Class for running CellAssign algorithm"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run_cell_typist(self):
        """Runs CellTypist"""
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
        adata.obs

        lib_size = adata.X.sum(1)
        adata.obs["size_factor"] = lib_size / np.mean(lib_size)
        adata.obs["lib_size"] = lib_size

        #  normalize adata to the log1p normalised format (to 10,000 counts per cell)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # download celltypist model and predict cell type
        if ".pkl" not in self.cell_typist_model:
            self.cell_typist_model = self.cell_typist_model + ".pkl"
        models.download_models(model=self.cell_typist_model)
        predictions = celltypist.annotate(adata, model=self.cell_typist_model)
        adata = predictions.to_adata(
            insert_labels=True, insert_conf=True, insert_prob=True
        )

        adata.obs.rename(columns={"predicted_labels": "cell_type"}, inplace=True)
        adata.obs[adata.obsm["spatial"].columns] = adata.obsm["spatial"]
        adata.obs[adata.obsm["stats"].columns] = adata.obsm["stats"]
        adata.obs["chunk_name"] = cell_lookup_df["chunk_name"]
        results_df = adata.obs.drop(columns=adata.obs["cell_type"].unique().tolist())
        results_df.to_csv(
            os.path.join(self.cellannotation_results_dir, "merged_results.csv")
        )


if __name__ == "__main__":
    # Creating CellAssignPipeline object
    cell_typist = CellTypistPipeline(configs_path="config/configs.yaml")
    cell_typist.run_cell_typist()
