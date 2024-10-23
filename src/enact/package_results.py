"""Class for defining methods to package pipeline outputs into AnnData objects
"""

import os

import anndata
import pandas as pd

# import squidpy as sq

from .pipeline import ENACT


class PackageResults(ENACT):
    """Class for packaging ENACT pipeline outputs"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.files_to_ignore = [
            "merged_results.csv",
            "merged_results_old.csv",
            "cells_adata.h5",
            ".ipynb_checkpoints",
        ]

    def merge_cellassign_output_files(self):
        """Merges the CellAssign results with gene counts

        Returns:
            _type_: _description_
        """
        chunks = os.listdir(self.bin_assign_dir)
        cell_by_gene_list = []
        for chunk_name in chunks:
            if chunk_name in self.files_to_ignore:
                continue
            index_lookup = pd.read_csv(
                os.path.join(self.cell_ix_lookup_dir, chunk_name)
            )
            trancript_counts = pd.read_csv(
                os.path.join(self.bin_assign_dir, chunk_name)
            ).drop(columns=["Unnamed: 0"])
            cell_by_gene_chunk = pd.concat(
                [index_lookup["id"], trancript_counts], axis=1
            )
            cell_by_gene_list.append(cell_by_gene_chunk)
        cell_by_gene_df = pd.concat(cell_by_gene_list, axis=0)
        return cell_by_gene_df

    def merge_sargent_output_files(self):
        """Merges the Sargent chunk results into a single results file

        Returns:
            _type_: _description_
        """
        os.makedirs(self.sargent_results_dir, exist_ok=True)
        # Merge the sargent_results_chunks data and gene_to_cell_assignment_chunks_ix_lookup
        chunks = os.listdir(self.sargent_results_dir)
        sargent_results_list = []
        cell_by_gene_list = []
        for chunk_name in chunks:
            if chunk_name in self.files_to_ignore:
                continue
            cell_labels = pd.read_csv(
                os.path.join(self.sargent_results_dir, chunk_name)
            )
            index_lookup = pd.read_csv(
                os.path.join(self.cell_ix_lookup_dir, chunk_name)
            )
            trancript_counts = pd.read_csv(
                os.path.join(self.bin_assign_dir, chunk_name)
            ).drop(columns=["Unnamed: 0"])

            sargent_result_chunk = pd.concat([index_lookup, cell_labels["x"]], axis=1)
            cell_by_gene_chunk = pd.concat(
                [index_lookup["id"], trancript_counts], axis=1
            )
            sargent_result_chunk.drop("Unnamed: 0", axis=1, inplace=True)
            sargent_results_list.append(sargent_result_chunk)
            cell_by_gene_list.append(cell_by_gene_chunk)
        sargent_results_df = pd.concat(sargent_results_list, axis=0)
        sargent_results_df = sargent_results_df.rename(columns={"x": "cell_type"})
        cell_by_gene_df = pd.concat(cell_by_gene_list, axis=0)
        sargent_results_df.to_csv(
            os.path.join(self.sargent_results_dir, "merged_results.csv"), index=False
        )
        return sargent_results_df, cell_by_gene_df

    def df_to_adata(self, results_df, cell_by_gene_df):
        """Converts pd.DataFrame object with pipeline results to AnnData

        Args:
            results_df (_type_): _description_

        Returns:
            anndata.AnnData: Anndata with pipeline outputs
        """
        file_columns = results_df.columns
        spatial_cols = ["cell_x", "cell_y"]
        stat_columns = ["num_shared_bins", "num_unique_bins", "num_transcripts"]
        results_df.loc[:, "id"] = results_df["id"].astype(str)
        results_df = results_df.set_index("id")
        results_df["num_transcripts"] = results_df["num_transcripts"].fillna(0)
        results_df["cell_type"] = results_df["cell_type"].str.lower()
        adata = anndata.AnnData(cell_by_gene_df.set_index("id").astype(int))

        # adata = anndata.AnnData(results_df[stat_columns].astype(int))
        adata.obsm["spatial"] = results_df[spatial_cols].astype(int)
        adata.obsm["stats"] = results_df[stat_columns].astype(int)
        # This column is the output of cell type inference pipeline
        adata.obs["cell_type"] = results_df[["cell_type"]].astype("category")
        adata.obs["patch_id"] = results_df[["chunk_name"]]
        return adata

    # def run_neighborhood_enrichment(self, adata):
    #     """Sample function to run Squidpy operations on AnnData object

    #     Args:
    #         adata (_type_): _description_

    #     Returns:
    #         _type_: _description_
    #     """
    #     sq.gr.spatial_neighbors(adata)
    #     sq.gr.nhood_enrichment(adata, cluster_key="cell_type")
    #     return adata

    def save_adata(self, adata):
        """Save the anndata object to disk

        Args:
            adata (_type_): _description_
        """
        adata.write(
            os.path.join(self.cellannotation_results_dir, "cells_adata.h5"),
            compression="gzip",
        )


if __name__ == "__main__":
    # Creating ENACT object
    so_hd = PackageResults(configs_path="config/configs.yaml")
    results_df, cell_by_gene_df = so_hd.merge_sargent_output_files()
    adata = so_hd.df_to_adata(results_df, cell_by_gene_df)
    # adata = so_hd.run_neighborhood_enrichment(adata) # Example integration with SquiPy
    so_hd.save_adata(adata)
