# Weighted by area method
import anndata
import numpy as np
import pandas as pd
from scipy import sparse
from tqdm import tqdm
from sklearn.cluster import KMeans


def apply_weights_to_adata_counts(expanded_adata, weights_df):
    """Applies the weights to the counts matrix

    Args:
        adata (AnnData): Counts AnnData

    Returns:
        AnnData: Weighted-adjusted AnnData
    """
    if weights_df.empty:
        return expanded_adata
    # Applying the weighting
    mask = (expanded_adata.obs_names.isin(weights_df.index)) & (
        expanded_adata.obs["id"].isin(weights_df["id"])
    )
    indices = np.where(mask)[0]
    # Apply weights to the entries in the expression matrix
    weights_matrix = np.ones(expanded_adata.shape)
    for idx in tqdm(indices, total=len(indices)):
        bin_id = expanded_adata.obs.iloc[idx]["index"]
        cell_id = expanded_adata.obs.iloc[idx]["id"]
        bin_rows = weights_df.loc[bin_id]
        weights = bin_rows[bin_rows["id"] == cell_id][expanded_adata.var_names]
        weights_matrix[idx] = weights.iloc[0].tolist()
    weighted_counts = expanded_adata.X.multiply(weights_matrix)
    # convert back to sparse
    expanded_adata.X = sparse.csr_matrix(weighted_counts)
    return expanded_adata


def weight_by_gene_assignment(
    result_spatial_join, expanded_adata, unique_cell_by_gene_adata
):
    # Getting the gene counts of the cells (unique signature for each cell)
    gene_counts_non_overlap = (
        pd.DataFrame(
            unique_cell_by_gene_adata.X.toarray(),
            index=unique_cell_by_gene_adata.obs_names,
            columns=unique_cell_by_gene_adata.var_names,
        )
        .groupby(unique_cell_by_gene_adata.obs["id"])
        .sum()
        .reset_index()
    )

    # Getting the bins that overlap with multiple cells
    overlapping_bins = result_spatial_join[~result_spatial_join["unique_bin"]]

    # Getting a table of bins with the parent cell and the parent cell's gene content
    overlap_merge = pd.merge(
        overlapping_bins[["index", "id"]], gene_counts_non_overlap, on="id", how="left"
    )
    overlap_merge.set_index("index", inplace=True)

    # Grouping the bins by the bin id
    grouped_overlap = overlap_merge.groupby("index")

    # Initialize progress bar for processing overlapping bins
    pbar = tqdm(grouped_overlap, desc="Processing overlapping bins", unit="bin")
    gene_columns = overlap_merge.columns.drop(["id"]).tolist()
    weights_list = []
    # Looping through the bins and splitting the counts
    for bin_index, group_rows in pbar:
        # getting total gene counts from the cells that share a bin
        gene_total = group_rows[gene_columns].sum(axis=0)
        # Dividing the cells gene counts by the total gene counts to get the weight
        gene_weights = group_rows[gene_columns].div(gene_total, axis=1).fillna(0)
        gene_weights["id"] = group_rows["id"]
        weights_list.append(gene_weights)
    # Getting a weights dataframe
    if weights_list:
        weights_df = pd.concat(weights_list, axis=0)
    else:
        weights_df = pd.DataFrame()
    pbar.close()
    expanded_adata = apply_weights_to_adata_counts(expanded_adata, weights_df)
    return result_spatial_join, expanded_adata


def weight_by_cluster_assignment(
    result_spatial_join, expanded_adata, unique_cell_by_gene_adata, n_clusters=4
):
    # Getting the gene counts of the cells (unique signature for each cell)
    gene_counts_non_overlap = (
        pd.DataFrame(
            unique_cell_by_gene_adata.X.toarray(),
            index=unique_cell_by_gene_adata.obs_names,
            columns=unique_cell_by_gene_adata.var_names,
        )
        .groupby(unique_cell_by_gene_adata.obs["id"])
        .sum()
        .reset_index()
    )

    # Getting the bins that overlap with multiple cells
    overlapping_bins = result_spatial_join[~result_spatial_join["unique_bin"]]

    gene_columns = gene_counts_non_overlap.columns.drop(["id"]).tolist()

    # clustering on gene counts from non-overlapping bins
    n_clusters = np.min([n_clusters, len(gene_counts_non_overlap)])
    kmeans = KMeans(n_clusters=n_clusters, random_state=0)
    clusters = kmeans.fit_predict(gene_counts_non_overlap[gene_columns])
    gene_counts_non_overlap["cluster"] = clusters
    cluster_means = gene_counts_non_overlap.groupby("cluster")[gene_columns].mean()

    # Getting a table of bins with the parent cell and the parent cell's gene content
    overlap_merge = pd.merge(
        overlapping_bins[["index", "id"]], gene_counts_non_overlap, on="id", how="left"
    )
    # merge cluster mean gene counts with overlapping bins -
    # using cluster gene counts instead of the bins's gene counts
    overlap_merge = pd.merge(
        overlap_merge[["index", "id", "cluster"]],
        cluster_means,
        left_on="cluster",
        right_index=True,
        how="left",
    )
    overlap_merge.set_index("index", inplace=True)

    grouped_overlap = overlap_merge.groupby("index")

    # Initialize progress bar for processing overlapping bins
    pbar = tqdm(grouped_overlap, desc="Processing overlapping bins", unit="bin")
    weights_list = []
    # Looping through the bins and splitting the counts
    for bin_index, group_rows in pbar:
        # getting total gene counts from the cells that share a bin
        gene_total = group_rows[gene_columns].sum(axis=0)
        # Dividing the cells gene counts by the total gene counts to get the weight
        gene_weights = group_rows[gene_columns].div(gene_total, axis=1).fillna(0)
        gene_weights["id"] = group_rows["id"]
        weights_list.append(gene_weights)
    # Getting a weights dataframe
    if weights_list:
        weights_df = pd.concat(weights_list, axis=0)
    else:
        weights_df = pd.DataFrame()
    pbar.close()
    expanded_adata = apply_weights_to_adata_counts(expanded_adata, weights_df)
    return result_spatial_join, expanded_adata
