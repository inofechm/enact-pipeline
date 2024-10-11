# Weighted by area method
import anndata
import numpy as np
from scipy import sparse


def apply_weights_to_adata_counts(adata):
    """Applies the weights to the counts matrix

    Args:
        adata (AnnData): Counts AnnData

    Returns:
        AnnData: Weighted-adjusted AnnData
    """
    weight = adata.obs["weight"]
    # Reshape weights to (130000, 1) for broadcasting
    weight = np.array(weight)
    weight = weight[:, np.newaxis]

    # OPTIMIZATION
    # Perform element-wise multiplication
    weighted_counts = adata.X.multiply(weight)

    # convert back to sparse
    adata.X = sparse.csr_matrix(weighted_counts)
    return adata


def weight_by_area_assignment(result_spatial_join, expanded_adata, cell_gdf_chunk):
    # Calculate overlapping area between cell and bin
    result_spatial_join["area"] = result_spatial_join.apply(
        lambda row: row["geometry"]
        .intersection(cell_gdf_chunk.loc[row["index_right"], "geometry"])
        .area,
        axis=1,
    )
    bin_area = result_spatial_join.iloc[0]["geometry"].area
    result_spatial_join["weight"] = result_spatial_join["area"] / bin_area
    result_spatial_join.loc[
        result_spatial_join["unique_bin"],
        "weight",
    ] = 1
    expanded_adata.obs["weight"] = result_spatial_join["weight"].tolist()
    expanded_adata = apply_weights_to_adata_counts(expanded_adata)
    return result_spatial_join, expanded_adata
