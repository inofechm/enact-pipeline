# Script runs the evaluation to compare ENACT cell annotations versus pathologist cell annotations

from shapely.geometry import shape
import plotly.express as px
import geopandas as gpd
import json
from shapely.geometry import Polygon, Point
from shapely import wkt
import pandas as pd
from sklearn.metrics import precision_recall_fscore_support, accuracy_score
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
import os

# from src.pipelines.enact_pipeline import ENACT

# so_hd = ENACT(configs_path="config/configs.yaml")

geojson_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/Visium_HD_Human_Colon_Cancer-wsi-40598_0_65263_22706.geojson"
segmentation_df_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/colon/cells_df.csv"
predictions_df_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/colon/chunks/weighted_by_cluster/cellassign_results/merged_results.csv"
# predictions_df_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/colon/chunks/weighted_by_cluster/sargent_results/merged_results.csv"
# predictions_df_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/colon/chunks/weighted_by_area/cellassign_results/merged_results.csv"
# predictions_df_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/colon/chunks/weighted_by_area/sargent_results/merged_results.csv"
# predictions_df_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/colon/chunks/weighted_by_gene/sargent_results/merged_results.csv"
# predictions_df_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/colon/chunks/naive/sargent_results/merged_results.csv"
# predictions_df_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/colon/chunks/weighted_by_gene/cellassign_results/merged_results.csv"
# predictions_df_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/colon/chunks/naive/cellassign_results/merged_results.csv"
# predictions_df_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/colon/chunks/naive/celltypist_results/merged_results.csv"
# predictions_df_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/colon/chunks/weighted_by_area/celltypist_results/merged_results.csv"
# predictions_df_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/colon/chunks/weighted_by_gene/celltypist_results/merged_results.csv"
# predictions_df_path = "/home/oneai/oneai-dda-spatialtr-visiumhd_analysis/cache/colon/chunks/weighted_by_cluster/celltypist_results/merged_results.csv"


results_eval_dir = os.path.join("/".join(predictions_df_path.split("/")[:-1]), "eval")
os.makedirs(results_eval_dir, exist_ok=True)


name_map = {
    'unclassified': "no label",
    'Immune': "immune cells",
    'Crypt cells': "epithelial cells",
    'Enterocytes': "epithelial cells",
    'Epithelial': "epithelial cells",
    'Smooth muscle cell': "stromal cells",
    'Fibroblast': "stromal cells",
    'Endothelial': "stromal cells",
    'Paneth cells': "epithelial cells",
    'Enteroendocrine cells': "epithelial cells",
    'Goblet cells': "epithelial cells",
    'Neuronal': "stromal cells",
    'ephitelial cells': "epithelial cells",
    'no label': "no label",
    "Ignore*": "no label",
    "B cells": "immune cells",
    "T cells": "immune cells",
    "NK cells": "immune cells",
    "Macrophages": "immune cells",
    "Neutrophils": "immune cells",
    "Eosinophils": "immune cells",
    'CD19+CD20+ B': "immune cells",               # B cells are immune cells
    'CD4+ T cells': "immune cells",               # CD4+ T cells are immune cells
    'CD8+ T cells': "immune cells",               # CD8+ T cells are immune cells
    'CMS1': "epithelial cells",                   # CMS (Consensus Molecular Subtypes) refer to tumor/epithelial cells
    'CMS2': "epithelial cells",                   # Same as above
    'CMS3': "epithelial cells",                   # Same as above
    'CMS4': "epithelial cells",                   # Same as above
    'Enteric glial cells': "stromal cells",       # Glial cells are part of the stromal tissue
    'Goblet cells': "epithelial cells",           # Goblet cells are epithelial cells
    'IgA+ Plasma': "immune cells",                # Plasma cells are immune cells (B-cell derivatives)
    'IgG+ Plasma': "immune cells",                # Same as above
    'Intermediate': "no label",                   # Ambiguous, no clear label
    'Lymphatic ECs': "stromal cells",             # Endothelial cells are considered stromal
    'Mast cells': "immune cells",                 # Mast cells are immune cells
    'Mature Enterocytes type 1': "epithelial cells", # Enterocytes are epithelial cells
    'Mature Enterocytes type 2': "epithelial cells", # Same as above
    'Myofibroblasts': "stromal cells",            # Fibroblasts are stromal cells
    'NK cells': "immune cells",                   # NK cells are immune cells
    'Pericytes': "stromal cells",                 # Pericytes are part of the vasculature (stromal)
    'Pro-inflammatory': "immune cells",           # Inflammation implies immune function
    'Proliferating': "no label",                  # Too vague to classify, no label
    'Proliferative ECs': "stromal cells",         # Endothelial cells are stromal
    'Regulatory T cells': "immune cells",         # T cells are immune cells
    'SPP1+': "no label",                          # Ambiguous, no clear label
    'Smooth muscle cells': "stromal cells",       # Smooth muscle cells are stromal cells
    'Stalk-like ECs': "stromal cells",            # Endothelial cells are stromal
    'Stem-like/TA': "epithelial cells",           # Stem cells in this context are usually epithelial
    'Stromal 1': "stromal cells",                 # Explicitly stromal
    'Stromal 2': "stromal cells",                 # Same as above
    'Stromal 3': "stromal cells",                 # Same as above
    'T follicular helper cells': "immune cells",  # T cells are immune cells
    'T helper 17 cells': "immune cells",          # Same as above
    'Tip-like ECs': "stromal cells",              # Endothelial cells are stromal
    'Unknown': "no label",                        # No clear label
    'cDC': "immune cells",                        # Conventional dendritic cells are immune cells
    'gamma delta T cells': "immune cells"         # T cells are immune cells
}


segmentation_df = pd.read_csv(segmentation_df_path)
predictions_df = pd.read_csv(predictions_df_path)
predictions_df = predictions_df.merge(segmentation_df[["id", "geometry"]], how="left", on="id")
predictions_df["geometry"] = predictions_df["geometry"].apply(wkt.loads)
pred_gpd = gpd.GeoDataFrame(predictions_df,geometry="geometry")

def load_path_annotations():
    annotation_names = []
    annotation_geometries = []
    with open(geojson_path) as f:
        regions = json.load(f)
    for region in regions["features"]:
        ann_type = region["properties"]["objectType"]
        if ann_type == "annotation":
            annotation_name = region["properties"]["classification"]["name"]
            if annotation_name in ["Region*"]:
                continue
            annotation_geometries.append(shape(region["geometry"]))
            annotation_names.append(annotation_name)
    annotations_gpd = gpd.GeoDataFrame({"geometry": annotation_geometries, "gt_label": annotation_names})
    annotations_gpd["ann_ix"] = [f"ID_{i}" for i in range(len(annotations_gpd))]
    return annotations_gpd

def get_gt_annotations(annotations_gpd):
    try:
        cells_within_ann_gpd = gpd.sjoin(annotations_gpd, pred_gpd[["cell_type", "cell_x", "cell_y", "geometry", "id"]], how='left', predicate='intersects')
    except:
        cells_within_ann_gpd = gpd.sjoin(annotations_gpd, pred_gpd[["cell_assign_results", "cell_x", "cell_y", "geometry", "id"]], how='left', predicate='intersects')
    cells_within_ann_gpd = cells_within_ann_gpd.drop_duplicates("ann_ix")
    try:
        cells_within_ann_gpd["cell_type"] = cells_within_ann_gpd["cell_type"].fillna("unclassified")
    except:
        cells_within_ann_gpd["cell_assign_results"] = cells_within_ann_gpd["cell_assign_results"].fillna("unclassified")
    return cells_within_ann_gpd

def validate_labels(cells_within_ann_gpd):
    try:
        cell_types_in_pred = set(cells_within_ann_gpd.cell_type.unique())
    except:
        cell_types_in_pred = set(cells_within_ann_gpd.cell_assign_results.unique())
    print(f"Cells in pred dataset: {cell_types_in_pred}")
    print (f"All cells are in the mapping!: {cell_types_in_pred.issubset(set(name_map.keys()))}")
    
def relabel_cells(cells_within_ann_gpd):
    # Renaming cell types
    for granular_name, generic_name in name_map.items():
        cells_within_ann_gpd.loc[cells_within_ann_gpd.gt_label == granular_name, "gt_label"] = generic_name
        try:
            cells_within_ann_gpd.loc[cells_within_ann_gpd.cell_type == granular_name, "pred_label_clean"] = generic_name
        except:
            cells_within_ann_gpd.loc[cells_within_ann_gpd.cell_assign_results == granular_name, "pred_label_clean"] = generic_name
    return cells_within_ann_gpd

def eval_annotations(results_table):
    cell_types = sorted(set(results_table.gt_label.unique().tolist() + results_table.pred_label_clean.unique().tolist()))
    cm = confusion_matrix(
        results_table.gt_label, 
        results_table.pred_label_clean, 
        labels=cell_types
        )
    cm_plot = ConfusionMatrixDisplay(
        confusion_matrix=cm,
        display_labels=cell_types
        )
    cm_plot.plot()

    averaging_methods = ["micro", "macro", "weighted"]
    eval_dict = {}
    for method in averaging_methods:
        eval_metrics = precision_recall_fscore_support(results_table.gt_label, results_table.pred_label_clean, average=method)
        precision, recall, fbeta_score, support = eval_metrics
        eval_dict[method] = eval_metrics
    num_correct_samples = accuracy_score(results_table.gt_label, results_table.pred_label_clean, normalize=False)
    accuracy = accuracy_score(results_table.gt_label, results_table.pred_label_clean, normalize=True)
    print(f"Experiment name: {predictions_df_path}")
    print (f"Number of GT annotations: {len(results_table)}\nNumber of correct predictions: {num_correct_samples}\nAccuracy: {accuracy}")
    print("__________")
    try:
        print(pd.DataFrame(results_table.cell_type.value_counts()))
    except:
        print(pd.DataFrame(results_table.cell_assign_results.value_counts()))
    print("__________")
    print(pd.DataFrame(results_table.pred_label_clean.value_counts()))
    print("__________")
    metrics_df = pd.DataFrame(eval_dict, index=["Precision", "Recall", "F-Score", "Support"])
    results_table.to_csv(os.path.join(results_eval_dir, "cell_annotation_eval.csv"), index=False)
    metrics_df.to_csv(os.path.join(results_eval_dir, "cell_annotation_eval_metrics.csv"), index=True)
    cm_plot.figure_.savefig(os.path.join(results_eval_dir, "confusion_matrix.png"),dpi=300)
    print (metrics_df)
    return results_table, metrics_df

if __name__ == "__main__":
    annotations_gpd =  load_path_annotations()
    cells_within_ann_gpd = get_gt_annotations(annotations_gpd)
    validate_labels(cells_within_ann_gpd)
    cells_within_ann_gpd = relabel_cells(cells_within_ann_gpd)
    results_table = cells_within_ann_gpd[(cells_within_ann_gpd["gt_label"] != "no label")]
    results_table, metrics_df = eval_annotations(results_table)