"""Contains all code for downloading hydrofabric data"""

import logging
from pathlib import Path
from typing import Any, cast

import geopandas as gpd
import pandas as pd
import polars as pl

from reference_builds.configs import PRVI
from reference_builds.graph import _validate_and_fix_geometries

logger = logging.getLogger(__name__)


def _load_and_concat_layers(gpkg_files: list[Path], layer_name: str) -> gpd.GeoDataFrame:
    """Load a specific layer from all gpkg files and concatenate."""
    gdfs = []
    for gpkg_path in gpkg_files:
        gdf = gpd.read_file(gpkg_path, layer=layer_name)
        gdfs.append(gdf)
    return pd.concat(gdfs, ignore_index=True)


def download_nhd_data(**context: dict[str, Any]) -> dict[str, pl.DataFrame]:
    """Opens local / downloads for the reference-build process

    Parameters
    ----------
    **context : dict
        Airflow-compatible context containing:
        - ti : TaskInstance for XCom operations
        - config : HFConfig with pipeline configuration
        - task_id : str identifier for this task
        - run_id : str identifier for this pipeline run
        - ds : str execution date
        - execution_date : datetime object

    Returns
    -------
    dict[str, gpd.GeoDataFrame]
        The reference flowpath and divides references in memory
    """
    cfg = cast(PRVI, context["config"])

    # find the gpkg files from the ScienceBase NHD folders
    matching_folders = list(cfg.output_dir.glob(cfg.input_file_regex))
    gpkg_files: list[Path] = []
    for folder in matching_folders:
        if folder.is_dir():
            gpkg_files.extend(folder.glob("*.gpkg"))

    # load layers
    layers = [
        "NHDFlowline",
        "NHDPlusCatchment",
        "NHDPlusFlowlineVAA",
    ]
    data = {layer: _load_and_concat_layers(gpkg_files, layer) for layer in layers}

    # filter/validate layers
    _flowpaths = _validate_and_fix_geometries(data["NHDFlowline"], geom_type="flowpaths")
    catchments = _validate_and_fix_geometries(data["NHDPlusCatchment"], geom_type="divides")
    _flowpaths_with_catchments = _flowpaths[_flowpaths["NHDPlusID"].isin(catchments["NHDPlusID"])]
    flowpaths = _flowpaths_with_catchments[
        _flowpaths_with_catchments["fcode_description"].isin(cfg.permitted_fcodes)
    ]

    return {
        "nhd_flowpaths": pl.from_pandas(flowpaths.to_wkb()),
        "nhd_divides": pl.from_pandas(catchments.to_wkb()),
        "nhd_connectivity": pl.from_pandas(data["NHDPlusFlowlineVAA"]),
    }
