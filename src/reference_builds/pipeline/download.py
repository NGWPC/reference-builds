"""Contains all code for downloading hydrofabric data"""

import logging
from pathlib import Path
from typing import Any, cast

import geopandas as gpd
import pandas as pd
import polars as pl

from reference_builds.configs import ReferenceConfig
from reference_builds.utils import _validate_and_fix_geometries

logger = logging.getLogger(__name__)


def _load_and_concat_layers(gpkg_files: list[Path], layer_name: str | None) -> gpd.GeoDataFrame:
    """Load a specific layer from all gpkg files and concatenate."""
    gdfs = []
    for gpkg_path in gpkg_files:
        if layer_name is None:
            gdf = gpd.read_file(gpkg_path, driver="GPKG")
        else:
            gdf = gpd.read_file(gpkg_path, layer=layer_name)
        gdfs.append(gdf)
    return pd.concat(gdfs, ignore_index=True)


def _load_and_concat_parquet(parquet_files: list[Path]) -> gpd.GeoDataFrame:
    """Load a specific layer from all parquet files and concatenate."""
    gdfs = []
    for parquet_path in parquet_files:
        gdf = gpd.read_parquet(parquet_path)
        gdfs.append(gdf)
    return pd.concat(gdfs, ignore_index=True)


def _merge_flowpaths_without_catchments(
    flowpaths: gpd.GeoDataFrame,
    catchments: gpd.GeoDataFrame,
    connectivity: pd.DataFrame,
) -> tuple[gpd.GeoDataFrame, pd.DataFrame]:
    """Merge flowpaths without catchments to its downstream neighbor

    Parameters
    ----------
    flowpaths : gpd.GeoDataFrame
        The flowpath geodataframe, must contain 'DnHydroSeq' and 'HydroSeq' columns
        for upstream/downstream connectivity and 'NHDPlusID' to match with its catchment.
    catchments : gpd.GeoDataFrame
        The catchment geodataframe, must contain 'NHDPlusID' column

    Returns
    -------
    tuple[gpd.GeoDataFrame, pd.DataFrame]
        The updated flowpath geodataframe and connectivity dataframe after merging flowpaths without catchments
    """
    _flowpaths = flowpaths.copy()
    _rename_mapping = {key: f"{key}_VAA" for key in connectivity.columns if key != "NHDPlusID"}
    _connectivity = connectivity.copy().rename(columns=_rename_mapping)

    # merge flowpaths with connectivity info
    n_connect = _flowpaths["NHDPlusID"].isin(_connectivity["NHDPlusID"]).sum()
    n_flowpaths = len(_flowpaths.index)
    if n_connect != n_flowpaths:
        n_missing = n_flowpaths - n_connect
        logger.warning(
            f"{n_missing}/{n_flowpaths} flowpaths are missing connectivity info. These flowpaths will be dropped from the reference build."
        )

    _flowpaths = _flowpaths.merge(_connectivity, on="NHDPlusID", how="inner")

    # identify flowpaths without catchments
    _has_catchment = _flowpaths["NHDPlusID"].isin(catchments["NHDPlusID"])
    _flowpaths_without_catchments = _flowpaths[~_has_catchment]

    logger.info(f"Merging {len(_flowpaths_without_catchments)} flowpaths without catchments")

    for _, row in _flowpaths_without_catchments.iterrows():
        hydroseq = row["HydroSeq_VAA"]
        dnhydroseq = row["DnHydroSeq_VAA"]
        downstream_flowpath = _flowpaths[_flowpaths["HydroSeq_VAA"] == dnhydroseq]

        if len(downstream_flowpath) == 0 or row["TerminalFl_VAA"] == 1:
            continue
        elif len(downstream_flowpath) > 1:
            logger.warning(f"Multiple downstream flowpaths found for {row['NHDPlusID']}")
            continue
        else:
            # merge geometry
            _gdf = gpd.GeoDataFrame([row, downstream_flowpath.iloc[0]], geometry="geometry")
            merged_geom = _gdf.geometry.union_all()
            # update geometry of downstream flowpath
            dn_idx = downstream_flowpath.index[0]
            _flowpaths.at[dn_idx, "geometry"] = merged_geom
            # update connectivity
            upstream_flowpaths = _flowpaths[_flowpaths["DnHydroSeq_VAA"] == hydroseq]
            if len(upstream_flowpaths) >= 1:
                _flowpaths.loc[upstream_flowpaths.index, "DnHydroSeq_VAA"] = dnhydroseq
                _flowpaths.loc[dn_idx, "FromNode_VAA"] = upstream_flowpaths.iloc[0]["ToNode_VAA"]
            # TODO: check if other attributes need to be updated??

    # drop flowpaths without catchments after merging
    _flowpaths = _flowpaths[_has_catchment]

    # derive updated connectivity from merged flowpaths
    connectivity = _flowpaths[list(_rename_mapping.values())].rename(
        columns={v: k for k, v in _rename_mapping.items()}
    )
    connectivity["NHDPlusID"] = _flowpaths["NHDPlusID"]
    flowpaths = _flowpaths.drop(columns=list(_rename_mapping.values()))

    return flowpaths, connectivity


def download_geoglows_data(**context: dict[str, Any]) -> dict[str, pl.DataFrame]:
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
    dict[str, pl.DataFrame]
        The reference flowpath and divides references in memory
    """
    cfg = cast(ReferenceConfig, context["config"])

    # find the gpkg files from the ScienceBase NHD folders
    gpkg_files = list(cfg.output_dir.glob(cfg.input_file_regex))

    assert cfg.geoglows_catchment_regex is not None, "Need to specify where the catchment parquet files are"
    parquet_files = list(cfg.output_dir.glob(cfg.geoglows_catchment_regex))
    # load layers
    __flowpaths = _load_and_concat_layers(gpkg_files, layer_name=None)
    __catchments = _load_and_concat_parquet(parquet_files)
    # filter/validate layers
    flowpaths = _validate_and_fix_geometries(__flowpaths, geom_type="flowpaths")
    catchments = _validate_and_fix_geometries(__catchments, geom_type="divides")

    return {
        "geoglows_flowpaths": pl.from_pandas(flowpaths.to_wkb()),
        "geoglows_divides": pl.from_pandas(catchments.to_wkb()),
    }


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
    dict[str, pl.DataFrame]
        The reference flowpath and divides references in memory
    """
    cfg = cast(ReferenceConfig, context["config"])

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
    _flowpaths = _flowpaths[_flowpaths["fcode_description"].isin(cfg.permitted_fcodes)]

    catchments = _validate_and_fix_geometries(data["NHDPlusCatchment"], geom_type="divides")

    return {
        "nhd_flowpaths": pl.from_pandas(_flowpaths.to_wkb()),
        "nhd_divides": pl.from_pandas(catchments.to_wkb()),
        "nhd_connectivity": pl.from_pandas(data["NHDPlusFlowlineVAA"]),
    }
