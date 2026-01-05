"""Contains all code for processing nhd data"""

import logging

import polars as pl

logger = logging.getLogger(__name__)


def _build_geoglows_graph(flowpaths: pl.DataFrame) -> dict[str, list[str]]:
    """Build a graph of upstream flowpath connections from GeoGLOWS data.

    Parameters
    ----------
    flowpaths : pl.DataFrame
        The GeoGLOWS flowpaths with LINKNO and DSLINKNO columns

    Returns
    -------
    dict[str, list[str]]
        The upstream dictionary containing upstream and downstream connections
        Key is the downstream flowpath ID, values are the upstream flowpath IDs
    """
    # Filter out terminal links (DSLINKNO == -1) for building upstream connections
    connectivity = flowpaths.select(
        [
            pl.col("LINKNO").cast(pl.Int64),
            pl.col("DSLINKNO").cast(pl.Int64),
        ]
    ).filter(pl.col("DSLINKNO") != -1)

    # Build upstream network: group by downstream link to get all upstream links
    upstream_network_df = connectivity.group_by(pl.col("DSLINKNO").cast(pl.Utf8).alias("downstream_id")).agg(
        pl.col("LINKNO").cast(pl.Utf8).alias("upstream_list")
    )

    upstream_dict: dict[str, list[str]] = dict(
        zip(
            upstream_network_df["downstream_id"].to_list(),
            upstream_network_df["upstream_list"].to_list(),
            strict=False,
        )
    )

    # Ensure all flowpath IDs are in the dictionary (even those with no upstream)
    all_flowpath_ids = flowpaths.select(pl.col("LINKNO").cast(pl.Utf8))["LINKNO"].to_list()

    for fp_id in all_flowpath_ids:
        if fp_id not in upstream_dict:
            upstream_dict[fp_id] = []

    return upstream_dict
