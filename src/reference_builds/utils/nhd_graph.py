"""Contains all code for processing nhd data"""

import logging

import polars as pl

logger = logging.getLogger(__name__)


def _build_graph(connectivity: pl.DataFrame, flowpaths: pl.DataFrame) -> dict[str, list[str]]:
    """Build a graph of upstream flowpath connections.

    Parameters
    ----------
    connectivity : pl.DataFrame
        The connectivity/flow table with FromNode and ToNode columns
    flowpaths : pl.DataFrame
        The reference flowpaths to filter to

    Returns
    -------
    dict[str, list[str]]
        The upstream dictionary containing upstream and downstream connections
        Key is the downstream flowpath ID, values are the upstream flowpath IDs
    """
    valid_ids = flowpaths.select(pl.col("NHDPlusID").cast(pl.Int64))["NHDPlusID"]

    filtered_connectivity = connectivity.select(
        [
            pl.col("NHDPlusID").cast(pl.Int64),
            pl.col("FromNode").cast(pl.Int64),
            pl.col("ToNode").cast(pl.Int64),
        ]
    ).filter(pl.col("NHDPlusID").is_in(valid_ids))

    tonode_lookup = filtered_connectivity.select(
        [
            pl.col("ToNode"),
            pl.col("NHDPlusID").cast(pl.Utf8).alias("upstream_id"),
        ]
    )

    fromnode_lookup = filtered_connectivity.select(
        [
            pl.col("FromNode"),
            pl.col("NHDPlusID").cast(pl.Utf8).alias("downstream_id"),
        ]
    )

    merged = tonode_lookup.join(fromnode_lookup, left_on="ToNode", right_on="FromNode", how="inner").select(
        ["upstream_id", "downstream_id"]
    )

    upstream_network_df = merged.group_by("downstream_id").agg(pl.col("upstream_id").alias("upstream_list"))

    upstream_dict: dict[str, list[str]] = dict(
        zip(
            upstream_network_df["downstream_id"].to_list(),
            upstream_network_df["upstream_list"].to_list(),
            strict=False,
        )
    )

    all_flowpath_ids = flowpaths.select(pl.col("NHDPlusID").cast(pl.Int64).cast(pl.Utf8))[
        "NHDPlusID"
    ].to_list()

    for fp_id in all_flowpath_ids:
        if fp_id not in upstream_dict:
            upstream_dict[fp_id] = []

    return upstream_dict
