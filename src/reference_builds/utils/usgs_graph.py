"""Contains all code for processing USGS Reference Hydrofabric data"""

import logging

import polars as pl

logger = logging.getLogger(__name__)


def _build_usgs_hf_graph(flowpaths: pl.DataFrame) -> dict[str, list[str]]:
    """Build a graph of upstream flowpath connections from USGS Reference Hydrofabric data.

    Parameters
    ----------
    flowpaths : pl.DataFrame
        The USGS Reference Hydrofabric flowpaths with hydroseq and dnhydroseq columns

    Returns
    -------
    dict[str, list[str]]
        The upstream dictionary containing upstream and downstream connections
        Key is the downstream flowpath ID, values are the upstream flowpath IDs
    """
    # Filter out terminal links (dnhydroseq == 0) for building upstream connections
    connectivity = flowpaths.select(
        [
            pl.col("hydroseq").cast(pl.Int64),
            pl.col("dnhydroseq").cast(pl.Int64),
        ]
    ).filter(pl.col("dnhydroseq") != 0)

    # Build upstream network: group by downstream link to get all upstream links
    upstream_network_df = connectivity.group_by(
        pl.col("dnhydroseq").cast(pl.Utf8).alias("downstream_id")
    ).agg(pl.col("hydroseq").cast(pl.Utf8).alias("upstream_list"))

    upstream_dict: dict[str, list[str]] = dict(
        zip(
            upstream_network_df["downstream_id"].to_list(),
            upstream_network_df["upstream_list"].to_list(),
            strict=False,
        )
    )

    # Ensure all flowpath IDs are in the dictionary (even those with no upstream)
    all_flowpath_ids = flowpaths.select(pl.col("hydroseq").cast(pl.Utf8))["hydroseq"].to_list()

    for fp_id in all_flowpath_ids:
        if fp_id not in upstream_dict:
            upstream_dict[fp_id] = []

    return upstream_dict
