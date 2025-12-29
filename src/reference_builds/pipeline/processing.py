"""Contains all code for downloading hydrofabric data"""

import logging
from typing import Any, cast

import polars as pl
import rustworkx as rx

from reference_builds.task_instance import TaskInstance

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

    return upstream_dict


def _build_rustworkx_object(
    upstream_network: dict[str, list[str]] | dict[int, list[int]],
) -> tuple[rx.PyDiGraph, dict[str, int] | dict[int, int]]:
    """Build a RustWorkX directed graph from upstream network dictionary.

    Parameters
    ----------
    upstream_network : dict[str, list[str]] | dict[int, list[int]]
        Dictionary mapping downstream flowpath IDs to lists of upstream flowpath IDs

    Returns
    -------
    tuple[rx.PyDiGraph, dict[str, int] | dict[int, int]]
        The flowpaths object in graph form and node indices for each object in the graph
    """
    graph = rx.PyDiGraph(check_cycle=True)
    node_indices: dict[Any, int] = {}
    for to_edge in sorted(upstream_network.keys()):
        from_edges = upstream_network[to_edge]  # type: ignore
        if to_edge not in node_indices:
            node_indices[to_edge] = graph.add_node(to_edge)
        for from_edge in from_edges:
            if from_edge not in node_indices:
                node_indices[from_edge] = graph.add_node(from_edge)
    for to_edge, from_edges in upstream_network.items():
        for from_edge in from_edges:
            graph.add_edge(node_indices[from_edge], node_indices[to_edge], None)
    return graph, node_indices


def build_graphs(**context: dict[str, Any]) -> dict[str, Any]:
    """Builds and processes graphs from NHD data

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
    dict[str, Any]
        The rustworkx graph object and node_indices for the NHD
    """
    ti = cast(TaskInstance, context["ti"])
    flowpaths: pl.DataFrame = ti.xcom_pull(task_id="download", key="nhd_flowpaths")
    connectivity: pl.DataFrame = ti.xcom_pull(task_id="download", key="nhd_connectivity")
    upstream_network = _build_graph(connectivity, flowpaths)
    graph, node_indices = _build_rustworkx_object(upstream_network)
    return {"graph": graph, "node_indices": node_indices}
