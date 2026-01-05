"""Contains all pipeline code for processing reference data"""

import logging
from typing import Any, cast

import polars as pl
import rustworkx as rx

from reference_builds.task_instance import TaskInstance
from reference_builds.utils.geoglows_graph import _build_geoglows_graph
from reference_builds.utils.nhd_graph import _build_graph

logger = logging.getLogger(__name__)


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
    if None in upstream_network:
        upstream_network.pop(None)  # type: ignore
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


def build_nhd_graphs(**context: dict[str, Any]) -> dict[str, Any]:
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


def build_geoglows_graphs(**context: dict[str, Any]) -> dict[str, Any]:
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
    flowpaths: pl.DataFrame = ti.xcom_pull(task_id="download", key="geoglows_flowpaths")
    upstream_network = _build_geoglows_graph(flowpaths)
    graph, node_indices = _build_rustworkx_object(upstream_network)
    return {"graph": graph, "node_indices": node_indices}
