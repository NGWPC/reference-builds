"""
Builds a graph object from the v2.2 Hydrofabric

Kudos to Nels Fraizer for assistance in making the graph building code:
https://github.com/DeepGroundwater/ddr/blob/ab4c3962c2c119e6a9182a77f2a9faceec19f2e0/engine/adjacency.py
"""

import sqlite3
from pathlib import Path

import polars as pl
import rustworkx as rx
from tqdm import tqdm


def create_matrix(fp: pl.LazyFrame, network: pl.LazyFrame) -> tuple[rx.PyDiGraph, dict[str, int]]:
    """
    Create a directed graph from flowpaths and network dataframes.

    Parameters
    ----------
    fp : pl.LazyFrame
        Flowpaths dataframe with 'toid' column indicating downstream nexus IDs.
    network : pl.LazyFrame
        Network dataframe with 'toid' column indicating downstream flowpath IDs.

    Returns
    -------
    tuple[rx.PyDiGraph, dict[str, int]]
        tuple[0]: A rustworkx directed graph
        tuple[1]: Mapping of flowpath IDs to graph node indices
    """
    fp = fp.with_row_index(name="idx").collect()
    network = network.collect().unique(subset=["id"])
    _values = zip(fp["idx"], fp["toid"], strict=False)
    fp = dict(zip(fp["id"], _values, strict=True))
    network = dict(zip(network["id"], network["toid"], strict=True))

    graph = rx.PyDiGraph(check_cycle=False, node_count_hint=len(fp), edge_count_hint=len(fp))
    gidx = graph.add_nodes_from(fp.keys())

    # Create mapping from flowpath ID to graph node index
    node_index = {graph.get_node_data(idx): idx for idx in gidx}

    for idx in tqdm(gidx, desc="Building network graph"):
        id = graph.get_node_data(idx)
        nex = fp[id][1]  # the downstream nexus id
        ds_wb = network.get(nex)
        if ds_wb is not None:
            graph.add_edge(idx, node_index[ds_wb], nex)

    return graph, node_index


def build_v22_graph(file_path: Path) -> tuple[rx.PyDiGraph, dict[str, int]]:
    """Builds a graph from the v2.2 hydrofabric

    Parameters
    ----------
    file_path : Path
        The path to the v2.2 geopackage

    Returns
    -------
    tuple[rx.PyDiGraph, dict[str, int]]
        tuple[0]: A rustworkx directed graph
        tuple[1]: Mapping of flowpath IDs to graph node indices
    """
    # Read hydrofabric geopackage using sqlite
    uri = "sqlite://" + str(file_path)
    query = "SELECT id,toid FROM flowpaths"
    conn = sqlite3.connect(uri)
    fp = pl.read_database(query=query, connection=conn)
    fp = fp.extend(pl.DataFrame({"id": ["wb-0"], "toid": [None]})).lazy()
    query = "SELECT id,toid FROM network"
    network = pl.read_database(query=query, connection=conn).lazy()
    network = network.filter(pl.col("id").str.starts_with("wb-").not_())
    return create_matrix(fp, network)
