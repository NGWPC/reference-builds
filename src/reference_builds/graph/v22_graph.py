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


def _find_outlets_by_hydroseq(reference_flowpaths: pl.DataFrame) -> list[str]:
    """Find outlets for the river using hydroseq.

    Parameters
    ----------
    reference_flowpaths : pl.DataFrame
        The flowpath reference

    Returns
    -------
    list[str]
        All outlets from the reference
    """
    df_pl = reference_flowpaths.select(pl.col(["flowpath_id", "hydroseq", "dnhydroseq", "totdasqkm"]))

    df_with_str_id = df_pl.with_columns(
        pl.col("flowpath_id").cast(pl.Float64).cast(pl.Int64).cast(pl.Utf8).alias("flowpath_id_str")
    )

    hydroseq_set: set[Any] = set(df_pl["hydroseq"].to_list())

    outlets_df = df_with_str_id.filter(
        (pl.col("dnhydroseq") == 0) | ~pl.col("dnhydroseq").is_in(hydroseq_set)
    ).sort("flowpath_id_str")  # dnhydroseq is 0, or doesn't exist in hydroseq

    # outlets_sorted = outlets_df.sort("totdasqkm", descending=True) # Commenting out until production
    outlets: list[str] = outlets_df["flowpath_id_str"].to_list()

    return outlets


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
    query = "SELECT id,toid FROM flowpaths"
    conn = sqlite3.connect(file_path)
    fp = pl.read_database(query=query, connection=conn)
    fp = fp.extend(pl.DataFrame({"id": ["wb-0"], "toid": [None]})).lazy()
    query = "SELECT id,toid FROM network"
    network = pl.read_database(query=query, connection=conn).lazy()
    network = network.filter(pl.col("id").str.starts_with("wb-").not_())
    return create_matrix(fp, network)
