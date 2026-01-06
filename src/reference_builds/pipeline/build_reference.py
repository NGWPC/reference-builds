"""Contains all code for building a reference fabric from the reference graph object"""

import logging
from typing import Any, cast

import geopandas as gpd
import pandas as pd
import polars as pl
import rustworkx as rx

from reference_builds.configs import ReferenceConfig
from reference_builds.task_instance import TaskInstance
from reference_builds.utils.geometries import _orient_flowpath_downstream

logger = logging.getLogger(__name__)


def _trace_attributes(
    graph: rx.PyDiGraph,
    node_indices: dict[str, int],
    flowpaths: gpd.GeoDataFrame,
    divides: gpd.GeoDataFrame,
    vpu_id: str,
) -> gpd.GeoDataFrame:
    """Trace flowpath attributes for the entire graph.

    Parameters
    ----------
    graph : rx.PyDiGraph
        The rustworkx directed graph (may contain multiple disconnected subgraphs)
    node_indices : dict[str, int]
        Mapping from NHDPlusID (as string) to node index
    flowpaths : gpd.GeoDataFrame
        The flowpaths GeoDataFrame with LengthKM
    divides : gpd.GeoDataFrame
        The divides GeoDataFrame with AreaSqKm

    Returns
    -------
    pl.DataFrame
        Traced attributes: totdasqkm, mainstemlp, pathlength, dnhydroseq, hydroseq, stream_order
    """
    flowpaths_lookup = flowpaths.set_index("NHDPlusID")["LengthKM"].to_dict()
    divides_lookup = divides.set_index("NHDPlusID")["AreaSqKm"].to_dict()
    fp_geom_lookup = flowpaths.set_index("NHDPlusID")["geometry"].to_dict()
    fp_fcode_lookup = flowpaths.set_index("NHDPlusID")["fcode_description"].to_dict()
    for node_idx in graph.node_indices():
        flowpath_id = str(graph[node_idx])
        nhd_id = int(flowpath_id)

        graph[node_idx] = {
            "flowpath_id": flowpath_id,
            "areasqkm": divides_lookup.get(nhd_id, 0.0),
            "lengthkm": flowpaths_lookup.get(nhd_id, 0.0),
            "totdasqkm": 0.0,
            "mainstemlp": None,
            "pathlength": 0.0,
            "dnhydroseq": None,
            "hydroseq": None,
            "streamorder": None,
            "fcode_description": fp_fcode_lookup[nhd_id],
            "geometry": fp_geom_lookup[nhd_id],
        }

    # Find all outlets (nodes with no downstream connections)
    outlets = [idx for idx in graph.node_indices() if graph.out_degree(idx) == 0]
    logger.info(f"build_nhd_reference task: Found {len(outlets)} outlets (disconnected subgraphs)")

    # Get topological order for entire graph
    try:
        topo_order = rx.topological_sort(graph)
    except rx.DAGHasCycle as e:
        raise AssertionError("Graph contains cycles") from e

    # PASS 1: Calculate pathlength and hydroseq (reverse topo order - upstream from outlets)
    current_hydroseq = 1

    # Initialize outlets
    for outlet_idx in outlets:
        graph[outlet_idx]["pathlength"] = 0.0
        graph[outlet_idx]["dnhydroseq"] = 0

    # Traverse in reverse topo order
    for node_idx in reversed(topo_order):
        # Assign hydroseq
        graph[node_idx]["hydroseq"] = current_hydroseq
        current_hydroseq += 1

        # Calculate pathlength based on downstream node
        out_edges = graph.out_edges(node_idx)
        if out_edges:
            downstream_nodes = [tgt_idx for _, tgt_idx, _ in out_edges]
            if downstream_nodes:
                downstream_idx = max(downstream_nodes, key=lambda idx: graph[idx]["pathlength"])
                graph[node_idx]["pathlength"] = (
                    graph[downstream_idx]["pathlength"] + graph[downstream_idx]["lengthkm"]
                )

    # Trace mainstems for each outlet's basin
    current_mainstem_id = 1
    processed: set[int] = set()

    for outlet_idx in outlets:
        # Trace main mainstem (longest path from outlet to headwater)
        current_idx = outlet_idx
        mainstem_nodes = []

        while current_idx not in processed:
            mainstem_nodes.append(current_idx)
            graph[current_idx]["mainstemlp"] = current_mainstem_id
            processed.add(current_idx)

            in_edges = list(graph.in_edges(current_idx))
            if not in_edges:
                break

            upstream_candidates = [src_idx for src_idx, _, _ in in_edges if src_idx not in processed]
            if not upstream_candidates:
                break

            current_idx = max(
                upstream_candidates,
                key=lambda idx: (graph[idx]["pathlength"], graph[idx]["totdasqkm"]),
            )

        current_mainstem_id += 1

    # Assign tributary mainstems for remaining nodes
    for node_idx in graph.node_indices():
        if node_idx not in processed:
            tributary_id = current_mainstem_id
            current_mainstem_id += 1

            trib_current = node_idx
            while trib_current not in processed:
                graph[trib_current]["mainstemlp"] = tributary_id
                processed.add(trib_current)

                in_edges = list(graph.in_edges(trib_current))
                upstream_in_basin = [src_idx for src_idx, _, _ in in_edges if src_idx not in processed]

                if not upstream_in_basin:
                    break

                trib_current = max(
                    upstream_in_basin,
                    key=lambda idx: (graph[idx]["pathlength"], graph[idx]["totdasqkm"]),
                )

    # Assign dnhydroseq and flowpath_toid based on graph edges
    for node_idx in graph.node_indices():
        out_edges = graph.out_edges(node_idx)
        downstream_nodes = [tgt_idx for _, tgt_idx, _ in out_edges]

        if downstream_nodes:
            downstream_idx = downstream_nodes[0]
            graph[node_idx]["dnhydroseq"] = graph[downstream_idx]["hydroseq"]
            graph[node_idx]["flowpath_toid"] = graph[downstream_idx]["flowpath_id"]
        else:
            graph[node_idx]["dnhydroseq"] = 0
            graph[node_idx]["flowpath_toid"] = "0"

    # PASS 2: Calculate totdasqkm and stream_order (forward topo order - downstream from headwaters)
    for node_idx in topo_order:
        in_edges = list(graph.in_edges(node_idx))

        # Accumulate upstream drainage area
        upstream_total = sum(graph[src_idx]["totdasqkm"] for src_idx, _, _ in in_edges)
        graph[node_idx]["totdasqkm"] = upstream_total + graph[node_idx]["areasqkm"]

        # Calculate Strahler stream order
        if not in_edges:
            graph[node_idx]["streamorder"] = 1
        else:
            upstream_orders = [graph[src_idx]["streamorder"] for src_idx, _, _ in in_edges]
            max_order = max(upstream_orders)
            count_max = upstream_orders.count(max_order)

            if count_max >= 2:
                graph[node_idx]["streamorder"] = max_order + 1
            else:
                graph[node_idx]["streamorder"] = max_order

    # PASS 3: Orient all flowpath geometries so they flow upstream -> downstream
    for node_idx in graph.node_indices():
        in_edges = graph.in_edges(node_idx)
        upstream_nodes = [src_idx for src_idx, _, _ in in_edges]

        # Check if this is an outlet (flowpath_toid == "0" means no downstream)
        is_outlet = graph[node_idx]["flowpath_toid"] == "0"

        if is_outlet and upstream_nodes:
            # Outlet: use upstream geometry
            upstream_idx = upstream_nodes[0]
            us_geom = graph[upstream_idx]["geometry"]
            graph[node_idx]["geometry"] = _orient_flowpath_downstream(
                graph[node_idx]["geometry"], ds_geom=None, us_geom=us_geom
            )
        else:
            # Normal case: use downstream geometry
            out_edges = graph.out_edges(node_idx)
            downstream_nodes = [tgt_idx for _, tgt_idx, _ in out_edges]

            ds_geom = None
            if downstream_nodes:
                downstream_idx = downstream_nodes[0]
                ds_geom = graph[downstream_idx]["geometry"]

            graph[node_idx]["geometry"] = _orient_flowpath_downstream(
                graph[node_idx]["geometry"], ds_geom=ds_geom, us_geom=None
            )

    # Extract results
    flowpath_ids = []
    flowpath_toids = []
    vpu_ids = []
    das = []
    lengthkms = []
    total_das = []
    mainstems = []
    pathlengths = []
    dnhydroseqs = []
    hydroseqs = []
    streamorders = []
    fcodes = []
    geometries = []

    for node_idx in graph.node_indices():
        node_data = graph[node_idx]
        flowpath_ids.append(node_data["flowpath_id"])
        flowpath_toids.append(node_data["flowpath_toid"])
        vpu_ids.append(vpu_id)
        das.append(node_data["areasqkm"])
        lengthkms.append(node_data["lengthkm"])
        total_das.append(node_data["totdasqkm"])
        mainstems.append(node_data["mainstemlp"])
        pathlengths.append(node_data["pathlength"])
        dnhydroseqs.append(node_data["dnhydroseq"])
        hydroseqs.append(node_data["hydroseq"])
        streamorders.append(node_data["streamorder"])
        fcodes.append(node_data["fcode_description"])
        geometries.append(node_data["geometry"])

    return gpd.GeoDataFrame(
        {
            "flowpath_id": flowpath_ids,
            "flowpath_toid": flowpath_toids,
            "VPUID": vpu_ids,
            "lengthkm": lengthkms,
            "areasqkm": das,
            "totdasqkm": total_das,
            "mainstemlp": mainstems,
            "pathlength": pathlengths,
            "dnhydroseq": dnhydroseqs,
            "hydroseq": hydroseqs,
            "streamorder": streamorders,
            "fcode_description": fcodes,
        },
        geometry=geometries,
        crs="EPSG:4269",
    )


def _trace_geoglows_attributes(
    graph: rx.PyDiGraph,
    node_indices: dict[str, int],
    flowpaths: gpd.GeoDataFrame,
    catchments: gpd.GeoDataFrame,
    vpu_id: str,
) -> gpd.GeoDataFrame:
    """Trace flowpath attributes for the entire GeoGLOWS graph.

    Parameters
    ----------
    graph : rx.PyDiGraph
        The rustworkx directed graph (may contain multiple disconnected subgraphs)
    node_indices : dict[str, int]
        Mapping from LINKNO (as string) to node index
    flowpaths : gpd.GeoDataFrame
        The GeoGLOWS flowpaths GeoDataFrame with LengthKM, strmOrder
    catchments : gpd.GeoDataFrame
        The GeoGLOWS catchments GeoDataFrame with linkno, areasqkm, and geometry
    vpu_id : str
        The VPUID for the domain

    Returns
    -------
    gpd.GeoDataFrame
        Traced attributes: totdasqkm, mainstemlp, pathlength, dnhydroseq, hydroseq, stream_order
    """
    length_lookup = flowpaths.set_index("LINKNO")["LengthKM"].to_dict()
    order_lookup = flowpaths.set_index("LINKNO")["strmOrder"].to_dict()
    fp_geom_lookup = flowpaths.set_index("LINKNO")["geometry"].to_dict()

    # Build catchment lookups (area and geometry keyed by linkno)
    catchment_area_lookup = catchments.set_index("linkno")["areasqkm"].to_dict()
    catchment_geom_lookup = catchments.set_index("linkno")["geometry"].to_dict()

    for node_idx in graph.node_indices():
        flowpath_id = str(graph[node_idx])
        link_id = int(flowpath_id)

        # Get length in km (already calculated from geometry)
        length_km = length_lookup.get(link_id, 0.0)

        # Get local catchment area (from catchments, keyed by linkno)
        local_areasqkm = catchment_area_lookup.get(link_id, 0.0)

        # Get geometries
        catchment_geometry = catchment_geom_lookup.get(link_id)
        flowpath_geometry = fp_geom_lookup.get(link_id)

        graph[node_idx] = {
            "flowpath_id": flowpath_id,
            "areasqkm": local_areasqkm,
            "lengthkm": length_km,
            "totdasqkm": 0.0,  # Will be accumulated in PASS 2
            "mainstemlp": None,
            "pathlength": 0.0,
            "dnhydroseq": None,
            "hydroseq": None,
            "streamorder": order_lookup.get(link_id, 1),
            "flowpath_geometry": flowpath_geometry,
            "catchment_geometry": catchment_geometry,
        }

    # Find all outlets (nodes with no downstream connections)
    outlets = [idx for idx in graph.node_indices() if graph.out_degree(idx) == 0]
    logger.info(f"build_geoglows_reference task: Found {len(outlets)} outlets (disconnected subgraphs)")

    # Get topological order for entire graph
    try:
        topo_order = rx.topological_sort(graph)
    except rx.DAGHasCycle as e:
        raise AssertionError("Graph contains cycles") from e

    # PASS 1: Calculate pathlength and hydroseq (reverse topo order - upstream from outlets)
    current_hydroseq = 1

    # Initialize outlets
    for outlet_idx in outlets:
        graph[outlet_idx]["pathlength"] = 0.0
        graph[outlet_idx]["dnhydroseq"] = 0

    # Traverse in reverse topo order
    for node_idx in reversed(topo_order):
        # Assign hydroseq
        graph[node_idx]["hydroseq"] = current_hydroseq
        current_hydroseq += 1

        # Calculate pathlength based on downstream node
        out_edges = graph.out_edges(node_idx)
        if out_edges:
            downstream_nodes = [tgt_idx for _, tgt_idx, _ in out_edges]
            if downstream_nodes:
                downstream_idx = max(downstream_nodes, key=lambda idx: graph[idx]["pathlength"])
                graph[node_idx]["pathlength"] = (
                    graph[downstream_idx]["pathlength"] + graph[downstream_idx]["lengthkm"]
                )

    # PASS 2: Calculate totdasqkm and stream_order (forward topo order - downstream from headwaters)
    for node_idx in topo_order:
        in_edges = list(graph.in_edges(node_idx))

        # Accumulate upstream drainage area
        upstream_total = sum(graph[src_idx]["totdasqkm"] for src_idx, _, _ in in_edges)
        graph[node_idx]["totdasqkm"] = upstream_total + graph[node_idx]["areasqkm"]

    # Trace mainstems for each outlet's basin
    current_mainstem_id = 1
    processed: set[int] = set()

    for outlet_idx in outlets:
        # Trace main mainstem (longest path from outlet to headwater)
        current_idx = outlet_idx

        while current_idx not in processed:
            graph[current_idx]["mainstemlp"] = current_mainstem_id
            processed.add(current_idx)

            in_edges = list(graph.in_edges(current_idx))
            if not in_edges:
                break

            upstream_candidates = [src_idx for src_idx, _, _ in in_edges if src_idx not in processed]
            if not upstream_candidates:
                break

            current_idx = max(
                upstream_candidates,
                key=lambda idx: (graph[idx]["pathlength"], graph[idx]["totdasqkm"]),
            )

        current_mainstem_id += 1

    # Assign tributary mainstems for remaining nodes
    for node_idx in graph.node_indices():
        if node_idx not in processed:
            tributary_id = current_mainstem_id
            current_mainstem_id += 1

            trib_current = node_idx
            while trib_current not in processed:
                graph[trib_current]["mainstemlp"] = tributary_id
                processed.add(trib_current)

                in_edges = list(graph.in_edges(trib_current))
                upstream_in_basin = [src_idx for src_idx, _, _ in in_edges if src_idx not in processed]

                if not upstream_in_basin:
                    break

                trib_current = max(
                    upstream_in_basin,
                    key=lambda idx: (graph[idx]["pathlength"], graph[idx]["totdasqkm"]),
                )

    # Assign dnhydroseq and flowpath_toid based on graph edges
    for node_idx in graph.node_indices():
        out_edges = graph.out_edges(node_idx)
        downstream_nodes = [tgt_idx for _, tgt_idx, _ in out_edges]

        if downstream_nodes:
            downstream_idx = downstream_nodes[0]
            graph[node_idx]["dnhydroseq"] = graph[downstream_idx]["hydroseq"]
            graph[node_idx]["flowpath_toid"] = graph[downstream_idx]["flowpath_id"]
        else:
            graph[node_idx]["dnhydroseq"] = 0
            graph[node_idx]["flowpath_toid"] = "0"

    # PASS 3: Orient all flowpath geometries so they flow upstream -> downstream
    for node_idx in graph.node_indices():
        in_edges = graph.in_edges(node_idx)
        upstream_nodes = [src_idx for src_idx, _, _ in in_edges]

        # Check if this is an outlet (flowpath_toid == "0" means no downstream)
        is_outlet = graph[node_idx]["flowpath_toid"] == "0"

        if is_outlet and upstream_nodes:
            # Outlet: use upstream geometry
            upstream_idx = upstream_nodes[0]
            us_geom = graph[upstream_idx]["flowpath_geometry"]
            graph[node_idx]["flowpath_geometry"] = _orient_flowpath_downstream(
                graph[node_idx]["flowpath_geometry"], ds_geom=None, us_geom=us_geom
            )
        else:
            # Normal case: use downstream geometry
            out_edges = graph.out_edges(node_idx)
            downstream_nodes = [tgt_idx for _, tgt_idx, _ in out_edges]

            ds_geom = None
            if downstream_nodes:
                downstream_idx = downstream_nodes[0]
                ds_geom = graph[downstream_idx]["flowpath_geometry"]

            graph[node_idx]["flowpath_geometry"] = _orient_flowpath_downstream(
                graph[node_idx]["flowpath_geometry"], ds_geom=ds_geom, us_geom=None
            )

    # Extract results for flowpaths
    flowpath_ids = []
    flowpath_toids = []
    vpu_ids = []
    das = []
    lengthkms = []
    total_das = []
    mainstems = []
    pathlengths = []
    dnhydroseqs = []
    hydroseqs = []
    streamorders = []
    flowpath_geometries = []

    for node_idx in graph.node_indices():
        node_data = graph[node_idx]
        flowpath_ids.append(node_data["flowpath_id"])
        flowpath_toids.append(node_data["flowpath_toid"])
        vpu_ids.append(vpu_id)
        das.append(node_data["areasqkm"])
        lengthkms.append(node_data["lengthkm"])
        total_das.append(node_data["totdasqkm"])
        mainstems.append(node_data["mainstemlp"])
        pathlengths.append(node_data["pathlength"])
        dnhydroseqs.append(node_data["dnhydroseq"])
        hydroseqs.append(node_data["hydroseq"])
        streamorders.append(node_data["streamorder"])
        flowpath_geometries.append(node_data["flowpath_geometry"])

    return gpd.GeoDataFrame(
        {
            "flowpath_id": flowpath_ids,
            "flowpath_toid": flowpath_toids,
            "VPUID": vpu_ids,
            "lengthkm": lengthkms,
            "areasqkm": das,
            "totdasqkm": total_das,
            "mainstemlp": mainstems,
            "pathlength": pathlengths,
            "dnhydroseq": dnhydroseqs,
            "hydroseq": hydroseqs,
            "streamorder": streamorders,
        },
        geometry=flowpath_geometries,
        crs="EPSG:3857",
    )


def _create_reference_divides(
    divides_df: gpd.GeoDataFrame, reference_flowpaths: gpd.GeoDataFrame, vpu_id: str
) -> gpd.GeoDataFrame:
    """A function to create the reference divides table

    Parameters
    ----------
    divides_df : gpd.GeoDataFrame
        the NHDCatchments table
    reference_flowpaths : gpd.GeoDataFrame
        The reference flowpaths
    vpu_id : str
        the VPUID we're working in

    Returns
    -------
    gpd.GeoDataFrame
        the outputted reference_divides
    """
    reference_divides = divides_df.rename(
        columns={"NHDPlusID": "divide_id", "VPUID": "vpuid", "AreaSqKm": "areasqkm"}
    )
    reference_divides["divide_id"] = reference_divides["divide_id"].astype(int).astype(str)
    reference_divides["vpuid"] = vpu_id
    mask = reference_divides["divide_id"].isin(reference_flowpaths["flowpath_id"])
    reference_divides["has_flowpath"] = mask
    reference_divides["flowpath_id"] = pd.NA
    reference_divides.loc[mask, "flowpath_id"] = reference_divides.loc[mask, "divide_id"]
    return reference_divides


def _create_geoglows_reference_divides(
    catchments_df: gpd.GeoDataFrame, reference_flowpaths: gpd.GeoDataFrame, vpu_id: str
) -> gpd.GeoDataFrame:
    """A function to create the reference divides table from GeoGLOWS catchments

    Parameters
    ----------
    catchments_df : gpd.GeoDataFrame
        The GeoGLOWS catchments table with linkno, areasqkm, and geometry
    reference_flowpaths : gpd.GeoDataFrame
        The reference flowpaths
    vpu_id : str
        The VPUID we're working in

    Returns
    -------
    gpd.GeoDataFrame
        The outputted reference_divides with catchment geometries
    """
    reference_divides = catchments_df.copy()
    reference_divides = reference_divides.rename(columns={"linkno": "divide_id"})
    reference_divides["divide_id"] = reference_divides["divide_id"].astype(int).astype(str)
    reference_divides["vpuid"] = vpu_id

    # Filter to only include catchments that have a corresponding flowpath
    mask = reference_divides["divide_id"].isin(reference_flowpaths["flowpath_id"])
    reference_divides["has_flowpath"] = mask
    reference_divides["flowpath_id"] = pd.NA
    reference_divides.loc[mask, "flowpath_id"] = reference_divides.loc[mask, "divide_id"]

    return reference_divides


def build_nhd_reference(**context: dict[str, Any]) -> dict[str, Any]:
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
    dict[str, Any]
        The reference flowpath and divides references in memory
    """
    ti = cast(TaskInstance, context["ti"])
    cfg = cast(ReferenceConfig, context["config"])
    graph: rx.PyDiGraph = ti.xcom_pull(task_id="build_nhd_graphs", key="graph")
    node_indices: dict[str, int] = ti.xcom_pull(task_id="build_nhd_graphs", key="node_indices")
    _flowpaths: pl.DataFrame = ti.xcom_pull(task_id="download", key="nhd_flowpaths")
    _divides: pl.DataFrame = ti.xcom_pull(task_id="download", key="nhd_divides")
    cycles_iter = rx.simple_cycles(graph)
    cycles: list[list[str]] = []
    cycle_ids: set[str] = set()
    for cycle in cycles_iter:
        _ids: list[Any] = [graph.get_node_data(node_idx) for node_idx in cycle]
        cycles.append(_ids)
        cycle_ids.update(_ids)
    if cycle_ids:
        raise NotImplementedError("Cycle Detected. Please create method for removing")

    _flowpaths_df = gpd.GeoDataFrame(
        _flowpaths.select(
            [
                pl.col("NHDPlusID"),
                pl.col("VPUID"),
                pl.col("LengthKM"),
                pl.col("fcode_description"),
            ]
        ).to_pandas(),
        geometry=gpd.GeoSeries.from_wkb(_flowpaths["geometry"]),
        crs="EPSG:4269",
    )

    _divides_df = gpd.GeoDataFrame(
        _divides.select(
            [
                pl.col("NHDPlusID"),
                pl.col("VPUID"),
                pl.col("AreaSqKm"),
            ]
        ).to_pandas(),
        geometry=gpd.GeoSeries.from_wkb(_divides["geometry"]),
        crs="EPSG:4269",
    )
    reference_flowpaths = _trace_attributes(graph, node_indices, _flowpaths_df, _divides_df, cfg.vpu_id)
    reference_divides = _create_reference_divides(_divides_df, reference_flowpaths, cfg.vpu_id)

    return {"reference_flowpaths": reference_flowpaths, "reference_divides": reference_divides}


def build_geoglows_reference(**context: dict[str, Any]) -> dict[str, Any]:
    """Builds reference fabric from GeoGLOWS data

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
        The reference flowpaths and divides in memory
    """
    ti = cast(TaskInstance, context["ti"])
    cfg = cast(ReferenceConfig, context["config"])
    graph: rx.PyDiGraph = ti.xcom_pull(task_id="build_geoglows_graphs", key="graph")
    node_indices: dict[str, int] = ti.xcom_pull(task_id="build_geoglows_graphs", key="node_indices")
    _flowpaths: pl.DataFrame = ti.xcom_pull(task_id="download", key="geoglows_flowpaths")
    _catchments: pl.DataFrame = ti.xcom_pull(task_id="download", key="geoglows_divides")

    # Check for cycles
    cycles_iter = rx.simple_cycles(graph)
    cycles: list[list[str]] = []
    cycle_ids: set[str] = set()
    for cycle in cycles_iter:
        _ids: list[Any] = [graph.get_node_data(node_idx) for node_idx in cycle]
        cycles.append(_ids)
        cycle_ids.update(_ids)
    if cycle_ids:
        raise NotImplementedError("Cycle Detected. Please create method for removing")

    _flowpaths_df = gpd.GeoDataFrame(
        _flowpaths.select(
            [
                pl.col("LINKNO"),
                pl.col("DSLINKNO"),
                pl.col("strmOrder"),
            ]
        ).to_pandas(),
        geometry=gpd.GeoSeries.from_wkb(_flowpaths["geometry"]),
        crs="EPSG:3857",
    )

    _flowpaths_df_projected = _flowpaths_df.to_crs(cfg.crs)
    _flowpaths_df["LengthKM"] = _flowpaths_df_projected.geometry.length / 1000

    _catchments_df = gpd.GeoDataFrame(
        _catchments.select([pl.col("linkno")]).to_pandas(),
        geometry=gpd.GeoSeries.from_wkb(_catchments["geometry"]),
        crs="EPSG:4326",
    )
    _catchments_df_projected = _catchments_df.to_crs(cfg.crs)
    _catchments_df["areasqkm"] = _catchments_df_projected.geometry.area / 1e6

    # Log any flowpaths without matching catchments
    flowpath_linkno_set = set(_flowpaths_df["LINKNO"].tolist())
    catchment_linkno_set = set(_catchments_df["linkno"].tolist())
    missing_catchments = flowpath_linkno_set - catchment_linkno_set
    if missing_catchments:
        logger.warning(
            f"build_geoglows_reference task: {len(missing_catchments)} flowpaths have no matching catchment"
        )

    reference_flowpaths = _trace_geoglows_attributes(
        graph, node_indices, _flowpaths_df, _catchments_df, cfg.vpu_id
    )
    reference_divides = _create_geoglows_reference_divides(_catchments_df, reference_flowpaths, cfg.vpu_id)

    return {"reference_flowpaths": reference_flowpaths, "reference_divides": reference_divides}
