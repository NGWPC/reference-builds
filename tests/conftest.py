"""Conftests for Pytest Suite"""

import geopandas as gpd
import pytest
import rustworkx as rx
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon


@pytest.fixture
def sample_graph() -> tuple[rx.PyDiGraph, dict[str, int]]:
    """Create a simple directed graph for testing.

    Graph structure:
        1 -> 3
        2 -> 3
        3 -> 4
        4 -> 5 (outlet)
    """
    graph = rx.PyDiGraph()

    # Add nodes (using NHDPlusID as string)
    node_data = [
        "85000100000001",
        "85000100000002",
        "85000100000003",
        "85000100000004",
        "85000100000005",
    ]

    node_indices = {}
    for fp_id in node_data:
        idx = graph.add_node(fp_id)
        node_indices[fp_id] = idx

    # Add edges (upstream -> downstream)
    graph.add_edge(node_indices["85000100000001"], node_indices["85000100000003"], None)
    graph.add_edge(node_indices["85000100000002"], node_indices["85000100000003"], None)
    graph.add_edge(node_indices["85000100000003"], node_indices["85000100000004"], None)
    graph.add_edge(node_indices["85000100000004"], node_indices["85000100000005"], None)

    return graph, node_indices


@pytest.fixture
def sample_flowpaths() -> gpd.GeoDataFrame:
    """Create sample flowpaths GeoDataFrame."""
    data = {
        "NHDPlusID": [
            85000100000001,
            85000100000002,
            85000100000003,
            85000100000004,
            85000100000005,
        ],
        "VPUID": ["2101", "2101", "2101", "2101", "2101"],
        "LengthKM": [0.333, 0.5, 1.201, 0.182, 0.387],
        "fcode_description": [
            "Stream/River: Hydrographic Category = Intermittent",
            "Artificial Path",
            "Stream/River: Hydrographic Category = Intermittent",
            "Artificial Path",
            "Stream/River: Hydrographic Category = Intermittent",
        ],
    }

    # Create simple linestring geometries
    geometries = [MultiLineString([LineString([(0, i), (1, i)])]) for i in range(5)]

    return gpd.GeoDataFrame(data, geometry=geometries, crs="EPSG:4269")


@pytest.fixture
def sample_divides() -> gpd.GeoDataFrame:
    """Create sample divides GeoDataFrame."""
    data = {
        "NHDPlusID": [
            85000100000001,
            85000100000002,
            85000100000003,
            85000100000004,
            85000100000005,
        ],
        "VPUID": ["2101", "2101", "2101", "2101", "2101"],
        "AreaSqKm": [0.1602, 0.6949, 0.0248, 0.1413, 0.9395],
    }

    # Create simple polygon geometries
    geometries = [MultiPolygon([Polygon([(i, 0), (i + 1, 0), (i + 1, 1), (i, 1)])]) for i in range(5)]

    return gpd.GeoDataFrame(data, geometry=geometries, crs="EPSG:4269")


@pytest.fixture
def disconnected_graph() -> tuple[rx.PyDiGraph, dict[str, int]]:
    """Create a graph with multiple disconnected subgraphs (multiple outlets).

    Subgraph 1:
        1 -> 2 (outlet)

    Subgraph 2:
        3 -> 4 -> 5 (outlet)
    """
    graph = rx.PyDiGraph()

    node_data = [
        "85000100000001",
        "85000100000002",
        "85000100000003",
        "85000100000004",
        "85000100000005",
    ]

    node_indices = {}
    for fp_id in node_data:
        idx = graph.add_node(fp_id)
        node_indices[fp_id] = idx

    # Subgraph 1
    graph.add_edge(node_indices["85000100000001"], node_indices["85000100000002"], None)

    # Subgraph 2
    graph.add_edge(node_indices["85000100000003"], node_indices["85000100000004"], None)
    graph.add_edge(node_indices["85000100000004"], node_indices["85000100000005"], None)

    return graph, node_indices
