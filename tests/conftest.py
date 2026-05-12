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


@pytest.fixture
def sample_geoglows_graph() -> tuple[rx.PyDiGraph, dict[str, int]]:
    """Create a sample GeoGLOWS graph for testing.

    Network topology:
        1 ──┐
            ├──> 3 ──> 4 ──> 5 (outlet)
        2 ──┘

    Where:
        - 1, 2 are headwaters (order 1)
        - 3 is confluence (order 2)
        - 4, 5 continue downstream (order 2)
    """
    graph = rx.PyDiGraph()

    node_data = [
        "810000001",
        "810000002",
        "810000003",
        "810000004",
        "810000005",
    ]

    node_indices = {}
    for fp_id in node_data:
        idx = graph.add_node(fp_id)
        node_indices[fp_id] = idx

    # Add edges (upstream -> downstream)
    graph.add_edge(node_indices["810000001"], node_indices["810000003"], None)
    graph.add_edge(node_indices["810000002"], node_indices["810000003"], None)
    graph.add_edge(node_indices["810000003"], node_indices["810000004"], None)
    graph.add_edge(node_indices["810000004"], node_indices["810000005"], None)

    return graph, node_indices


@pytest.fixture
def sample_geoglows_flowpaths() -> gpd.GeoDataFrame:
    """Create sample GeoGLOWS flowpaths GeoDataFrame."""
    return gpd.GeoDataFrame(
        {
            "LINKNO": [810000001, 810000002, 810000003, 810000004, 810000005],
            "DSLINKNO": [810000003, 810000003, 810000004, 810000005, -1],
            "strmOrder": [1, 1, 2, 2, 2],
            "LengthKM": [0.5, 0.6, 0.8, 0.7, 0.4],
        },
        geometry=[
            LineString([(0, 0), (1, 1)]),
            LineString([(2, 0), (1, 1)]),
            LineString([(1, 1), (1, 2)]),
            LineString([(1, 2), (1, 3)]),
            LineString([(1, 3), (1, 4)]),
        ],
        crs="EPSG:3857",
    )


@pytest.fixture
def sample_geoglows_catchments() -> gpd.GeoDataFrame:
    """Create sample GeoGLOWS catchments GeoDataFrame."""
    return gpd.GeoDataFrame(
        {
            "linkno": [810000001, 810000002, 810000003, 810000004, 810000005],
            "areasqkm": [1.0, 1.2, 0.8, 0.9, 0.6],
        },
        geometry=[
            MultiPolygon([Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]),
            MultiPolygon([Polygon([(1, 0), (2, 0), (2, 1), (1, 1)])]),
            MultiPolygon([Polygon([(0, 1), (2, 1), (2, 2), (0, 2)])]),
            MultiPolygon([Polygon([(0, 2), (2, 2), (2, 3), (0, 3)])]),
            MultiPolygon([Polygon([(0, 3), (2, 3), (2, 4), (0, 4)])]),
        ],
        crs="EPSG:3857",
    )


@pytest.fixture
def disconnected_geoglows_graph() -> tuple[rx.PyDiGraph, dict[str, int]]:
    """Create a disconnected GeoGLOWS graph with two separate networks.

    Network 1: 1 -> 3 (outlet)
    Network 2: 2 -> 4 -> 5 (outlet)
    """
    graph = rx.PyDiGraph()

    node_data = [
        "810000001",
        "810000002",
        "810000003",
        "810000004",
        "810000005",
    ]

    node_indices = {}
    for fp_id in node_data:
        idx = graph.add_node(fp_id)
        node_indices[fp_id] = idx

    # Network 1
    graph.add_edge(node_indices["810000001"], node_indices["810000003"], None)
    # Network 2
    graph.add_edge(node_indices["810000002"], node_indices["810000004"], None)
    graph.add_edge(node_indices["810000004"], node_indices["810000005"], None)

    return graph, node_indices


@pytest.fixture
def reversed_flowpaths() -> gpd.GeoDataFrame:
    """Create sample flowpaths with reversed digitization (downstream to upstream)."""
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
            "Stream/River",
            "Artificial Path",
            "Stream/River",
            "Artificial Path",
            "Stream/River",
        ],
    }

    # Create geometries that connect but are digitized downstream-to-upstream
    # Network: 1,2 -> 3 -> 4 -> 5 (outlet)
    geometries = [
        MultiLineString([LineString([(1, 3), (0, 4)])]),  # 1: connects to 3 at (1,3)
        MultiLineString([LineString([(1, 3), (2, 4)])]),  # 2: connects to 3 at (1,3)
        MultiLineString([LineString([(1, 2), (1, 3)])]),  # 3: connects to 4 at (1,2)
        MultiLineString([LineString([(1, 1), (1, 2)])]),  # 4: connects to 5 at (1,1)
        MultiLineString([LineString([(1, 0), (1, 1)])]),  # 5: outlet at (1,0)
    ]

    return gpd.GeoDataFrame(data, geometry=geometries, crs="EPSG:4269")


@pytest.fixture
def correctly_oriented_flowpaths() -> gpd.GeoDataFrame:
    """Create sample flowpaths with correct digitization (upstream to downstream)."""
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
            "Stream/River",
            "Artificial Path",
            "Stream/River",
            "Artificial Path",
            "Stream/River",
        ],
    }

    # Create geometries digitized upstream-to-downstream
    # Network: 1,2 -> 3 -> 4 -> 5 (outlet)
    geometries = [
        MultiLineString([LineString([(0, 4), (1, 3)])]),  # 1: ends at confluence (1,3)
        MultiLineString([LineString([(2, 4), (1, 3)])]),  # 2: ends at confluence (1,3)
        MultiLineString([LineString([(1, 3), (1, 2)])]),  # 3: from confluence to (1,2)
        MultiLineString([LineString([(1, 2), (1, 1)])]),  # 4: continues downstream
        MultiLineString([LineString([(1, 1), (1, 0)])]),  # 5: outlet ends at (1,0)
    ]

    return gpd.GeoDataFrame(data, geometry=geometries, crs="EPSG:4269")
