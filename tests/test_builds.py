"""Tests for build_reference module"""

import geopandas as gpd
import numpy as np
import pytest
import rustworkx as rx
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Polygon

from reference_builds.pipeline.build_reference import (
    _create_reference_divides,
    _trace_attributes,
)


class TestTraceAttributes:
    """Tests for _trace_attributes function."""

    def test_output_columns(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that output has expected columns."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, sample_flowpaths, sample_divides, "21")

        expected_columns = {
            "flowpath_id",
            "VPUID",
            "areasqkm",
            "totdasqkm",
            "mainstemlp",
            "pathlength",
            "dnhydroseq",
            "hydroseq",
            "stream_order",
            "fcode_description",
            "geometry",
        }

        assert set(result.columns) == expected_columns

    def test_all_flowpaths_traced(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that all flowpaths in graph are traced."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, sample_flowpaths, sample_divides, "21")

        assert len(result) == graph.num_nodes()

    def test_hydroseq_unique(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that hydroseq values are unique."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, sample_flowpaths, sample_divides, "21")

        assert result["hydroseq"].nunique() == len(result)

    def test_outlet_has_zero_dnhydroseq(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that outlet node has dnhydroseq of 0."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, sample_flowpaths, sample_divides, "21")

        # Outlet is 85000100000005
        outlet_row = result[result["flowpath_id"] == "85000100000005"]
        assert outlet_row["dnhydroseq"].iloc[0] == 0

    def test_outlet_has_zero_pathlength(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that outlet node has pathlength of 0."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, sample_flowpaths, sample_divides, "21")

        outlet_row = result[result["flowpath_id"] == "85000100000005"]
        assert outlet_row["pathlength"].iloc[0] == 0.0

    def test_headwaters_have_stream_order_1(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that headwater nodes have stream order 1."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, sample_flowpaths, sample_divides, "21")

        # Headwaters are 85000100000001 and 85000100000002
        headwater_rows = result[result["flowpath_id"].isin(["85000100000001", "85000100000002"])]
        assert (headwater_rows["stream_order"] == 1).all()

    def test_strahler_order_increases_at_confluence(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that stream order increases when two streams of same order meet."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, sample_flowpaths, sample_divides, "21")

        # Node 3 is confluence of two order-1 streams, should be order 2
        confluence_row = result[result["flowpath_id"] == "85000100000003"]
        assert confluence_row["stream_order"].iloc[0] == 2

    def test_totdasqkm_accumulates_downstream(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that total drainage area accumulates downstream."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, sample_flowpaths, sample_divides, "21")

        # Outlet should have largest totdasqkm
        outlet_row = result[result["flowpath_id"] == "85000100000005"]
        max_da = result["totdasqkm"].max()

        assert outlet_row["totdasqkm"].iloc[0] == max_da

    def test_totdasqkm_equals_sum_of_all_areas_at_outlet(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that outlet totdasqkm equals sum of all upstream areas."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, sample_flowpaths, sample_divides, "21")

        outlet_row = result[result["flowpath_id"] == "85000100000005"]
        total_area = sample_divides["AreaSqKm"].sum()

        assert np.isclose(outlet_row["totdasqkm"].iloc[0], total_area)

    def test_pathlength_increases_upstream(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that pathlength increases going upstream."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, sample_flowpaths, sample_divides, "21")

        # Headwaters should have larger pathlength than outlet
        headwater_row = result[result["flowpath_id"] == "85000100000001"]
        outlet_row = result[result["flowpath_id"] == "85000100000005"]

        assert headwater_row["pathlength"].iloc[0] > outlet_row["pathlength"].iloc[0]

    def test_mainstemlp_assigned_to_all_nodes(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that all nodes have a mainstem level path assigned."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, sample_flowpaths, sample_divides, "21")

        assert result["mainstemlp"].notna().all()

    def test_multiple_outlets_handled(
        self,
        disconnected_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that disconnected subgraphs with multiple outlets are handled."""
        graph, node_indices = disconnected_graph
        result = _trace_attributes(graph, node_indices, sample_flowpaths, sample_divides, "21")

        # Should have two outlets with dnhydroseq = 0
        outlets = result[result["dnhydroseq"] == 0]
        assert len(outlets) == 2


class TestCreateReferenceDivides:
    """Tests for _create_reference_divides function."""

    @pytest.fixture
    def sample_reference_flowpaths(self, sample_flowpaths: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        """Create sample reference flowpaths (subset of divides)."""
        # Only include some flowpaths (simulating filtered network)
        return gpd.GeoDataFrame(
            {
                "flowpath_id": ["85000100000001", "85000100000003", "85000100000005"],
                "VPUID": ["21", "21", "21"],
            },
            geometry=sample_flowpaths.geometry.iloc[:3].values,
            crs="EPSG:4269",
        )

    def test_returns_geodataframe(
        self,
        sample_divides: gpd.GeoDataFrame,
        sample_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that _create_reference_divides returns a GeoDataFrame."""
        result = _create_reference_divides(sample_divides, sample_reference_flowpaths, "21")

        assert isinstance(result, gpd.GeoDataFrame)

    def test_all_divides_included(
        self,
        sample_divides: gpd.GeoDataFrame,
        sample_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that all divides are included in output."""
        result = _create_reference_divides(sample_divides, sample_reference_flowpaths, "21")

        assert len(result) == len(sample_divides)

    def test_has_flowpath_flag_correct(
        self,
        sample_divides: gpd.GeoDataFrame,
        sample_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that has_flowpath flag is set correctly."""
        result = _create_reference_divides(sample_divides, sample_reference_flowpaths, "21")

        # 3 divides should have flowpaths
        assert result["has_flowpath"].sum() == 3

    def test_flowpath_id_assigned_when_has_flowpath(
        self,
        sample_divides: gpd.GeoDataFrame,
        sample_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that flowpath_id is assigned when has_flowpath is True."""
        result = _create_reference_divides(sample_divides, sample_reference_flowpaths, "21")

        with_flowpath = result[result["has_flowpath"]]
        assert with_flowpath["flowpath_id"].notna().all()
        assert (with_flowpath["flowpath_id"] == with_flowpath["divide_id"]).all()

    def test_flowpath_id_na_when_no_flowpath(
        self,
        sample_divides: gpd.GeoDataFrame,
        sample_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that flowpath_id is NA when has_flowpath is False."""
        result = _create_reference_divides(sample_divides, sample_reference_flowpaths, "21")

        without_flowpath = result[~result["has_flowpath"]]
        assert without_flowpath["flowpath_id"].isna().all()


class TestGraphWithCycles:
    """Tests for handling graphs with cycles."""

    @pytest.fixture
    def cyclic_graph(self) -> tuple[rx.PyDiGraph, dict[str, int]]:
        """Create a graph with a cycle.

        1 -> 2 -> 3 -> 1 (cycle)
        """
        graph = rx.PyDiGraph()

        node_data = ["85000100000001", "85000100000002", "85000100000003"]

        node_indices = {}
        for fp_id in node_data:
            idx = graph.add_node(fp_id)
            node_indices[fp_id] = idx

        graph.add_edge(node_indices["85000100000001"], node_indices["85000100000002"], None)
        graph.add_edge(node_indices["85000100000002"], node_indices["85000100000003"], None)
        graph.add_edge(node_indices["85000100000003"], node_indices["85000100000001"], None)

        return graph, node_indices

    def test_trace_attributes_raises_on_cycle(
        self,
        cyclic_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that _trace_attributes raises AssertionError on cyclic graph."""
        graph, node_indices = cyclic_graph

        with pytest.raises(AssertionError, match="Graph contains cycles"):
            _trace_attributes(graph, node_indices, sample_flowpaths, sample_divides, "21")


class TestEdgeCases:
    """Tests for edge cases."""

    @pytest.fixture
    def single_node_graph(self) -> tuple[rx.PyDiGraph, dict[str, int]]:
        """Create a graph with a single node (headwater that is also outlet)."""
        graph = rx.PyDiGraph()
        idx = graph.add_node("85000100000001")
        return graph, {"85000100000001": idx}

    @pytest.fixture
    def single_flowpath(self) -> gpd.GeoDataFrame:
        """Create single flowpath GeoDataFrame."""
        return gpd.GeoDataFrame(
            {
                "NHDPlusID": [85000100000001],
                "VPUID": ["2101"],
                "LengthKM": [0.5],
                "fcode_description": ["Stream/River"],
            },
            geometry=[MultiLineString([LineString([(0, 0), (1, 1)])])],
            crs="EPSG:4269",
        )

    @pytest.fixture
    def single_divide(self) -> gpd.GeoDataFrame:
        """Create single divide GeoDataFrame."""
        return gpd.GeoDataFrame(
            {
                "NHDPlusID": [85000100000001],
                "VPUID": ["2101"],
                "AreaSqKm": [1.0],
            },
            geometry=[MultiPolygon([Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])])],
            crs="EPSG:4269",
        )

    def test_single_node_graph(
        self,
        single_node_graph: tuple[rx.PyDiGraph, dict[str, int]],
        single_flowpath: gpd.GeoDataFrame,
        single_divide: gpd.GeoDataFrame,
    ) -> None:
        """Test that single node graph is handled correctly."""
        graph, node_indices = single_node_graph
        result = _trace_attributes(graph, node_indices, single_flowpath, single_divide, "21")

        assert len(result) == 1
        assert result["stream_order"].iloc[0] == 1
        assert result["dnhydroseq"].iloc[0] == 0
        assert result["pathlength"].iloc[0] == 0.0
        assert result["totdasqkm"].iloc[0] == 1.0

    def test_missing_divide_defaults_to_zero_area(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that missing divides default to zero area."""
        graph, node_indices = sample_graph

        # Create divides missing one entry
        partial_divides = gpd.GeoDataFrame(
            {
                "NHDPlusID": [85000100000001, 85000100000002],
                "VPUID": ["2101", "2101"],
                "AreaSqKm": [0.5, 0.5],
            },
            geometry=[
                MultiPolygon([Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])]),
                MultiPolygon([Polygon([(1, 0), (2, 0), (2, 1), (1, 1)])]),
            ],
            crs="EPSG:4269",
        )

        result = _trace_attributes(graph, node_indices, sample_flowpaths, partial_divides, "21")

        # Nodes without divides should have areasqkm = 0
        missing_divide_rows = result[~result["flowpath_id"].isin(["85000100000001", "85000100000002"])]
        assert (missing_divide_rows["areasqkm"] == 0.0).all()
