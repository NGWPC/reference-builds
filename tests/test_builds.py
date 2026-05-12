"""Tests for build_nhd_reference module"""

import geopandas as gpd
import numpy as np
import pytest
import rustworkx as rx
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Point, Polygon

from reference_builds.pipeline.build_reference import (
    _create_geoglows_reference_divides,
    _create_reference_divides,
    _trace_attributes,
    _trace_geoglows_attributes,
)


class TestTraceGeoglowsAttributes:
    """Tests for _trace_geoglows_attributes function."""

    def test_output_columns(
        self,
        sample_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that output has expected columns."""
        graph, node_indices = sample_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, sample_geoglows_flowpaths, sample_geoglows_catchments, "701"
        )

        expected_columns = {
            "flowpath_id",
            "flowpath_toid",
            "VPUID",
            "lengthkm",
            "areasqkm",
            "totdasqkm",
            "mainstemlp",
            "pathlength",
            "dnhydroseq",
            "hydroseq",
            "streamorder",
            "geometry",
        }

        assert set(result.columns) == expected_columns

    def test_all_flowpaths_traced(
        self,
        sample_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that all flowpaths in graph are traced."""
        graph, node_indices = sample_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, sample_geoglows_flowpaths, sample_geoglows_catchments, "701"
        )

        assert len(result) == graph.num_nodes()

    def test_hydroseq_unique(
        self,
        sample_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that hydroseq values are unique."""
        graph, node_indices = sample_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, sample_geoglows_flowpaths, sample_geoglows_catchments, "701"
        )

        assert result["hydroseq"].nunique() == len(result)

    def test_totdasqkm_accumulates_downstream(
        self,
        sample_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that total drainage area accumulates downstream."""
        graph, node_indices = sample_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, sample_geoglows_flowpaths, sample_geoglows_catchments, "701"
        )

        outlet_row = result[result["flowpath_id"] == "810000005"]
        max_da = result["totdasqkm"].max()

        assert outlet_row["totdasqkm"].iloc[0] == max_da

    def test_totdasqkm_equals_sum_of_all_areas_at_outlet(
        self,
        sample_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that outlet totdasqkm equals sum of all upstream areas."""
        graph, node_indices = sample_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, sample_geoglows_flowpaths, sample_geoglows_catchments, "701"
        )

        outlet_row = result[result["flowpath_id"] == "810000005"]
        total_area = sample_geoglows_catchments["areasqkm"].sum()

        assert np.isclose(outlet_row["totdasqkm"].iloc[0], total_area)

    def test_pathlength_increases_upstream(
        self,
        sample_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that pathlength increases going upstream."""
        graph, node_indices = sample_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, sample_geoglows_flowpaths, sample_geoglows_catchments, "701"
        )

        headwater_row = result[result["flowpath_id"] == "810000001"]
        outlet_row = result[result["flowpath_id"] == "810000005"]

        assert headwater_row["pathlength"].iloc[0] > outlet_row["pathlength"].iloc[0]

    def test_pathlength_calculation_correct(
        self,
        sample_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that pathlength is calculated correctly as sum of downstream lengths."""
        graph, node_indices = sample_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, sample_geoglows_flowpaths, sample_geoglows_catchments, "701"
        )

        # Headwater 1's pathlength should be sum of downstream lengths (5->4->3)
        # 0.4 + 0.7 + 0.8 = 1.9
        headwater_row = result[result["flowpath_id"] == "810000001"]
        expected_pathlength = 0.4 + 0.7 + 0.8

        assert np.isclose(headwater_row["pathlength"].iloc[0], expected_pathlength)

    def test_mainstemlp_assigned_to_all_nodes(
        self,
        sample_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that all nodes have a mainstem level path assigned."""
        graph, node_indices = sample_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, sample_geoglows_flowpaths, sample_geoglows_catchments, "701"
        )

        assert result["mainstemlp"].notna().all()

    def test_multiple_outlets_handled(
        self,
        disconnected_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that disconnected subgraphs with multiple outlets are handled."""
        graph, node_indices = disconnected_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, sample_geoglows_flowpaths, sample_geoglows_catchments, "701"
        )

        outlets = result[result["dnhydroseq"] == 0]
        assert len(outlets) == 2

    def test_vpuid_assigned_correctly(
        self,
        sample_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that VPUID is assigned correctly to all rows."""
        graph, node_indices = sample_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, sample_geoglows_flowpaths, sample_geoglows_catchments, "701"
        )

        assert (result["VPUID"] == "701").all()

    def test_lengthkm_preserved_from_input(
        self,
        sample_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that lengthkm values are preserved from input."""
        graph, node_indices = sample_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, sample_geoglows_flowpaths, sample_geoglows_catchments, "701"
        )

        for _, row in sample_geoglows_flowpaths.iterrows():
            linkno = str(row["LINKNO"])
            expected_length = row["LengthKM"]
            actual_length = result[result["flowpath_id"] == linkno]["lengthkm"].iloc[0]
            assert np.isclose(actual_length, expected_length)

    def test_areasqkm_from_catchments(
        self,
        sample_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that areasqkm values come from catchments."""
        graph, node_indices = sample_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, sample_geoglows_flowpaths, sample_geoglows_catchments, "701"
        )

        for _, row in sample_geoglows_catchments.iterrows():
            linkno = str(row["linkno"])
            expected_area = row["areasqkm"]
            actual_area = result[result["flowpath_id"] == linkno]["areasqkm"].iloc[0]
            assert np.isclose(actual_area, expected_area)


class TestCreateGeoglowsReferenceDivides:
    """Tests for _create_geoglows_reference_divides function."""

    @pytest.fixture
    def sample_geoglows_reference_flowpaths(
        self, sample_geoglows_flowpaths: gpd.GeoDataFrame
    ) -> gpd.GeoDataFrame:
        """Create sample reference flowpaths (subset of catchments)."""
        return gpd.GeoDataFrame(
            {
                "flowpath_id": ["810000001", "810000003", "810000005"],
                "VPUID": ["701", "701", "701"],
            },
            geometry=sample_geoglows_flowpaths.geometry.iloc[:3].values,
            crs="EPSG:3857",
        )

    def test_returns_geodataframe(
        self,
        sample_geoglows_catchments: gpd.GeoDataFrame,
        sample_geoglows_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that _create_geoglows_reference_divides returns a GeoDataFrame."""
        result = _create_geoglows_reference_divides(
            sample_geoglows_catchments, sample_geoglows_reference_flowpaths, "701"
        )

        assert isinstance(result, gpd.GeoDataFrame)

    def test_all_catchments_included(
        self,
        sample_geoglows_catchments: gpd.GeoDataFrame,
        sample_geoglows_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that all catchments are included in output."""
        result = _create_geoglows_reference_divides(
            sample_geoglows_catchments, sample_geoglows_reference_flowpaths, "701"
        )

        assert len(result) == len(sample_geoglows_catchments)

    def test_has_flowpath_flag_correct(
        self,
        sample_geoglows_catchments: gpd.GeoDataFrame,
        sample_geoglows_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that has_flowpath flag is set correctly."""
        result = _create_geoglows_reference_divides(
            sample_geoglows_catchments, sample_geoglows_reference_flowpaths, "701"
        )

        assert result["has_flowpath"].sum() == 3

    def test_flowpath_id_assigned_when_has_flowpath(
        self,
        sample_geoglows_catchments: gpd.GeoDataFrame,
        sample_geoglows_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that flowpath_id is assigned when has_flowpath is True."""
        result = _create_geoglows_reference_divides(
            sample_geoglows_catchments, sample_geoglows_reference_flowpaths, "701"
        )

        with_flowpath = result[result["has_flowpath"]]
        assert with_flowpath["flowpath_id"].notna().all()
        assert (with_flowpath["flowpath_id"] == with_flowpath["divide_id"]).all()

    def test_flowpath_id_na_when_no_flowpath(
        self,
        sample_geoglows_catchments: gpd.GeoDataFrame,
        sample_geoglows_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that flowpath_id is NA when has_flowpath is False."""
        result = _create_geoglows_reference_divides(
            sample_geoglows_catchments, sample_geoglows_reference_flowpaths, "701"
        )

        without_flowpath = result[~result["has_flowpath"]]
        assert without_flowpath["flowpath_id"].isna().all()

    def test_vpuid_assigned_correctly(
        self,
        sample_geoglows_catchments: gpd.GeoDataFrame,
        sample_geoglows_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that vpuid is assigned correctly to all rows."""
        result = _create_geoglows_reference_divides(
            sample_geoglows_catchments, sample_geoglows_reference_flowpaths, "701"
        )

        assert (result["vpuid"] == "701").all()

    def test_divide_id_is_string(
        self,
        sample_geoglows_catchments: gpd.GeoDataFrame,
        sample_geoglows_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that divide_id is converted to string."""
        result = _create_geoglows_reference_divides(
            sample_geoglows_catchments, sample_geoglows_reference_flowpaths, "701"
        )

        assert result["divide_id"].dtype == object

    def test_areasqkm_preserved(
        self,
        sample_geoglows_catchments: gpd.GeoDataFrame,
        sample_geoglows_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that areasqkm values are preserved."""
        result = _create_geoglows_reference_divides(
            sample_geoglows_catchments, sample_geoglows_reference_flowpaths, "701"
        )

        assert "areasqkm" in result.columns
        original_areas = set(sample_geoglows_catchments["areasqkm"].tolist())
        result_areas = set(result["areasqkm"].tolist())
        assert original_areas == result_areas

    def test_geometry_preserved(
        self,
        sample_geoglows_catchments: gpd.GeoDataFrame,
        sample_geoglows_reference_flowpaths: gpd.GeoDataFrame,
    ) -> None:
        """Test that catchment geometries are preserved."""
        result = _create_geoglows_reference_divides(
            sample_geoglows_catchments, sample_geoglows_reference_flowpaths, "701"
        )

        assert result.geometry is not None
        assert len(result.geometry) == len(sample_geoglows_catchments)


# =============================================================================
# Tests for GeoGLOWS Edge Cases
# =============================================================================


class TestGeoglowsEdgeCases:
    """Tests for GeoGLOWS edge cases."""

    @pytest.fixture
    def single_node_geoglows_graph(self) -> tuple[rx.PyDiGraph, dict[str, int]]:
        """Create a graph with a single node (headwater that is also outlet)."""
        graph = rx.PyDiGraph()
        idx = graph.add_node("810000001")
        return graph, {"810000001": idx}

    @pytest.fixture
    def single_geoglows_flowpath(self) -> gpd.GeoDataFrame:
        """Create single GeoGLOWS flowpath GeoDataFrame."""
        return gpd.GeoDataFrame(
            {
                "LINKNO": [810000001],
                "DSLINKNO": [-1],
                "strmOrder": [1],
                "LengthKM": [0.5],
            },
            geometry=[LineString([(0, 0), (1, 1)])],
            crs="EPSG:3857",
        )

    @pytest.fixture
    def single_geoglows_catchment(self) -> gpd.GeoDataFrame:
        """Create single GeoGLOWS catchment GeoDataFrame."""
        return gpd.GeoDataFrame(
            {
                "linkno": [810000001],
                "areasqkm": [1.0],
            },
            geometry=[MultiPolygon([Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])])],
            crs="EPSG:3857",
        )

    def test_single_node_graph(
        self,
        single_node_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        single_geoglows_flowpath: gpd.GeoDataFrame,
        single_geoglows_catchment: gpd.GeoDataFrame,
    ) -> None:
        """Test that single node graph is handled correctly."""
        graph, node_indices = single_node_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, single_geoglows_flowpath, single_geoglows_catchment, "701"
        )

        assert len(result) == 1
        assert result["streamorder"].iloc[0] == 1
        assert result["dnhydroseq"].iloc[0] == 0
        assert result["pathlength"].iloc[0] == 0.0
        assert result["totdasqkm"].iloc[0] == 1.0


class TestGeoglowsGraphWithCycles:
    """Tests for handling GeoGLOWS graphs with cycles."""

    @pytest.fixture
    def cyclic_geoglows_graph(self) -> tuple[rx.PyDiGraph, dict[str, int]]:
        """Create a GeoGLOWS graph with a cycle."""
        graph = rx.PyDiGraph()

        node_data = ["810000001", "810000002", "810000003"]

        node_indices = {}
        for fp_id in node_data:
            idx = graph.add_node(fp_id)
            node_indices[fp_id] = idx

        graph.add_edge(node_indices["810000001"], node_indices["810000002"], None)
        graph.add_edge(node_indices["810000002"], node_indices["810000003"], None)
        graph.add_edge(node_indices["810000003"], node_indices["810000001"], None)

        return graph, node_indices

    def test_trace_geoglows_attributes_raises_on_cycle(
        self,
        cyclic_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        sample_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that _trace_geoglows_attributes raises AssertionError on cyclic graph."""
        graph, node_indices = cyclic_geoglows_graph

        with pytest.raises(AssertionError, match="Graph contains cycles"):
            _trace_geoglows_attributes(
                graph, node_indices, sample_geoglows_flowpaths, sample_geoglows_catchments, "701"
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
            "flowpath_toid",
            "VPUID",
            "lengthkm",
            "areasqkm",
            "totdasqkm",
            "mainstemlp",
            "pathlength",
            "dnhydroseq",
            "hydroseq",
            "streamorder",
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
        assert (headwater_rows["streamorder"] == 1).all()

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
        assert confluence_row["streamorder"].iloc[0] == 2

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
        assert result["streamorder"].iloc[0] == 1
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


class TestFlowpathOrientation:
    """Tests for flowpath geometry orientation in _trace_attributes."""

    def test_reversed_flowpaths_are_corrected(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        reversed_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that reversed flowpaths are oriented correctly."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, reversed_flowpaths, sample_divides, "21")

        # Check that flowpath 4 ends at the confluence with 5
        fp4 = result[result["flowpath_id"] == "85000100000004"].iloc[0]
        fp5 = result[result["flowpath_id"] == "85000100000005"].iloc[0]

        # Get endpoints
        from reference_builds.utils.geometries import _get_endpoints

        _, fp4_end = _get_endpoints(fp4.geometry)
        fp5_start, _ = _get_endpoints(fp5.geometry)

        # fp4 should end where fp5 starts (or very close)
        assert Point(fp4_end).distance(Point(fp5_start)) < 0.001

    def test_correctly_oriented_flowpaths_unchanged(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        correctly_oriented_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that correctly oriented flowpaths remain unchanged."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, correctly_oriented_flowpaths, sample_divides, "21")

        from reference_builds.utils.geometries import _get_endpoints

        # Check original and result have same orientation
        for _, row in correctly_oriented_flowpaths.iterrows():
            fp_id = str(int(row["NHDPlusID"]))
            original_start, original_end = _get_endpoints(row.geometry)

            result_row = result[result["flowpath_id"] == fp_id].iloc[0]
            result_start, result_end = _get_endpoints(result_row.geometry)

            assert original_start == result_start
            assert original_end == result_end

    def test_outlet_oriented_correctly(
        self,
        sample_graph: tuple[rx.PyDiGraph, dict[str, int]],
        reversed_flowpaths: gpd.GeoDataFrame,
        sample_divides: gpd.GeoDataFrame,
    ) -> None:
        """Test that outlet flowpath is oriented correctly using upstream."""
        graph, node_indices = sample_graph
        result = _trace_attributes(graph, node_indices, reversed_flowpaths, sample_divides, "21")

        from reference_builds.utils.geometries import _get_endpoints

        # Outlet is 85000100000005
        outlet = result[result["flowpath_id"] == "85000100000005"].iloc[0]
        upstream = result[result["flowpath_id"] == "85000100000004"].iloc[0]

        _, upstream_end = _get_endpoints(upstream.geometry)
        outlet_start, _ = _get_endpoints(outlet.geometry)

        # Outlet should start where upstream ends
        assert Point(outlet_start).distance(Point(upstream_end)) < 0.001


class TestGeoglowsFlowpathOrientation:
    """Tests for flowpath orientation in _trace_geoglows_attributes."""

    @pytest.fixture
    def reversed_geoglows_flowpaths(self) -> gpd.GeoDataFrame:
        """Create GeoGLOWS flowpaths with reversed digitization."""
        return gpd.GeoDataFrame(
            {
                "LINKNO": [810000001, 810000002, 810000003, 810000004, 810000005],
                "DSLINKNO": [810000003, 810000003, 810000004, 810000005, -1],
                "strmOrder": [1, 1, 2, 2, 2],
                "LengthKM": [0.5, 0.6, 0.8, 0.7, 0.4],
            },
            geometry=[
                LineString([(1, 1), (0, 0)]),  # Reversed: should be (0,0)->(1,1)
                LineString([(1, 1), (2, 0)]),  # Reversed: should be (2,0)->(1,1)
                LineString([(1, 2), (1, 1)]),  # Reversed: should be (1,1)->(1,2)
                LineString([(1, 3), (1, 2)]),  # Reversed: should be (1,2)->(1,3)
                LineString([(1, 4), (1, 3)]),  # Reversed: should be (1,3)->(1,4)
            ],
            crs="EPSG:3857",
        )

    def test_reversed_geoglows_flowpaths_are_corrected(
        self,
        sample_geoglows_graph: tuple[rx.PyDiGraph, dict[str, int]],
        reversed_geoglows_flowpaths: gpd.GeoDataFrame,
        sample_geoglows_catchments: gpd.GeoDataFrame,
    ) -> None:
        """Test that reversed GeoGLOWS flowpaths are oriented correctly."""
        graph, node_indices = sample_geoglows_graph
        result = _trace_geoglows_attributes(
            graph, node_indices, reversed_geoglows_flowpaths, sample_geoglows_catchments, "701"
        )

        from reference_builds.utils.geometries import _get_endpoints

        # Check that flowpath 4 ends where flowpath 5 starts
        fp4 = result[result["flowpath_id"] == "810000004"].iloc[0]
        fp5 = result[result["flowpath_id"] == "810000005"].iloc[0]

        _, fp4_end = _get_endpoints(fp4.geometry)
        fp5_start, _ = _get_endpoints(fp5.geometry)

        assert Point(fp4_end).distance(Point(fp5_start)) < 0.001
