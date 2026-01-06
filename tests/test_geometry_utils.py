"""Tests for geometry utility functions"""

from shapely import wkb
from shapely.geometry import LineString, MultiLineString

from reference_builds.utils.geometries import (
    _ensure_geometry,
    _get_endpoints,
    _orient_flowpath_downstream,
    _reverse_line,
)


class TestEnsureGeometry:
    """Tests for _ensure_geometry function."""

    def test_returns_geometry_unchanged(self) -> None:
        """Test that Shapely geometry is returned unchanged."""
        line = LineString([(0, 0), (1, 1)])
        result = _ensure_geometry(line)
        assert result == line

    def test_converts_wkb_bytes_to_geometry(self) -> None:
        """Test that WKB bytes are converted to Shapely geometry."""
        line = LineString([(0, 0), (1, 1)])
        wkb_bytes = wkb.dumps(line)
        result = _ensure_geometry(wkb_bytes)
        assert result.equals(line)

    def test_handles_multilinestring(self) -> None:
        """Test that MultiLineString is handled correctly."""
        multi = MultiLineString([[(0, 0), (1, 1)], [(2, 2), (3, 3)]])
        result = _ensure_geometry(multi)
        assert result == multi

    def test_handles_multilinestring_wkb(self) -> None:
        """Test that MultiLineString WKB is converted correctly."""
        multi = MultiLineString([[(0, 0), (1, 1)], [(2, 2), (3, 3)]])
        wkb_bytes = wkb.dumps(multi)
        result = _ensure_geometry(wkb_bytes)
        assert result.equals(multi)


class TestGetEndpoints:
    """Tests for _get_endpoints function."""

    def test_linestring_endpoints(self) -> None:
        """Test getting endpoints from LineString."""
        line = LineString([(0, 0), (1, 1), (2, 2)])
        start, end = _get_endpoints(line)
        assert start == (0, 0)
        assert end == (2, 2)

    def test_multilinestring_endpoints(self) -> None:
        """Test getting endpoints from MultiLineString."""
        multi = MultiLineString([[(0, 0), (1, 1)], [(2, 2), (3, 3)]])
        start, end = _get_endpoints(multi)
        assert start == (0, 0)
        assert end == (3, 3)

    def test_single_segment_linestring(self) -> None:
        """Test LineString with only two points."""
        line = LineString([(5, 5), (10, 10)])
        start, end = _get_endpoints(line)
        assert start == (5, 5)
        assert end == (10, 10)

    def test_handles_wkb_input(self) -> None:
        """Test that WKB input is handled."""
        line = LineString([(0, 0), (1, 1)])
        wkb_bytes = wkb.dumps(line)
        start, end = _get_endpoints(wkb_bytes)
        assert start == (0, 0)
        assert end == (1, 1)


class TestReverseLine:
    """Tests for _reverse_line function."""

    def test_reverse_linestring(self) -> None:
        """Test reversing a LineString."""
        line = LineString([(0, 0), (1, 1), (2, 2)])
        result = _reverse_line(line)
        assert list(result.coords) == [(2, 2), (1, 1), (0, 0)]

    def test_reverse_multilinestring(self) -> None:
        """Test reversing a MultiLineString."""
        multi = MultiLineString([[(0, 0), (1, 1)], [(2, 2), (3, 3)]])
        result = _reverse_line(multi)

        # Should reverse order of parts and coords within each part
        parts = list(result.geoms)
        assert len(parts) == 2
        assert list(parts[0].coords) == [(3, 3), (2, 2)]
        assert list(parts[1].coords) == [(1, 1), (0, 0)]

    def test_reverse_preserves_geometry_type(self) -> None:
        """Test that reversed geometry has same type."""
        line = LineString([(0, 0), (1, 1)])
        result = _reverse_line(line)
        assert result.geom_type == "LineString"

        multi = MultiLineString([[(0, 0), (1, 1)]])
        result = _reverse_line(multi)
        assert result.geom_type == "MultiLineString"

    def test_handles_wkb_input(self) -> None:
        """Test that WKB input is handled."""
        line = LineString([(0, 0), (1, 1)])
        wkb_bytes = wkb.dumps(line)
        result = _reverse_line(wkb_bytes)
        assert list(result.coords) == [(1, 1), (0, 0)]


class TestOrientFlowpathDownstream:
    """Tests for _orient_flowpath_downstream function."""

    def test_already_correct_orientation_with_downstream(self) -> None:
        """Test that correctly oriented line is unchanged."""
        # Line goes from (0,0) to (1,1), downstream is at (1,1) to (1,2)
        line = LineString([(0, 0), (1, 1)])
        ds_geom = LineString([(1, 1), (1, 2)])

        result = _orient_flowpath_downstream(line, ds_geom=ds_geom)

        # End point (1,1) is closer to downstream, so no change
        assert list(result.coords) == [(0, 0), (1, 1)]

    def test_reversed_orientation_with_downstream(self) -> None:
        """Test that reversed line is corrected."""
        # Line goes from (1,1) to (0,0), but downstream is at (1,1)
        line = LineString([(1, 1), (0, 0)])
        ds_geom = LineString([(1, 1), (1, 2)])

        result = _orient_flowpath_downstream(line, ds_geom=ds_geom)

        # Start point (1,1) is closer to downstream, so should reverse
        assert list(result.coords) == [(0, 0), (1, 1)]

    def test_outlet_with_upstream_geometry_correct(self) -> None:
        """Test outlet orientation using upstream geometry (already correct)."""
        # Outlet line goes from (1,1) to (1,2), upstream ends at (1,1)
        line = LineString([(1, 1), (1, 2)])
        us_geom = LineString([(0, 0), (1, 1)])

        result = _orient_flowpath_downstream(line, ds_geom=None, us_geom=us_geom)

        # Start (1,1) is closer to upstream, which is correct
        assert list(result.coords) == [(1, 1), (1, 2)]

    def test_outlet_with_upstream_geometry_reversed(self) -> None:
        """Test outlet orientation using upstream geometry (needs reversal)."""
        # Outlet line goes from (1,2) to (1,1), but upstream ends at (1,1)
        line = LineString([(1, 2), (1, 1)])
        us_geom = LineString([(0, 0), (1, 1)])

        result = _orient_flowpath_downstream(line, ds_geom=None, us_geom=us_geom)

        # End (1,1) is closer to upstream, so should reverse
        assert list(result.coords) == [(1, 1), (1, 2)]

    def test_no_reference_returns_unchanged(self) -> None:
        """Test that line is unchanged when no reference geometry provided."""
        line = LineString([(0, 0), (1, 1)])
        result = _orient_flowpath_downstream(line, ds_geom=None, us_geom=None)
        assert list(result.coords) == [(0, 0), (1, 1)]

    def test_multilinestring_orientation(self) -> None:
        """Test orientation of MultiLineString."""
        # Multi goes from (0,0) to (2,2), downstream starts at (2,2)
        multi = MultiLineString([[(0, 0), (1, 1)], [(1, 1), (2, 2)]])
        ds_geom = LineString([(2, 2), (3, 3)])

        result = _orient_flowpath_downstream(multi, ds_geom=ds_geom)

        # End (2,2) is closer to downstream, so no change
        start, end = _get_endpoints(result)
        assert start == (0, 0)
        assert end == (2, 2)

    def test_multilinestring_reversed(self) -> None:
        """Test reversal of MultiLineString."""
        # Multi goes from (2,2) to (0,0), but downstream is at (2,2)
        multi = MultiLineString([[(2, 2), (1, 1)], [(1, 1), (0, 0)]])
        ds_geom = LineString([(2, 2), (3, 3)])

        result = _orient_flowpath_downstream(multi, ds_geom=ds_geom)

        # Start (2,2) is closer to downstream, so should reverse
        start, end = _get_endpoints(result)
        assert start == (0, 0)
        assert end == (2, 2)

    def test_handles_wkb_input(self) -> None:
        """Test that WKB input is handled for both geometries."""
        line = LineString([(0, 0), (1, 1)])
        ds_geom = LineString([(1, 1), (1, 2)])

        line_wkb = wkb.dumps(line)
        ds_wkb = wkb.dumps(ds_geom)

        result = _orient_flowpath_downstream(line_wkb, ds_geom=ds_wkb)
        assert list(result.coords) == [(0, 0), (1, 1)]

    def test_downstream_takes_priority_over_upstream(self) -> None:
        """Test that downstream geometry is used when both are provided."""
        line = LineString([(1, 1), (0, 0)])
        ds_geom = LineString([(1, 1), (1, 2)])  # Would cause reversal
        us_geom = LineString([(0, 0), (0, -1)])  # Would not cause reversal

        result = _orient_flowpath_downstream(line, ds_geom=ds_geom, us_geom=us_geom)

        # Should use downstream, so should reverse
        assert list(result.coords) == [(0, 0), (1, 1)]


class TestOrientFlowpathDownstreamEdgeCases:
    """Edge case tests for _orient_flowpath_downstream."""

    def test_identical_endpoints(self) -> None:
        """Test handling when both endpoints are equidistant from downstream."""
        # This is a degenerate case - line perpendicular to downstream
        line = LineString([(0, 0), (2, 0)])
        ds_geom = LineString([(1, 0), (1, 1)])

        # Should return unchanged (or reversed, but consistently)
        result = _orient_flowpath_downstream(line, ds_geom=ds_geom)
        assert result.geom_type == "LineString"

    def test_touching_downstream(self) -> None:
        """Test when line endpoint touches downstream geometry."""
        line = LineString([(0, 0), (1, 1)])
        ds_geom = LineString([(1, 1), (2, 2)])  # Starts exactly at line end

        result = _orient_flowpath_downstream(line, ds_geom=ds_geom)

        # End touches downstream (distance 0), so should not reverse
        assert list(result.coords) == [(0, 0), (1, 1)]

    def test_touching_upstream(self) -> None:
        """Test outlet when line start touches upstream geometry."""
        line = LineString([(1, 1), (2, 2)])
        us_geom = LineString([(0, 0), (1, 1)])  # Ends exactly at line start

        result = _orient_flowpath_downstream(line, ds_geom=None, us_geom=us_geom)

        # Start touches upstream (distance 0), so should not reverse
        assert list(result.coords) == [(1, 1), (2, 2)]

    def test_long_multilinestring(self) -> None:
        """Test with MultiLineString with many segments."""
        multi = MultiLineString(
            [
                [(0, 0), (1, 0)],
                [(1, 0), (2, 0)],
                [(2, 0), (3, 0)],
                [(3, 0), (4, 0)],
            ]
        )
        ds_geom = LineString([(4, 0), (5, 0)])

        result = _orient_flowpath_downstream(multi, ds_geom=ds_geom)

        start, end = _get_endpoints(result)
        assert start == (0, 0)
        assert end == (4, 0)
