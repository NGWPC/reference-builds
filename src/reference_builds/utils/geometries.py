"""A file for all geometry related internal functions"""

import logging

import geopandas as gpd
import pandas as pd
from shapely import Geometry, wkb
from shapely.geometry import LineString, MultiLineString, MultiPolygon, Point

logger = logging.getLogger(__name__)


def _ensure_geometry(geom):  # type: ignore[no-untyped-def]
    """Convert bytes to Shapely geometry if needed."""
    if isinstance(geom, bytes):
        return wkb.loads(geom)
    return geom


def _get_endpoints(geom):  # type: ignore[no-untyped-def]
    """Get start and end points of a line or multiline geometry."""
    geom = _ensure_geometry(geom)
    if geom.geom_type == "MultiLineString":
        start_coord = list(geom.geoms)[0].coords[0]
        end_coord = list(geom.geoms)[-1].coords[-1]
    else:
        start_coord = geom.coords[0]
        end_coord = geom.coords[-1]
    return start_coord, end_coord


def _reverse_line(geom):  # type: ignore[no-untyped-def]
    """Reverse a LineString or MultiLineString."""
    geom = _ensure_geometry(geom)
    if geom.geom_type == "MultiLineString":
        # Reverse each component and reverse the order of components
        reversed_parts = [LineString(part.coords[::-1]) for part in reversed(geom.geoms)]
        return MultiLineString(reversed_parts)
    else:
        return LineString(geom.coords[::-1])


def _orient_flowpath_downstream(geom, ds_geom=None, us_geom=None):  # type: ignore[no-untyped-def]
    """
    Orient a flowpath so coords go from upstream to downstream.

    Parameters
    ----------
    geom : geometry
        The flowpath geometry
    ds_geom : geometry, optional
        The downstream flowpath geometry for direction detection
    us_geom : geometry, optional
        The upstream flowpath geometry (used for outlets when ds_geom is None)

    Returns
    -------
    geometry
        The oriented geometry (start = upstream, end = downstream)
    """
    geom = _ensure_geometry(geom)
    start_coord, end_coord = _get_endpoints(geom)
    start_pt = Point(start_coord)
    end_pt = Point(end_coord)

    # Primary: use downstream geometry
    if ds_geom is not None:
        ds_geom = _ensure_geometry(ds_geom)
        dist_start = start_pt.distance(ds_geom)
        dist_end = end_pt.distance(ds_geom)

        # If start is closer to downstream, the line is reversed
        if dist_start < dist_end:
            return _reverse_line(geom)
        return geom

    # Fallback for outlets: use upstream geometry
    if us_geom is not None:
        us_geom = _ensure_geometry(us_geom)
        dist_start = start_pt.distance(us_geom)
        dist_end = end_pt.distance(us_geom)

        # Upstream end should be CLOSER to upstream geometry
        # So if end is closer to upstream, line is reversed
        if dist_end < dist_start:
            return _reverse_line(geom)
        return geom

    # No reference at all - return as-is
    return geom


def _drop_exclaves(geom: Geometry) -> Geometry:
    """Find and destroy non-contiguous parts of MultiPolygons"""
    if geom.geom_type != "MultiPolygon":
        return geom
    main_part = geom.geoms[0]
    for part in geom.geoms:
        if part.area > main_part.area:
            main_part = part
    main_part_buffered = main_part.buffer(1.0)
    for part in geom.geoms:
        if part.intersects(main_part_buffered):
            main_part = main_part.union(part)
    return main_part


def _find_exclaves(geom: Geometry) -> pd.Series:
    """Find and exclude non-contiguous parts of MultiPolygons, appending them to a list to be resolved later"""
    exclaves = []
    if geom.geom_type != "MultiPolygon":
        return pd.Series(data={"geometry": geom, "exclaves": exclaves}, index=["geometry", "exclaves"])

    main_part = geom.geoms[0]
    for part in geom.geoms:
        if part.area > main_part.area:
            main_part = part
    included_parts = [main_part]
    main_part_buffered = main_part.buffer(1.0)
    for part in geom.geoms:
        if part.intersects(main_part_buffered):
            included_parts.append(part)
            main_part = MultiPolygon(included_parts)
        else:
            exclaves.append(part)
    return pd.Series(data={"geometry": main_part, "exclaves": exclaves}, index=["geometry", "exclaves"])


def _fix_divide_exclaves(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Remove exclaves from catchment geometries, merging them in with the neighbor it intersects the most with"""
    exclaves = gdf["geometry"].apply(_find_exclaves)

    fixed_count = 0

    gdf["geometry"] = exclaves["geometry"]
    for idx, e in exclaves["exclaves"].items():
        for exclave in e:
            buffer = exclave.buffer(1.0)
            intersection_areas = gdf.intersection(buffer).area
            intersection_areas[idx] = 0.0
            best_idx = intersection_areas.argmax()
            gdf.loc[best_idx, "geometry"] = gdf["geometry"][best_idx].union(exclave)
            fixed_count += 1

    logger.info(f"fix_divide_exclaves: fixed/merged {fixed_count} polygons")
    return gdf


def _validate_and_fix_geometries(gdf: gpd.GeoDataFrame, geom_type: str) -> gpd.GeoDataFrame:
    """Validate and fix invalid geometries in a GeoDataFrame.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoDataFrame to validate
    geom_type : str
        Description for logging (e.g., "flowpaths", "divides")

    Returns
    -------
    gpd.GeoDataFrame
        GeoDataFrame with fixed geometries

    Raises
    ------
    ValueError
        If geometries cannot be fixed or invalid geometries remain
    """
    invalid_mask = ~gdf.geometry.is_valid
    invalid_count = invalid_mask.sum()

    if invalid_count == 0:
        return gdf  # No invalid geometries

    geometries = gdf[invalid_mask].geometry
    gdf.loc[invalid_mask, "geometry"] = geometries.make_valid(method="structure")

    if len(gdf[~gdf.geometry.is_valid]) > 0:
        raise ValueError(f"Could not fix invalid geometries in {geom_type}")

    still_invalid = (~gdf.geometry.is_valid).sum()
    if still_invalid > 0:
        raise ValueError(f"Invalid Geometries remain: {gdf[~gdf.geometry.is_valid]}")

    return gdf
