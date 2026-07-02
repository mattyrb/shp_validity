#!/usr/bin/env python3
"""
shp_validity.py - Shapefile geometry validator and cleaner

Validates and repairs polygon geometries in an ESRI shapefile, separating
valid features from rejected ones. Designed for cleaning vector products
derived from Landsat imagery (e.g., classification polygons, water masks,
field delineations) into flat, valid polygons fit for spatial analysis and
modeling of evapotranspiration and applied irrigation water.

USAGE
-----
Run from a directory containing one or more .shp files. The script lists
shapefiles in the current working directory and prompts you to pick one.

OUTPUTS (written to a 'cleaned/' subdirectory next to the input)
---------------------------------------------------------------
    <name>_valid.shp            Valid polygon / multipolygon features.
                                Optionally includes boolean 'repaired' and
                                'sliver' flags and overlap fields -- decline
                                the flag-fields prompt to keep the output
                                attributes identical to the input.
    <name>_rejected_poly.shp    Unrepairable polygons (only if any).
    <name>_rejected_line.shp    Filtered (multi)line features (only if any).
    <name>_rejected_point.shp   Filtered (multi)point features (only if any).
    <name>_slivers.shp          Thin sliver features for review (only if
                                sliver flagging is on and any are found).
    <name>_overlaps.csv         Overlapping feature pairs with overlap area
                                (only if overlap checking is on).
    <name>_overlap_zones.shp    The overlap geometries themselves, for map
                                review (only if overlap checking is on).
    <name>_report.csv           Per-feature audit log with row index, the
                                user-chosen identifier column, the original
                                geometry type, validity, repair status,
                                optional sliver / overlap metrics, and
                                Shapely's validity diagnostic message.

REPAIR PIPELINE
---------------
    1. make_valid()  - Shapely's modern geometry repair.
    2. Extract polygon parts from any resulting GeometryCollection.
    3. buffer(0) fallback for any remaining cases.
    4. Reject if still bad.

QUALITY CHECKS FOR ET / APPLIED-WATER WORK
------------------------------------------
    - CRS: warns if the layer is geographic (degrees), since area, sliver
      width, and overlap area are then not meaningful for water calculations.
    - Sliver flagging: flags thin features by characteristic width
      (2 * area / perimeter), so artifacts that are thin but still cover a
      meaningful area are caught. A Polsby-Popper thinness ratio is also
      reported. Long but legitimately wide fields are not flagged.
    - Overlap flagging: reports overlapping polygons (which double-count ET
      and applied water), with per-feature overlap area, a pair list, and the
      overlap geometries for review.

ROBUSTNESS HANDLING
-------------------
    - Warns and prompts before continuing if the input has no .prj / CRS.
    - Z and M dimensions are stripped automatically (with a note).
    - If the input already has a flag column, the flag is written under a
      non-conflicting name.
    - Falls back to CP1252 encoding if the default read fails.
    - All output geometries are normalized to multipart.

REQUIREMENTS
------------
    geopandas, shapely >= 2.0

Author: written for Matt Bromley (DRI), GIS / remote sensing workflows.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import sys
import time
from collections import Counter
from pathlib import Path

import geopandas as gpd
import shapely
from shapely.geometry import (
    GeometryCollection,
    LineString,
    MultiLineString,
    MultiPoint,
    MultiPolygon,
    Point,
    Polygon,
)
from shapely.geometry.base import BaseGeometry
from shapely.validation import explain_validity, make_valid


def _is_geometry(geom) -> bool:
    """True if ``geom`` is a real shapely geometry object."""
    return isinstance(geom, BaseGeometry)


def _to_multipart(geom):
    """Promote singleton geometries to their Multi* equivalents."""
    if not _is_geometry(geom):
        return geom
    if isinstance(geom, Polygon):
        return MultiPolygon([geom])
    if isinstance(geom, LineString):
        return MultiLineString([geom])
    if isinstance(geom, Point):
        return MultiPoint([geom])
    return geom


def read_shapefile_robust(shp_path: Path):
    """Read a shapefile, falling back to CP1252 encoding if the default fails."""
    try:
        return gpd.read_file(shp_path)
    except UnicodeDecodeError:
        pass
    except Exception as e:
        msg = str(e).lower()
        if "decode" not in msg and "encoding" not in msg and "utf" not in msg:
            raise

    print("  Note: default encoding failed; retrying with CP1252.")
    return gpd.read_file(shp_path, encoding="cp1252")


def strip_z_dimensions(gdf):
    """Strip Z (and any M) dimensions from all geometries, warning if any."""
    z_count = sum(1 for g in gdf.geometry if _is_geometry(g) and g.has_z)
    if z_count > 0:
        print(f"  Note: {z_count:,} geometr{'y' if z_count == 1 else 'ies'} "
              "had Z/M dimensions -- stripping to 2D.")
        gdf = gdf.copy()
        gdf["geometry"] = gdf.geometry.apply(
            lambda g: shapely.force_2d(g) if _is_geometry(g) else g
        )
    return gdf


def pick_unique_field_name(gdf, preferred: str, alternatives=None) -> str:
    """Return a field name that does not already exist in ``gdf.columns``."""
    if alternatives is None:
        alternatives = ["was_fixed", "fixed", "is_fixed", "repair_flg"]

    if preferred not in gdf.columns:
        return preferred

    for alt in alternatives:
        if alt not in gdf.columns:
            return alt

    for i in range(2, 100):
        candidate = f"{preferred}{i}"[:10]
        if candidate not in gdf.columns:
            return candidate

    raise RuntimeError("Could not find an unused field name for a flag column.")


# Geometry-type buckets ------------------------------------------------------
POLY_TYPES = ("Polygon", "MultiPolygon")
LINE_TYPES = ("LineString", "MultiLineString")
POINT_TYPES = ("Point", "MultiPoint")


# Status tags written in full to the CSV audit report.
REPAIRED_STATUSES = {
    "repaired_make_valid",
    "repaired_collection_extracted",
    "repaired_buffer0",
}


# ---------------------------------------------------------------------------
# Interactive prompts
# ---------------------------------------------------------------------------

def list_shapefiles(directory: Path) -> list[Path]:
    """Return a sorted list of .shp files in the given directory."""
    return sorted(directory.glob("*.shp"))


def prompt_shapefile_choice(shapefiles: list[Path]) -> Path:
    """Prompt the user to pick a shapefile from the list."""
    print("\nShapefiles available in this directory:")
    for i, shp in enumerate(shapefiles, 1):
        print(f"  {i}. {shp.name}")

    while True:
        choice = input(f"\nSelect a shapefile [1-{len(shapefiles)}]: ").strip()
        try:
            idx = int(choice) - 1
            if 0 <= idx < len(shapefiles):
                return shapefiles[idx]
        except ValueError:
            pass
        print(f"  Please enter a number between 1 and {len(shapefiles)}.")


def prompt_yes_no(question: str, default: bool = True) -> bool:
    """Prompt for yes/no input, returning a bool."""
    suffix = " [Y/n]: " if default else " [y/N]: "
    while True:
        ans = input(question + suffix).strip().lower()
        if ans == "":
            return default
        if ans in ("y", "yes"):
            return True
        if ans in ("n", "no"):
            return False
        print("  Please answer 'y' or 'n'.")


def prompt_float(question: str, default: float) -> float:
    """Prompt for a numeric value, returning a float. Enter accepts default."""
    while True:
        ans = input(f"{question} [{default}]: ").strip()
        if ans == "":
            return float(default)
        try:
            return float(ans)
        except ValueError:
            print("  Please enter a number (or press Enter for the default).")


def prompt_id_field(gdf) -> str | None:
    """Prompt the user to choose an attribute column as a feature identifier."""
    geom_col = gdf.geometry.name
    candidates = [c for c in gdf.columns if c != geom_col]

    if not candidates:
        print("\n  No attribute columns found -- using row index as identifier.")
        return None

    print("\nAttribute columns available for use as an identifier:")
    for i, col in enumerate(candidates, 1):
        sample = gdf[col].iloc[0] if len(gdf) else ""
        sample_str = str(sample)
        if len(sample_str) > 30:
            sample_str = sample_str[:27] + "..."
        print(f"  {i:>2}. {col:<15s}  e.g. {sample_str}")

    while True:
        choice = input(
            f"\nSelect an identifier column [1-{len(candidates)}], "
            "or press Enter to use row index only: "
        ).strip()
        if choice == "":
            return None
        try:
            idx = int(choice) - 1
            if 0 <= idx < len(candidates):
                return candidates[idx]
        except ValueError:
            pass
        print(f"  Please enter a number between 1 and {len(candidates)}, "
              "or press Enter to skip.")


# ---------------------------------------------------------------------------
# Geometry handling
# ---------------------------------------------------------------------------

def extract_polygons(geom):
    """Pull (Multi)Polygon parts out of a geometry, discarding line / point
    fragments. Returns a Polygon, MultiPolygon, or None."""
    if geom is None or geom.is_empty:
        return None

    if isinstance(geom, Polygon):
        return geom if not geom.is_empty else None

    if isinstance(geom, MultiPolygon):
        parts = [p for p in geom.geoms if not p.is_empty]
        if not parts:
            return None
        return parts[0] if len(parts) == 1 else MultiPolygon(parts)

    if isinstance(geom, GeometryCollection):
        polys = []
        for g in geom.geoms:
            if isinstance(g, Polygon) and not g.is_empty:
                polys.append(g)
            elif isinstance(g, MultiPolygon):
                polys.extend(p for p in g.geoms if not p.is_empty)
        if not polys:
            return None
        return polys[0] if len(polys) == 1 else MultiPolygon(polys)

    return None


def repair_geometry(geom):
    """Try to coerce ``geom`` into a valid (Multi)Polygon.

    Returns ``(repaired_geom_or_None, status_tag)``.
    """
    if not _is_geometry(geom):
        return None, "null_geom"
    if geom.is_empty:
        return None, "empty_geom"

    geom_type = geom.geom_type

    if geom_type in LINE_TYPES or geom_type in POINT_TYPES:
        return None, f"wrong_type_{geom_type.lower()}"
    if geom_type not in POLY_TYPES and geom_type != "GeometryCollection":
        return None, f"wrong_type_{geom_type.lower()}"

    if geom_type in POLY_TYPES and geom.is_valid:
        return geom, "valid_original"

    try:
        repaired = make_valid(geom)
        polys = extract_polygons(repaired)
        if polys is not None and polys.is_valid and not polys.is_empty:
            tag = (
                "repaired_collection_extracted"
                if isinstance(repaired, GeometryCollection)
                else "repaired_make_valid"
            )
            return polys, tag
    except Exception:
        pass

    try:
        buffered = geom.buffer(0)
        polys = extract_polygons(buffered)
        if polys is not None and polys.is_valid and not polys.is_empty:
            return polys, "repaired_buffer0"
    except Exception:
        pass

    return None, "unrepairable"


def geometry_metrics(geom):
    """Return ``(area, perimeter, width, thinness)`` for a polygonal geometry.

    width    = 2 * area / perimeter, a mean-width estimate. For a long thin
               strip this approaches the true width regardless of length, so
               it flags slivers that still cover a meaningful area. A
               legitimately wide-but-long field is not flagged, which a pure
               area threshold cannot distinguish.
    thinness = 4 * pi * area / perimeter**2 (Polsby-Popper compactness): 1.0
               for a circle, approaching 0 for thin or ragged shapes. Reported
               for reference so you can apply your own cutoff if you prefer it.
    """
    if not _is_geometry(geom) or geom.is_empty:
        return 0.0, 0.0, 0.0, 0.0
    area = float(geom.area)
    perim = float(geom.length)
    if perim <= 0:
        return area, perim, 0.0, 0.0
    width = 2.0 * area / perim
    thinness = 4.0 * math.pi * area / (perim * perim)
    return area, perim, width, thinness


def find_overlaps(gdf, min_area: float = 0.0):
    """Find overlapping polygons within ``gdf`` using its spatial index.

    Returns ``(overlap_area, overlap_count, pairs)`` where overlap_area and
    overlap_count are dicts keyed by the gdf index, and pairs is a list of
    ``(label_a, label_b, area, intersection_geom)`` for each overlapping pair
    (a < b by position). Pure edge touches (zero-area intersections) are
    ignored; ``min_area`` raises the threshold further if you want it.
    """
    geoms = gdf.geometry
    labels = list(gdf.index)
    sindex = gdf.sindex

    overlap_area = {lab: 0.0 for lab in labels}
    overlap_count = {lab: 0 for lab in labels}
    pairs = []
    touch_floor = 1e-9
    floor = max(float(min_area), touch_floor)

    for pos_i, lab_i in enumerate(labels):
        g = geoms.iloc[pos_i]
        if not _is_geometry(g) or g.is_empty:
            continue
        for pos_j in sindex.query(g, predicate="intersects"):
            if pos_j <= pos_i:
                continue
            h = geoms.iloc[pos_j]
            if not _is_geometry(h) or h.is_empty:
                continue
            try:
                inter = g.intersection(h)
            except Exception:
                continue
            a = float(inter.area)
            if a >= floor:
                lab_j = labels[pos_j]
                overlap_area[lab_i] += a
                overlap_area[lab_j] += a
                overlap_count[lab_i] += 1
                overlap_count[lab_j] += 1
                pairs.append((lab_i, lab_j, a, inter))

    return overlap_area, overlap_count, pairs


# Parts smaller than this (CRS units squared) are treated as degenerate
# "orphaned vertex" debris and dropped regardless of the user threshold.
DEGENERATE_AREA_EPS = 1e-6


def explode_polygons(geom):
    """Return a flat list of Polygon parts from any polygonal geometry."""
    if not _is_geometry(geom) or geom.is_empty:
        return []
    if isinstance(geom, Polygon):
        return [geom]
    if isinstance(geom, MultiPolygon):
        return [g for g in geom.geoms if not g.is_empty]
    if isinstance(geom, GeometryCollection):
        out = []
        for g in geom.geoms:
            out.extend(explode_polygons(g))
        return out
    return []


def clean_multipart_parts(geom, min_part_area: float = 0.0):
    """Explode a (multi)polygon, repair each part, and drop degenerate or
    undersized parts. Returns ``(cleaned_geom_or_None, dropped)`` where
    dropped is a list of ``(part_geom, area)`` for each removed part.

    A part is dropped if it cannot be repaired to a valid polygon (the
    "couple of orphaned vertices" case), or if its area is below the
    degenerate epsilon or the user's ``min_part_area``. Surviving parts are
    reassembled into a Polygon (one part) or MultiPolygon (several). This
    automates the manual explode / identify-sliver / erase workflow.
    """
    parts = explode_polygons(geom)
    if not parts:
        return None, []

    threshold = max(DEGENERATE_AREA_EPS, float(min_part_area))
    keep = []
    dropped = []
    for part in parts:
        repaired = part
        if not part.is_valid:
            try:
                repaired = extract_polygons(make_valid(part))
            except Exception:
                repaired = None
        if repaired is None or repaired.is_empty:
            dropped.append((part, float(part.area)))
            continue
        if repaired.area < threshold:
            dropped.append((repaired, float(repaired.area)))
            continue
        keep.append(repaired)

    if not keep:
        return None, dropped

    flat = []
    for k in keep:
        if isinstance(k, MultiPolygon):
            flat.extend(k.geoms)
        else:
            flat.append(k)
    cleaned = flat[0] if len(flat) == 1 else MultiPolygon(flat)
    return cleaned, dropped


# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------

def _bucket_for_geom_type(geom_type: str) -> str | None:
    """Map an original geometry type to a rejected-file bucket, or None to skip."""
    if geom_type in POLY_TYPES:
        return "poly"
    if geom_type in LINE_TYPES:
        return "line"
    if geom_type in POINT_TYPES:
        return "point"
    return None


def write_rejected_files(rejected_records, cleaned_dir: Path, base: str, crs):
    """Split rejected features by geometry type and write each as its own .shp."""
    buckets: dict[str, list] = {"poly": [], "line": [], "point": []}

    for rec, status, geom_type in rejected_records:
        bucket = _bucket_for_geom_type(geom_type)
        if bucket is None:
            continue
        buckets[bucket].append(dict(rec))

    written = []
    for bucket, rows in buckets.items():
        if not rows:
            continue
        gdf = gpd.GeoDataFrame(rows, crs=crs)
        gdf["geometry"] = gdf.geometry.apply(_to_multipart)
        out_path = cleaned_dir / f"{base}_rejected_{bucket}.shp"
        gdf.to_file(out_path)
        written.append((out_path, len(gdf)))
    return written


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------

def process_shapefile(shp_path: Path, write_report: bool = True,
                      add_repaired_field: bool = True,
                      flag_slivers: bool = False, sliver_width: float = 0.0,
                      flag_overlaps: bool = False, min_overlap_area: float = 0.0,
                      clean_parts: bool = False, min_part_area: float = 0.0,
                      id_field="__prompt__", interactive: bool = True,
                      out_dir: Path = None, in_place: bool = False,
                      backup_dir: Path = None):
    """Run the full validate / repair / write pipeline on one shapefile.

    If ``add_repaired_field`` is False, the valid output keeps the input
    attributes exactly as they were -- no flag fields are added. Repair,
    sliver, and overlap information is still recorded in the CSV report and
    the side files, so nothing is lost; it just lives outside the shapefile.
    """
    print(f"\nReading: {shp_path.name}")
    try:
        gdf = read_shapefile_robust(shp_path)
    except Exception as e:
        print(f"  ERROR reading shapefile: {e}")
        if interactive:
            sys.exit(1)
        raise

    print(f"  Features: {len(gdf):,}")
    print(f"  CRS:      {gdf.crs}")

    # --- CRS classification + sanity checks --------------------------------
    crs_kind = "none"
    crs_units = "unknown"
    if gdf.crs is not None:
        crs_kind = "geographic" if gdf.crs.is_geographic else "projected"
        try:
            crs_units = gdf.crs.axis_info[0].unit_name
        except Exception:
            crs_units = "unknown"

    if gdf.crs is None:
        print("\n  WARNING: this shapefile has no CRS (the .prj file is")
        print("  missing or unreadable). Output shapefiles will also lack")
        print("  a .prj, and downstream tools may misinterpret coordinates.")
        if interactive and not prompt_yes_no("  Proceed anyway?", default=False):
            print("  Aborting.")
            sys.exit(1)
    elif gdf.crs.is_geographic:
        print(f"\n  WARNING: this layer is in a geographic CRS (units: {crs_units}).")
        print("  Area, sliver width, and overlap area are computed in those")
        print("  units and are not meaningful for ET or applied-water work.")
        print("  Reproject to a projected CRS (for example UTM) before relying")
        print("  on any area-based output.")
        if (interactive and (flag_slivers or flag_overlaps)
                and not prompt_yes_no(
                    "  Continue with area-based checks anyway?", default=False)):
            print("  Aborting.")
            sys.exit(1)

    # --- Strip Z / M dimensions -------------------------------------------
    gdf = strip_z_dimensions(gdf)

    # --- Pick non-colliding names for any flag fields ----------------------
    # Only relevant if we are tagging the output. When flag fields are off,
    # the valid output attributes are left identical to the input.
    repaired_field = None
    sliver_field = None
    ov_area_field = None
    ov_n_field = None
    parts_drop_field = None
    if add_repaired_field:
        repaired_field = pick_unique_field_name(gdf, "repaired")
        if repaired_field != "repaired":
            print(f"  Note: input already has a 'repaired' column; "
                  f"using '{repaired_field}' for the repair flag instead.")
        if clean_parts:
            parts_drop_field = pick_unique_field_name(
                gdf, "parts_drop", ["pts_drop", "ndropped"])
        if flag_slivers:
            sliver_field = pick_unique_field_name(
                gdf, "sliver", ["is_sliver", "sliver_flg", "slivr"])
        if flag_overlaps:
            ov_area_field = pick_unique_field_name(
                gdf, "ov_area", ["ovlap_area", "ov_area2"])
            ov_n_field = pick_unique_field_name(
                gdf, "ov_n", ["ov_count", "ovlap_n"])

    # Determine the feature-identifier column for the report.
    if id_field == "__prompt__":
        id_field = prompt_id_field(gdf) if interactive else None
    elif id_field is not None:
        if id_field == gdf.geometry.name or id_field not in gdf.columns:
            print(f"  Note: identifier column '{id_field}' not found; "
                  "using row index.")
            id_field = None

    valid_rows = []
    rejected_records = []
    sliver_rows = []
    removed_part_geoms = []
    removed_part_attrs = []
    n_part_features = 0
    dropped_area_total = 0.0
    report_rows = []
    report_by_idx = {}

    for row_idx, row in gdf.iterrows():
        geom = row.geometry
        has_geom = _is_geometry(geom)
        original_type = geom.geom_type if has_geom else "None"
        was_valid = bool(geom.is_valid) if has_geom else False

        validity_msg = ""
        if has_geom and not was_valid:
            try:
                validity_msg = explain_validity(geom)
            except Exception as e:
                validity_msg = f"explain_validity error: {e}"

        repaired_geom, status = repair_geometry(geom)
        was_repaired = status in REPAIRED_STATUSES

        # --- Clean degenerate / sliver parts out of multiparts -------------
        parts_dropped = 0
        part_drop_area = 0.0
        if clean_parts and repaired_geom is not None:
            cleaned, dropped = clean_multipart_parts(repaired_geom, min_part_area)
            if dropped:
                parts_dropped = len(dropped)
                part_drop_area = sum(a for _g, a in dropped)
                for dg, da in dropped:
                    removed_part_geoms.append(_to_multipart(dg))
                    removed_part_attrs.append(
                        {"parent_idx": row_idx, "part_area": round(da, 6)})
                was_repaired = True
                n_part_features += 1
                dropped_area_total += part_drop_area
            repaired_geom = cleaned
            if repaired_geom is None:
                status = "parts_all_dropped"

        report_row = {"row_idx": row_idx}
        if id_field is not None:
            report_row[id_field] = row[id_field]
        report_row["original_type"] = original_type
        report_row["was_valid"] = was_valid
        report_row["status"] = status
        report_row["repaired"] = was_repaired
        if clean_parts:
            report_row["parts_drop"] = parts_dropped
            report_row["drop_area"] = round(part_drop_area, 4)

        # --- Sliver metrics on the output geometry -------------------------
        sliver_flag = False
        width = thinness = 0.0
        if flag_slivers:
            if repaired_geom is not None:
                _a, _p, width, thinness = geometry_metrics(repaired_geom)
                sliver_flag = width > 0.0 and width < sliver_width
                report_row["width"] = round(width, 4)
                report_row["thinness"] = round(thinness, 5)
                report_row["sliver"] = sliver_flag
            else:
                report_row["width"] = ""
                report_row["thinness"] = ""
                report_row["sliver"] = ""

        if flag_overlaps:
            report_row["ov_area"] = 0.0 if repaired_geom is not None else ""
            report_row["ov_n"] = 0 if repaired_geom is not None else ""

        report_row["validity_msg"] = validity_msg
        report_rows.append(report_row)
        report_by_idx[row_idx] = report_row

        if repaired_geom is not None:
            new_row = row.copy()
            new_row.geometry = repaired_geom
            if add_repaired_field and repaired_field is not None:
                new_row[repaired_field] = was_repaired
            if add_repaired_field and sliver_field is not None:
                new_row[sliver_field] = sliver_flag
            if add_repaired_field and parts_drop_field is not None:
                new_row[parts_drop_field] = parts_dropped
            valid_rows.append(new_row)
            if flag_slivers and sliver_flag:
                srow = row.copy()
                srow.geometry = repaired_geom
                srow["width"] = round(width, 4)
                srow["thinness"] = round(thinness, 5)
                sliver_rows.append(srow)
        else:
            rejected_records.append((row.to_dict(), status, original_type))

    # Output directory for reports / rejected / side files.
    cleaned_dir = out_dir if out_dir is not None else (shp_path.parent / "cleaned")
    cleaned_dir.mkdir(parents=True, exist_ok=True)
    base = shp_path.stem

    n_slivers = len(sliver_rows) if flag_slivers else None
    n_overlapping = None
    total_overlap_area = 0.0
    overlap_pairs = []

    # --- Valid output ------------------------------------------------------
    if valid_rows:
        valid_gdf = gpd.GeoDataFrame(valid_rows, crs=gdf.crs)
        valid_gdf["geometry"] = valid_gdf.geometry.apply(_to_multipart)

        # --- Overlap detection among valid features ------------------------
        if flag_overlaps:
            overlap_area, overlap_count, overlap_pairs = find_overlaps(
                valid_gdf, min_area=min_overlap_area)
            n_overlapping = sum(1 for v in overlap_count.values() if v > 0)
            total_overlap_area = sum(p[2] for p in overlap_pairs)
            for idx, a in overlap_area.items():
                rr = report_by_idx.get(idx)
                if rr is not None:
                    rr["ov_area"] = round(a, 4)
                    rr["ov_n"] = overlap_count.get(idx, 0)
            if add_repaired_field and ov_area_field is not None:
                valid_gdf[ov_area_field] = [
                    round(overlap_area.get(i, 0.0), 4) for i in valid_gdf.index]
                valid_gdf[ov_n_field] = [
                    overlap_count.get(i, 0) for i in valid_gdf.index]

        if in_place:
            if backup_dir is not None:
                n_bk = backup_shapefile(shp_path, backup_dir)
                print(f"\n  Backed up original ({n_bk} file(s)) to:")
                print(f"    {backup_dir}")
            _overwrite_in_place(valid_gdf, shp_path)
            print(f"  Fixed in place ({len(valid_gdf):,} valid features): "
                  f"{shp_path.name}")
        else:
            valid_path = cleaned_dir / f"{base}_valid.shp"
            valid_gdf.to_file(valid_path)
            print(f"\n  Wrote {len(valid_gdf):,} valid features:")
            print(f"    {valid_path}")
    else:
        if in_place:
            print("\n  No valid features -- left the original file untouched.")
        else:
            print("\n  No valid features to write.")

    # --- Rejected output (split by geometry type) --------------------------
    if rejected_records:
        written = write_rejected_files(rejected_records, cleaned_dir, base, gdf.crs)
        if written:
            print(f"  Wrote rejected features to:")
            for path, n in written:
                print(f"    {path.name}  ({n:,} features)")
        spatial_count = sum(n for _, n in written)
        unspatial = len(rejected_records) - spatial_count
        if unspatial:
            print(f"    ({unspatial} rejected feature(s) had no usable geometry "
                  f"-- see CSV report)")

    # --- Sliver output -----------------------------------------------------
    if flag_slivers and sliver_rows:
        sliver_gdf = gpd.GeoDataFrame(sliver_rows, crs=gdf.crs)
        sliver_gdf["geometry"] = sliver_gdf.geometry.apply(_to_multipart)
        sliver_path = cleaned_dir / f"{base}_slivers.shp"
        sliver_gdf.to_file(sliver_path)
        print(f"  Wrote {len(sliver_gdf):,} thin sliver feature(s) for review:")
        print(f"    {sliver_path.name}")

    # --- Removed multipart parts ------------------------------------------
    if clean_parts and removed_part_geoms:
        removed_gdf = gpd.GeoDataFrame(
            removed_part_attrs, geometry=removed_part_geoms, crs=gdf.crs)
        removed_path = cleaned_dir / f"{base}_removed_parts.shp"
        removed_gdf.to_file(removed_path)
        print(f"  Removed {len(removed_gdf):,} degenerate/sliver part(s) "
              f"from {n_part_features:,} multipart feature(s):")
        print(f"    {removed_path.name}")

    # --- Overlap output ----------------------------------------------------
    if flag_overlaps and overlap_pairs:
        id_lookup = gdf[id_field].to_dict() if id_field is not None else None
        ov_csv = cleaned_dir / f"{base}_overlaps.csv"
        fld = ["row_idx_a"]
        if id_field is not None:
            fld.append(f"{id_field}_a")
        fld.append("row_idx_b")
        if id_field is not None:
            fld.append(f"{id_field}_b")
        fld.append("overlap_area")
        with open(ov_csv, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fld)
            writer.writeheader()
            for a_idx, b_idx, area, _g in overlap_pairs:
                rowd = {"row_idx_a": a_idx, "row_idx_b": b_idx,
                        "overlap_area": round(area, 4)}
                if id_field is not None:
                    rowd[f"{id_field}_a"] = id_lookup.get(a_idx)
                    rowd[f"{id_field}_b"] = id_lookup.get(b_idx)
                writer.writerow(rowd)
        print(f"  Wrote {len(overlap_pairs):,} overlapping pair(s):")
        print(f"    {ov_csv.name}")

        zone_geoms = []
        zone_attr = []
        for a_idx, b_idx, area, g in overlap_pairs:
            poly = extract_polygons(g)
            if poly is not None and not poly.is_empty:
                zone_geoms.append(_to_multipart(poly))
                zone_attr.append({"row_a": a_idx, "row_b": b_idx,
                                  "ov_area": round(area, 4)})
        if zone_geoms:
            zones = gpd.GeoDataFrame(zone_attr, geometry=zone_geoms, crs=gdf.crs)
            zone_path = cleaned_dir / f"{base}_overlap_zones.shp"
            zones.to_file(zone_path)
            print(f"    {zone_path.name}  ({len(zones):,} overlap zones)")

    # --- Audit report ------------------------------------------------------
    if write_report:
        report_path = cleaned_dir / f"{base}_report.csv"
        fieldnames = ["row_idx"]
        if id_field is not None:
            fieldnames.append(id_field)
        fieldnames += ["original_type", "was_valid", "status", "repaired"]
        if clean_parts:
            fieldnames += ["parts_drop", "drop_area"]
        if flag_slivers:
            fieldnames += ["width", "thinness", "sliver"]
        if flag_overlaps:
            fieldnames += ["ov_area", "ov_n"]
        fieldnames += ["validity_msg"]

        with open(report_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(report_rows)
        print(f"  Wrote audit report:")
        print(f"    {report_path}")

    meta = {
        "crs": str(gdf.crs),
        "crs_kind": crs_kind,
        "crs_units": crs_units,
        "n_features": len(gdf),
        "n_slivers": n_slivers,
        "sliver_width": sliver_width if flag_slivers else None,
        "n_overlapping": n_overlapping,
        "total_overlap_area": total_overlap_area,
        "n_part_features": n_part_features if clean_parts else None,
        "dropped_area_total": dropped_area_total,
        "n_valid": len(valid_rows),
        "n_rejected": len(rejected_records),
        "base": base,
        "in_place": in_place,
    }
    return report_rows, meta


def print_summary(report_rows, meta=None):
    """Print a tidy summary of what happened to each feature."""
    counts = Counter(r["status"] for r in report_rows)
    total = len(report_rows)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    if meta:
        print(f"CRS: {meta.get('crs')}  "
              f"({meta.get('crs_kind')}, units: {meta.get('crs_units')})")
        if meta.get("crs_kind") == "geographic":
            print("  NOTE: geographic CRS -- area, width, and overlap values")
            print("        are in degrees and not suitable for ET / water work.")
    print(f"Total features processed: {total:,}\n")
    print(f"  {'Status':<35s} {'Count':>8s} {'%':>7s}")
    print(f"  {'-' * 35} {'-' * 8} {'-' * 7}")
    for status, n in sorted(counts.items()):
        pct = 100 * n / total if total else 0
        print(f"  {status:<35s} {n:>8,d} {pct:>6.1f}%")

    if meta:
        units = meta.get("crs_units", "units")
        if meta.get("n_slivers") is not None:
            print(f"\n  Thin sliver features flagged: {meta['n_slivers']:,} "
                  f"(mean width < {meta.get('sliver_width')} {units})")
        if meta.get("n_overlapping") is not None:
            print(f"  Features with overlaps:       {meta['n_overlapping']:,}")
            print(f"  Total overlap area:           "
                  f"{meta.get('total_overlap_area', 0):,.2f} {units}^2")
        if meta.get("n_part_features") is not None:
            print(f"  Multipart features cleaned:   {meta['n_part_features']:,} "
                  f"(area removed: {meta.get('dropped_area_total', 0):,.4f} {units}^2)")
    print()


# ---------------------------------------------------------------------------
# Bulk / in-place processing
# ---------------------------------------------------------------------------

# Suffixes appended to output shapefiles, so batch mode can skip its own output.
GENERATED_SUFFIXES = (
    "_valid", "_rejected_poly", "_rejected_line", "_rejected_point",
    "_slivers", "_overlap_zones", "_removed_parts",
)
# Subdirectories batch mode creates and should never treat as input.
SKIP_DIRS = {"cleaned", "_review"}
# Stale sidecars an overwrite would leave behind (spatial indexes, metadata).
STALE_SIDECARS = (".sbn", ".sbx", ".qix", ".fix", ".shp.xml", ".aih", ".ain")


def backup_shapefile(shp_path: Path, backup_dir: Path) -> int:
    """Copy a shapefile and all its sidecars into ``backup_dir``. Returns the
    number of files copied."""
    backup_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for f in shp_path.parent.glob(shp_path.stem + ".*"):
        shutil.copy2(f, backup_dir / f.name)
        n += 1
    return n


def _overwrite_in_place(valid_gdf, shp_path: Path):
    """Overwrite the original shapefile with the cleaned valid features, then
    remove stale sidecars (spatial index, metadata) that no longer match."""
    valid_gdf.to_file(shp_path)
    for ext in STALE_SIDECARS:
        stale = (shp_path.parent / (shp_path.name + ".xml")
                 if ext == ".shp.xml" else shp_path.with_suffix(ext))
        if stale.exists():
            try:
                stale.unlink()
            except OSError:
                pass


def _is_generated(shp_path: Path) -> bool:
    return any(shp_path.stem.endswith(suf) for suf in GENERATED_SUFFIXES)


def list_batch_shapefiles(directory: Path, recursive: bool = False) -> list:
    """List input shapefiles in a directory, skipping generated output files
    and the directories this tool creates."""
    it = directory.rglob("*.shp") if recursive else directory.glob("*.shp")
    result = []
    for p in sorted(it):
        rel_dirs = p.relative_to(directory).parts[:-1]
        if any(d in SKIP_DIRS or d.startswith("_backup") for d in rel_dirs):
            continue
        if _is_generated(p):
            continue
        result.append(p)
    return result


def print_batch_summary(directory, metas, errors, in_place, backup_dir, review_dir):
    print("\n" + "=" * 60)
    print("BATCH SUMMARY")
    print("=" * 60)
    total_feat = sum(m["n_features"] for _, m in metas)
    total_valid = sum(m["n_valid"] for _, m in metas)
    total_rej = sum(m["n_rejected"] for _, m in metas)
    total_parts = sum((m.get("n_part_features") or 0) for _, m in metas)
    geographic = [p.name for p, m in metas if m.get("crs_kind") == "geographic"]
    print(f"Files processed:      {len(metas):,}")
    if errors:
        print(f"Files with errors:    {len(errors):,}")
    print(f"Total features:       {total_feat:,}")
    print(f"Valid features kept:  {total_valid:,}")
    print(f"Rejected features:    {total_rej:,}")
    print(f"Multipart features cleaned: {total_parts:,}")
    if geographic:
        print(f"\n  WARNING: {len(geographic)} file(s) were in a geographic CRS "
              "(area-based\n  values are not meaningful): "
              + ", ".join(geographic[:5])
              + (" ..." if len(geographic) > 5 else ""))
    if in_place:
        print("\nOne-shot fix: originals were overwritten in place.")
        if backup_dir is not None:
            print(f"  Backups:              {backup_dir}")
        if review_dir is not None:
            print(f"  Reports and rejected: {review_dir}")
    for p, e in errors:
        print(f"  ERROR {p.name}: {e}")
    print()


def process_directory(directory: Path, *, recursive=False, in_place=False,
                      backup=True, backup_dir=None, dry_run=False, yes=False,
                      **kw):
    """Process every input shapefile in ``directory``."""
    directory = Path(directory)
    if not directory.is_dir():
        print(f"Not a directory: {directory}")
        sys.exit(1)

    shps = list_batch_shapefiles(directory, recursive)
    if not shps:
        print(f"No shapefiles found in {directory}"
              + (" (recursive)" if recursive else "") + ".")
        return

    print(f"Found {len(shps)} shapefile(s) in {directory}"
          + (" (recursive)" if recursive else "") + ":")
    for p in shps:
        print(f"  - {p.relative_to(directory)}")

    if dry_run:
        mode = "overwrite in place" if in_place else "write cleaned/ outputs for"
        print(f"\nDry run: no files will be modified. Would {mode} the "
              f"{len(shps)} file(s) above.")
        if in_place and backup:
            print("Originals would be backed up first.")
        return

    resolved_backup = None
    review_dir = None
    if in_place:
        review_dir = directory / "_review"
        if backup:
            resolved_backup = (Path(backup_dir) if backup_dir is not None
                               else directory / f"_backup_{time.strftime('%Y%m%d_%H%M%S')}")

    if in_place and not yes:
        msg = f"\nOne-shot fix: about to OVERWRITE {len(shps)} shapefile(s) in place"
        if resolved_backup is not None:
            msg += f" (originals backed up to {resolved_backup.name})"
        if not prompt_yes_no(msg + ". Continue?", default=False):
            print("Aborting; no files changed.")
            return

    metas, errors = [], []
    for p in shps:
        try:
            _, meta = process_shapefile(
                p, interactive=False, in_place=in_place,
                out_dir=(review_dir if in_place else None),
                backup_dir=resolved_backup, **kw)
            metas.append((p, meta))
        except Exception as e:  # keep going on a bad file
            print(f"  ERROR processing {p.name}: {e}")
            errors.append((p, str(e)))

    print_batch_summary(directory, metas, errors, in_place,
                        resolved_backup, review_dir)


def _run_interactive():
    cwd = Path.cwd()
    print(f"Working directory: {cwd}")

    shapefiles = list_shapefiles(cwd)
    if not shapefiles:
        print("\nNo .shp files found in this directory.")
        sys.exit(1)

    chosen = prompt_shapefile_choice(shapefiles)
    write_report = prompt_yes_no("Write per-feature CSV audit report?", default=True)
    add_repaired_field = prompt_yes_no(
        "Add flag fields (repaired, and sliver/overlap if enabled) to the valid output?\n"
        "  (No keeps the output attributes identical to the input; all flags\n"
        "  are still recorded in the CSV report and side files.)",
        default=True,
    )

    flag_slivers = prompt_yes_no(
        "Flag thin 'sliver' polygons (artifacts of prior overlays)?", default=False)
    sliver_width = 0.0
    if flag_slivers:
        sliver_width = prompt_float(
            "  Maximum sliver width in CRS units (mean width = 2*area/perimeter)",
            default=5.0)

    flag_overlaps = prompt_yes_no(
        "Check for overlapping polygons (they double-count area)?", default=False)
    min_overlap_area = 0.0
    if flag_overlaps:
        min_overlap_area = prompt_float(
            "  Minimum overlap area to flag in CRS units squared (0 = any overlap)",
            default=0.0)

    clean_parts = prompt_yes_no(
        "Clean degenerate / sliver parts out of multipart polygons?", default=False)
    min_part_area = 0.0
    if clean_parts:
        min_part_area = prompt_float(
            "  Minimum part area to keep in CRS units squared\n"
            "  (0 = drop only near-zero-area orphan parts)",
            default=0.0)

    report_rows, meta = process_shapefile(
        chosen,
        write_report=write_report,
        add_repaired_field=add_repaired_field,
        flag_slivers=flag_slivers,
        sliver_width=sliver_width,
        flag_overlaps=flag_overlaps,
        min_overlap_area=min_overlap_area,
        clean_parts=clean_parts,
        min_part_area=min_part_area,
    )
    print_summary(report_rows, meta)


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description="Validate, repair, and clean shapefile polygon geometries. "
                    "Run with no arguments for the interactive single-file mode.")
    p.add_argument("path", nargs="?",
                   help="A .shp file to process non-interactively. Omit "
                        "(and omit --batch) for interactive mode.")
    p.add_argument("--batch", metavar="DIR",
                   help="Process every top-level .shp in DIR non-interactively.")
    p.add_argument("--recursive", action="store_true",
                   help="With --batch, also recurse into subdirectories.")
    p.add_argument("--one-shot", "--in-place", dest="in_place",
                   action="store_true",
                   help="One-shot fix: overwrite each input with its cleaned "
                        "valid output (originals are backed up first unless "
                        "--no-backup). --in-place is an accepted alias.")
    p.add_argument("--no-backup", action="store_true",
                   help="Do not back up originals before overwriting (dangerous).")
    p.add_argument("--backup-dir", metavar="DIR",
                   help="Directory to hold backups of the originals.")
    p.add_argument("--out-dir", metavar="DIR",
                   help="Directory for cleaned/report outputs when not in place.")
    p.add_argument("--id-field", metavar="NAME",
                   help="Attribute column to use as the identifier in reports.")
    p.add_argument("--no-report", action="store_true",
                   help="Do not write the per-feature CSV audit report.")
    p.add_argument("--flag-fields", action="store_true",
                   help="Add repaired/sliver/overlap flag fields to the output "
                        "(off by default in batch to preserve the schema).")
    p.add_argument("--no-clean-parts", action="store_true",
                   help="Disable multipart part cleanup (on by default in batch).")
    p.add_argument("--min-part-area", type=float, default=0.0,
                   help="Minimum multipart part area to keep, CRS units squared.")
    p.add_argument("--flag-slivers", action="store_true",
                   help="Flag thin sliver features (sidecar report only).")
    p.add_argument("--sliver-width", type=float, default=5.0,
                   help="Maximum sliver width in CRS units (default 5).")
    p.add_argument("--flag-overlaps", action="store_true",
                   help="Check for overlapping polygons (sidecar report only).")
    p.add_argument("--min-overlap-area", type=float, default=0.0,
                   help="Minimum overlap area to flag, CRS units squared.")
    p.add_argument("--dry-run", action="store_true",
                   help="With --batch, list what would happen without writing.")
    p.add_argument("-y", "--yes", action="store_true",
                   help="Skip the in-place overwrite confirmation prompt.")
    return p


def main():
    args = _build_arg_parser().parse_args()

    # No target given -> classic interactive mode.
    if not args.batch and not args.path:
        _run_interactive()
        return

    common = dict(
        write_report=not args.no_report,
        add_repaired_field=args.flag_fields,
        flag_slivers=args.flag_slivers,
        sliver_width=args.sliver_width,
        flag_overlaps=args.flag_overlaps,
        min_overlap_area=args.min_overlap_area,
        clean_parts=not args.no_clean_parts,
        min_part_area=args.min_part_area,
        id_field=args.id_field if args.id_field else None,
    )

    if args.batch:
        process_directory(
            Path(args.batch), recursive=args.recursive, in_place=args.in_place,
            backup=not args.no_backup, backup_dir=args.backup_dir,
            dry_run=args.dry_run, yes=args.yes, **common)
        return

    # Single file, non-interactive.
    shp = Path(args.path)
    if not shp.exists():
        print(f"File not found: {shp}")
        sys.exit(1)
    backup_dir = None
    if args.in_place and not args.no_backup:
        backup_dir = (Path(args.backup_dir) if args.backup_dir
                      else shp.parent / f"_backup_{time.strftime('%Y%m%d_%H%M%S')}")
    if args.in_place and not args.yes:
        if not prompt_yes_no(f"Overwrite {shp.name} in place?", default=False):
            print("Aborting; no files changed.")
            return
    report_rows, meta = process_shapefile(
        shp, interactive=False, in_place=args.in_place,
        out_dir=(Path(args.out_dir) if args.out_dir else None),
        backup_dir=backup_dir, **common)
    print_summary(report_rows, meta)


if __name__ == "__main__":
    main()
