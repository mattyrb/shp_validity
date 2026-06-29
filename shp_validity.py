#!/usr/bin/env python3
"""
shp_validity.py - Shapefile geometry validator and cleaner

Validates and repairs polygon geometries in an ESRI shapefile, separating
valid features from rejected ones. Designed for cleaning vector products
derived from Landsat imagery (e.g., classification polygons, water masks,
field delineations).

USAGE
-----
Run from a directory containing one or more .shp files. The script lists
shapefiles in the current working directory and prompts you to pick one.

OUTPUTS (written to a 'cleaned/' subdirectory next to the input)
---------------------------------------------------------------
    <name>_valid.shp            Valid polygon / multipolygon features.
                                Optionally includes a boolean 'repaired'
                                field -- filter on it in GIS to inspect
                                repairs. Decline the prompt to keep the
                                output attributes identical to the input.
    <name>_rejected_poly.shp    Unrepairable polygons (only if any).
    <name>_rejected_line.shp    Filtered (multi)line features (only if any).
    <name>_rejected_point.shp   Filtered (multi)point features (only if any).
    <name>_report.csv           Per-feature audit log with row index, the
                                user-chosen identifier column, the original
                                geometry type, validity, repair status, and
                                Shapely's validity diagnostic message. Use
                                this to look up why a given feature ended
                                up in a rejected file.

REPAIR PIPELINE
---------------
    1. make_valid()  - Shapely's modern geometry repair.
    2. Extract polygon parts from any resulting GeometryCollection.
       (Line / point fragments produced by repair are discarded; logged
       in the report as 'repaired_collection_extracted'.)
    3. buffer(0) fallback for any remaining cases.
    4. Reject if still bad.

ROBUSTNESS HANDLING
-------------------
    - Warns and prompts before continuing if the input has no .prj / CRS.
    - Z and M dimensions are stripped automatically (with a note).
    - If the input already has a 'repaired' column, the repair flag is
      written under a non-conflicting name (was_fixed, fixed, etc.).
    - Falls back to CP1252 encoding if the default read fails (handles
      older DBF files with non-UTF-8 attribute values).
    - All output geometries are normalized to multipart for predictable
      shapefile output.

REQUIREMENTS
------------
    geopandas, shapely >= 2.0

Author: written for Matt Bromley (DRI), GIS / remote sensing workflows.
"""

from __future__ import annotations

import csv
import sys
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
    """
    True if ``geom`` is a real shapely geometry object.

    Shapefiles with missing geometries can produce ``None`` *or* float NaN
    once they're loaded by geopandas, depending on the driver. A simple
    ``is not None`` check misses the NaN case and blows up on attribute
    access, so we check the actual type instead.
    """
    return isinstance(geom, BaseGeometry)


def _to_multipart(geom):
    """
    Promote singleton geometries to their Multi* equivalents.

    Mixed singlepart / multipart geometries in one shapefile occasionally
    confuse the pyogrio and fiona drivers. Normalizing to multipart before
    writing produces predictable output. (ESRI shapefile format itself does
    not distinguish single vs. multi, so re-reading the file gives back the
    same Polygon / MultiPolygon mix users expect.)
    """
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
    """
    Read a shapefile, falling back to CP1252 encoding if the default fails.

    Older shapefiles -- especially those with non-ASCII attribute values
    (Spanish place names, accented characters) -- are often encoded in
    CP1252 / Windows-1252 rather than UTF-8. The default read may either
    raise an encoding error or silently produce mojibake. We try the
    default first, then retry once with CP1252.
    """
    try:
        return gpd.read_file(shp_path)
    except UnicodeDecodeError:
        pass
    except Exception as e:
        # pyogrio / fiona may wrap encoding errors in their own exceptions.
        msg = str(e).lower()
        if "decode" not in msg and "encoding" not in msg and "utf" not in msg:
            raise

    print("  Note: default encoding failed; retrying with CP1252.")
    return gpd.read_file(shp_path, encoding="cp1252")


def strip_z_dimensions(gdf):
    """
    Strip Z (and any M) dimensions from all geometries, warning if any
    were present. The user does not want Z/M values in these datasets.
    """
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
    """
    Return a field name that does not already exist in ``gdf.columns``.
    Tries the preferred name first, then a list of alternatives, then
    numbered fallbacks. Shapefile field names are bound to 10 chars, so
    we keep all candidates within that limit.
    """
    if alternatives is None:
        alternatives = ["was_fixed", "fixed", "is_fixed", "repair_flg"]

    if preferred not in gdf.columns:
        return preferred

    for alt in alternatives:
        if alt not in gdf.columns:
            return alt

    # Numeric fallbacks, keeping within shapefile's 10-char field-name limit.
    for i in range(2, 100):
        candidate = f"repaired{i}"[:10]
        if candidate not in gdf.columns:
            return candidate

    raise RuntimeError("Could not find an unused field name for the repair flag.")


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


def prompt_id_field(gdf) -> str | None:
    """
    Prompt the user to choose an attribute column to use as a feature
    identifier in the audit report. Returns the column name, or None if
    the user opts to use the row index alone.
    """
    # Exclude the geometry column from the candidate list.
    geom_col = gdf.geometry.name
    candidates = [c for c in gdf.columns if c != geom_col]

    if not candidates:
        print("\n  No attribute columns found -- using row index as identifier.")
        return None

    print("\nAttribute columns available for use as an identifier:")
    for i, col in enumerate(candidates, 1):
        # Show a sample value to help the user recognize the column.
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
    """
    Pull (Multi)Polygon parts out of a geometry, discarding line / point
    fragments. Returns a Polygon, MultiPolygon, or None if there are no
    polygon parts to keep.
    """
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

    # LineString, Point, etc. -- nothing polygonal to recover.
    return None


def repair_geometry(geom):
    """
    Try to coerce ``geom`` into a valid (Multi)Polygon.

    Returns a tuple ``(repaired_geom_or_None, status_tag)``. The status_tag
    is a short string describing what happened, suitable for logging in the
    audit report.
    """
    if not _is_geometry(geom):
        return None, "null_geom"
    if geom.is_empty:
        return None, "empty_geom"

    geom_type = geom.geom_type

    # Drop non-polygon source geometries outright. They keep their original
    # geometry on the way to the rejected file for spatial QA.
    if geom_type in LINE_TYPES or geom_type in POINT_TYPES:
        return None, f"wrong_type_{geom_type.lower()}"
    if geom_type not in POLY_TYPES and geom_type != "GeometryCollection":
        return None, f"wrong_type_{geom_type.lower()}"

    # Already-valid (Multi)Polygon -- nothing to do.
    if geom_type in POLY_TYPES and geom.is_valid:
        return geom, "valid_original"

    # Step 1: Shapely's make_valid() handles most real-world breakage.
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

    # Step 2: buffer(0) fallback -- older trick that occasionally salvages
    # cases make_valid does not.
    try:
        buffered = geom.buffer(0)
        polys = extract_polygons(buffered)
        if polys is not None and polys.is_valid and not polys.is_empty:
            return polys, "repaired_buffer0"
    except Exception:
        pass

    return None, "unrepairable"


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
    return None  # null, empty, or GeometryCollection -> CSV only


def write_rejected_files(rejected_records, cleaned_dir: Path, base: str, crs):
    """
    Split rejected features by geometry type and write each as its own .shp.

    Features with no usable geometry (null, empty, or GeometryCollection)
    are skipped here -- they live only in the CSV report.

    Returns a list of (path, feature_count) tuples for what was written.
    """
    buckets: dict[str, list] = {"poly": [], "line": [], "point": []}

    for rec, status, geom_type in rejected_records:
        bucket = _bucket_for_geom_type(geom_type)
        if bucket is None:
            continue
        # Original attributes pass through unchanged; the rejection reason
        # for each feature lives in the CSV audit report (matched by
        # row_idx and the user's chosen identifier column).
        buckets[bucket].append(dict(rec))

    written = []
    for bucket, rows in buckets.items():
        if not rows:
            continue
        gdf = gpd.GeoDataFrame(rows, crs=crs)
        # Normalize to multipart so each rejected file has a single,
        # predictable geometry type.
        gdf["geometry"] = gdf.geometry.apply(_to_multipart)
        out_path = cleaned_dir / f"{base}_rejected_{bucket}.shp"
        gdf.to_file(out_path)
        written.append((out_path, len(gdf)))
    return written


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------

def process_shapefile(shp_path: Path, write_report: bool = True,
                      add_repaired_field: bool = True):
    """Run the full validate / repair / write pipeline on one shapefile.

    If ``add_repaired_field`` is False, the valid output keeps the input
    attributes exactly as they were -- no extra column is added. The
    per-feature repair status is still recorded in the CSV report, so you
    do not lose the information; it just lives outside the shapefile.
    """
    print(f"\nReading: {shp_path.name}")
    try:
        gdf = read_shapefile_robust(shp_path)
    except Exception as e:
        print(f"  ERROR reading shapefile: {e}")
        sys.exit(1)

    print(f"  Features: {len(gdf):,}")
    print(f"  CRS:      {gdf.crs}")

    # --- CRS sanity check --------------------------------------------------
    if gdf.crs is None:
        print("\n  WARNING: this shapefile has no CRS (the .prj file is")
        print("  missing or unreadable). Output shapefiles will also lack")
        print("  a .prj, and downstream tools may misinterpret coordinates.")
        if not prompt_yes_no("  Proceed anyway?", default=False):
            print("  Aborting.")
            sys.exit(1)

    # --- Strip Z / M dimensions -------------------------------------------
    gdf = strip_z_dimensions(gdf)

    # --- Pick a non-colliding field name for the repaired flag ------------
    # Only relevant if we are tagging the output. When the flag is off, the
    # valid output attributes are left identical to the input.
    if add_repaired_field:
        repaired_field = pick_unique_field_name(gdf, "repaired")
        if repaired_field != "repaired":
            print(f"  Note: input already has a 'repaired' column; "
                  f"using '{repaired_field}' for the repair flag instead.")
    else:
        repaired_field = None

    # Let the user pick an attribute column to use as a feature identifier
    # in the CSV report. The row index is always also included.
    id_field = prompt_id_field(gdf)

    valid_rows = []
    rejected_records = []  # tuples of (row_dict, status, original_geom_type)
    report_rows = []

    for row_idx, row in gdf.iterrows():
        geom = row.geometry
        has_geom = _is_geometry(geom)
        original_type = geom.geom_type if has_geom else "None"
        was_valid = bool(geom.is_valid) if has_geom else False

        # explain_validity is only meaningful for invalid geometries.
        validity_msg = ""
        if has_geom and not was_valid:
            try:
                validity_msg = explain_validity(geom)
            except Exception as e:
                validity_msg = f"explain_validity error: {e}"

        repaired_geom, status = repair_geometry(geom)
        was_repaired = status in REPAIRED_STATUSES

        report_row = {"row_idx": row_idx}
        if id_field is not None:
            report_row[id_field] = row[id_field]
        report_row.update({
            "original_type": original_type,
            "was_valid": was_valid,
            "status": status,
            # The report always records repair status under a fixed column,
            # independent of whether the shapefile gets a flag field.
            "repaired": was_repaired,
            "validity_msg": validity_msg,
        })
        report_rows.append(report_row)

        if repaired_geom is not None:
            new_row = row.copy()
            new_row.geometry = repaired_geom
            # Optional bool flag in the output shapefile so you can filter
            # repaired features in QGIS / ArcGIS. Skipped when the user wants
            # the output attributes to match the input exactly.
            if repaired_field is not None:
                new_row[repaired_field] = was_repaired
            valid_rows.append(new_row)
        else:
            rejected_records.append((row.to_dict(), status, original_type))

    # Output directory lives next to the input shapefile.
    cleaned_dir = shp_path.parent / "cleaned"
    cleaned_dir.mkdir(exist_ok=True)
    base = shp_path.stem

    # --- Valid output ------------------------------------------------------
    if valid_rows:
        valid_gdf = gpd.GeoDataFrame(valid_rows, crs=gdf.crs)
        # Normalize to multipart for predictable shapefile output.
        valid_gdf["geometry"] = valid_gdf.geometry.apply(_to_multipart)
        valid_path = cleaned_dir / f"{base}_valid.shp"
        valid_gdf.to_file(valid_path)
        print(f"\n  Wrote {len(valid_gdf):,} valid features:")
        print(f"    {valid_path}")
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

    # --- Audit report ------------------------------------------------------
    if write_report:
        report_path = cleaned_dir / f"{base}_report.csv"
        # Build the fieldnames in a stable order: row_idx, optional user ID,
        # then the standard descriptor columns.
        fieldnames = ["row_idx"]
        if id_field is not None:
            fieldnames.append(id_field)
        fieldnames.extend(["original_type", "was_valid", "status",
                           "repaired", "validity_msg"])

        with open(report_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(report_rows)
        print(f"  Wrote audit report:")
        print(f"    {report_path}")

    return report_rows


def print_summary(report_rows):
    """Print a tidy summary of what happened to each feature."""
    counts = Counter(r["status"] for r in report_rows)
    total = len(report_rows)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total features processed: {total:,}\n")
    print(f"  {'Status':<35s} {'Count':>8s} {'%':>7s}")
    print(f"  {'-' * 35} {'-' * 8} {'-' * 7}")
    for status, n in sorted(counts.items()):
        pct = 100 * n / total if total else 0
        print(f"  {status:<35s} {n:>8,d} {pct:>6.1f}%")
    print()


def main():
    cwd = Path.cwd()
    print(f"Working directory: {cwd}")

    shapefiles = list_shapefiles(cwd)
    if not shapefiles:
        print("\nNo .shp files found in this directory.")
        sys.exit(1)

    chosen = prompt_shapefile_choice(shapefiles)
    write_report = prompt_yes_no("Write per-feature CSV audit report?", default=True)
    add_repaired_field = prompt_yes_no(
        "Add a 'repaired' flag field to the valid output?\n"
        "  (No keeps the output attributes identical to the input; repair\n"
        "  status is still recorded in the CSV report.)",
        default=True,
    )

    report_rows = process_shapefile(
        chosen,
        write_report=write_report,
        add_repaired_field=add_repaired_field,
    )
    print_summary(report_rows)


if __name__ == "__main__":
    main()
