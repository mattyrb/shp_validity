# shp_validity

A Python utility for validating, repairing, and cleaning geometries in
ESRI shapefiles. Designed for GIS and remote sensing workflows where
vector products — classification polygons, water masks, field boundaries —
often arrive with messy topology, mixed geometry types, and the
occasional null geometry.

The script is interactive: run it from a directory containing one or
more shapefiles, pick the one to clean, and it produces a tidy
`cleaned/` subdirectory with the valid features, any rejected features
(split by geometry type), and a per-feature CSV audit log.

## Features

- Interactive prompts — no command-line flags to remember.
- Validates every feature with Shapely's `is_valid` and `explain_validity`.
- Two-stage geometry repair: `make_valid()` first, then `buffer(0)` as a fallback.
- Drops non-polygon geometries (points, lines, etc.) and routes them to a separate rejected file.
- Extracts polygon parts from any `GeometryCollection` produced during repair.
- Optionally adds a boolean `repaired` field to the valid output so you can filter / symbolize repaired features in QGIS or ArcGIS. Decline the prompt to keep the output attributes identical to the input (repair status is still recorded in the CSV report).
- Flags thin "sliver" polygons by characteristic width (`2 * area / perimeter`), catching overlay artifacts that are thin but still cover real area, while leaving legitimately long, wide fields alone. Writes a separate slivers file for review.
- Flags overlapping polygons (which double-count ET and applied water), with per-feature overlap area, a pair list, and the overlap geometries as a layer you can open on a map.
- Warns when the layer is in a geographic CRS (degrees), since area, sliver width, and overlap area are then not meaningful for water calculations, and logs the CRS in the run summary.
- Cleans degenerate / sliver parts out of multipart polygons (the classic "one good part plus a couple of orphaned vertices" artifact): explodes each multipart, repairs each part, drops the junk, reassembles, and writes the removed parts to a side file for review.
- Runs interactively on one shapefile, or non-interactively in bulk over a directory with `--batch`, including a one-shot `--one-shot` fix that overwrites the originals after backing them up.
- Per-feature CSV audit log includes a user-selected identifier column for cross-reference.
- Robustness handling for missing CRS, Z/M dimensions, field-name collisions, encoding issues, and mixed singlepart/multipart output.

## Installation

```bash
pip install -r requirements.txt
```

Requires Python 3.9 or later and Shapely 2.0+ (for `force_2d` and the
modern `make_valid` behavior).

## Try it out

Drop a small shapefile into `sample_data/` and run the script there to
see the whole pipeline without touching your real working data:

```bash
cd sample_data
python ../shp_validity.py
```

Results land in `sample_data/cleaned/`, which is git-ignored, so demo
runs never get committed. See [`sample_data/README.md`](sample_data/README.md).

## Usage

Open a terminal in the directory containing your shapefile and run:

```bash
python shp_validity.py
```

The script will:

1. List the `.shp` files in the current working directory and ask which one to clean.
2. Ask whether to write the CSV audit report.
3. Ask whether to add flag fields (`repaired`, plus `sliver` and overlap fields when those checks are on). Decline to keep the output attributes identical to the input.
4. Ask whether to flag thin sliver polygons, and if so the maximum sliver width in CRS units.
5. Ask whether to check for overlapping polygons, and if so the minimum overlap area to flag.
6. Ask whether to clean degenerate / sliver parts out of multipart polygons, and if so the minimum part area to keep.
7. Warn you if the input has no `.prj`, or if the CRS is geographic, and ask whether to proceed.
8. Strip any Z / M dimensions automatically (with a note).
9. Show the input's attribute columns and ask which one to use as a feature identifier in the CSV report.
10. Run the pipeline and print a summary, including the CRS and any sliver / overlap / part-cleanup counts.

If you add flag fields and your input already has a column of that name, the script picks a non-conflicting name (`was_fixed`, `is_sliver`, etc.) and tells you which it used.

See [example_usage.md](example_usage.md) for a full walkthrough with
sample console output.

## Command line and bulk mode

Run with no arguments for the interactive single-file mode above. Pass a path
or `--batch` to run non-interactively:

```bash
# one file, cleaned output written to a cleaned/ folder
python shp_validity.py path/to/fields.shp

# every shapefile in a directory, one-shot fix in place (originals backed up first)
python shp_validity.py --batch path/to/dir --one-shot

# preview what a one-shot run would touch, without writing anything
python shp_validity.py --batch path/to/dir --one-shot --dry-run
```

The **one-shot** approach (`--one-shot`, also accepted as `--in-place`)
overwrites each input shapefile with its cleaned,
valid features. Before overwriting, the originals and all their sidecars are
copied to a timestamped `_backup_...` folder (use `--no-backup` to skip, or
`--backup-dir` to choose the location). Features that cannot be made into valid
polygons are dropped from the fixed file and written, together with the
per-feature report, to a `_review/` folder, so nothing is lost silently. By
default the attribute schema is preserved (no flag fields are added) and
multipart part cleanup runs; use `--flag-fields`, `--no-clean-parts`, or
`--min-part-area` to change that, and `--yes` to skip the overwrite prompt in
scripted runs.

Other flags: `--recursive`, `--id-field NAME`, `--no-report`, `--out-dir DIR`,
`--flag-slivers` / `--sliver-width`, and `--flag-overlaps` / `--min-overlap-area`.
Run `python shp_validity.py --help` for the full list.

## Outputs

All outputs land in a `cleaned/` subdirectory next to the input shapefile:

| File | Contents |
| --- | --- |
| `<name>_valid.shp` | Valid `(Multi)Polygon` features, including features that were successfully repaired. Has a boolean `repaired` field unless you decline that prompt. |
| `<name>_rejected_poly.shp` | Polygons that could not be repaired (only created if any exist). |
| `<name>_rejected_line.shp` | Filtered `(Multi)LineString` features (only created if any exist). |
| `<name>_rejected_point.shp` | Filtered `(Multi)Point` features (only created if any exist). |
| `<name>_slivers.shp` | Thin sliver features flagged for review (only when sliver flagging is on and any are found). |
| `<name>_overlaps.csv` | Overlapping feature pairs with their overlap area (only when overlap checking is on). |
| `<name>_overlap_zones.shp` | The overlap geometries themselves, for map review (only when overlap checking is on). |
| `<name>_removed_parts.shp` | Degenerate / sliver parts removed from multipart features, with the parent row index and part area (only when part cleanup is on and any are removed). |
| `<name>_report.csv` | Per-feature audit log — the place to look up exactly why a feature ended up where it did. |

The CSV report contains one row per input feature with the columns:
`row_idx`, `<your_id_field>` (if you picked one), `original_type`,
`was_valid`, `status`, `repaired`, then `parts_drop`, `drop_area` (when part
cleanup is on), `width`, `thinness`, `sliver` (when sliver flagging is on),
`ov_area`, `ov_n` (when overlap checking is on), and finally `validity_msg`.

The `status` column uses these descriptive tags:

| Status | Meaning |
| --- | --- |
| `valid_original` | Feature was already valid; untouched. |
| `repaired_make_valid` | Repaired by `shapely.validation.make_valid()`. |
| `repaired_collection_extracted` | `make_valid()` produced a `GeometryCollection`; only the polygon parts were kept (line / point fragments discarded). |
| `repaired_buffer0` | Repaired by `geom.buffer(0)` after `make_valid()` failed. |
| `unrepairable` | Started as a polygon but neither repair method produced a valid result. |
| `null_geom` | Geometry was null / NaN. |
| `empty_geom` | Geometry was empty (no coordinates). |
| `wrong_type_point` etc. | Original geometry was not a polygon (point, line, etc.). |

## Quality checks for ET and applied-water work

These optional checks target problems that distort zonal area, and therefore
distort evapotranspiration and applied-water totals.

**CRS check.** ET depth times area, and applied-water volumes, are only
meaningful in a projected CRS. If the layer is geographic (degrees) the script
warns, and if you have area-based checks on it asks before continuing. The CRS
and its units are logged in the run summary.

**Sliver flagging.** Overlay operations leave thin sliver polygons that can
still cover a meaningful area, so an area threshold alone does not find them.
The script flags features by characteristic width, `2 * area / perimeter`,
which approaches the true width of an elongated shape regardless of its length.
A long but legitimately wide field is not flagged, while a 2 m wide artifact is,
even when both have similar area. A Polsby-Popper thinness ratio
(`4 * pi * area / perimeter**2`) is reported alongside so you can apply your own
cutoff. Flagged features are written to `<name>_slivers.shp` for review.

**Overlap flagging.** Overlapping polygons double-count ET and applied water.
The script finds overlaps with a spatial index, records each feature's total
overlap area and neighbor count, writes the overlapping pairs to
`<name>_overlaps.csv`, and writes the overlap geometries to
`<name>_overlap_zones.shp` so you can see exactly where they are on a map. This
flags overlaps for review; it does not resolve them, so you decide how to assign
or erase the shared area.

**Multipart part cleanup.** A multipart polygon can carry a perfectly good part
alongside a degenerate one, such as a couple of orphaned vertices left by an
earlier overlay. When the bad part is invalid it usually makes the whole feature
invalid, and the repair step above already drops it. When the bad part is a tiny
but technically valid sliver, the whole feature reads as valid and slips through
untouched. Part cleanup handles both: it explodes each multipart, repairs each
part on its own, drops any part that will not repair to a valid polygon or that
falls below the degenerate-area floor (or your `min_part_area`), then reassembles
the survivors. Removed parts are written to `<name>_removed_parts.shp` with their
parent feature index and area, so you can confirm what was erased. This automates
the manual explode / find-sliver / erase workflow.

## Repair pipeline

For each feature with a polygon-like input:

1. **`make_valid()`** — Shapely's modern geometry repair. Handles most real-world breakage (self-intersections, bowties, repeated vertices).
2. **Extract polygon parts** — if `make_valid()` returns a `GeometryCollection` (common when repairing self-intersections), keep only the `(Multi)Polygon` components and discard line / point fragments.
3. **`buffer(0)` fallback** — the older trick, kept around for rare cases `make_valid()` cannot handle.
4. **Reject** — if all repair attempts fail, send the feature to the rejected output with its original geometry preserved.

## Robustness handling

- **Missing `.prj`** — script warns clearly and prompts for confirmation before proceeding.
- **Z and M dimensions** — stripped automatically with a console note (this script is intended for 2D output).
- **`repaired` field name collision** — if the input already has a `repaired` column, the script picks an alternate name (`was_fixed`, `fixed`, etc.) and tells you which one it used.
- **DBF encoding** — falls back to CP1252 if the default UTF-8 read fails (handles older shapefiles with non-ASCII attribute values).
- **Mixed singlepart / multipart output** — all output geometries are normalized to `Multi*` before writing, sidestepping pyogrio / fiona quirks. ESRI shapefile format doesn't distinguish single vs. multi at the format level, so re-reading the output returns `Polygon` / `MultiPolygon` exactly as expected.

## Tips for using the output

- **Inspect repaired features in QGIS / ArcGIS** by filtering `"repaired" = 1` (or `'true'` depending on driver) on the valid output. Categorical symbology on the same field also works.
- **Trace a rejected feature back to its origin** by opening the CSV report and matching on `row_idx` or your chosen identifier column.
- **Re-running the script** silently overwrites previous outputs in `cleaned/`.

## Requirements

- Python 3.9+
- [`geopandas`](https://geopandas.org/) >= 0.13
- [`shapely`](https://shapely.readthedocs.io/) >= 2.0

See [`requirements.txt`](requirements.txt) for the exact constraints.

## License

MIT — see [LICENSE](LICENSE).

## Author

Matt Bromley, Desert Research Institute
