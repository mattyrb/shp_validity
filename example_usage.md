# Example usage

An illustrative walkthrough of running `shp_validity.py`, showing the
prompts and the files it produces. The dataset and numbers below are made
up for illustration; your run will reflect your own data.

## Setup

```bash
# From the root of this repo
pip install -r requirements.txt
```

## Running the script

The script reads from the current working directory, so `cd` to the data
first. To try it, drop a small shapefile in `sample_data/` and run:

```bash
cd sample_data
python ../shp_validity.py
```

The script reads from the current working directory, so it does not
matter where the `.py` lives. `cd` to the data first.

## Sample interactive session

The example below uses a fictional `field_boundaries.shp` with 3,500
features: some invalid polygons, a few stray lines, and a couple of null
geometries.

```
Working directory: /data/fields

Shapefiles available in this directory:
  1. field_boundaries.shp
  2. irrigation_districts.shp

Select a shapefile [1-2]: 1
Write per-feature CSV audit report? [Y/n]: y
Add a 'repaired' flag field to the valid output?
  (No keeps the output attributes identical to the input; repair
  status is still recorded in the CSV report.) [Y/n]: y

Reading: field_boundaries.shp
  Features: 3,500
  CRS:      EPSG:32613
  Note: 8 geometries had Z/M dimensions -- stripping to 2D.

Attribute columns available for use as an identifier:
   1. OBJECTID         e.g. 1
   2. FIELD_ID         e.g. NM-BERN-00142
   3. OWNER            e.g. Sandia Pueblo
   4. AREA_ACRES       e.g. 124.3
   5. CROP_2024        e.g. alfalfa

Select an identifier column [1-5], or press Enter to use row index only: 2

  Wrote 3,421 valid features:
    /data/fields/cleaned/field_boundaries_valid.shp
  Wrote rejected features to:
    field_boundaries_rejected_poly.shp  (51 features)
    field_boundaries_rejected_line.shp  (24 features)
    (4 rejected feature(s) had no usable geometry -- see CSV report)
  Wrote audit report:
    /data/fields/cleaned/field_boundaries_report.csv

============================================================
SUMMARY
============================================================
Total features processed: 3,500

  Status                                Count       %
  ----------------------------------- -------- -------
  empty_geom                                 2    0.1%
  null_geom                                  2    0.1%
  repaired_buffer0                           3    0.1%
  repaired_collection_extracted              7    0.2%
  repaired_make_valid                       38    1.1%
  unrepairable                              51    1.5%
  valid_original                         3,373   96.4%
  wrong_type_linestring                     24    0.7%
```

## Keeping the output attributes unchanged

Answer **n** to the `repaired` flag prompt when you need the valid output
to carry exactly the same fields as the input, for example when the file
feeds a downstream model or join that expects a fixed schema. The repair
status for every feature is still written to the CSV report, so nothing
is lost. It just lives outside the shapefile.

## What gets produced

After the run above, the `cleaned/` directory holds the valid layer, any
rejected layers (split by geometry type), and the CSV report:

```
cleaned/
├── field_boundaries_valid.shp          (+ .shx, .dbf, .prj, .cpg)
├── field_boundaries_rejected_poly.shp  (+ .shx, .dbf, .prj, .cpg)
├── field_boundaries_rejected_line.shp  (+ .shx, .dbf, .prj, .cpg)
└── field_boundaries_report.csv
```

Rejected shapefiles only appear when there are features of that type, so
cleaner datasets often produce just the `_valid` files and the CSV.

## The audit report

`cleaned/<name>_report.csv` has one row per input feature. A few
representative rows from the example above:

| row_idx | FIELD_ID | original_type | was_valid | status | repaired | validity_msg |
| --- | --- | --- | --- | --- | --- | --- |
| 12 | NM-BERN-00013 | Polygon | False | repaired_make_valid | True | Self-intersection[-106.612 35.084] |
| 217 | NM-DONA-00041 | Polygon | True | valid_original | False | |
| 884 | NM-CHAV-00322 | LineString | False | wrong_type_linestring | False | |
| 1402 | NM-RIO-00109 | Polygon | False | unrepairable | False | Ring Self-intersection[-105.872 34.420] |
| 2901 | NM-OTERO-00057 | None | False | null_geom | False | |

To trace a rejected feature back to its origin, match on `row_idx` or
your chosen identifier column. Feature `NM-RIO-00109` lives in
`field_boundaries_rejected_poly.shp`; opening it in QGIS or ArcGIS shows
the original (broken) geometry so you can inspect what went wrong.

## Filtering repaired features in QGIS

If you added the `repaired` flag, you can see only the features repaired
during the run:

1. Open `<name>_valid.shp` in QGIS.
2. Right-click the layer → **Filter…**.
3. Enter `"repaired" = 1` (or `"repaired" = 'true'`, depending on which
   driver wrote the field).
4. Apply.

The same field works for categorical symbology
(Properties → Symbology → Categorized → `repaired`) so repaired features
stand out against features that were already valid.

## Sliver and overlap checks

Two optional prompts target artifacts that distort zonal area:

```
Flag thin 'sliver' polygons (artifacts of prior overlays)? [y/N]: y
  Maximum sliver width in CRS units (mean width = 2*area/perimeter) [5.0]: 5
Check for overlapping polygons (they double-count area)? [y/N]: y
  Minimum overlap area to flag in CRS units squared (0 = any overlap) [0.0]: 0
```

With these on, the run adds `width`, `thinness`, and `sliver` (and `ov_area`,
`ov_n`) columns to the report, and the summary ends with counts:

```
  Thin sliver features flagged: 1 (mean width < 5.0 metre)
  Features with overlaps:       2
  Total overlap area:           5,000.00 metre^2
```

It also writes `<name>_slivers.shp` (the flagged thin features),
`<name>_overlaps.csv` (the overlapping pairs and their area), and
`<name>_overlap_zones.shp` (the overlap polygons themselves, for map review).
The width metric flags a 2 m wide artifact while leaving a legitimately long,
40 m wide field alone, so length by itself never triggers a false positive.
Run on a geographic (degree) layer, the script warns first that area, width,
and overlap values will not be meaningful.

## Cleaning multipart parts

If your data has multipart polygons where one part is good and another is a
couple of orphaned vertices, answer yes to the part-cleanup prompt:

```
Clean degenerate / sliver parts out of multipart polygons? [y/N]: y
  Minimum part area to keep in CRS units squared
  (0 = drop only near-zero-area orphan parts) [0.0]: 0
```

The script explodes each multipart, repairs each part, drops the junk, and
reassembles the survivors. The removed parts are written to
`<name>_removed_parts.shp` (with the parent row index and area) so you can
confirm what was erased, and the report gains `parts_drop` and `drop_area`
columns. Set a `min_part_area` above 0 to also drop small-but-real sliver
parts; a feature whose parts are all dropped is sent to the rejected output.

## Re-running

Re-running the script silently overwrites the contents of `cleaned/`.
To keep a previous run for comparison, rename or move the `cleaned/`
folder before re-running.
