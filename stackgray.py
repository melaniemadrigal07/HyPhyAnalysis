#!/usr/bin/env python3

from pathlib import Path
import argparse
import sys
import re
from collections import defaultdict

import numpy as np

try:
    # tifffile supports 16-bit and ImageJ stacks
    import tifffile
except ImportError:
    print("This script requires 'tifffile'. Install with: pip install tifffile", file=sys.stderr)
    sys.exit(1)


# ---- Config defaults
DEFAULT_ROOT = r"/Volumes/MMMHD/8.11.25/Plates1-8_FinalProtocl_WithDataExport_MM_Nooutput_1/GWAS_CT "
DEFAULT_PLATES = ["Plate_1", "Plate_2", "Plate_4", "Plate_5", "Plate_6", "Plate_7", "Plate_8"]


# ---- Helpers ----
def is_tiff(p: Path) -> bool:
    return p.suffix.lower() in {".tif", ".tiff"}


def parse_name_flexible(name: str):
    """
    Attempt to parse (well, colorNum, frame) from filenames like:
      A1_xxx_2_xxx_xxx_0001.tif
    Primary path: split by underscores and pick parts[0], parts[2], and stem(parts[5]).
    Fallback: regex.
    """
    parts = name.split("_")
    if len(parts) >= 6:
        well = parts[0]
        colorNum = parts[2]
        frame = Path(parts[5]).stem  # "0001" from "0001.tif"
        return well, colorNum, frame

    m = re.match(
        r"^([A-H]\d{1,2})_[^_]*_(\d)_(?:[^_]*_)??(?:[^_]*_)??(\d+)\.(?:tif|tiff)$",
        name,
        flags=re.IGNORECASE,
    )
    if m:
        return m.group(1), m.group(2), m.group(3)

    return None, None, None


def load_gray(path: Path) -> np.ndarray:
    """Return 2D grayscale array; preserves dtype."""
    arr = tifffile.imread(str(path))
    if arr.ndim == 3:
        # handle (H,W,1) or squeezed
        if arr.shape[-1] == 1:
            arr = arr[..., 0]
        else:
            arr = np.squeeze(arr)
    if arr.ndim != 2:
        raise ValueError(f"Expected 2D grayscale for {path.name}, got {arr.shape}")
    return arr


def save_grayscale_stack(ch_paths: dict, out_path: Path, order=("R", "G", "B"), imagej=True, dry_run=False):
    """
    Write a single multi-page grayscale TIFF containing R/G/B as separate pages.
    Keeps only channels that match the first valid channel's shape & dtype.
    Returns number of pages written (0 if none).
    """
    if dry_run:
        # Simulate reading to report shapes/dtypes succinctly
        labels_present = [k for k in order if k in ch_paths]
        print(f"   -> would stack pages: {','.join(labels_present)} → {out_path.name}")
        return len(labels_present)

    planes = []
    labels = []
    ref_shape = None
    ref_dtype = None

    # Load in requested order; establish shape/dtype from first valid
    for key in order:
        if key not in ch_paths:
            continue
        try:
            arr = load_gray(ch_paths[key])
        except Exception as e:
            print(f"     Skipping {key} ({ch_paths[key].name}): load error: {e}")
            continue

        if ref_shape is None:
            ref_shape = arr.shape
            ref_dtype = arr.dtype
            planes.append(arr)
            labels.append(key)
        else:
            if arr.shape != ref_shape:
                print(f"     Skipping {key}: shape mismatch {arr.shape} vs {ref_shape}")
                continue
            if arr.dtype != ref_dtype:
                print(f"     Skipping {key}: dtype mismatch {arr.dtype} vs {ref_dtype}")
                continue
            planes.append(arr)
            labels.append(key)

    if not planes:
        return 0

    data = np.stack(planes, axis=0)  # (pages, H, W)
    # Ensure grayscale; ImageJ=True for nice stack behavior
    tifffile.imwrite(
        str(out_path),
        data,
        imagej=imagej,
        photometric="minisblack",
        metadata={"axes": "ZYX", "Labels": labels},
    )
    return data.shape[0]


# Main
def main():
    ap = argparse.ArgumentParser(description="Stack 1/2/3 channels into a multi-page grayscale TIFF by (well, frame).")
    ap.add_argument("--root", type=str, default=DEFAULT_ROOT, help="Root directory containing plate folders")
    ap.add_argument("--plates", type=str, nargs="*", default=DEFAULT_PLATES, help="Plate folder names (relative to root)")
    ap.add_argument("--dry-run", action="store_true", help="List actions without writing files")
    args = ap.parse_args()

    root = Path(args.root).expanduser().resolve()
    print(f" Running on: {sys.platform}")
    print(f" Root: {root}")

    if not root.exists():
        print(f"ERROR: Root directory does not exist: {root}", file=sys.stderr)
        sys.exit(1)

    for plate in args.plates:
        input_dir = root / plate
        if not input_dir.exists():
            print(f"⚠️  Skipping missing plate folder: {input_dir}")
            continue

        output_dir = input_dir / "stacked"
        if not args.dry_run:
            output_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n📁 Processing: {input_dir}")

        # Map: group_id -> {'R': Path, 'G': Path, 'B': Path}
        groups = defaultdict(dict)

        files = [p for p in input_dir.iterdir() if p.is_file() and is_tiff(p)]
        for f in files:
            well, colorNum, frame = parse_name_flexible(f.name)
            if well is None or colorNum is None or frame is None:
                print(f"     Skipping (unexpected name): {f.name}")
                continue

            # Map color number to channel letter
            if colorNum == "1":
                chan = "R"
            elif colorNum == "2":
                chan = "G"
            elif colorNum == "3":
                chan = "B"
            else:
                print(f"     Unknown color number in: {f.name}")
                continue

            group_id = f"{well}_{frame}"
            groups[group_id][chan] = f

        # Process each (well, frame) group
        total_groups = len(groups)
        done, empty = 0, 0

        for group_id, ch in sorted(groups.items()):
            present = sorted(ch.keys())
            if not present:
                print(f"   (no channels) {group_id}")
                empty += 1
                continue

            out_path = output_dir / f"{group_id}.tif"
            n_pages = save_grayscale_stack(ch, out_path, order=("R", "G", "B"), imagej=True, dry_run=args.dry_run)

            if n_pages == 0:
                print(f"    Could not create stack for {group_id} (no valid channels)")
            else:
                print(f"    Saved grayscale stack ({n_pages} page{'s' if n_pages!=1 else ''}): {group_id}")
                done += 1

        print(f"📊 Summary for {plate}: {done}/{total_groups} stacks written, {empty} empty groups")

    print("\n✅ Finished all plates.")


if __name__ == "__main__":
    main()
