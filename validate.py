#!/usr/bin/env python3
"""Check a BIDS output directory against its cohort manifest.

Usage:
    python3 validate.py <cohort_manifest.csv> <bids_dir>

Prints how many selected scans were converted and lists any that are missing.
"""
import csv
import glob
import os
import sys


def main(manifest_path, bids_dir):
    # scans we expected to convert (one row per selected T1)
    expected = list(csv.DictReader(open(manifest_path)))

    # T1w files actually produced, by BIDS name (without the .nii.gz extension)
    produced = set()
    for path in glob.glob(f"{bids_dir}/sub-*/ses-*/anat/*_T1w.nii.gz"):
        produced.add(os.path.basename(path).replace(".nii.gz", ""))

    missing = [row for row in expected if row["nifti"] not in produced]

    print(f"expected:  {len(expected)}")
    print(f"converted: {len(expected) - len(missing)}")
    print(f"missing:   {len(missing)}")
    for row in missing:
        print(f"  - {row['nifti']}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: python3 validate.py <cohort_manifest.csv> <bids_dir>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
