#!/usr/bin/env python3
"""ADNI3 QC-selected T1 cohort: reproducible selection -> manifest -> dcm2bids (BIDS).

Selection (CSV-only): MRIQC `T1w` rows with `MRIProtocolPhase=ADNI3` (all 3T), grouped by
(subject, visit). Keep the MAYO-flagged scan (`SERIES_SELECTED=TRUE`, joined on
LONI_IMAGE == MRIQC image_id); tie / none-selected -> lowest image_id. Then locate DICOMs
in the ADNI3 raw zips and, with --bids, run dcm2bids (acq-3T).

Usage:
  python3 adni3_build.py                 # manifest + download list (no raws needed)
  python3 adni3_build.py --bids          # convert (after ADNI3 raws downloaded)
  python3 adni3_build.py --bids -n 3
"""
import argparse
import csv
import json
import re
import shutil
import subprocess
import tempfile
import zipfile
from collections import defaultdict
from glob import glob
from pathlib import Path

ROOT = Path("/scratch/moncia/ADNI")
MRIQC = ROOT / "MRIQC_01Jun2026.csv"
MAYO = ROOT / "MAYOADIRL_MRI_QUALITY_ADNI3_01Jun2026.csv"
# ADNI3 raw zips land here once downloaded (drop anywhere matching these globs)
ZIPS = sorted(set(glob(str(ROOT / "3" / "*.zip")) + glob(str(ROOT / "*adni3*.zip")) +
                  glob(str(ROOT / "*ADNI3*.zip")) + glob(str(ROOT / "*adni_3*.zip"))))
OUT_MANIFEST = ROOT / "adni3_cohort_manifest.csv"
OUT_IDS = ROOT / "adni3_image_ids.txt"
BIDS_DIR = ROOT / "bids_adni3"

DCM_IN_FOLDER = re.compile(r"/I(\d+)/[^/]+\.dcm$")
# guard against MRIQC SeriesType mislabels (e.g. a FLAIR tagged T1w): drop non-T1 descriptions
NON_T1 = re.compile(r"flair|t2|dti|dwi|\basl\b|swi|perfus|localiz|scout|field.?map", re.I)


def select():
    """One per (subject, visit): MAYO SERIES_SELECTED preferred, else lowest image_id."""
    sel, qual = set(), {}
    for r in csv.DictReader(open(MAYO)):
        if r["SERIES_SELECTED"] == "TRUE":
            sel.add(r["LONI_IMAGE"])
        qual[r["LONI_IMAGE"]] = r.get("SERIES_QUALITY", "")
    mq = [r for r in csv.DictReader(open(MRIQC))
          if r["MRIProtocolPhase"] == "ADNI3" and r["SeriesType"] == "T1w"
          and not NON_T1.search(r["SeriesDescription"]) and r["VISCODE2"].strip()]
    grp = defaultdict(list)
    for r in mq:
        grp[(r["ParticipantID"], r["VISCODE2"], r["MagneticFieldStrength"])].append(r)
    out = []
    for v in grp.values():
        flagged = [r for r in v if r["image_id"] in sel]
        pool = flagged if flagged else v
        w = min(pool, key=lambda r: int(r["image_id"]))
        w["_reason"] = "mayo_selected" if flagged else "unselected_fallback"
        w["_qual"] = qual.get(w["image_id"], "")
        out.append(w)
    return out


def locate():
    loc = {}
    for z in ZIPS:
        with zipfile.ZipFile(z) as zf:
            members = defaultdict(list)
            for n in zf.namelist():
                m = DCM_IN_FOLDER.search(n)
                if m:
                    members[m.group(1)].append(n)
        for iid, names in members.items():
            loc.setdefault(iid, (z, names))
    return loc


def resolve_run_splits(anat, base):
    """A Philips multi-ImageType series makes dcm2bids write run-01/run-02 instead of <base>.
    Rename run-01 -> canonical <base>, run-NN (N>=2) -> <...rec-NN..> (keeps every variant)."""
    if (anat / f"{base}.nii.gz").exists():
        return
    prefix = base[:-4]  # strip '_T1w'
    for nii in sorted(anat.glob(f"{prefix}_run-*_T1w.nii.gz")):
        m = re.search(r"_run-(\d+)_T1w\.nii\.gz$", nii.name)
        if not m:
            continue
        n = int(m.group(1))
        new = base if n == 1 else f"{prefix}_rec-{n}_T1w"
        for ext in ("nii.gz", "json"):
            src = anat / nii.name.replace("nii.gz", ext)
            if src.exists():
                src.rename(anat / f"{new}.{ext}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bids", action="store_true")
    ap.add_argument("-n", type=int, default=0)
    ap.add_argument("--subject", default=None)
    args = ap.parse_args()

    sel = select()
    loc = locate()
    rows = []
    for r in sel:
        iid = r["image_id"]
        l = loc.get(iid)
        subj, vis = r["ParticipantID"], r["VISCODE2"]
        rows.append({
            "image_id": iid, "subject": subj, "visit": vis, "field_strength": "3.0",
            "series": r["LONISeries"], "description": r["SeriesDescription"],
            "series_quality": r["_qual"], "reason": r["_reason"],
            "zip": Path(l[0]).name if l else "MISSING", "n_dcm": len(l[1]) if l else 0,
            "nifti": f"sub-{subj.replace('_', '')}_ses-{vis}_acq-3T_T1w",
        })
    with open(OUT_MANIFEST, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    with open(OUT_IDS, "w") as f:
        f.write(", ".join(sorted("I" + r["image_id"] for r in rows)))
    on = sum(1 for r in rows if r["zip"] != "MISSING")
    fb = sum(1 for r in rows if r["reason"] == "unselected_fallback")
    print(f"manifest: {len(rows)} scans -> {OUT_MANIFEST}")
    print(f"  unselected_fallback (no MAYO flag in session): {fb}")
    print(f"  on disk: {on} | to download: {len(rows) - on} -> {OUT_IDS}")

    if not args.bids:
        return
    if not shutil.which("dcm2bids"):
        print("ERROR: dcm2bids not on PATH (module load dcm2bids).")
        return
    if not ZIPS:
        print("ERROR: no ADNI3 zips found — download raws first (see adni3_image_ids.txt).")
        return
    todo = [r for r in rows if r["zip"] != "MISSING"]
    if args.subject:
        todo = [r for r in todo if r["subject"] == args.subject]
    if args.n:
        todo = todo[:args.n]
    BIDS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"dcm2bids: {len(todo)} scans -> {BIDS_DIR}")
    for i, r in enumerate(todo, 1):
        z, names = loc[r["image_id"]]
        tmp = Path(tempfile.mkdtemp())
        try:
            with zipfile.ZipFile(z) as zf:
                for n in names:
                    zf.extract(n, tmp)
            # two rules: dcm2niix's equalized "_Eq" twin -> rec-eq, everything else -> canonical T1w
            cfg = tmp / "cfg.json"
            cfg.write_text(json.dumps({"descriptions": [
                {"datatype": "anat", "suffix": "T1w", "custom_entities": ["acq-3T", "rec-eq"],
                 "criteria": {"SidecarFilename": "*_Eq*"}},
                {"datatype": "anat", "suffix": "T1w", "custom_entities": ["acq-3T"],
                 "criteria": {"Modality": "MR"}}]}))
            subprocess.run(["dcm2bids", "-d", str(tmp),
                            "-p", r["subject"].replace("_", ""), "-s", r["visit"],
                            "-c", str(cfg), "-o", str(BIDS_DIR), "--clobber"],
                           check=False, capture_output=True)
            anat = BIDS_DIR / f"sub-{r['subject'].replace('_', '')}" / f"ses-{r['visit']}" / "anat"
            resolve_run_splits(anat, r["nifti"])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        if i % 200 == 0:
            print(f"  {i}/{len(todo)}")
    print("done ->", BIDS_DIR)


if __name__ == "__main__":
    main()
