#!/usr/bin/env python3
"""ADNI1 QC-selected T1 cohort: reproducible selection -> manifest -> dcm2bids (BIDS).

Selection (CSV-only, reproducible): for each (subject, visit, field-strength) keep the
best MRIMPRANK-ranked T1 (MRIQC `T1w` rows joined to MRIMPRANK on series number;
tie -> lowest image_id). Then locate each chosen image's DICOMs in the raw zips and,
with --bids, run dcm2bids once per (subject, session) -> a valid BIDS dataset.
Field strength routes to acq-15T / acq-3T via adni_t1_dcm2bids_config.json.

Usage:
  python3 adni1_build.py                 # write manifest only
  python3 adni1_build.py --bids          # also dcm2bids all chosen scans
  python3 adni1_build.py --bids -n 3     # convert first 3 sessions (smoke test)

Requires dcm2bids (+ dcm2niix) on PATH only for --bids (module load first).
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
MRIMPRANK = ROOT / "MRIMPRANK_01Jun2026.csv"
MRIQC = ROOT / "MRIQC_01Jun2026.csv"
# all zips that may contain ADNI1 raw T1 (original 5 + any added drops)
ZIPS = sorted(glob(str(ROOT / "1" / "adni_1_rage_t1_3d*.zip"))) + \
       sorted(glob(str(ROOT / "*.zip")))
OUT_MANIFEST = ROOT / "adni1_cohort_manifest.csv"
BIDS_DIR = ROOT / "bids_adni1"

DCM_IN_FOLDER = re.compile(r"/I(\d+)/[^/]+\.dcm$")


def fs_tag(fs):
    return {"1.5": "15T", "3.0": "3T"}.get(fs, fs.replace(".", "p") + "T")


def select():
    """One best-ranked T1 per (subject, visit, field strength)."""
    rank = {}
    for r in csv.DictReader(open(MRIMPRANK)):
        s = r["LONIUID"].lstrip("S")
        v = int(r["RANK"]); v = 999 if v < 0 else v
        if s not in rank or v < rank[s]:
            rank[s] = v
    mq = [r for r in csv.DictReader(open(MRIQC))
          if r["MRIProtocolPhase"] == "ADNI1" and r["SeriesType"] == "T1w"]
    grp = defaultdict(list)
    for r in mq:
        r["_rk"] = rank.get(r["LONISeries"], 999)
        grp[(r["ParticipantID"], r["VISCODE2"], r["MagneticFieldStrength"])].append(r)
    return [min(v, key=lambda r: (r["_rk"], int(r["image_id"]))) for v in grp.values()]


def locate():
    """image_id -> (zip_path, [member names of its DICOM folder])."""
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
    ap.add_argument("-n", type=int, default=0, help="limit to first N scans (0 = all)")
    ap.add_argument("--subject", default=None, help="only this subject (e.g. XXX_S_XXXX), for testing")
    args = ap.parse_args()

    sel = select()
    loc = locate()

    rows = []
    for r in sel:
        iid = r["image_id"]
        l = loc.get(iid)
        subj, vis = r["ParticipantID"], r["VISCODE2"]
        rows.append({
            "image_id": iid, "subject": subj, "visit": vis,
            "field_strength": r["MagneticFieldStrength"], "series": r["LONISeries"],
            "rank": r["_rk"] if r["_rk"] != 999 else "",
            "description": r["SeriesDescription"],
            "zip": Path(l[0]).name if l else "MISSING",
            "n_dcm": len(l[1]) if l else 0,
            "nifti": f"sub-{subj.replace('_', '')}_ses-{vis}_acq-{fs_tag(r['MagneticFieldStrength'])}_T1w",
        })
    with open(OUT_MANIFEST, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    miss = sum(1 for r in rows if r["zip"] == "MISSING")
    print(f"manifest: {len(rows)} scans -> {OUT_MANIFEST}")
    print(f"  on disk: {len(rows) - miss} | missing from zips: {miss}")

    if not args.bids:
        return
    if not shutil.which("dcm2bids"):
        print("ERROR: dcm2bids not on PATH (module load dcm2bids first).")
        return

    # one dcm2bids call per scan; field strength picks the config (acq-15T / acq-3T).
    # dual-field-strength sessions get two calls into the same ses- (different acq).
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
            # acq (field strength) is ours, not the sidecar (MagneticFieldStrength unreliable:
            # 1.5 / 1.494 / 15000 Gauss). Two rules: equalized "_Eq" twin -> rec-eq, else canonical.
            acq = fs_tag(r["field_strength"])
            cfg = tmp / "cfg.json"
            cfg.write_text(json.dumps({"descriptions": [
                {"datatype": "anat", "suffix": "T1w", "custom_entities": [f"acq-{acq}", "rec-eq"],
                 "criteria": {"SidecarFilename": "*_Eq*"}},
                {"datatype": "anat", "suffix": "T1w", "custom_entities": [f"acq-{acq}"],
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
    print("validate: any unmatched scans land in", BIDS_DIR / "tmp_dcm2bids")


if __name__ == "__main__":
    main()
