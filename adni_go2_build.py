#!/usr/bin/env python3
"""ADNIGO/2 QC-selected T1 cohort: reproducible selection -> manifest -> dcm2bids (BIDS).

Keeps BOTH the accelerated and non-accelerated T1 per session (ADNI grades them equally,
so we don't choose). Selection (CSV-only): MRIQC `T1w` rows with
`MRIProtocolPhase=ADNIGO/ADNI2`, dropping non-T1 mislabels (FLAIR/T2/...) and grade-4
failures (`MAYOADIRL_MRI_IMAGEQC.series_quality`). One best per
(subject, visit, field-strength, accel-status); accel from IMAGEQC `T1_ACCELERATED`,
falling back to `MRIQC.Acceleration` for scans newer than the 2015 QC file.
BIDS acq = field-strength + accel, e.g. `acq-3Taccel`, `acq-3Tunaccel`.

Usage:
  python3 adni_go2_build.py               # manifest + download list
  python3 adni_go2_build.py --bids        # convert (after raws downloaded)
  python3 adni_go2_build.py --bids -n 3
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
IMAGEQC = ROOT / "MAYOADIRL_MRI_IMAGEQC_05_07_15_02Jun2026.csv"
ZIPS = sorted(set(glob(str(ROOT / "go" / "*.zip")) + glob(str(ROOT / "2" / "*.zip")) +
                  glob(str(ROOT / "*adni_go*.zip")) + glob(str(ROOT / "*adnigo*.zip")) +
                  glob(str(ROOT / "*go_2*.zip")) + glob(str(ROOT / "*go2*.zip")) +
                  glob(str(ROOT / "*adni_2_qc*.zip"))))
OUT_MANIFEST = ROOT / "adni_go2_cohort_manifest.csv"
OUT_IDS = ROOT / "adni_go2_image_ids.txt"
BIDS_DIR = ROOT / "bids_adni_go2"

DCM_IN_FOLDER = re.compile(r"/I(\d+)/[^/]+\.dcm$")
NON_T1 = re.compile(r"flair|t2|dti|dwi|\basl\b|swi|perfus|localiz|scout|field.?map", re.I)
FS = {"1.5": "15T", "3.0": "3T"}
# Some scmri-visit scans IDA will not return by Original search (real T1s per MRIQC but
# unservable), so they are excluded from the cohort. The excluded LONI image IDs are kept
# out of version control (ADNI-derived) in adni_go2_exclude.txt; missing file -> no exclusions.
EXCLUDE_FILE = ROOT / "adni_go2_exclude.txt"
EXCLUDE = ({l.split("#")[0].strip() for l in EXCLUDE_FILE.read_text().splitlines()} - {""}
           if EXCLUDE_FILE.exists() else set())


def select():
    qc = {r["loni_image"].lstrip("I"): r for r in csv.DictReader(open(IMAGEQC))}

    def grade(i):
        g = qc.get(i, {}).get("series_quality", "")
        return int(g) if g not in ("", "-1") else 99

    def accel(i, mq_acc):
        a = qc.get(i, {}).get("T1_ACCELERATED", "")
        if a == "1":
            return "accel"
        if a == "0":
            return "unaccel"
        return {"Accelerated": "accel", "Unaccelerated": "unaccel"}.get(mq_acc, "unaccel")

    mq = [r for r in csv.DictReader(open(MRIQC))
          if r["MRIProtocolPhase"] == "ADNIGO/ADNI2" and r["SeriesType"] == "T1w"
          and not NON_T1.search(r["SeriesDescription"]) and grade(r["image_id"]) != 4
          and r["image_id"] not in EXCLUDE]
    grp = defaultdict(list)
    for r in mq:
        fs = FS[r["MagneticFieldStrength"]]
        a = accel(r["image_id"], r["Acceleration"])
        r["_acq"], r["_grade"], r["_accel"], r["_fs"] = fs + a, grade(r["image_id"]), a, fs
        grp[(r["ParticipantID"], r["VISCODE2"], fs, a)].append(r)
    return [min(v, key=lambda r: (r["_grade"], int(r["image_id"]))) for v in grp.values()]


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
            "image_id": iid, "subject": subj, "visit": vis,
            "field_strength": r["_fs"], "acceleration": r["_accel"], "acq": r["_acq"],
            "series": r["LONISeries"], "description": r["SeriesDescription"],
            "grade": "" if r["_grade"] == 99 else r["_grade"],
            "zip": Path(l[0]).name if l else "MISSING", "n_dcm": len(l[1]) if l else 0,
            "nifti": f"sub-{subj.replace('_', '')}_ses-{vis}_acq-{r['_acq']}_T1w",
        })
    with open(OUT_MANIFEST, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    with open(OUT_IDS, "w") as f:
        f.write(", ".join(sorted("I" + r["image_id"] for r in rows)))
    on = sum(1 for r in rows if r["zip"] != "MISSING")
    from collections import Counter
    print(f"manifest: {len(rows)} scans -> {OUT_MANIFEST}")
    print("  by acq:", dict(Counter(r["acq"] for r in rows)))
    print(f"  on disk: {on} | to download: {len(rows) - on} -> {OUT_IDS}")

    if not args.bids:
        return
    if not shutil.which("dcm2bids"):
        print("ERROR: dcm2bids not on PATH (module load dcm2bids).")
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
            # per-scan config (sidecar field strength unreliable, so acq is ours).
            # two rules: dcm2niix's equalized "_Eq" twin -> rec-eq, everything else -> canonical.
            cfg = tmp / "cfg.json"
            cfg.write_text(json.dumps({"descriptions": [
                {"datatype": "anat", "suffix": "T1w",
                 "custom_entities": [f"acq-{r['acq']}", "rec-eq"],
                 "criteria": {"SidecarFilename": "*_Eq*"}},
                {"datatype": "anat", "suffix": "T1w",
                 "custom_entities": [f"acq-{r['acq']}"],
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
