#!/usr/bin/env python3
"""ADNI4 T1 -> BIDS, two sets:
  bids_adni4_all/  — every ADNI4 T1w (no QC filter)
  bids_adni4_qc/   — symlinked subset that passed in-lab cobralab QC (grade 0 or 1)

ADNI4's Mayo MRI QC isn't released, so selection uses the in-lab `ADNI4_QC_t1_cobralab.tsv`
(qc_olivier: 0=good, 1=some artifacts, 2=fail), joined per-scan on subject+date+sequence.
Spine = MRIQC ADNI4 T1w (all 3T). acq = acceleration (Accelerated/Unaccelerated/Ultrafast=CS);
run- only where a session has >1 of the same acceleration (~4 cases).

Usage:
  python3 adni4_build.py                 # manifest + download list
  python3 adni4_build.py --bids          # convert all -> bids_adni4_all/ + qc symlinks
  python3 adni4_build.py --bids -n 3
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
COBRALAB = ROOT / "ADNI4_QC_t1_cobralab.tsv"
ZIPS = sorted(set(glob(str(ROOT / "4" / "*.zip")) + glob(str(ROOT / "*adni4*.zip")) +
                  glob(str(ROOT / "*adni_4*.zip")) + glob(str(ROOT / "*ADNI4*.zip"))))
OUT_MANIFEST = ROOT / "adni4_cohort_manifest.csv"
OUT_IDS = ROOT / "adni4_image_ids.txt"
BIDS_ALL = ROOT / "bids_adni4_all"
BIDS_QC = ROOT / "bids_adni4_qc"

DCM_IN_FOLDER = re.compile(r"/I(\d+)/[^/]+\.dcm$")
NON_T1 = re.compile(r"flair|t2|dti|dwi|\basl\b|swi|perfus|localiz|scout|field.?map", re.I)
ACQ = {"Accelerated": "3Taccel", "Unaccelerated": "3Tunaccel", "Ultrafast": "3Tultrafast"}
PASS = {"0", "1"}  # cobralab grades kept in the QC set


def norm(s):
    return re.sub(r"[^a-z0-9]", "", s.lower())


def load_grades():
    """ (subject, date, normseq) -> grade ; plus (subject, date) -> set(grades) fallback. """
    exact, byday = {}, defaultdict(set)
    for r in csv.DictReader(open(COBRALAB), delimiter="\t"):
        exact[(r["ID"], r["date"], norm(r["mri_seq"]))] = r["qc_olivier"]
        byday[(r["ID"], r["date"])].add(r["qc_olivier"])
    return exact, byday


def select():
    exact, byday = load_grades()
    mq = [r for r in csv.DictReader(open(MRIQC))
          if r["MRIProtocolPhase"] == "ADNI4" and r["SeriesType"] == "T1w"
          and not NON_T1.search(r["SeriesDescription"]) and r["VISCODE2"].strip()]
    # attach acq + grade
    for r in mq:
        r["_acq"] = ACQ.get(r["Acceleration"], "3T" + r["Acceleration"].lower())
        d = r["StudyDate"].split()[0]
        g = exact.get((r["ParticipantID"], d, norm(r["SeriesDescription"])))
        if g is None:
            gs = byday.get((r["ParticipantID"], d), set())
            g = gs.pop() if len(gs) == 1 else ""   # fall back only if the day is unambiguous
        r["_grade"] = g
    # run- only where a (subject, visit, acq) group has >1
    grp = defaultdict(list)
    for r in mq:
        grp[(r["ParticipantID"], r["VISCODE2"], r["_acq"])].append(r)
    out = []
    for v in grp.values():
        multi = len(v) > 1
        for i, r in enumerate(sorted(v, key=lambda x: int(x["image_id"])), 1):
            r["_run"] = f"{i:02d}" if multi else ""
            out.append(r)
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
    """Philips multi-ImageType split -> run-01/run-02 instead of <base>: run-01 -> canonical,
    run-NN -> rec-NN. (For scans already carrying a run- entity this is a no-op.)"""
    if (anat / f"{base}.nii.gz").exists():
        return
    prefix = base[:-4]
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


def convert(todo, loc, out_dir, label):
    """dcm2bids-convert the given scans into out_dir, independently (each: extract -> dcm2bids
    -> resolve Philips run-splits). Used for both the all-set and the QC subset."""
    if not shutil.which("dcm2bids"):
        print("ERROR: dcm2bids not on PATH (module load dcm2bids).")
        return
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"dcm2bids: {len(todo)} scans -> {out_dir}  ({label})")
    for i, r in enumerate(todo, 1):
        z, names = loc[r["image_id"]]
        tmp = Path(tempfile.mkdtemp())
        try:
            with zipfile.ZipFile(z) as zf:
                for n in names:
                    zf.extract(n, tmp)
            ce = [f"acq-{r['acq']}"] + ([f"run-{r['run']}"] if r["run"] else [])
            cfg = tmp / "cfg.json"
            cfg.write_text(json.dumps({"descriptions": [
                {"datatype": "anat", "suffix": "T1w", "custom_entities": ce,
                 "criteria": {"Modality": "MR"}}]}))
            subprocess.run(["dcm2bids", "-d", str(tmp),
                            "-p", r["subject"].replace("_", ""), "-s", r["visit"],
                            "-c", str(cfg), "-o", str(out_dir), "--clobber"],
                           check=False, capture_output=True)
            anat = out_dir / f"sub-{r['subject'].replace('_', '')}" / f"ses-{r['visit']}" / "anat"
            resolve_run_splits(anat, r["nifti"])
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        if i % 200 == 0:
            print(f"  {i}/{len(todo)}")
    print(f"done -> {out_dir}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bids", action="store_true", help="convert all T1w -> bids_adni4_all/")
    ap.add_argument("--qc", action="store_true", help="copy grade-0/1 subset -> bids_adni4_qc/")
    ap.add_argument("-n", type=int, default=0)
    ap.add_argument("--subject", default=None)
    args = ap.parse_args()

    sel = select()
    loc = locate()
    rows = []
    for r in sel:
        iid = r["image_id"]
        l = loc.get(iid)
        subj, vis, acq, run = r["ParticipantID"], r["VISCODE2"], r["_acq"], r["_run"]
        ents = f"_run-{run}" if run else ""
        rows.append({
            "image_id": iid, "subject": subj, "visit": vis, "acq": acq, "run": run,
            "description": r["SeriesDescription"], "qc_grade": r["_grade"],
            "qc_pass": r["_grade"] in PASS,
            "zip": Path(l[0]).name if l else "MISSING", "n_dcm": len(l[1]) if l else 0,
            "nifti": f"sub-{subj.replace('_', '')}_ses-{vis}_acq-{acq}{ents}_T1w",
        })
    with open(OUT_MANIFEST, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    with open(OUT_IDS, "w") as f:
        f.write(", ".join(sorted("I" + r["image_id"] for r in rows)))
    from collections import Counter
    on = sum(1 for r in rows if r["zip"] != "MISSING")
    print(f"manifest: {len(rows)} ADNI4 T1w -> {OUT_MANIFEST}")
    print("  grade:", dict(Counter(r["qc_grade"] or "ungraded" for r in rows)))
    print(f"  QC-pass (grade 0/1): {sum(1 for r in rows if r['qc_pass'])}")
    print(f"  on disk: {on} | to download: {len(rows) - on} -> {OUT_IDS}")

    if not (args.bids or args.qc):
        return
    todo = [r for r in rows if r["zip"] != "MISSING"]
    if args.subject:
        todo = [r for r in todo if r["subject"] == args.subject]
    if args.n:
        todo = todo[:args.n]
    # each flag converts INDEPENDENTLY into its own dir, straight from the zips
    if args.bids:
        convert(todo, loc, BIDS_ALL, "all")
    if args.qc:
        convert([r for r in todo if r["qc_pass"]], loc, BIDS_QC, "QC grade-0/1")


if __name__ == "__main__":
    main()
