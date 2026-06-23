# ADNI T1w to BIDS

Reproducible selection and BIDSification of T1ws of ADNI1, 2/GO/3, and 4.
ADNI 1,2/GO,3 QC determined by ADNIs QC standards.
ADNI 4 QC has not been released.
Two version will be available: all ADNI4 T1ws, and a set QC'd by CoBra Lab member Olivier Parent.
Raw DICOMS were obtained from IDA and converted to bids with [dcm2bids](https://unfmontreal.github.io/Dcm2Bids/).

ADNI-derived data: subject to the ADNI Data Use Agreement, not for redistribution.

## Selection

| Phase | Script              | QC source                         | Rule                                                       | n           |
| ----- | ------------------- | --------------------------------- | ---------------------------------------------------------- | ----------- |
| ADNI1 | `adni1_build.py`    | `MRIMPRANK`                       | min `RANK` per (subject, visit, field strength)            | 3,535       |
| GO/2  | `adni_go2_build.py` | `MAYOADIRL_MRI_IMAGEQC`           | keep accelerated + unaccelerated per session; drop grade 4 | 8,861       |
| ADNI3 | `adni3_build.py`    | `MAYOADIRL_MRI_QUALITY_ADNI3`     | `SERIES_SELECTED=TRUE` per session                         | 2,329       |
| ADNI4 | `adni4_build.py`    | in-lab `ADNI4_QC_t1_cobralab.tsv` | `--bids`: all T1w; `--qc`: grade 0/1                       | 1,470 / 398 |

ADNI4 Mayo QC is unreleased; `qc_olivier` = 0 good, 1 artifacts, 2 fail.

The scripts build the list of images that pass QC.
The output is then pasted into the IDA Advanced Search **by Image-ID**.
This will give a list of images that can be added to collection.
The images and the corresponding <phase>\_cohort_manifest.csv can then be downloaded.

## Requirements

- Python 3.12.3  
- dcm2niix/v1.0.20250506  
- dcm2bids/3.1.1  

## Inputs

QC/metadata sheets (IDA -> Study Data, newest dated version):

| File                                   | Phase | IDA name                                      |
| -------------------------------------- | ----- | --------------------------------------------- |
| `MRIQC_*.csv`                          | all   | MAYO ADIR LAB MRI Quality [ADNI1,GO,2,3,4]    |
| `MRIMPRANK_*.csv`                      | ADNI1 | MRI MPRAGE Ranking [ADNI1]                    |
| `MAYOADIRL_MRI_QUALITY_ADNI3_*.csv`    | ADNI3 | Mayo ADNI3 MRI QC                             |
| `MAYOADIRL_MRI_IMAGEQC_05_07_15_*.csv` | GO/2  | Mayo (Jack Lab) ADNI GO/2 MRI QC (`_Archive`) |
| `ADNI4_QC_t1_cobralab.tsv`             | ADNI4 | in-lab                                        |

Set `ROOT` (top of each script) to the directory holding the sheets and zips.

## Usage

```bash
python3 adni1_build.py            # -> <phase>_cohort_manifest.csv + Image-ID download list
module load dcm2bids
python3 adni1_build.py --bids     # -> bids_adni1/
```

ADNIGO/2 and 3 are run in the same way. ADNI4 has two independent targets:

```bash
python3 adni4_build.py --bids     # bids_adni4_all/  (all T1w)
python3 adni4_build.py --qc       # bids_adni4_qc/   (grade 0/1; standalone)
```

Validate (converted vs missing, per phase):

```bash
python3 validate.py adni1_cohort_manifest.csv bids_adni1
```

## Output

One BIDS dataset per phase (`bids_adni1/`, `bids_adni3/`, `bids_adni_go2/`; ADNI4 yields two,
`bids_adni4_all/` and `bids_adni4_qc/`), each laid out as:

```
bids_<phase>/
├── dataset_description.json
└── sub-<ID>/
    └── ses-<VISIT>/
        └── anat/
            └── sub-<ID>_ses-<VISIT>_acq-<ACQ>[_run-N][_rec-N]_T1w.nii.gz   (+ .json)
```

Entities:
- `acq` — field strength / acceleration: `15T`, `3T`, `3Taccel`, `3Tunaccel`, `3Tultrafast`
- `run-N` — present only when a session has >1 scan of the same `acq`
- `rec-N` — Philips reconstruction variant (see Notes)

## Notes

- `acq` (field strength / acceleration) is taken from `MRIQC`, not the DICOM sidecar
  (`MagneticFieldStrength` is inconsistent: 1.5 / 1.494 / 15000 G); injected per-scan via a
  `Modality=MR` dcm2bids config.
- dcm2niix may emit >1 NIfTI/scan: Philips multi-`ImageType` → `run-01/02`, resolved to
  canonical + `rec-NN` (`resolve_run_splits`); `_Eq`/sidecar-less outputs are left under
  `tmp_dcm2bids/`. The canonical `_T1w` is always written.
- ADNI1 m48 visits scanned in the GO/2 era are filed under the GO/2 project in IDA (search
  unscoped to retrieve); they remain ADNI1 by `MRIProtocolPhase`.

## Acknowledgements

Pipeline code was developed with assistance from large language models, including
locally-hosted (self-contained) LLMs and Claude (Anthropic).
