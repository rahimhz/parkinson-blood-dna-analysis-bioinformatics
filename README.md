# Parkinson's Disease DNA Methylation Analysis

Beginner-friendly analysis of the processed DNA methylation M-value matrix for GEO accession **GSE165081**. The project validates the GEO metadata, compares Parkinson's disease (PD) and Control samples when explicitly requested, and writes per-run outputs.

The analysis does not establish causality, diagnostic performance, or clinical utility.

## Data and scope

The pipeline downloads the official GEO supplementary files:

- `GSE165081_M.csv.gz`: processed methylation M-value measurements.
- `GSE165081_Sample_Sheet.PD-CN.csv.gz`: sample metadata.

`GSE165081_RAW.tar` contains raw IDAT data. **This project does not download, extract, or process raw IDAT files.** The current workflow starts from processed M-values; raw IDAT preprocessing would require a separate platform-specific workflow.

No sample counts, headers, sample IDs, group labels, or pair IDs are assumed. The code inspects downloaded files and stops with an error when it cannot identify or validate the actual structure.

## Statistical design

The GEO study description mentions matched controls, but that wording alone does not verify one-to-one pairs.

- Default `auto`: validates inputs and sample metadata, writes `validation_summary.json`, archives it, and exits successfully without statistical testing or statistical result files.
- `-p PAIR_COLUMN`: requests a paired t-test only after validating exactly one PD and one Control per pair ID.
- `-d independent`: explicitly requests an exploratory Welch two-sample t-test, which does not model reported matching. Use only after reviewing the study design.
- `-p` and `-d independent` conflict and are rejected.

Statistical modes apply Benjamini-Hochberg FDR correction across tested features. The simple tests do not adjust for covariates such as age, sex, or batch.

## Project layout

```text
.
├── README.md
├── SOURCES.md
├── LIEN_GIT.txt
├── run.sh
├── cron.example
├── environment.yml
├── data/{raw,processed}/
├── scripts/{analyze.py,run_options.sh}
├── results/runs/
├── logs/
├── report/
└── tests/
```

## Five-member team task allocation

Replace placeholders with the actual team members and revise tasks to reflect real contributions.

| Member | Suggested responsibility |
|---|---|
| `[Member 1 name]` | GEO record, data provenance, assignment requirements |
| `[Member 2 name]` | Inspect and document sample metadata and labels |
| `[Member 3 name]` | Review matrix structure and sample validation |
| `[Member 4 name]` | Review statistical design, tests, and interpretation |
| `[Member 5 name]` | Coordinate reproducibility, figures, and final report |

## Ubuntu installation

```bash
sudo apt update
sudo apt install -y bash wget gzip coreutils tar
```

Install Miniconda or Miniforge if needed, then from the repository directory:

```bash
conda env create -f environment.yml
conda activate parkinson-methylation
chmod +x run.sh
```

## Offline tests

Tests use small synthetic data only and do not download GEO files:

```bash
python -m unittest discover -s tests -v
bash -n run.sh scripts/run_options.sh
```

## Run the pipeline

```bash
./run.sh -h
```

Default validation-only run:

```bash
./run.sh
```

The default run validates gzip integrity, checksums, metadata labels and sample matching, inspects matrix orientation, creates `validation_summary.json`, archives the validation output, and exits without statistical testing.

Paired analysis requires verification of a genuine pair-ID column first:

```bash
./run.sh -p VERIFIED_PAIR_COLUMN
```

Explicit exploratory independent comparison:

```bash
./run.sh -d independent
```

If automatic metadata-column detection is ambiguous, inspect the logged sample-sheet preview and supply actual headers:

```bash
./run.sh -d independent -s ACTUAL_SAMPLE_ID_HEADER -g ACTUAL_GROUP_HEADER
```

Options: `-h` help, `-f` force re-download, `-d auto|independent`, `-s` sample-ID column, `-g` group column, and `-p` verified pair-ID column. A valid `-p` resolves `auto` to `paired`.

Existing valid downloads are reused. `-f` forces a fresh download. Each invocation gets its own unique results directory and archive.

## Outputs

Files are stored under `results/runs/RUN_ID/` and are excluded from Git.

- `validation_summary.json`: generated in auto mode; contains detected structure, validated sample IDs/counts, and validation status.
- `differential_methylation.csv`, `analysis_summary.json`, `figures/methylation_volcano.png`: generated only for statistical modes.
- `results_RUN_ID.tar.gz`: archive of regular output files from this run only.
- `logs/run_RUN_ID.log`: timestamped log.
- `logs/SHA256SUMS_RUN_ID.txt`: checksums of the downloaded compressed files.

A local checksum manifest verifies file consistency against the digest generated during the same run. It is not independent proof of authenticity without a trusted provider digest.

## Cron scheduling

`cron.example` is disabled and validation-only. It does not choose a statistical test. To schedule it, replace the path, uncomment the line, and install it with `crontab -e`. Check with `crontab -l`; remove the line with `crontab -e` to disable. Ensure cron uses the intended Python environment.

## Report deliverable

The assignment requires a `report.pdf` of no more than four pages. It is a future deliverable and is not included or claimed complete. See `report/README.md`. Do not invent results or a commit hash.

## Limitations

The workflow analyzes submitted processed M-values and does not process raw IDAT. Statistical tests do not model covariates; independent testing does not account for matching. No real data results are asserted by this starter project.
