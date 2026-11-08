# Parkinson’s Disease DNA Methylation Analysis

A reproducible, validation-first bioinformatics project investigating DNA methylation in whole-blood samples from people with Parkinson’s disease (PD) and control individuals.

**GEO accession:** [GSE165081](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165081)

**Research question:** Quels sites de méthylation de l’ADN présentent des différences entre les échantillons sanguins de personnes atteintes de la maladie de Parkinson (PD) et ceux des témoins (Control) ?

This project is intended to be exploratory and reproducible. It does not claim causal inference, clinical diagnosis, or treatment relevance from the methylation signal alone.

> Current status: the workflow defaults to metadata and input validation. It does not claim that the real GEO dataset has been analyzed, that statistical results are final, or that the final PDF report is complete.

## 1. Dataset and scope

This project uses the publicly available GEO accession **GSE165081** and focuses on the processed methylation matrix rather than raw IDAT data.

### Files used

| File | Description | Purpose |
|---|---|---|
| `GSE165081_M.csv.gz` | Processed methylation M-value matrix | Main analysis input |
| `GSE165081_Sample_Sheet.PD-CN.csv.gz` | Sample metadata table | Sample IDs, group labels, and matching information |
| `GSE165081_RAW.tar` | Raw IDAT archive | Not used in this workflow |

The workflow downloads the processed GEO supplementary files from the official NCBI GEO FTP site when needed. It does not download, extract, or process the raw IDAT archive; raw-data preprocessing would require a separate platform-specific workflow.

The downloaded files are stored locally in `data/raw/` and `data/processed/`. The repository’s `.gitignore` is intended to keep dataset files and local generated data out of Git. Do not force-add GEO files or other large data artefacts.

## 2. Analysis design

The pipeline supports a validation-first mode and optional statistical analysis modes.

### Default mode: validation only

Run:

```bash
./run.sh
```

The default execution is intended to:

1. Download or reuse the processed matrix and sample sheet.
2. Validate gzip integrity and generate local checksums.
3. Inspect sample metadata and matrix structure.
4. Verify that sample IDs and group labels are consistent with the matrix.
5. Write a validation summary and archive the run outputs.
6. Exit without claiming statistical findings.

The project does not assume sample counts, headers, labels, or pairing information. If the actual file structure is ambiguous or invalid, the pipeline is designed to stop with an error and require the metadata to be resolved before analysis proceeds.

### Paired analysis

The GEO study description may mention matched controls, but that wording alone does not guarantee valid one-to-one pairing. A pair-ID column must be verified from the actual metadata before running a paired test.

Run paired analysis only after confirming a true pair identifier:

```bash
./run.sh -p VERIFIED_PAIR_COLUMN
```

The paired test should be used only when each pair contains exactly one PD sample and one Control sample, as confirmed by the downloaded sample sheet.

### Independent exploratory analysis

If valid paired information is not confirmed, an exploratory unpaired comparison can be requested explicitly:

```bash
./run.sh -d independent
```

This mode is intended for exploratory testing only. It does not model the reported matching in the study design and must be interpreted with that limitation in mind.

### Multiple-testing correction

When statistical mode is enabled, the analysis is intended to apply the Benjamini–Hochberg false-discovery-rate (FDR) adjustment across tested features.

The current statistical workflow does not account for potential confounding variables such as age, sex, or batch effects.

## 3. Command-line options

```text
-h  Display help and exit.
-f  Force re-download of the GEO supplementary files.
-d  Select design mode: auto (default) or independent.
-s  Specify the sample-ID column when metadata detection is ambiguous.
-g  Specify the group-label column when detection is ambiguous.
-p  Specify a verified one-to-one pair-ID column for paired analysis.
```

`-p` and `-d independent` are mutually exclusive and should not be combined.

## 4. Installation on Ubuntu

Install the required system packages:

```bash
sudo apt update
sudo apt install -y git bash wget gzip coreutils tar
```

Clone the project and create the conda environment:

```bash
git clone https://github.com/rahimhz/parkinson-blood-dna-analysis-bioinformatics.git
cd parkinson-blood-dna-analysis-bioinformatics

conda env create -f environment.yml
conda activate parkinson-methylation

chmod +x run.sh
```

The pipeline requires internet access to download the GEO supplementary files when they are not already present locally.

## 5. Running tests

The repository includes synthetic tests designed to validate analysis logic without downloading the full GEO dataset.

```bash
python -m unittest discover -s tests -v
bash -n run.sh scripts/run_options.sh
```

A passing synthetic test run does not validate the structure of the real GEO files or replace a full pipeline run on the actual dataset.

## 6. Repository layout

This describes the current repository structure:

```text
parkinson-blood-dna-analysis-bioinformatics/
├── .gitignore
├── DATA_POLICY.md
├── LIEN_GIT.txt
├── README.md
├── SOURCES.md
├── cron.example
├── environment.yml
├── data/
│   ├── raw/.gitkeep
│   └── processed/.gitkeep
├── logs/.gitkeep
├── report/
│   └── README.md
├── results/
│   └── runs/.gitkeep
├── run.sh
├── scripts/
│   ├── __init__.py
│   ├── analyze.py
│   └── run_options.sh
└── tests/
    ├── __init__.py
    ├── test_analysis.py
    └── test_run_options.py
```

The `.gitkeep` files are placeholders. Data downloaded from GEO, cleaned metadata, logs, and per-run analysis outputs are created locally and are intended to remain outside Git.

## 7. Generated outputs

Each execution creates a unique directory under:

```text
results/runs/RUN_ID/
```

In validation mode, expected outputs include:

- `validation_summary.json`: validation status and detected metadata/matrix information.
- `results_RUN_ID.tar.gz`: archive of the run outputs.

In statistical modes, additional outputs may include:

- `differential_methylation.csv`: per-feature statistical results and adjusted p-values.
- `analysis_summary.json`: summary of the statistical design and findings.
- `figures/methylation_volcano.png`: volcano plot for methylation effect versus significance.

Execution logs and checksum manifests are stored under `logs/`. These are local run artefacts and should not be published to GitHub.

## 8. Data protection

The project is intended to download the required GEO files automatically when `./run.sh` is executed on Ubuntu. The downloaded files are local-only and should not be added to the repository.

This is enforced through `.gitignore` and documented in `DATA_POLICY.md`.

Do not run `git add -f` on GEO data or other large dataset files.

## 9. Cron scheduling

`cron.example` is a disabled validation-only template. It does not perform statistical testing and should only be enabled after verifying the local path, Python environment, and intended execution mode.

Before activating it, check the repository path and your conda environment. Do not schedule a statistical comparison until the actual metadata and design have been reviewed.

## 10. Team task distribution

The project is divided into five complementary responsibilities. Each member is responsible for a specific part of the bioinformatics workflow, including implementation, verification, documentation, and presentation of their contribution.

| **Member** | **Main responsibility** | **Key files / areas** | **Expected deliverables** |
|---|---|---|---|
| **Person 1** | Dataset acquisition and scientific background: study the GEO record, identify the supplementary files, implement or verify data downloading, document data sources, and explain the biological context of Parkinson's disease and DNA methylation. | `SOURCES.md`, `DATA_POLICY.md`, `data/`, GEO documentation | Verified dataset sources, documented download procedure, checksums, and scientific background |
| **Person 2** | Data inspection and validation: inspect the downloaded matrix and sample metadata, identify sample and group columns, verify PD/Control labels, check sample matching, and document data-quality issues. | `scripts/analyze.py`, `data/`, validation outputs | Verified metadata structure, sample matching report, and validation results |
| **Person 3** | Bash automation and pipeline management: implement and maintain download orchestration, command-line options, error handling, logging, file permissions, run directories, and result archiving. | `run.sh`, `scripts/run_options.sh`, `cron.example`, `.gitignore` | Functional and reproducible execution pipeline with logging and automated file management |
| **Person 4** | Statistical analysis and visualization: implement and review the PD/Control comparison, justify the statistical test based on the verified study design, apply multiple-testing correction, and generate figures and numerical results. | `scripts/analyze.py`, `results/runs/`, statistical methods | Statistical analysis, adjusted p-values, generated figures, and scientific interpretation |
| **Person 5** | Testing, reproducibility, and final reporting: develop and execute synthetic tests, perform end-to-end pipeline checks, verify generated outputs, organize the final report and presentation, and document reproducibility. | `tests/`, `report/`, `results/runs/`, project documentation | Automated test results, pipeline verification, final scientific report, and presentation materials |

### Shared responsibilities

- Every member must understand the complete workflow and be able to explain the code during the oral defense.
- Each member must document their contribution and help review another member's work.
- The team must jointly review the final results, scientific interpretation, and submitted files.
- Responsibilities may overlap during integration and debugging, but each member remains accountable for their assigned deliverables.
## 11. Final academic deliverables

The assignment requires a final scientific report in PDF format of no more than four pages. `report/README.md` describes the intended contents; it is not itself a completed report.

Do not claim completion of the final PDF report, ZIP archive, or final commit hash until each has been generated and checked.

The final submission may additionally require a ZIP archive containing the appropriate project files while excluding the original GEO dataset and other large local data files.

`LIEN_GIT.txt` is intended to hold the repository URL and the final commit hash. Fill in the commit hash only after the final commit actually exists.

## 12. Limitations

- The workflow operates on processed M-values and does not process raw IDAT data.
- Actual GEO headers, sample IDs, group labels, and pair information must be confirmed from the downloaded metadata.
- A paired test is appropriate only when the metadata and study design verify one-to-one pairs.
- An independent Welch test does not model the matching implied by the study description.
- The current statistical framework does not account for age, sex, batch, or other confounding variables.
- Results are exploratory associations and should not be interpreted as causal findings or clinical diagnostic claims.

## 13. Sources

- [GSE165081 – NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165081)
- [Official GEO supplementary directory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165081/suppl/)

