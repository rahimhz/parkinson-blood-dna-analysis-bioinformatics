# Local GEO data policy

The project downloads the processed GEO inputs automatically when `./run.sh`
is executed from a fresh clone. The source files are stored locally in
`data/raw/`; the cleaned sample-sheet copy is stored in `data/processed/`.

These dataset and derived input files are intentionally excluded from Git by
`.gitignore`, including compressed CSV files, TAR archives, IDAT files, and all
other files placed in `data/raw/` or `data/processed/`. Only the `.gitkeep`
placeholders in those folders are tracked. Do not force-add GEO data files.

The current pipeline uses `GSE165081_M.csv.gz` and
`GSE165081_Sample_Sheet.PD-CN.csv.gz` from the official NCBI GEO supplementary
site. It does not download the raw `GSE165081_RAW.tar` archive.
