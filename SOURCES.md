# Sources and data notes

## GEO accession

- GSE165081 record: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165081
- Official supplementary directory: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165081/suppl/

## Files used

- `GSE165081_M.csv.gz` — processed methylation M-value matrix.
- `GSE165081_Sample_Sheet.PD-CN.csv.gz` — sample metadata.

The pipeline does not assume actual headers, sample counts, IDs, labels, or pair IDs. It inspects the data at run time and stops on ambiguous structure or invalid sample matching.

## Processed and raw data

The workflow uses the processed M-value matrix above. `GSE165081_RAW.tar` contains raw IDAT data and is deliberately not downloaded, extracted, or processed by this implementation.

## Design and interpretation

The GEO description mentions matched controls. This does not establish that the sample metadata contains verified one-to-one pairs. Paired analysis is permitted only with a supplied pair-ID column that passes validation. Independent analysis is explicit and exploratory. Neither analysis establishes causality or diagnostic utility.

## Citation

Use the citation and linked publication on the GEO record for the academic report. Bibliographic details should be verified from GEO; none are invented here.
